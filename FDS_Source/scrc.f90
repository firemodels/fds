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
PUBLIC SCARC_STORAGE  
PUBLIC SCARC_STENCIL  
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
CHARACTER(10) :: SCARC_METHOD     = 'null'                 ! requested solver method (KRYLOV/MULTIGRID)
CHARACTER(10) :: SCARC_STORAGE    = 'BANDWISE'             ! matrix storage technique (BANDWISE/COMPACT)
CHARACTER(10) :: SCARC_STENCIL    = 'NPOINT'               ! underlying matrix stencil (NPOINT/AGGRESSIVE)
CHARACTER(10) :: SCARC_INITIAL    = 'NONE'                 ! initial solution (currently only default is used)

!!! General iteration parameters
REAL (EB)     :: SCARC_RESIDUAL             = -1.0_EB      ! residual of global selected solver
INTEGER       :: SCARC_ITERATIONS           =  0           ! number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_CAPPA                =  1.0_EB      ! convergence rate of selected ScarC solver
REAL (EB)     :: SCARC_ACCURACY_DIVERGENCE  = 1.E+6_EB     ! divergence epsilon for all solvers
REAL (EB)     :: SCARC_ACCURACY_RELATIVE    = 1.E-2_EB     ! minimum relative accuracy for all solvers
CHARACTER(10) :: SCARC_ACCURACY             = 'ABSOLUTE'   ! accuracy type (ABSOLUTE/RELATIVE)

!!! Parameters for multigrid-type methods
CHARACTER(10) :: SCARC_MULTIGRID            = 'GEOMETRIC'  ! type of MG-method (GEOMETRIC/ALGEBRAIC)
INTEGER       :: SCARC_MULTIGRID_LEVEL      = -1           ! User defined number of MG-levels (optionally)
CHARACTER(1)  :: SCARC_MULTIGRID_CYCLE      = 'V'          ! Cycling type  (F/V/W)
CHARACTER(10) :: SCARC_MULTIGRID_COARSENING = 'STANDARD'   ! Coarsening strategy  (STANDARD/AGGRESSIVE)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS = 1000         ! max number of iterations
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY   = 1.E-15_EB    ! requested accuracy for convergence

!!! Parameters for Krylov-type methods
CHARACTER(10) :: SCARC_KRYLOV            = 'CG'            ! type of Krylov-method (CG/BICG)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 1000            ! max number of iterations 
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-15_EB       ! requested accuracy for convergence

!!! Parameters for smoothing method (used in multigrids-methods)
CHARACTER(10) :: SCARC_SMOOTH            = 'JACOBI'        ! smoother for MG (JACOBI/SSOR/GSTRIX)
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 1000            ! max number of iterations 
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-15_EB       ! requested accuracy for convergence
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.90E+0_EB      ! relaxation parameter 

!!! Parameters for preconditioning method (used in Krylov-methods)
CHARACTER(10) :: SCARC_PRECON            = 'JACOBI'        ! preconditioner for CG/BICG (JACOBI/SSOR/GSTRIX/MG)
INTEGER       :: SCARC_PRECON_ITERATIONS = 1000            ! max number of iterations 
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-15_EB       ! requested accuracy for convergence
REAL (EB)     :: SCARC_PRECON_OMEGA      = 0.90E+0_EB      ! relaxation parameter 

!!! Parameters for coarse grid method
CHARACTER(10) :: SCARC_COARSE            = 'CG'            ! coarse grid solver (CG/Gaussian elimination)
INTEGER       :: SCARC_COARSE_ITERATIONS = 1000            ! max number of iterations 
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-15_EB       ! requested accuracy for convergence
REAL (EB)     :: SCARC_COARSE_OMEGA      = 1.5E+0_EB       ! relaxation parameter 
CHARACTER(10) :: SCARC_COARSE_PRECON     = 'SSOR'          ! preconditioner
 
!!! debugging parameters
CHARACTER(10) :: SCARC_DEBUG = 'NONE'                      ! debugging level (NONE/LESS/MEDIUM/MUCH)
CHARACTER(40) :: SCARC_FN                                  ! file name for ScaRC debug messages
INTEGER       :: SCARC_LU                                  ! unit number for ScaRC debug file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!----------------------------------------------------------------------------------------------------
!!! Global constants
!!!----------------------------------------------------------------------------------------------------
INTEGER, PARAMETER :: NSCARC_SCOPE_NONE            = -1, &         
                      NSCARC_SCOPE_MAIN            =  1, &      ! method used as main solver
                      NSCARC_SCOPE_SMOOTH          =  2, &      ! method used as smoother
                      NSCARC_SCOPE_PRECON          =  3, &      ! method used as preconditiner
                      NSCARC_SCOPE_COARSE          =  4     ! method used as coarse grid solver

INTEGER, PARAMETER :: NSCARC_METHOD_NONE           = -1, &
                      NSCARC_METHOD_KRYLOV         =  1, &      ! Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID      =  2         ! multigrid-method as global solver

INTEGER, PARAMETER :: NSCARC_KRYLOV_NONE           = -1, &
                      NSCARC_KRYLOV_CG             =  1, &      ! CG   as Krylov solver
                      NSCARC_KRYLOV_BICG           =  2         ! BICG as Krylov solver

INTEGER, PARAMETER :: NSCARC_MULTIGRID_NONE        = -1, &
                      NSCARC_MULTIGRID_GEOMETRIC   =  1, &      ! geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC   =  2         ! algebraic multigrid

INTEGER, PARAMETER :: NSCARC_STORAGE_NONE          = -1, &
                      NSCARC_STORAGE_BANDWISE      =  1, &      ! bandwise storage technique
                      NSCARC_STORAGE_COMPACT       =  2         ! compact  storage technique

INTEGER, PARAMETER :: NSCARC_STENCIL_NONE          = -1, &
                      NSCARC_STENCIL_NPOINT        =  1, &      ! standard n-point stencil (2D: n=5, 3D: n=7)
                      NSCARC_STENCIL_AGGRESSIVE    =  2         ! stencil resulting from aggressive coarsening

INTEGER, PARAMETER :: NSCARC_EXCHANGE_NONE         = -1, &
                      NSCARC_EXCHANGE_INIT         =  1, &      ! initialize communication
                      NSCARC_EXCHANGE_MATVEC       =  2, &      ! matrix-vector communication 
                      NSCARC_EXCHANGE_BDRY         =  3         ! exchange internal boundaries

INTEGER, PARAMETER :: NSCARC_SMOOTH_NONE           = -1, &
                      NSCARC_SMOOTH_JACOBI         =  1, &      ! smoothing by JACOBI-method
                      NSCARC_SMOOTH_SSOR           =  2, &      ! smoothing by SSOR-method
                      NSCARC_SMOOTH_GSTRIX         =  3, &      ! smoothing by GSTRIX-method
                      NSCARC_SMOOTH_FFT            =  4         ! smoothing by FFT-method

INTEGER, PARAMETER :: NSCARC_PRECON_NONE           = -1, &
                      NSCARC_PRECON_JACOBI         =  1, &      ! preconditioning by JACOBI-method
                      NSCARC_PRECON_SSOR           =  2, &      ! preconditioning by SSOR-method
                      NSCARC_PRECON_GSTRIX         =  3, &      ! preconditioning by GSTRIX-method
                      NSCARC_PRECON_FFT            =  4, &      ! preconditioning by FFT-method
                      NSCARC_PRECON_MG             =  5         ! preconditioning by MG-method

INTEGER, PARAMETER :: NSCARC_CYCLE_NONE            = -1, &
                      NSCARC_CYCLE_F               =  0, &      ! F-cycle for mg-method
                      NSCARC_CYCLE_V               =  1, &      ! V-cycle for mg-method
                      NSCARC_CYCLE_W               =  2, &      ! W-cycle for mg-method
                      NSCARC_CYCLE_INIT            =  3, &      ! initialize cycle counts
                      NSCARC_CYCLE_RESET           =  4, &      ! reset cycle counts
                      NSCARC_CYCLE_PROCEED         =  5, &      ! proceed cycle counts
                      NSCARC_CYCLE_PRESMOOTH       =  6, &      ! presmoothing cycle
                      NSCARC_CYCLE_POSTSMOOTH      =  7, &      ! postsmoothing cycle
                      NSCARC_CYCLE_NEXT            =  8, &      ! perform next cycling loop
                      NSCARC_CYCLE_EXIT            =  9         ! exit cycling loop

INTEGER, PARAMETER :: NSCARC_LOOP_PROCEED          =  0, &      ! no convergence and no divergence
                      NSCARC_LOOP_CONV             =  1, &      ! convergence
                      NSCARC_LOOP_DIVG             =  2         ! divergence

INTEGER, PARAMETER :: NSCARC_DEBUG_NONE            = -1, &      ! no debugging requested
                      NSCARC_DEBUG_LESS            =  1, &      ! low    level of debugging requested
                      NSCARC_DEBUG_MEDIUM          =  2, &      ! medium level of debugging requested
                      NSCARC_DEBUG_MUCH            =  3         ! strong level of debugging requested

INTEGER, PARAMETER :: NSCARC_COARSENING_NONE       = -1, &
                      NSCARC_COARSENING_STANDARD   =  1, &      ! standard   coarsening strategy
                      NSCARC_COARSENING_AGGRESSIVE =  2         ! aggressive coarsening strategy

INTEGER, PARAMETER :: NSCARC_COARSE_NONE           = -1, &
                      NSCARC_COARSE_CG             =  1, &      ! coarse grid solution by cg-method
                      NSCARC_COARSE_GE             =  2         ! coarse grid solution by Gaussian elimination

INTEGER, PARAMETER :: NSCARC_VECTOR_NONE           = -1, &
                      NSCARC_VECTOR_X              =  1, &      ! selection parameter for vector X
                      NSCARC_VECTOR_Y              =  2, &      ! selection parameter for vector Y
                      NSCARC_VECTOR_G              =  3, &      ! selection parameter for vector G
                      NSCARC_VECTOR_R              =  4, &      ! selection parameter for vector R
                      NSCARC_VECTOR_D              =  5, &      ! selection parameter for vector D
                      NSCARC_VECTOR_Z              =  6, &      ! selection parameter for vector Z
                      NSCARC_VECTOR_X2             =  7, &      ! selection parameter for vector X2
                      NSCARC_VECTOR_D2             =  8, &      ! selection parameter for vector D2
                      NSCARC_VECTOR_R2             =  9, &      ! selection parameter for vector R2
                      NSCARC_VECTOR_Y2             = 10         ! selection parameter for vector Y2

INTEGER, PARAMETER :: NSCARC_ACCURACY_NONE         = -1, &
                      NSCARC_ACCURACY_ABSOLUTE     =  1, &      ! absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE     =  2         ! relative accuracy must be reached

INTEGER, PARAMETER :: NSCARC_TIME_NONE             = -1, &
                      NSCARC_TIME_COMPLETE         =  1, &      ! time for complete ScaRC part of FDS
                      NSCARC_TIME_SETUP            =  2, &      ! time for setup phase
                      NSCARC_TIME_SOLVER           =  3, &      ! time for ScaRC solver 
                      NSCARC_TIME_KRYLOV           =  4, &      ! time for Krylov solver
                      NSCARC_TIME_MULTIGRID        =  5, &      ! time for multigrid solver
                      NSCARC_TIME_PRECON           =  6, &      ! time for preconditioner
                      NSCARC_TIME_SMOOTH           =  7, &      ! time for smoother
                      NSCARC_TIME_COARSE           =  8, &      ! time for coarse grid solver
                      NSCARC_TIME_LOCAL_MATVEC     =  9, &      ! time for local matrix-vector product
                      NSCARC_TIME_GLOBAL_MATVEC    = 10, &      ! time for global matrix-vector product
                      NSCARC_TIME_LOCAL_SCALPROD   = 11, &      ! time for local scalar product
                      NSCARC_TIME_GLOBAL_SCALPROD  = 12, &      ! time for global scalar product
                      NSCARC_TIME_GLOBAL_L2NORM    = 13, &      ! time for global scalar product
                      NSCARC_TIME_EXCHANGE_BDRY    = 14         ! time for exchange of internal boundary

INTEGER, PARAMETER :: NSCARC_INITIAL_NONE          = -1         ! no predefined initial solution used

INTEGER, PARAMETER :: NSCARC_LEVEL_MAX             =  8         ! max number of levels currently allowed
 
INTEGER, PARAMETER :: NSCARC_DUMMY                 = -1         ! dummy variable (needed at several places)


!!!----------------------------------------------------------------------------------------------------
!!! Global variables 
!!!----------------------------------------------------------------------------------------------------
!!! use integer types for the user defined input data (based on SCARC_TYPE_... variables)
INTEGER :: TYPE_SCOPE      = NSCARC_SCOPE_NONE
INTEGER :: TYPE_METHOD     = NSCARC_METHOD_NONE
INTEGER :: TYPE_KRYLOV     = NSCARC_KRYLOV_NONE
INTEGER :: TYPE_MULTIGRID  = NSCARC_MULTIGRID_NONE
INTEGER :: TYPE_STORAGE    = NSCARC_STORAGE_NONE
INTEGER :: TYPE_STENCIL    = NSCARC_STENCIL_NONE
INTEGER :: TYPE_ACCURACY   = NSCARC_ACCURACY_NONE
INTEGER :: TYPE_SMOOTH     = NSCARC_SMOOTH_NONE
INTEGER :: TYPE_PRECON     = NSCARC_PRECON_NONE
INTEGER :: TYPE_CYCLE      = NSCARC_CYCLE_NONE
INTEGER :: TYPE_COARSENING = NSCARC_COARSENING_NONE
INTEGER :: TYPE_COARSE     = NSCARC_COARSE_NONE
INTEGER :: TYPE_DEBUG      = NSCARC_DEBUG_NONE   
INTEGER :: TYPE_INITIAL    = NSCARC_INITIAL_NONE

!!! range of meshes which must be processed for MYID
INTEGER :: NMESHES_MIN, NMESHES_MAX                 
 
!!! total, minimum and maximum number of multigrid levels
INTEGER :: NLEVEL, NLEVEL_MIN, NLEVEL_MAX                              

!!! additional arrays for data exchange
INTEGER :: NREQ_SCARC, NCOM_SCARC, TAG_SCARC, SNODE, RNODE, STATUS2_SCARC(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION (:)    :: REQ_SCARC

!!! time measurements with ScaRC
INTEGER, PARAMETER :: N_TIMERS_SCARC=14         
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TUSED_SCARC

!!! auxiliary variables for global and local scalproducts and number of cells
REAL(EB), ALLOCATABLE, DIMENSION (:) :: SP_LOCAL
INTEGER,  ALLOCATABLE, DIMENSION (:) :: NC_GLOBAL, NC_LOCAL

!!! number of couplings in given matrix stencil and pointer to indices in matrix stencil on given level 
INTEGER, POINTER :: NCPL, ID, ILX, ILY, ILZ, IUX, IUY, IUZ

 
!!! Private type declarations

!!!----------------------------------------------------------------------------------------------------
!!! SCARC type for information on a single grid level of own mesh
!!!----------------------------------------------------------------------------------------------------
TYPE SCARC_TYPE_LEVEL
 
!!! local numbers of cells (also per direction), wall cells, matrix entries
INTEGER :: IBAR, JBAR, KBAR
INTEGER :: N_CELLS, N_MATRIX_ENTRIES, N_EXTERNAL_WALL_CELLS

!!! mesh coordinates and step sizes
REAL (EB) :: DX,    DY,    DZ
REAL (EB) :: DXI,   DYI,   DZI
REAL (EB) :: DXI2,  DYI2,  DZI2
REAL (EB) :: DXMIN, DYMIN, DZMIN
REAL (EB) :: DX_MIN, DY_MIN, DZ_MIN
REAL (EB) :: DX_MAX, DY_MAX, DZ_MAX
REAL (EB) :: XX_MIN, YY_MIN, ZZ_MIN
REAL (EB) :: XX_MAX, YY_MAX, ZZ_MAX
REAL (EB), POINTER, DIMENSION (:) :: XX, YY, ZZ

!!! neighbourship structures 
INTEGER, POINTER, DIMENSION (:, :) :: IJKW
INTEGER, POINTER, DIMENSION (:)    :: PRESSURE_BC_INDEX

!!! matrix information: stencils, number of couplings, positions in matrix stencil, subdiagonal entries
INTEGER, POINTER, DIMENSION (:) :: STENCIL_STD, STENCIL_AGG
INTEGER :: N_COUPLED, N_COUPLED_STD, N_COUPLED_AGG
INTEGER :: ND_STD, NLX_STD, NLY_STD, NLZ_STD, NUX_STD, NUY_STD, NUZ_STD 
INTEGER :: ND_AGG, NLX_AGG, NLY_AGG, NLZ_AGG, NUX_AGG, NUY_AGG, NUZ_AGG 
REAL (EB) :: ASUB(3)

!!! algebraic multigrid information
INTEGER , POINTER, DIMENSION (:)   :: CPOINTS, FPOINTS
REAL(EB), POINTER, DIMENSION (:,:) :: INTERPOL
REAL(EB), POINTER, DIMENSION (:)   :: WEIGHTS2
INTEGER,  POINTER, DIMENSION (:)   :: WEIGHTS2_ROW, WEIGHTS2_COL
LOGICAL , POINTER, DIMENSION (:,:) :: COUPLED
INTEGER:: N_CPOINTS, N_FPOINTS

!!! matrices and iteration vectors in case of bandwise storage technique
REAL (EB), POINTER, DIMENSION (:, :)    :: BA
REAL (EB), POINTER, DIMENSION (:, :, :) :: BX , BF , BD , BY , BG , BR , BZ
REAL (EB), POINTER, DIMENSION (:, :, :) :: BX2, BF2, BD2, BY2, BG2, BR2, BZ2

!!! matrices and iteration vectors in case of compact storage technique
REAL (EB), POINTER, DIMENSION (:) :: CA
INTEGER,   POINTER, DIMENSION (:) :: ROW, COL
REAL (EB), POINTER, DIMENSION (:) :: CX , CF , CD , CY , CG , CR , CZ
REAL (EB), POINTER, DIMENSION (:) :: CX2, CF2, CD2, CY2, CG2, CR2, CZ2

!!! auxiliary workspace for FFT- and GSTRI/GSADI preconditioners
REAL (EB), POINTER, DIMENSION (:, :, :) :: FFT
REAL (EB), POINTER, DIMENSION (:)       :: MDX, LDX, UDX, MDY, MDZ, DWORK, PERIOD
 
END TYPE SCARC_TYPE_LEVEL
 
!!!----------------------------------------------------------------------------------------------------
!!! SCARC type for information on different grid levels of other mesh
!!!----------------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE_LEVEL
 
!!! local neighborship structure, dimensions and number of local grid and wall cells 
INTEGER, POINTER, DIMENSION (:, :) :: IJKW
INTEGER :: IBAR, JBAR, KBAR
INTEGER :: N_CELLS, N_EXTERNAL_WALL_CELLS

!!! communication information
INTEGER :: I_MIN_R=-10,I_MAX_R=-10,J_MIN_R=-10,J_MAX_R=-10,K_MIN_R=-10,K_MAX_R=-10,NIC_R=0, &
           I_MIN_S=-10,I_MAX_S=-10,J_MIN_S=-10,J_MAX_S=-10,K_MIN_S=-10,K_MAX_S=-10,NIC_S=0

END TYPE OSCARC_TYPE_LEVEL
 

!!!----------------------------------------------------------------------------------------------------
!!! General ScaRC type on own mesh with communication vectors, MG-cycling information and
!!!! neighboring ScaRC-structures
!!!----------------------------------------------------------------------------------------------------
TYPE SCARC_TYPE

REAL (EB), POINTER, DIMENSION (:) :: SEND_BUF, RECV_BUF
INTEGER :: CYCLE_COUNT(2,NSCARC_LEVEL_MAX) = 0

TYPE (OSCARC_TYPE),      POINTER, DIMENSION (:) :: OSCARC
TYPE (SCARC_TYPE_LEVEL), POINTER, DIMENSION (:) :: LEVEL

END TYPE SCARC_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! General ScaRC type on other mesh with communication vectors and original communication 
!!! information from main 
!!!----------------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE

REAL (EB), POINTER, DIMENSION (:) :: SEND_BUF, RECV_BUF
INTEGER :: NICMAX_R, NICMAX_S

TYPE (OSCARC_TYPE_LEVEL), POINTER, DIMENSION (:) :: LEVEL

END TYPE OSCARC_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! Auxiliary type for vector pointers to bandwise/compact storage technique structures
!!!   - NX, NY, NZ : correspond to IBAR, JBAR, KBAR
!!!   - NC, NA, NW : correspond to N_CELLS, N_MATRIX_ENTRIES, N_EXTERNAL_WALL_CELLS
!!!----------------------------------------------------------------------------------------------------
TYPE BANDWISE_SYSTEM_POINTER_TYPE
INTEGER,   POINTER :: NX, NY, NZ, NC, NA, NW
REAL (EB), POINTER, DIMENSION (:, :)    :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: X , F , D , Y , G , R, Z
REAL (EB), POINTER, DIMENSION (:, :, :) :: X2, F2, D2, Y2, G2, R2
END TYPE BANDWISE_SYSTEM_POINTER_TYPE

TYPE COMPACT_SYSTEM_POINTER_TYPE
INTEGER,   POINTER :: NX, NY, NZ, NC, NA, NW
INTEGER,   POINTER, DIMENSION (:) :: ROW, COL
REAL (EB), POINTER, DIMENSION (:) :: A
REAL (EB), POINTER, DIMENSION (:) :: X , F , D , Y , G , R, Z
REAL (EB), POINTER, DIMENSION (:) :: X2, F2, D2, Y2, G2, R2
END TYPE COMPACT_SYSTEM_POINTER_TYPE


!!!----------------------------------------------------------------------------------------------------
!!!  Type variables
!!!----------------------------------------------------------------------------------------------------
TYPE (SCARC_TYPE),  SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: SCARC
TYPE (OSCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: OSCARC

TYPE (SCARC_TYPE),  POINTER ::   S,  SNM,  SNOM
TYPE (OSCARC_TYPE), POINTER ::  OS, OSNM, OSNOM

TYPE (SCARC_TYPE_LEVEL),  POINTER ::  SL,  SLMAX,  SLHI,  SLLO,  SNML,  SNOML
TYPE (OSCARC_TYPE_LEVEL), POINTER :: OSL, OSLMAX, OSLHI, OSLLO, OSNML, OSNOML

TYPE (MESH_TYPE) , POINTER ::   M
TYPE (OMESH_TYPE), POINTER ::   OM


!!!----------------------------------------------------------------------------------------------------
!!!  Interface definitions
!!!----------------------------------------------------------------------------------------------------
INTERFACE SCARC_INITIALIZE_SOLUTION
   MODULE PROCEDURE SCARC_INITIALIZE_SOLUTION_BANDWISE, SCARC_INITIALIZE_SOLUTION_COMPACT
END INTERFACE SCARC_INITIALIZE_SOLUTION

INTERFACE SCARC_TRANSFER_DATA
   MODULE PROCEDURE SCARC_TRANSFER_SOLUTION_BANDWISE, SCARC_TRANSFER_SOLUTION_COMPACT
END INTERFACE SCARC_TRANSFER_DATA

INTERFACE SCARC_LOCAL_MATVEC
   MODULE PROCEDURE SCARC_LOCAL_MATVEC_BANDWISE, SCARC_LOCAL_MATVEC_COMPACT
END INTERFACE SCARC_LOCAL_MATVEC

INTERFACE SCARC_LOCAL_SCALPROD
   MODULE PROCEDURE SCARC_LOCAL_SCALPROD_BANDWISE, SCARC_LOCAL_SCALPROD_COMPACT
END INTERFACE SCARC_LOCAL_SCALPROD

INTERFACE SCARC_RESTRICTION
   MODULE PROCEDURE SCARC_RESTRICTION_BANDWISE, SCARC_RESTRICTION_COMPACT
END INTERFACE SCARC_RESTRICTION

INTERFACE SCARC_PROLONGATION
   MODULE PROCEDURE SCARC_PROLONGATION_BANDWISE, SCARC_PROLONGATION_COMPACT
END INTERFACE SCARC_PROLONGATION

INTERFACE SCARC_PRECONDITIONER
   MODULE PROCEDURE SCARC_PRECONDITIONER_BANDWISE, SCARC_PRECONDITIONER_COMPACT
END INTERFACE SCARC_PRECONDITIONER

INTERFACE SCARC_SMOOTHER
   MODULE PROCEDURE SCARC_SMOOTHER_BANDWISE, SCARC_SMOOTHER_COMPACT
END INTERFACE SCARC_SMOOTHER

INTERFACE SCARC_SHOW
   MODULE PROCEDURE SCARC_SHOW_BANDWISE, SCARC_SHOW_COMPACT
END INTERFACE SCARC_SHOW

INTERFACE SCARC_SHOW_LEVEL
   MODULE PROCEDURE SCARC_SHOW_LEVEL_BANDWISE, SCARC_SHOW_LEVEL_COMPACT
END INTERFACE SCARC_SHOW_LEVEL

CONTAINS
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC_SETUP : Initialize ScaRC structures 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP

INTEGER :: NM, IERR
REAL(EB):: TNOW_SETUP

TNOW_SETUP = SECOND()
IERR = 0

!!!----------------------------------------------------------------------------------------------------
!!! Allocate SCARC structure and initialize time measurement array 
!!!----------------------------------------------------------------------------------------------------
ALLOCATE (SCARC(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERR)
 
ALLOCATE(TUSED_SCARC(0:N_TIMERS_SCARC,NMESHES),STAT=IERR)
CALL ChkMemErr('SCARC_SETUP','TUSED_SCARC',IERR)

TUSED_SCARC = 0._EB
TUSED_SCARC(NSCARC_TIME_COMPLETE,:) = SECOND()


!!!----------------------------------------------------------------------------------------------------
!!! Call different setup routines
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_TYPES                       ! define types of used structures
CALL SCARC_SETUP_DEBUGGING                   ! open debug file if requested
CALL SCARC_SETUP_MESHES                      ! define set of meshes which must be processed on MYID
CALL SCARC_SETUP_LEVELS                      ! define number of necessary grid levels 

MESH_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   CALL SCARC_SETUP_GEOMETRY(NM)             ! set geometry information
   CALL SCARC_SETUP_EXCHANGE(NM)             ! set information for data exchange
   CALL SCARC_SETUP_WALLCELLS(NM)            ! set wall cell information
   CALL SCARC_SETUP_SYSTEM(NM)               ! assemble system matrix including boundary conditions
   CALL SCARC_SETUP_SOLVER(NM)               ! allocate workspace for different solvers
ENDDO MESH_LOOP

CALL SCARC_SETUP_GLOBALS                      ! define some global variables
CALL SCARC_SETUP_NEIGHBORS                    ! compute information about abutting neighbors

TUSED_SCARC(NSCARC_TIME_SETUP   ,:)=TUSED_SCARC(NSCARC_TIME_SETUP   ,:)+SECOND()-TNOW_SETUP
TUSED_SCARC(NSCARC_TIME_COMPLETE,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE,:)+SECOND()-TNOW_SETUP
END SUBROUTINE SCARC_SETUP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Define a set of logical variables (just to simplify the code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      SELECT CASE(TRIM(SCARC_PRECON))
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
      SELECT CASE(TRIM(SCARC_SMOOTH))                          ! use same parameters as for preconditioner
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
   SELECT CASE(TRIM(SCARC_MULTIGRID))

      CASE ('GEOMETRIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC

      CASE ('ALGEBRAIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC

         !!! set type of coarsening strategy (STANDARD/AGGRESSIVE)
         SELECT CASE(TRIM(SCARC_MULTIGRID_COARSENING))
            CASE ('STANDARD')
               TYPE_COARSENING = NSCARC_COARSENING_STANDARD
            CASE ('AGGRESSIVE')
               TYPE_COARSENING = NSCARC_COARSENING_AGGRESSIVE
            CASE DEFAULT
               WRITE(CMESSAGE,1002) 'coarsening',TRIM(SCARC_MULTIGRID_COARSENING),&
                                    'AMG-Krylov-preconditioner','STANDARD','AGGRESSIVE'
               CALL SCARC_SHUTDOWN(CMESSAGE)
         END SELECT

      CASE DEFAULT
         WRITE(CMESSAGE,1002) 'multigrid',TRIM(SCARC_MULTIGRID),&
                              'Krylov-preconditioner','GEOMETRIC','ALGEBRAIC'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

   !!! set type of coarse grid solver (CG/GE)
   SELECT CASE(TRIM(SCARC_COARSE))
      CASE ('CG')
         TYPE_COARSE = NSCARC_COARSE_CG
      CASE ('GE')
         TYPE_COARSE = NSCARC_COARSE_GE
      CASE DEFAULT
         WRITE(CMESSAGE,1002) 'coarse grid solver',TRIM(SCARC_COARSE),'multigrid','CG','GE'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

   !!! set type of cycling pattern (F/V/W)
   SELECT CASE(TRIM(SCARC_MULTIGRID_CYCLE))
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
!!! set storage type (BANDWISE/COMPACT)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE(TRIM(SCARC_STORAGE))
   CASE ('BANDWISE')
      TYPE_STORAGE = NSCARC_STORAGE_BANDWISE
   CASE ('COMPACT')
      TYPE_STORAGE = NSCARC_STORAGE_COMPACT
   CASE DEFAULT
      WRITE(CMESSAGE,1002) 'storage',TRIM(SCARC_STORAGE),'ScaRC','BANDWISE','COMPACT'
      CALL SCARC_SHUTDOWN(CMESSAGE)
END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! set type of stencil (NPOINT/AGGRESSIVE)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE(TRIM(SCARC_STENCIL))
   CASE ('NPOINT')
      TYPE_STENCIL = NSCARC_STENCIL_NPOINT
   CASE ('AGGRESSIVE')
      TYPE_STENCIL = NSCARC_STENCIL_AGGRESSIVE
   CASE DEFAULT
      WRITE(CMESSAGE,1002) 'stencil',TRIM(SCARC_STENCIL),'ScaRC','NPOINT','AGGRESSIVE'
      CALL SCARC_SHUTDOWN(CMESSAGE)
END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! set type of accuracy (ABSOLUTE/RELATIVE)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE(TRIM(SCARC_ACCURACY))
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
SELECT CASE(TRIM(SCARC_DEBUG))
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
SELECT CASE(TRIM(SCARC_INITIAL))
   CASE ('NONE')
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
             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/A',A,/,&
             'ScaRC: Aborting program ...')
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
SUBROUTINE SCARC_SETUP_MESHES
IF (USE_MPI) THEN
   NMESHES_MIN = MYID+1
   NMESHES_MAX = MYID+1
ELSE
   NMESHES_MIN = 1
   NMESHES_MAX = NMESHES
ENDIF
END SUBROUTINE SCARC_SETUP_MESHES

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine number of grid levels  (1 for CG/BICG-method, NLEVEL for MG-method)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_LEVELS
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM

SELECT CASE (TYPE_MULTIGRID)

   !!!------------------------------------------------------------------------------------------------------
   !!! predefined hierarchy of levels in case of geometric multigrid-method
   !!!------------------------------------------------------------------------------------------------------
   CASE (NSCARC_MULTIGRID_GEOMETRIC)

      NLEVEL = NSCARC_LEVEL_MAX
      KLEVEL = NSCARC_LEVEL_MAX

      DO NM=1,NMESHES

         KLEVEL(1)=SCARC_GET_MAXLEVEL(MESHES(NM)%IBAR,1)
         IF (.NOT.TWO_D) KLEVEL(2)=SCARC_GET_MAXLEVEL(MESHES(NM)%JBAR,2)
         KLEVEL(3)=SCARC_GET_MAXLEVEL(MESHES(NM)%KBAR,3)

         KLEVEL_MIN = MINVAL(KLEVEL)
         IF (KLEVEL_MIN<NLEVEL) NLEVEL=KLEVEL_MIN

      ENDDO 
      NLEVEL_MAX  = NLEVEL
      IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
         NLEVEL_MIN  = NLEVEL-SCARC_MULTIGRID_LEVEL+1
      ELSE
         NLEVEL_MIN  = 1
      ENDIF

   !!!------------------------------------------------------------------------------------------------------
   !!! first, only finest level is set, further levels are defined during coarsening process
   !!!------------------------------------------------------------------------------------------------------
   CASE (NSCARC_MULTIGRID_ALGEBRAIC)

      NLEVEL     = NSCARC_LEVEL_MAX
      NLEVEL_MAX = NSCARC_LEVEL_MAX
      NLEVEL_MIN = 1

   !!!------------------------------------------------------------------------------------------------------
   !!! no multigrid-hierachy needed in case of a Krylov-method: use only one level
   !!!------------------------------------------------------------------------------------------------------
   CASE DEFAULT

      NLEVEL     = 1
      NLEVEL_MAX = 1
      NLEVEL_MIN = 1

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
   SELECT CASE(IOR0)
      CASE(1)
         WRITE(*,1000) IBAR, NC
      CASE(2)
         WRITE(*,1000) JBAR, NC
      CASE(3)
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
SUBROUTINE SCARC_SETUP_GEOMETRY (NM)

INTEGER, INTENT(IN) :: NM
INTEGER :: IERR, NL
INTEGER :: IBAR0, JBAR0, KBAR0
INTEGER :: I, J, K

IERR=0
M => MESHES(NM)
S => SCARC (NM)

!!!
!!!  one fine grid on every MESH for CG-method
!!!  complete hierarchy of grids on every MESH for MG-method
!!!
ALLOCATE (S%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GEOMETRY', 'LEVEL', IERR)

IBAR0=M%IBAR
JBAR0=M%JBAR
KBAR0=M%KBAR

LEVEL_LOOP3D: DO NL=NLEVEL_MAX,NLEVEL_MIN,-1

   SL => S%LEVEL(NL)
   
   !!!-------------------------------------------------------------------------------------------------
   !!! numbers of cells in x-, y- and z-direction for level 'NL'
   !!!-------------------------------------------------------------------------------------------------
   SL%IBAR=IBAR0
   SL%JBAR=JBAR0
   SL%KBAR=KBAR0
   
   !!!-------------------------------------------------------------------------------------------------
   !!! step widths in x-, y- and z-direction for level 'NL'
   !!!-------------------------------------------------------------------------------------------------
   SL%DX=(M%XF-M%XS)/REAL(SL%IBAR,EB)
   SL%DY=(M%YF-M%YS)/REAL(SL%JBAR,EB)
   SL%DZ=(M%ZF-M%ZS)/REAL(SL%KBAR,EB)
   
   SL%DXI=1.0_EB/SL%DX
   SL%DYI=1.0_EB/SL%DY
   SL%DZI=1.0_EB/SL%DZ
   
   !!!-------------------------------------------------------------------------------------------------
   !!! x-, y-, and z-coordinates for level 'NL'
   !!! -> compute minimum and maximum coordinates in each direction 
   !!! -> compute minimum and maximum step sizes  in each direction
   !!!-------------------------------------------------------------------------------------------------
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
   
   ALLOCATE (SL%YY(0:SL%JBAR), STAT=IERR)
   CALL CHKMEMERR ('SCARC_INIT_MESHES3D', 'SL%YY', IERR)
   SL%YY_MIN= 1.0E+5_EB
   SL%YY_MAX=-1.0E+5_EB
   SL%DY_MIN= 1.0E+5_EB
   SL%DY_MAX=-1.0E+5_EB
   DO J=0,SL%JBAR
      SL%YY(J)=M%YS+J*SL%DY
      SL%YY_MIN=MIN(SL%YY_MIN,SL%YY(J))
      SL%YY_MAX=MAX(SL%YY_MAX,SL%YY(J))
      IF (J>0) THEN
        SL%DY_MIN=MIN(SL%DY_MIN,ABS(SL%YY(J)-SL%YY(J-1)))
        SL%DY_MAX=MAX(SL%DY_MAX,ABS(SL%YY(J)-SL%YY(J-1)))
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
        SL%DZ_MIN=MIN(SL%DZ_MIN,ABS(SL%ZZ(K)-SL%ZZ(K-1)))
        SL%DZ_MAX=MAX(SL%DZ_MAX,ABS(SL%ZZ(K)-SL%ZZ(K-1)))
      ENDIF
   ENDDO
   
   !!!-------------------------------------------------------------------------------------------------
   !!! Get global number of grid cells
   !!!-------------------------------------------------------------------------------------------------
   SL%N_CELLS = SL%IBAR * SL%JBAR * SL%KBAR 
    
   IBAR0=IBAR0/2
   IF (.NOT.TWO_D) JBAR0=JBAR0/2
   KBAR0=KBAR0/2

ENDDO LEVEL_LOOP3D
 
END SUBROUTINE SCARC_SETUP_GEOMETRY
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup communication structure for data exchange 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_EXCHANGE (NM)
 
INTEGER, INTENT(IN) :: NM
INTEGER :: NOM, NL, IREFINE
INTEGER :: IERR

IERR = 0
M => MESHES(NM)
S => SCARC(NM)
 
!!! Initialize communication counter for ScaRC, use same TAG for all communications
NCOM_SCARC = 0
TAG_SCARC  = 99

!!!----------------------------------------------------------------------------------------------------
!!! Initialize level structures on neighboring meshes
!!! Note that OSLMAX%IJKW corresponds to MESHES(NM)%OMESH(NOM)%IJKW
!!!----------------------------------------------------------------------------------------------------
ALLOCATE (S%OSCARC(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'OSCARC', IERR)
   
OLEVEL_LOOP: DO NOM = 1, NMESHES

   IF (NOM == NM) CYCLE OLEVEL_LOOP

   OM => M%OMESH(NOM)
   OS => S%OSCARC(NOM)

   OS%NICMAX_S = OM%NIC_S
   OS%NICMAX_R = OM%NIC_R

   IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0)  CYCLE OLEVEL_LOOP

   !!! Allocate OSCARC structure for mesh NM
   ALLOCATE (OS%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'OS%LEVEL', IERR)

   OSLMAX => OS%LEVEL(NLEVEL_MAX)

   OSLMAX%IBAR = MESHES(NOM)%IBAR 
   OSLMAX%JBAR = MESHES(NOM)%JBAR 
   OSLMAX%KBAR = MESHES(NOM)%KBAR 

   OSLMAX%N_EXTERNAL_WALL_CELLS = 2*OSLMAX%IBAR * OSLMAX%JBAR + &
                                  2*OSLMAX%IBAR * OSLMAX%KBAR + &
                                  2*OSLMAX%JBAR * OSLMAX%KBAR  

   !!! In case of multigrid with a grid hierarchy allocate corresponding level-structures
   IF (NLEVEL_MAX>NLEVEL_MIN) THEN                   

      IREFINE=1
      DO NL=NLEVEL_MAX-1,NLEVEL_MIN,-1

         IREFINE=IREFINE*2

         OSLLO => OS%LEVEL(NL)                            ! pointer to lower level
         OSLHI => OS%LEVEL(NL+1)                          ! pointer to higher level

         ! get number of internal cells and external wall cells on neighbor NOM for corresponding level
         OSLLO%IBAR=OSLHI%IBAR/IREFINE
         IF (TWO_D) THEN
            OSLLO%JBAR=1
         ELSE
            OSLLO%JBAR=OSLHI%JBAR/IREFINE
         ENDIF
         OSLLO%KBAR=OSLHI%KBAR/IREFINE

         OSLLO%N_EXTERNAL_WALL_CELLS= 2*OSLLO%IBAR * OSLLO%JBAR + &
                                      2*OSLLO%IBAR * OSLLO%KBAR + &
                                      2*OSLLO%JBAR * OSLLO%KBAR  

         ! allocate array ijkw for neighbor on corresponding level
         ALLOCATE (OSLLO%IJKW(15,OSLLO%N_EXTERNAL_WALL_CELLS), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'OSLLO%IJKW', IERR)
         OSLLO%IJKW = 0

      ENDDO
   ENDIF

ENDDO OLEVEL_LOOP


!!!----------------------------------------------------------------------------------------------------
!!! Initialize communication structures on every level
!!!----------------------------------------------------------------------------------------------------
EXCHANGE_LEVEL_LOOP: DO NL=NLEVEL_MAX,NLEVEL_MIN,-1

   NBR_LOOP: DO NOM = 1, NMESHES

      IF (NOM == NM) CYCLE NBR_LOOP

!!! NOT FINISHED YET
      OS => SCARC(NM)%OSCARC(NOM)
      IF (OS%NICMAX_S == 0 .AND. OS%NICMAX_R == 0) CYCLE NBR_LOOP

   ENDDO NBR_LOOP
 
ENDDO EXCHANGE_LEVEL_LOOP
 
END SUBROUTINE SCARC_SETUP_EXCHANGE

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup neighborship structures and boundary conditions on finest level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_WALLCELLS(NM)
 
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES

INTEGER, INTENT(IN) :: NM
INTEGER :: NL, IERR, NOM
INTEGER :: IOFFSET, I_LO, J_LO, K_LO, IW_LO, IREFINE
INTEGER :: IW

IERR=0
M => MESHES(NM)
S => SCARC(NM)
   
!!!----------------------------------------------------------------------------------------------------
!!! For all solver: 
!!! Determine array IJKW and PRESSURE_BC_INDEX on finest level
!!!----------------------------------------------------------------------------------------------------
SLMAX => S%LEVEL(NLEVEL_MAX)

SLMAX%IJKW              => M%IJKW
SLMAX%PRESSURE_BC_INDEX => M%PRESSURE_BC_INDEX

SLMAX%N_EXTERNAL_WALL_CELLS =  M%N_EXTERNAL_WALL_CELLS

DO IW = 1, M%N_EXTERNAL_WALL_CELLS
   IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
      M%PRESSURE_BC_INDEX(IW)=DIRICHLET
   ELSE IF (M%IJKW(9,IW)/=0) THEN
      M%PRESSURE_BC_INDEX(IW)=INTERNAL
   ELSE IF (M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
      M%PRESSURE_BC_INDEX(IW)=DIRICHLET
   ELSE
      M%PRESSURE_BC_INDEX(IW)=NEUMANN
   ENDIF
ENDDO

OTHER_MESH_LOOP: DO NOM = 1, NMESHES
   
   IF (NOM == NM) CYCLE OTHER_MESH_LOOP

   OSLMAX => S%OSCARC(NOM)%LEVEL(NLEVEL_MAX)

   OSLMAX%NIC_S   = M%OMESH(NOM)%NIC_S
   OSLMAX%I_MIN_S = M%OMESH(NOM)%I_MIN_S
   OSLMAX%I_MAX_S = M%OMESH(NOM)%I_MAX_S
   OSLMAX%J_MIN_S = M%OMESH(NOM)%J_MIN_S
   OSLMAX%J_MAX_S = M%OMESH(NOM)%J_MAX_S
   OSLMAX%K_MIN_S = M%OMESH(NOM)%K_MIN_S
   OSLMAX%K_MAX_S = M%OMESH(NOM)%K_MAX_S
    
   OSLMAX%NIC_R   = M%OMESH(NOM)%NIC_R
   OSLMAX%I_MIN_R = M%OMESH(NOM)%I_MIN_R
   OSLMAX%I_MAX_R = M%OMESH(NOM)%I_MAX_R
   OSLMAX%J_MIN_R = M%OMESH(NOM)%J_MIN_R
   OSLMAX%J_MAX_R = M%OMESH(NOM)%J_MAX_R
   OSLMAX%K_MIN_R = M%OMESH(NOM)%K_MIN_R
   OSLMAX%K_MAX_R = M%OMESH(NOM)%K_MAX_R

   NCOM_SCARC = NCOM_SCARC+1

ENDDO OTHER_MESH_LOOP


!!!----------------------------------------------------------------------------------------------------
!!! Only in case of MG-method:
!!! Determine arrays IJKW for coarser levels
!!!----------------------------------------------------------------------------------------------------
MG_IF: IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN

   IREFINE=1
   INIT_NBR_LEVEL: DO NL=NLEVEL_MAX-1,NLEVEL_MIN,-1
   
      SLLO => S%LEVEL(NL)
      SLHI => S%LEVEL(NL+1)

      IREFINE=IREFINE*2

      !!!
      !!! set number of external wall cells on lower level and allocate corresponding arrays
      !!!
      SLLO%N_EXTERNAL_WALL_CELLS = 2 * SLLO%IBAR * SLLO%JBAR + &
                                   2 * SLLO%IBAR * SLLO%KBAR + &
                                   2 * SLLO%JBAR * SLLO%KBAR

      ALLOCATE (SLLO%IJKW(15,SLLO%N_EXTERNAL_WALL_CELLS), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_WALLCELLS', 'SLLO%IJKW', IERR)
      SLLO%IJKW=0
   
      ALLOCATE (SLLO%PRESSURE_BC_INDEX(SLLO%N_EXTERNAL_WALL_CELLS), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_WALLCELLS', 'SLLO%PRESSURE_BC_INDEX', IERR)
      SLLO%PRESSURE_BC_INDEX=0
   
      !!! 
      !!! set wall cells for coarser grid and define corresponding IJKW
      !!! 
      IW_LO=1

      !!! wall cells along IOR=1
      IOFFSET = 0
      DO K_LO=1,SLLO%KBAR
         DO J_LO=1,SLLO%JBAR
            CALL SCARC_SETUP_IJKW(SLLO%IJKW, SLHI%IJKW, SLLO%PRESSURE_BC_INDEX, SLHI%PRESSURE_BC_INDEX, &
                                  IW_LO,  1, IOFFSET, IREFINE, 0, J_LO, K_LO, SLHI%JBAR)
         ENDDO
      ENDDO

      !!! wall cells along IOR=-1
      IOFFSET = SLHI%JBAR*SLHI%KBAR
      DO K_LO=1,SLLO%KBAR
         DO J_LO=1,SLLO%JBAR
            CALL SCARC_SETUP_IJKW(SLLO%IJKW, SLHI%IJKW, SLLO%PRESSURE_BC_INDEX, SLHI%PRESSURE_BC_INDEX, &
                                  IW_LO, -1, IOFFSET, IREFINE, SLLO%IBAR+1, J_LO, K_LO, SLHI%JBAR)
         ENDDO
      ENDDO
   
      !!! wall cells along IOR=2
      IOFFSET = 2*SLHI%JBAR*SLHI%KBAR
      DO K_LO=1,SLLO%KBAR
         DO I_LO=1,SLLO%IBAR
            CALL SCARC_SETUP_IJKW(SLLO%IJKW, SLHI%IJKW, SLLO%PRESSURE_BC_INDEX, SLHI%PRESSURE_BC_INDEX, &
                                  IW_LO,  2, IOFFSET, IREFINE, I_LO, 0, K_LO, SLHI%IBAR)
         ENDDO
      ENDDO

      !!! wall cells along IOR=-2
      IOFFSET = 2*SLHI%JBAR*SLHI%KBAR + SLHI%IBAR*SLHI%KBAR
      DO K_LO=1,SLLO%KBAR
         DO I_LO=1,SLLO%IBAR
            CALL SCARC_SETUP_IJKW(SLLO%IJKW, SLHI%IJKW, SLLO%PRESSURE_BC_INDEX, SLHI%PRESSURE_BC_INDEX, &
                                  IW_LO, -2, IOFFSET, IREFINE, I_LO, SLLO%JBAR+1, K_LO, SLHI%IBAR)
         ENDDO
      ENDDO

      !!! wall cells along IOR=3
      IOFFSET = 2*SLHI%JBAR*SLHI%KBAR + 2*SLHI%IBAR*SLHI%KBAR
      DO J_LO=1,SLLO%JBAR
         DO I_LO=1,SLLO%IBAR
            CALL SCARC_SETUP_IJKW(SLLO%IJKW, SLHI%IJKW, SLLO%PRESSURE_BC_INDEX, SLHI%PRESSURE_BC_INDEX, &
                                  IW_LO,  3, IOFFSET, IREFINE, I_LO, J_LO, 0, SLHI%IBAR)
         ENDDO
      ENDDO

      !!! wall cells along IOR=-3
      IOFFSET = 2*SLHI%JBAR*SLHI%KBAR + 2*SLHI%IBAR*SLHI%KBAR + SLHI%IBAR*SLHI%JBAR
      DO J_LO=1,SLLO%JBAR
         DO I_LO=1,SLLO%IBAR
            CALL SCARC_SETUP_IJKW(SLLO%IJKW, SLHI%IJKW, SLLO%PRESSURE_BC_INDEX, SLHI%PRESSURE_BC_INDEX, &
                                  IW_LO, -3, IOFFSET, IREFINE, I_LO, J_LO, SLLO%KBAR+1, SLHI%IBAR)
         ENDDO
      ENDDO

      !!!
      !!! compute analogues to NIC, I_MIN, I_MAX, K_MIN, K_MAX on level 'NL'
      !!!
      CALL SCARC_SETUP_COARSE_DIMENSIONS(IREFINE, NM, NL)

   ENDDO INIT_NBR_LEVEL
ENDIF MG_IF

END SUBROUTINE SCARC_SETUP_WALLCELLS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set wall cell information on coarse level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_IJKW(IJKW_LO, IJKW_HI, PBCI_LO, PBCI_HI, &
                            IW_LO, IOR0, IOFFSET, IREFINE, IP, JP, KP, ILEN_HI)
INTEGER, POINTER, DIMENSION(:,:), INTENT(OUT) :: IJKW_LO
INTEGER, POINTER, DIMENSION(:,:), INTENT(IN)  :: IJKW_HI
INTEGER, POINTER, DIMENSION(:),   INTENT(OUT) :: PBCI_LO
INTEGER, POINTER, DIMENSION(:),   INTENT(IN)  :: PBCI_HI
INTEGER, INTENT(INOUT) :: IW_LO
INTEGER, INTENT(IN) :: IOR0, IOFFSET, IREFINE
INTEGER, INTENT(IN) :: IP, JP, KP, ILEN_HI
INTEGER :: IW_HI(4) , IBC_HI(4), NOM_HI(4)
INTEGER :: I, I1, I2, J1, J2, K1, K2
INTEGER :: IDIFF, JDIFF, KDIFF

!!!
!!! Set IJKW_LO(1:8, IW_LO)
!!!
SELECT CASE(IOR0)
   CASE(1)
      CALL SCARC_SETUP_IJKW1(IJKW_LO, IW_LO, IOR0, IP, JP, KP, IP+1, JP  , KP  )
   CASE(-1)
      CALL SCARC_SETUP_IJKW1(IJKW_LO, IW_LO, IOR0, IP, JP, KP, IP-1, JP  , KP  )
   CASE(2)
      CALL SCARC_SETUP_IJKW1(IJKW_LO, IW_LO, IOR0, IP, JP, KP, IP  , JP+1, KP  )
   CASE(-2)
      CALL SCARC_SETUP_IJKW1(IJKW_LO, IW_LO, IOR0, IP, JP, KP, IP  , JP-1, KP  )
   CASE(3)
      CALL SCARC_SETUP_IJKW1(IJKW_LO, IW_LO, IOR0, IP, JP, KP, IP  , JP  , KP+1)
   CASE(-3)
      CALL SCARC_SETUP_IJKW1(IJKW_LO, IW_LO, IOR0, IP, JP, KP, IP  , JP  , KP-1)
END SELECT


!!!----------------------------------------------------------------------------------------------------
!!! 2D-version
!!!----------------------------------------------------------------------------------------------------
IF (TWO_D) THEN

   !!! set IJKW(1:8,IW_LO) for coarser grid IW_LO
   SELECT CASE(ABS(IOR0))
      CASE( 1)
         IW_HI(1) = IOFFSET + 2*KP-1
      CASE( 2)
         IW_HI(1) = IOFFSET + (2*KP-1)*ILEN_HI + 2*IP - 1
      CASE( 3)
         IW_HI(1) = IOFFSET + 2*IP-1
   END SELECT
   IW_HI(2) = IW_HI(1)+1

   !!! set neighbors IJKW(9,IW_LO) for coarser grid IW_LO
   NOM_HI(1) = IJKW_HI(9,IW_HI(1))
   NOM_HI(2) = IJKW_HI(9,IW_HI(2))
   IF (NOM_HI(1)/=NOM_HI(2)) THEN
      WRITE(*,*) 'SCARC_SETUP_IJKW: Inconsistent neighbors on IOR=', IOR0,' not allowed!'
      STOP
   ENDIF

   IJKW_LO(9,IW_LO)=NOM_HI(1) 

   !!! set corresponding pressure_bc_index on lower level
   IBC_HI(1) = PBCI_HI(IW_HI(1))
   IBC_HI(2) = PBCI_HI(IW_HI(2))
   IF (IBC_HI(1)==INTERNAL .OR. IBC_HI(2)==INTERNAL) THEN
      PBCI_LO(IW_LO)=INTERNAL
   ELSE IF (IBC_HI(1)==DIRICHLET .OR. IBC_HI(2)==DIRICHLET) THEN
      PBCI_LO(IW_LO)=DIRICHLET
   ELSE
      PBCI_LO(IW_LO)=NEUMANN
   ENDIF

   !!! in case of an internal boundary set IJKW(10:15,IW_LO)
   IF (NOM_HI(1) > 0) THEN   

      J1 = 1
      J2 = 1
      SELECT CASE(ABS(IOR0))
         CASE(1)
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
         CASE(3)
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

      SELECT CASE(IOR0)
         CASE(1)
            I1 = MESHES(NOM_HI(1))%IBAR/IREFINE
            I2 = I1
         CASE(-1)
            I1 = 1
            I2 = 1
         CASE(3)
            K1 = MESHES(NOM_HI(1))%KBAR/IREFINE
            K2 = K1
         CASE(-3)
            K1 = 1
            K2 = 1
      END SELECT

      !!!
      !!! Set IJKW_LO(10:15, IW_LO)
      !!!
      CALL SCARC_SETUP_IJKW2(IJKW_LO, IW_LO, I1, J1, K1, I2, J2, K2)

   ENDIF
      

!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE

   !!! set IJKW(1:8,IW_LO) for coarser grid IW_LO
   SELECT CASE(ABS(IOR0))
      CASE(1)
         IW_HI(1) = IOFFSET + (2*KP-2)*ILEN_HI + 2*JP - 1
         IW_HI(3) = IOFFSET + (2*KP-1)*ILEN_HI + 2*JP - 1
      CASE(2)
         IW_HI(1) = IOFFSET + (2*KP-2)*ILEN_HI + 2*IP - 1
         IW_HI(3) = IOFFSET + (2*KP-1)*ILEN_HI + 2*IP - 1
      CASE(3)
         IW_HI(1) = IOFFSET + (2*JP-2)*ILEN_HI + 2*IP - 1
         IW_HI(3) = IOFFSET + (2*JP-1)*ILEN_HI + 2*IP - 1
   END SELECT
   IW_HI(2) = IW_HI(1)+1
   IW_HI(4) = IW_HI(3)+1

   !!! set neighbors IJKW(9,IW_LO) for coarser grid IW_LO
   DO I=1,4
      NOM_HI(I) = IJKW_HI(9,IW_HI(I))
   ENDDO
      
   IF (NOM_HI(1)/=NOM_HI(2) .OR. NOM_HI(1)/=NOM_HI(3) .OR. NOM_HI(1)/=NOM_HI(4)) THEN
      WRITE(*,*) 'SCARC_SETUP_IJKW: Inconsistent neighbors on IOR=', IOR0,' not allowed!'
      STOP
   ENDIF
   IJKW_LO(9,IW_LO)=NOM_HI(1) 

   !!! set corresponding pressure_bc_index on lower level
   DO I=1,4
      IBC_HI(I) = PBCI_HI(IW_HI(I))
   ENDDO
   IF (IBC_HI(1)==INTERNAL.OR.IBC_HI(2)==INTERNAL.OR.IBC_HI(3)==INTERNAL.OR.IBC_HI(4)==INTERNAL) THEN
      PBCI_LO(IW_LO)=INTERNAL
   ELSE IF (IBC_HI(1)==DIRICHLET.OR.IBC_HI(2)==DIRICHLET.OR.IBC_HI(3)==DIRICHLET.OR.IBC_HI(4)==DIRICHLET) THEN
      PBCI_LO(IW_LO)=DIRICHLET
   ELSE
      PBCI_LO(IW_LO)=NEUMANN
   ENDIF

   !!! in case of an internal boundary set IJKW(10:15,IW_LO)
   IF (NOM_HI(1) > 0) THEN   

      SELECT CASE(ABS(IOR0))
         CASE(1)
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
         CASE(2)
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
         CASE(3)
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

      SELECT CASE(IOR0)
         CASE(1)
            I1 = MESHES(NOM_HI(1))%IBAR/IREFINE
            I2 = I1
         CASE(-1)
            I1 = 1
            I2 = I1
         CASE(2)
            J1 = MESHES(NOM_HI(1))%JBAR/IREFINE
            J2 = J1
         CASE(-2)
            J1 = 1
            J2 = J1
         CASE(3)
            K1 = MESHES(NOM_HI(1))%KBAR/IREFINE
            K2 = K1
         CASE(-3)
            K1 = 1
            K2 = K1
      END SELECT

      !!!
      !!! Set IJKW_LO(10:15, IW_LO)
      !!!
      CALL SCARC_SETUP_IJKW2(IJKW_LO, IW_LO, I1, J1, K1, I2, J2, K2)

   ENDIF

ENDIF

IW_LO = IW_LO + 1

END SUBROUTINE SCARC_SETUP_IJKW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! First part: Set correct IJKW-values in SCARC_INIT_NEIGHBOURS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_IJKW1(IJKW0, IW0, IOR0, IMIN, JMIN, KMIN, IMAX, JMAX, KMAX)
INTEGER, POINTER, DIMENSION(:,:), INTENT(OUT) :: IJKW0
INTEGER, INTENT(IN) :: IW0, IOR0
INTEGER, INTENT(IN) :: IMIN, JMIN, KMIN, IMAX, JMAX, KMAX

IJKW0(1,IW0)=IMIN           
IJKW0(2,IW0)=JMIN    
IJKW0(3,IW0)=KMIN 
IJKW0(4,IW0)=IOR0     
IJKW0(5,IW0)=0        
IJKW0(6,IW0)=IMAX           
IJKW0(7,IW0)=JMAX 
IJKW0(8,IW0)=KMAX 

END SUBROUTINE SCARC_SETUP_IJKW1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Second part: Set correct IJKW-values in SCARC_INIT_NEIGHBOURS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_IJKW2(IJKW0,IW0, IMIN, JMIN, KMIN, IMAX, JMAX, KMAX)
INTEGER, POINTER, DIMENSION(:,:), INTENT(OUT) :: IJKW0
INTEGER, INTENT(IN) :: IW0, IMIN, JMIN, KMIN, IMAX, JMAX, KMAX

IJKW0(10,IW0)=IMIN
IJKW0(11,IW0)=JMIN
IJKW0(12,IW0)=KMIN
IJKW0(13,IW0)=IMAX
IJKW0(14,IW0)=JMAX
IJKW0(15,IW0)=KMAX

END SUBROUTINE SCARC_SETUP_IJKW2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set analogues to I_MIN, IMAX, J_MIN, J_MAX, K_MIN, K_MAX on coarse level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_COARSE_DIMENSIONS(IREFINE, NM, NL)
INTEGER, INTENT(IN) :: IREFINE, NM, NL
INTEGER :: NOM, IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IERR, IW
LOGICAL :: FOUND

IERR = 0

OTHER_MESH_LOOP: DO NOM = 1, NMESHES

   IF (NOM == NM) CYCLE OTHER_MESH_LOOP
   OSNML => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
    
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
   
   OSNML%NIC_S = 0
   FOUND = .FALSE.

   SEARCH_LOOP: DO IW=1,SL%N_EXTERNAL_WALL_CELLS
   
      ! neighborship structure already known from finest level
      IF (MESHES(NM)%IJKW(9,IW)/=NOM) CYCLE SEARCH_LOOP
      OSNML%NIC_S = OSNML%NIC_S + 1
      FOUND = .TRUE.
   
      SELECT CASE(IJKW(4,IW))
         CASE( 1)
            IMIN=MAX(IMIN,IJKW(10,IW)-1)
         CASE(-1) 
            IMAX=MIN(IMAX,IJKW(13,IW))
         CASE( 2) 
            JMIN=MAX(JMIN,IJKW(11,IW)-1)
         CASE(-2) 
            JMAX=MIN(JMAX,IJKW(14,IW))
         CASE( 3) 
            KMIN=MAX(KMIN,IJKW(12,IW)-1)
         CASE(-3)
            KMAX=MIN(KMAX,IJKW(15,IW))
      END SELECT
   ENDDO SEARCH_LOOP
   
   IF (.NOT.FOUND) CYCLE OTHER_MESH_LOOP

   NCOM_SCARC = NCOM_SCARC+1

   OSNML%I_MIN_R = IMIN
   OSNML%I_MAX_R = IMAX
   OSNML%J_MIN_R = JMIN
   OSNML%J_MAX_R = JMAX
   OSNML%K_MIN_R = KMIN
   OSNML%K_MAX_R = KMAX
   
ENDDO OTHER_MESH_LOOP

END SUBROUTINE SCARC_SETUP_COARSE_DIMENSIONS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_SOLVER (NM)
 
INTEGER, INTENT(IN) :: NM
INTEGER :: IERR, NL, IBP1, JBP1, KBP1
 
IERR = 0
S => SCARC(NM)

SELECT_STORAGE: SELECT CASE(TYPE_STORAGE) 

   !!!-------------------------------------------------------------------------------------------------
   !!! Bandwise storage technique
   !!!-------------------------------------------------------------------------------------------------
   CASE(NSCARC_STORAGE_BANDWISE)
      
      LEVEL_LOOP_BANDWISE: DO NL=NLEVEL_MAX,NLEVEL_MIN,-1
      
         SL => SCARC(NM)%LEVEL(NL)
      
         IBP1 = SL%IBAR + 1
         JBP1 = SL%JBAR + 1
         KBP1 = SL%IBAR + 1

         SELECT_METHOD_BANDWISE: SELECT CASE(TYPE_METHOD)
      
            !!! working and auxiliary vectors for global CG/BICG-method
            CASE(NSCARC_METHOD_KRYLOV)
      
               ALLOCATE (SL%BX(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'X', IERR)
               SL%BX = 0.0_EB
          
               ALLOCATE (SL%BF(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'B', IERR)
               SL%BF = 0.0_EB
          
               ALLOCATE (SL%BD(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'D', IERR)
               SL%BD = 0.0_EB
      
               ALLOCATE (SL%BR(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'R', IERR)
               SL%BR = 0.0_EB
          
               ALLOCATE (SL%BG(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'G', IERR)
               SL%BG = 0.0_EB
          
               ALLOCATE (SL%BY(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'Y', IERR)
               SL%BY = 0.0_EB
      
               ALLOCATE (SL%BZ(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'Z', IERR)
               SL%BZ = 0.0_EB
      
      
               IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
            
                  ALLOCATE (SL%BX2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X2', IERR)
                  SL%BX2 = 0.0_EB
            
                  ALLOCATE (SL%BF2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'F2', IERR)
                  SL%BF2 = 0.0_EB
            
                  ALLOCATE (SL%BD2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D2', IERR)
                  SL%BD2 = 0.0_EB
      
                  IF (NL==NLEVEL_MIN) THEN
         
                     ALLOCATE (SL%BR2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'R2', IERR)
                     SL%BR2 = 0.0_EB
             
                     ALLOCATE (SL%BG2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'G2', IERR)
                     SL%BG2 = 0.0_EB
             
                     ALLOCATE (SL%BY2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'Y2', IERR)
                     SL%BY2 = 0.0_EB
         
                  ENDIF
                
               ENDIF
      
               IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN 
                  ALLOCATE (SL%FFT(1:IBP1, 1:JBP1, 1:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                  SL%FFT = 0.0_EB
               ENDIF
      
            !!! working and auxiliary vectors for global GMG/AMG-method
            CASE(NSCARC_METHOD_MULTIGRID)
      
               ALLOCATE (SL%BX(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'X', IERR)
               SL%BX = 0.0_EB
          
               ALLOCATE (SL%BF(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'B', IERR)
               SL%BF = 0.0_EB
          
               ALLOCATE (SL%BD(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'D', IERR)
               SL%BD = 0.0_EB
      
               IF (NL==NLEVEL_MIN) THEN
      
                  ALLOCATE (SL%BR(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'R', IERR)
                  SL%BR = 0.0_EB
          
                  ALLOCATE (SL%BG(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'G', IERR)
                  SL%BG = 0.0_EB
          
                  ALLOCATE (SL%BY(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Y', IERR)
                  SL%BY = 0.0_EB
      
               ENDIF
      
               IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
                  ALLOCATE (SL%FFT(1:IBP1, 1:JBP1, 1:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                  SL%FFT = 0.0_EB
               ENDIF
      
         END SELECT SELECT_METHOD_BANDWISE
      ENDDO LEVEL_LOOP_BANDWISE
      

   !!!-------------------------------------------------------------------------------------------------
   !!! Compact storage technique
   !!!-------------------------------------------------------------------------------------------------
   CASE(NSCARC_STORAGE_COMPACT)

      LEVEL_LOOP_COMPACT: DO NL=NLEVEL_MAX,NLEVEL_MIN,-1
      
         SL => SCARC(NM)%LEVEL(NL)
      
         SELECT_METHOD_COMPACT: SELECT CASE(TRIM(SCARC_METHOD))
      
            !!! working and auxiliary vectors for global CG/BICG-method
            CASE('KRYLOV')
      
               ALLOCATE (SL%CX(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'X', IERR)
               SL%CX = 0.0_EB
          
               ALLOCATE (SL%CF(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'B', IERR)
               SL%CF = 0.0_EB
          
               ALLOCATE (SL%CD(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'D', IERR)
               SL%CD = 0.0_EB
      
               ALLOCATE (SL%CR(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'R', IERR)
               SL%CR = 0.0_EB
          
               ALLOCATE (SL%CG(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'G', IERR)
               SL%CG = 0.0_EB
          
               ALLOCATE (SL%CY(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'Y', IERR)
               SL%CY = 0.0_EB
      
               ALLOCATE (SL%CZ(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'Z', IERR)
               SL%CZ = 0.0_EB
      
      
               IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
            
                  ALLOCATE (SL%CX2(SL%N_CELLS), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X2', IERR)
                  SL%CX2 = 0.0_EB
            
                  ALLOCATE (SL%CF2(SL%N_CELLS), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'F2', IERR)
                  SL%CF2 = 0.0_EB
            
                  ALLOCATE (SL%CD2(SL%N_CELLS), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D2', IERR)
                  SL%CD2 = 0.0_EB
      
                  IF (NL==NLEVEL_MIN) THEN
         
                     ALLOCATE (SL%CR2(SL%N_CELLS), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'R2', IERR)
                     SL%CR2 = 0.0_EB
             
                     ALLOCATE (SL%CG2(SL%N_CELLS), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'G2', IERR)
                     SL%CG2 = 0.0_EB
             
                     ALLOCATE (SL%CY2(SL%N_CELLS), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'Y2', IERR)
                     SL%CY2 = 0.0_EB
         
                  ENDIF
                
               ENDIF
      
               IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
                  ALLOCATE (SL%FFT(1:SL%IBAR+1, 1:SL%JBAR+1, 1:SL%KBAR+1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                  SL%FFT = 0.0_EB
               ENDIF
      
            !!! working and auxiliary vectors for global GMG/AMG-method
            CASE('MULTIGRID')
      
               ALLOCATE (SL%CX(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'X', IERR)
               SL%CX = 0.0_EB
          
               ALLOCATE (SL%CF(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'B', IERR)
               SL%CF = 0.0_EB
          
               ALLOCATE (SL%CD(SL%N_CELLS), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'D', IERR)
               SL%CD = 0.0_EB
      
               IF (NL==NLEVEL_MIN) THEN
      
                  ALLOCATE (SL%CR(SL%N_CELLS), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'R', IERR)
                  SL%CR = 0.0_EB
          
                  ALLOCATE (SL%CG(SL%N_CELLS), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'G', IERR)
                  SL%CG = 0.0_EB
          
                  ALLOCATE (SL%CY(SL%N_CELLS), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Y', IERR)
                  SL%CY = 0.0_EB
      
               ENDIF
      
               IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
                  ALLOCATE (SL%FFT(1:SL%IBAR+1, 1:SL%JBAR+1, 1:SL%KBAR+1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                  SL%FFT = 0.0_EB
               ENDIF
      
         END SELECT SELECT_METHOD_COMPACT
      ENDDO LEVEL_LOOP_COMPACT

END SELECT SELECT_STORAGE

END SUBROUTINE SCARC_SETUP_SOLVER

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate several global structures for data exchange 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NL, IERR
INTEGER :: IREFINE, IBAR1, JBAR1, KBAR1
REAL(EB) :: SCAL
IERR = 0

!!!
!!! Allocate arrays which are used for (global) communciations
!!!
ALLOCATE(NC_GLOBAL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'NC_GLOBAL', IERR)

ALLOCATE(NC_LOCAL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'NC_LOCAL', IERR)

ALLOCATE(SP_LOCAL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBAL', 'SP_LOCAL', IERR)

!!!
!!! Compute global number of cells for all levels
!!!
IREFINE=0
DO NL=NLEVEL_MAX,NLEVEL_MIN,-1

   DO NM = NMESHES_MIN, NMESHES_MAX
      NC_LOCAL(NM)=SCARC(NM)%LEVEL(NL)%N_CELLS
   ENDDO

   IF (NMESHES>1) THEN

      ! Initialize communication structures
      NREQ_SCARC=0
      CALL SCARC_RECEIVE (NSCARC_EXCHANGE_INIT, NL)    
      CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_INIT, NL, NSCARC_DUMMY, NSCARC_DUMMY)

      ! Determine global number of cells for all levels 
      IF (USE_MPI) THEN
         CALL MPI_ALLREDUCE(NC_LOCAL,NC_GLOBAL(NL),NMESHES,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         NC_GLOBAL(NL)=0
         DO NM=1,NMESHES
            SCAL=1.0_EB/REAL(2**IREFINE,EB)
            IBAR1=SCAL*MESHES(NM)%IBAR
            IF (TWO_D) THEN
               JBAR1=1
            ELSE
               JBAR1=SCAL*MESHES(NM)%JBAR
            ENDIF
            KBAR1=SCAL*MESHES(NM)%KBAR
            NC_GLOBAL(NL)=NC_GLOBAL(NL) + IBAR1*JBAR1*KBAR1
         ENDDO
         IREFINE=IREFINE+1
      ENDIF
   ELSE
      NC_GLOBAL(NL) = NC_LOCAL(1)
   ENDIF
ENDDO

END SUBROUTINE SCARC_SETUP_GLOBALS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate several global structures for data exchange 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_NEIGHBORS
INTEGER :: NM, NOM, NL, IERR
IERR = 0

!!!
!!! Exchange information about lengths of abutting faces
!!!
IF (NMESHES>1) THEN
   ALLOCATE (REQ_SCARC(NCOM_SCARC*40))
   CALL CHKMEMERR ('SCARC_SETUP_GLOBAL', 'REQ_SCARC', IERR)
   REQ_SCARC = MPI_REQUEST_NULL
ENDIF

IF (NMESHES < 0) THEN

   LEVEL_LOOP: DO NL = NLEVEL_MAX-1, NLEVEL_MIN, -1
   
      DO NM=NMESHES_MIN,NMESHES_MAX
         DO NOM=1,NMESHES
   
            OSNML  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            OSNOML => SCARC(NOM)%OSCARC(NM)%LEVEL(NL)
   
            IF (USE_MPI) THEN
               CALL MPI_SEND(OSNML%I_MIN_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OSNML%I_MAX_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OSNML%J_MIN_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OSNML%J_MAX_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OSNML%K_MIN_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OSNML%K_MAX_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OSNML%NIC_S,  1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
            ELSE
               OSNOML%I_MIN_S = OSNML%I_MIN_R
               OSNOML%I_MAX_S = OSNML%I_MAX_R
               OSNOML%J_MIN_S = OSNML%J_MIN_R
               OSNOML%J_MAX_S = OSNML%J_MAX_R
               OSNOML%K_MIN_S = OSNML%K_MIN_R
               OSNOML%K_MAX_S = OSNML%K_MAX_R
               OSNOML%NIC_R   = OSNML%NIC_S
            ENDIF
         ENDDO
      ENDDO
   
      DO NM=1,NMESHES
         DO NOM=NMESHES_MIN,NMESHES_MAX
   
            OSNOML => SCARC(NOM)%OSCARC(NM)%LEVEL(NL)
   
            IF (USE_MPI) THEN
               CALL MPI_RECV(OSNOML%I_MIN_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSNOML%I_MAX_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSNOML%J_MIN_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSNOML%J_MAX_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSNOML%K_MIN_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSNOML%K_MAX_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSNOML%NIC_R,  1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
            ENDIF
   
         ENDDO 
      ENDDO 
   
   ENDDO LEVEL_LOOP

ENDIF

END SUBROUTINE SCARC_SETUP_NEIGHBORS
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup system of equation:
!!! Define matrix stencils and initialize matrices and boundary conditions on all needed levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_SYSTEM (NM)
 
INTEGER, INTENT(IN) :: NM
INTEGER :: NL
 

!!! for all methods the standard n-point stencil is used on finest level
TYPE_STENCIL = NSCARC_STENCIL_NPOINT

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
               !!!    -  use bandwise storage technique on all levels
               !!!    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
               CASE (NSCARC_MULTIGRID_GEOMETRIC)

                  IF (SCARC_STORAGE == 'null') TYPE_STORAGE = NSCARC_STORAGE_BANDWISE

                  DO NL = NLEVEL_MAX, NLEVEL_MIN, -1
                     CALL SCARC_SETUP_STENCIL (NM, NL)
                     CALL SCARC_SETUP_MATRIX  (NM, NL)
                     CALL SCARC_SETUP_BOUNDARY(NM, NL)
                  ENDDO 

               !!! algebraic multigrid:
               !!!    -  use compact storage technique on all levels
               !!!    -  assemble standard n-point-matrix only on finest level 
               !!!    -  construct all coarser levels by requested coarsening strategy
               CASE (NSCARC_MULTIGRID_ALGEBRAIC)

                  IF (SCARC_STORAGE == 'null') TYPE_STORAGE = NSCARC_STORAGE_COMPACT

                  CALL SCARC_SETUP_STENCIL (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)

                  CALL SCARC_SETUP_AMG(NM)

            END SELECT SELECT_PRECON_MG

         !!!-------------------------------------------------------------------------------------------
         !!! in case of one-level preconditioners (JACOBI/SSOR/GSTRIX/FFT)
         !!!    -  use bandwise storage technique on finest level (no other levels needed)
         !!!    -  assemble standard n-point-matrix on finest level 
         !!!-------------------------------------------------------------------------------------------
         CASE DEFAULT

            IF (SCARC_STORAGE == 'null') TYPE_STORAGE = NSCARC_STORAGE_BANDWISE

            CALL SCARC_SETUP_STENCIL (NM, NLEVEL_MAX)
            CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MAX)
            CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)

      END SELECT SELECT_PRECON

   !!!-------------------------------------------------------------------------------------------------
   !!! Multigrid as main solver
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_MULTIGRID)

      SELECT_MG: SELECT CASE(TYPE_MULTIGRID)

         !!!-------------------------------------------------------------------------------------------
         !!! geometric multigrid:
         !!!    -  use bandwise storage technique on all levels
         !!!    -  assemble standard n-point-matrix hierarchy on all levels 
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            IF (SCARC_STORAGE == 'null') TYPE_STORAGE = NSCARC_STORAGE_BANDWISE

            DO NL = NLEVEL_MAX, NLEVEL_MIN, -1
               CALL SCARC_SETUP_STENCIL (NM, NL)
               CALL SCARC_SETUP_MATRIX  (NM, NL)
               CALL SCARC_SETUP_BOUNDARY(NM, NL)
            ENDDO 

         !!!-------------------------------------------------------------------------------------------
         !!! algebraic multigrid:
         !!!    -  use bandwise storage technique on finest level, compact storage technique on rest
         !!!    -  assemble standard n-point-matrix only on finest level
         !!!    -  construct all coarser levels by requested coarsening strategy
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)

            IF (SCARC_STORAGE == 'null') TYPE_STORAGE = NSCARC_STORAGE_COMPACT
   
            CALL SCARC_SETUP_STENCIL (NM, NLEVEL_MAX)
            CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MAX)
            CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
   
            CALL SCARC_SETUP_AMG(NM)

      END SELECT SELECT_MG

END SELECT SELECT_SOLVER

END SUBROUTINE SCARC_SETUP_SYSTEM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_MATRIX (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, IERR, IP

SL => SCARC(NM)%LEVEL(NL)

SELECT CASE(TYPE_STORAGE)

!!!----------------------------------------------------------------------------------------------------
!!! Bandwise storage technique:
!!!
!!! The matrix is stored "bandwise", namely diagonal for diagonal.
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
   CASE(NSCARC_STORAGE_BANDWISE)

      SL%N_MATRIX_ENTRIES = SL%N_CELLS * SL%N_COUPLED_STD                     

      ALLOCATE (SL%BA(1:SL%N_CELLS, 1:SL%N_COUPLED_STD), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%BA', IERR)
      SL%BA = 0.0_EB
      
      !!! 
      !!! 2D-version
      !!!
      IF (TWO_D) THEN
      
          SL%DXI2=1.0_EB/(SL%DX)**2
          SL%DYI2=1.0_EB
          SL%DZI2=1.0_EB/(SL%DZ)**2
       
          DO K = 1, SL%KBAR
             DO I = 1, SL%IBAR
           
                IC = (K-1) * SL%IBAR + I
           
                SL%BA(IC, ID) = - 2.0_EB*SL%DXI2 - 2.0_EB*SL%DZI2  ! main diagonal
       
                IF (I > 1)       SL%BA(IC,ILX) = SL%DXI2           ! lower subdiagonal in x-direction
                IF (K > 1)       SL%BA(IC,ILZ) = SL%DZI2           ! lower subdiagonal in z-direction
                IF (I < SL%IBAR) SL%BA(IC,IUX) = SL%DXI2           ! upper subdiagonal in x-direction
                IF (K < SL%KBAR) SL%BA(IC,IUZ) = SL%DZI2           ! upper subdiagonal in z-direction
           
             ENDDO
          ENDDO
      
      !!!
      !!! 3D-version
      !!!
      ELSE
      
         SL%DXI2=1.0_EB/(SL%DX)**2
         SL%DYI2=1.0_EB/(SL%DY)**2
         SL%DZI2=1.0_EB/(SL%DZ)**2
      
         DO K = 1, SL%KBAR
            DO J = 1, SL%JBAR
               DO I = 1, SL%IBAR
          
                  IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
          
                  SL%BA(IC, ID) = - 2.0_EB*SL%DXI2 - 2.0*SL%DYI2 - 2.0_EB*SL%DZI2  ! main diagonal
      
                  IF (I > 1)       SL%BA(IC,ILX) = SL%DXI2         ! lower subdiagonal in x-direction
                  IF (J > 1)       SL%BA(IC,ILY) = SL%DYI2         ! lower subdiagonal in y-direction
                  IF (K > 1)       SL%BA(IC,ILZ) = SL%DZI2         ! lower subdiagonal in z-direction
                  IF (I < SL%IBAR) SL%BA(IC,IUX) = SL%DXI2         ! upper subdiagonal in x-direction
                  IF (J < SL%JBAR) SL%BA(IC,IUY) = SL%DYI2         ! upper subdiagonal in y-direction
                  IF (K < SL%KBAR) SL%BA(IC,IUZ) = SL%DZI2         ! upper subdiagonal in z-direction
          
               ENDDO
            ENDDO
         ENDDO
             
      ENDIF
             
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
   CASE (NSCARC_STORAGE_COMPACT)

      ALLOCATE (SL%CA(SL%N_CELLS*SL%N_COUPLED_STD), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%CA', IERR)
      SL%CA = 0.0_EB
      
      ALLOCATE (SL%ROW(SL%N_CELLS+1), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%ROW', IERR)
      SL%ROW = 0
      
      ALLOCATE (SL%COL(SL%N_CELLS*SL%N_COUPLED_STD), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%COL', IERR)
      SL%COL = 0.0_EB
      
      !!! 
      !!! 2D-version
      !!!
      IF (TWO_D) THEN
      
         SL%DXI2=1.0_EB/(SL%DX)**2
         SL%DYI2=1.0_EB
         SL%DZI2=1.0_EB/(SL%DZ)**2
       
         IP = 1
         DO K = 1, SL%KBAR
            DO I = 1, SL%IBAR
          
               IC = (K-1) * SL%IBAR + I
          
               !!! main diagonal
               SL%CA(IP)   = - 2.0_EB*SL%DXI2 - 2.0_EB*SL%DZI2      
               SL%ROW(IC) = IP                                     
               SL%COL(IP) = IC                                     
               IP = IP + 1
      
               ! lower subdiagonal in z-direction
               IF (K > 1) THEN
                  SL%CA(IP)    = SL%DZI2           
                  SL%COL(IP) = IC + SL%STENCIL_STD(ILZ)
                  IP = IP + 1
               ENDIF

               !!! lower subdiagonal in x-direction
               IF (I > 1) THEN
                  SL%CA(IP)    =  SL%DXI2                            
                  SL%COL(IP) = IC + SL%STENCIL_STD(ILX)
                  IP = IP + 1
               ENDIF

               ! upper subdiagonal in x-direction
               IF (I < SL%IBAR) THEN
                  SL%CA(IP)    = SL%DXI2           
                  SL%COL(IP) = IC + SL%STENCIL_STD(IUX)
                  IP = IP + 1
               ENDIF

               ! upper subdiagonal in z-direction
               IF (K < SL%KBAR) THEN
                  SL%CA(IP)    = SL%DZI2           
                  SL%COL(IP) = IC + SL%STENCIL_STD(IUZ)
                  IP = IP + 1
               ENDIF

            ENDDO
         ENDDO
      
      !!!
      !!! 3D-version
      !!!
      ELSE
      
         SL%DXI2=1.0_EB/(SL%DX)**2
         SL%DYI2=1.0_EB/(SL%DY)**2
         SL%DZI2=1.0_EB/(SL%DZ)**2
      
         DO K = 1, SL%KBAR
            DO J = 1, SL%JBAR
               DO I = 1, SL%IBAR
          
                  IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
          
                  !!! main diagonal
                  SL%CA(IP)   = - 2.0_EB*SL%DXI2 - 2.0*SL%DYI2 - 2.0_EB*SL%DZI2  
                  SL%ROW(IC) = IP                                     
                  SL%COL(IP) = IC 
                  IP = IP + 1
         
                  ! lower subdiagonal in z-direction
                  IF (K > 1) THEN
                     SL%CA(IP)   = SL%DZI2           
                     SL%COL(IP) = IC + SL%STENCIL_STD(ILZ)
                     IP = IP + 1
                  ENDIF
  
                  !!! lower subdiagonal in y-direction
                  IF (J > 1) THEN
                     SL%CA(IP)   =  SL%DYI2                            
                     SL%COL(IP) = IC + SL%STENCIL_STD(ILY)
                     IP = IP + 1
                  ENDIF
  
                  !!! lower subdiagonal in x-direction
                  IF (I > 1) THEN
                     SL%CA(IP)   =  SL%DXI2                            
                     SL%COL(IP) = IC + SL%STENCIL_STD(ILX)
                     IP = IP + 1
                  ENDIF
  
                  ! upper subdiagonal in x-direction
                  IF (I < SL%IBAR) THEN
                     SL%CA(IP)   = SL%DXI2           
                     SL%COL(IP) = IC + SL%STENCIL_STD(IUX)
                     IP = IP + 1
                  ENDIF
  
                  ! upper subdiagonal in y-direction
                  IF (J < SL%IBAR) THEN
                     SL%CA(IP)   = SL%DYI2           
                     SL%COL(IP) = IC + SL%STENCIL_STD(IUY)
                     IP = IP + 1
                  ENDIF
  
                  ! upper subdiagonal in z-direction
                  IF (K < SL%KBAR) THEN
                     SL%CA(IP)   = SL%DZI2           
                     SL%COL(IP) = IC + SL%STENCIL_STD(IUZ)
                     IP = IP + 1
                  ENDIF

               ENDDO
            ENDDO
         ENDDO
          
      ENDIF
      SL%ROW(SL%N_CELLS+1) = IP
      SL%N_MATRIX_ENTRIES = IP - 1
             
END SELECT

SL%ASUB(1) = SL%DXI2
SL%ASUB(2) = SL%DYI2
SL%ASUB(3) = SL%DZI2

END SUBROUTINE SCARC_SETUP_MATRIX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set pointer for different structures on level NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER  :: I, J, K, IOR0, IW, IC, NOM, BC_INDEX, IP
REAL(EB) :: DBC

SL => SCARC(NM)%LEVEL(NL)

!!!----------------------------------------------------------------------------------------------------
!!! 2D-version
!!!----------------------------------------------------------------------------------------------------
IF (TWO_D) THEN

   WALL_CELL_LOOP2D: DO IW = 1, SL%N_EXTERNAL_WALL_CELLS
    
      IOR0 = SL%IJKW(4, IW)
      IF (ABS(IOR0) == 2) CYCLE            ! 2D: in case of y-boundary cycle

      I    = SL%IJKW(6, IW)
      K    = SL%IJKW(8, IW)
      NOM  = SL%IJKW(9, IW)

      BC_INDEX = SL%PRESSURE_BC_INDEX(IW)

      SELECT CASE(IOR0)
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

      SELECT_STORAGE2D: SELECT CASE (TYPE_STORAGE)
         !!!
         !!! bandwise storage technique
         !!!
         CASE (NSCARC_STORAGE_BANDWISE)

            SELECT CASE(BC_INDEX)
               CASE (DIRICHLET)                    ! set Dirichlet BC's along open boundary cells
                  SL%BA(IC,ID) = SL%BA(IC,ID) - DBC
               !CASE (INTERNAL)                    ! do nothing along internal boundaries (just only debug message)
               CASE (NEUMANN)                      ! set Neumann BC's at all other nodes
                  SL%BA(IC,ID) = SL%BA(IC,ID) + DBC
            END SELECT

         !!!
         !!! compact storage technique
         !!!
         CASE (NSCARC_STORAGE_COMPACT) 

            IP = SL%ROW(IC)
            SELECT CASE(BC_INDEX)
               CASE (DIRICHLET)                    ! set Dirichlet BC's along open boundary cells
                  SL%CA(IP) = SL%CA(IP) - DBC
               !CASE (INTERNAL)                    ! do nothing along internal boundaries (just only debug message)
               CASE (NEUMANN)                      ! set Neumann BC's at all other nodes
                  SL%CA(IP) = SL%CA(IP) + DBC
            END SELECT

      END SELECT SELECT_STORAGE2D
 
   ENDDO WALL_CELL_LOOP2D


!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE

   WALL_CELL_LOOP3D: DO IW = 1, SL%N_EXTERNAL_WALL_CELLS

      IOR0 = SL%IJKW(4, IW)
      I    = SL%IJKW(6, IW)
      J    = SL%IJKW(7, IW)
      K    = SL%IJKW(8, IW)
      NOM  = SL%IJKW(9, IW)

      BC_INDEX = SL%PRESSURE_BC_INDEX(IW)

      SELECT CASE(IOR0)
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

      SELECT_STORAGE3D: SELECT CASE (TYPE_STORAGE)
         !!!
         !!! bandwise storage technique
         !!!
         CASE (NSCARC_STORAGE_BANDWISE)

            SELECT CASE(BC_INDEX)
               CASE (DIRICHLET)                    ! set Dirichlet BC's at open and null boundary cells
                  SL%BA(IC,ID) = SL%BA(IC,ID) - DBC
               !CASE (INTERNAL)                    ! do nothing along internal boundaries (only debug message)
               CASE (NEUMANN)                      ! set Neumann BC's at all other cells
                  SL%BA(IC,ID) = SL%BA(IC,ID) + DBC
            END SELECT

         !!!
         !!! compact storage technique
         !!!
         CASE (NSCARC_STORAGE_COMPACT) 

            IP = SL%ROW(IC)
            SELECT CASE(BC_INDEX)
               CASE (DIRICHLET)                    ! set Dirichlet BC's at open and null boundary cells
                  SL%CA(IP) = SL%CA(IP) - DBC
               !CASE (INTERNAL)                    ! do nothing along internal boundaries (only debug message)
               CASE (NEUMANN)                      ! set Neumann BC's at all other cells
                  SL%CA(IP) = SL%CA(IP) + DBC
            END SELECT

      END SELECT SELECT_STORAGE3D

   ENDDO WALL_CELL_LOOP3D
     
ENDIF

END SUBROUTINE SCARC_SETUP_BOUNDARY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize global 3D-solver methods (cg/mg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_AMG (NM)
INTEGER, INTENT(IN) :: NM
INTEGER :: IERR, NMEASURE, NL
INTEGER, ALLOCATABLE, DIMENSION(:)   :: MEASURE_STD
INTEGER, ALLOCATABLE, DIMENSION(:)   :: CPOINTS_STD, CPOINTS_AGG
INTEGER, ALLOCATABLE, DIMENSION(:)   :: FPOINTS_STD, FPOINTS_AGG
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: COUPLED_STD, COUPLED_AGG

IERR = 0
NL   = NLEVEL_MAX 

!!! Determine number of multigrid levels
GET_NUMBER_OF_LEVELS_LOOP: DO 
   
   SL => SCARC(NM)%LEVEL(NL)

   TYPE_STENCIL    = NSCARC_STENCIL_NPOINT
   TYPE_COARSENING = NSCARC_COARSENING_STANDARD

   SL%N_COUPLED = SL%N_COUPLED_STD

   !!!-------------------------------------------------------------------------------------------------
   !!! Allocate and initialize some auxiliary arrays:
   !!! Define number of couplings (i.e. number of coupled neighboring cells) and offsets between
   !!! Define measure (i.e. total number) and indices of strongly coupled neighbors for all cells
   !!!-------------------------------------------------------------------------------------------------
   ALLOCATE (MEASURE_STD(1:SL%N_CELLS), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'MEASURE_STD', IERR)
   MEASURE_STD = 0
   
   ALLOCATE (COUPLED_STD(1:SL%N_CELLS, 2:SL%N_COUPLED), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'COUPLED', IERR)
   COUPLED_STD = .FALSE.
   
   CALL SCARC_AMG_MEASURE_STD(MEASURE_STD, COUPLED_STD, NM, NL)
   
   
   !!!-------------------------------------------------------------------------------------------------
   !!! Perform standard AMG-coarsening :
   !!! For each cell IC determine the set of all strongly coupled C- and F-points
   !!! strongly coupled C- and F-points are set to '1', otherwise '0'
   !!!-------------------------------------------------------------------------------------------------
   ALLOCATE (CPOINTS_STD(1:SL%N_CELLS), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'CPOINTS_STD', IERR)
   CPOINTS_STD = 0
   
   ALLOCATE (FPOINTS_STD(1:SL%N_CELLS), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'FPOINTS_STD', IERR)
   FPOINTS_STD = 0
   
   CALL SCARC_AMG_COARSENING_STD(COUPLED_STD, CPOINTS_STD, FPOINTS_STD, &
                                 MEASURE_STD, NMEASURE, NM, NL)
   
   !!!-------------------------------------------------------------------------------------------------
   !!! If requested, use additional aggressive coarsening strategy leading to less coarse grid cells :
   !!! Define offsets and indices of strongly coupled neighbors for aggressive coarsening 
   !!! Note1: strongly coupled F-points are set to '1', weakly coupled F-points are set to '2'
   !!! Note2: Aggressive coarsening is only used on finest level
   !!!-------------------------------------------------------------------------------------------------
   IF (NL == NLEVEL_MAX .AND. TYPE_COARSENING == NSCARC_COARSENING_AGGRESSIVE) THEN
   
      TYPE_COARSENING = NSCARC_COARSENING_AGGRESSIVE
      TYPE_STENCIL    = NSCARC_STENCIL_AGGRESSIVE

      SL%N_COUPLED = SL%N_COUPLED_AGG

      CALL SCARC_SETUP_STENCIL(NM, NL)

      ALLOCATE (COUPLED_AGG(1:SL%N_CELLS, 2:SL%N_COUPLED), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_AMG', 'COUPLED_AGG', IERR)
      COUPLED_AGG = .FALSE.
   
      ALLOCATE (CPOINTS_AGG(1:SL%N_CELLS), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_AMG', 'CPOINTS_AGG', IERR)
      CPOINTS_AGG = CPOINTS_STD
   
      ALLOCATE (FPOINTS_AGG(1:SL%N_CELLS), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_AMG', 'FPOINTS_AGG', IERR)
      FPOINTS_AGG = 2*FPOINTS_STD              ! Initially FPOINTS_AGG are thought to be weakly coupled
   
      CALL SCARC_AMG_COARSENING_AGG(COUPLED_STD, FPOINTS_STD, COUPLED_AGG, CPOINTS_AGG, FPOINTS_AGG, NM, NL)
   
   ENDIF
   
   
   !!!-------------------------------------------------------------------------------------------------
   !!! Save information of C- and F-point information and compute interpolation weights
   !!! Note, only {F,C}POINT1 or {F,C}POINT2 must be saved, depending on coarsening strategy 
   !!!-------------------------------------------------------------------------------------------------
   ALLOCATE (SL%COUPLED(1:SL%N_CPOINTS, 2:SL%N_COUPLED), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'SL%CPOINTS', IERR)
   SL%COUPLED = .FALSE.

   ALLOCATE (SL%CPOINTS(1:SL%N_CPOINTS), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'SL%CPOINTS', IERR)
   SL%CPOINTS = 0

   ALLOCATE (SL%FPOINTS(1:SL%N_FPOINTS), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'SL%FPOINTS', IERR)
   SL%FPOINTS = 0

   ALLOCATE (SL%INTERPOL(1:SL%N_CELLS, 1:SL%N_CPOINTS), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'SL%INTERPOL', IERR)
   SL%INTERPOL = 0.0_EB

   ALLOCATE (SL%WEIGHTS2(SL%N_CELLS*SL%N_COUPLED), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'SL%INTERPOL', IERR)
   SL%WEIGHTS2 = 0.0_EB

   ALLOCATE (SL%WEIGHTS2_ROW(SL%N_CELLS+1), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'SL%INTERPOL', IERR)
   SL%WEIGHTS2_ROW = 0.0_EB

   ALLOCATE (SL%WEIGHTS2_COL(SL%N_CELLS*SL%N_COUPLED), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_AMG', 'SL%INTERPOL', IERR)
   SL%WEIGHTS2_COL = 0.0_EB


   SELECT CASE(TYPE_COARSENING)
      CASE(NSCARC_COARSENING_STANDARD)
         CALL SCARC_AMG_INTERPOLATION_STD(COUPLED_STD, CPOINTS_STD, FPOINTS_STD, NM, NL)
      CASE(NSCARC_COARSENING_AGGRESSIVE)
         ! Not yet finished
         !CALL SCARC_AMG_INTERPOLATION_AGG(COUPLED_STD, COUPLED_AGG, CPOINTS_AGG, FPOINTS_AGG, NM, NL)
         DEALLOCATE(FPOINTS_AGG)
         DEALLOCATE(CPOINTS_AGG)
   END SELECT
   
   CALL SCARC_AMG_TRANSFER(COUPLED_STD, NM, NL)

   DEALLOCATE(FPOINTS_STD)
   DEALLOCATE(CPOINTS_STD)
   DEALLOCATE(MEASURE_STD)

WRITE(*,*) 'ACHTUNG, NUR VORRUEBERGEHEND'
WRITE(SCARC_LU,*) 'ACHTUNG, NUR VORRUEBERGEHEND'
   IF (NL > 0) EXIT GET_NUMBER_OF_LEVELS_LOOP
   NL = NL - 1


ENDDO GET_NUMBER_OF_LEVELS_LOOP


END SUBROUTINE SCARC_SETUP_AMG

!!!----------------------------------------------------------------------------------------------------
!!! Setup matrix stencil corresponding to chosen coarsening strategy
!!!----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STENCIL(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IERR

IERR = 0
SL => SCARC(NM)%LEVEL(NL)

SELECT CASE(TYPE_STENCIL)

   !!!-------------------------------------------------------------------------------------------------
   !!! Setup structure of 'standard' matrix stencil on finest level
   !!!   2D : 5-point stencil which couples 5 cells
   !!!   3D : 7-point stencil which couples 7 cells
   !!! N_COUPLED_STD is the number of couplings in the stencil
   !!! ID, ILZ, ILY, ILX, IUX, IUY, IUZ are pointers to offsets in diagonal,  lower x-, y-, z- 
   !!! and upper x-, y-, z-directions
   !!!-------------------------------------------------------------------------------------------------
   CASE(NSCARC_STENCIL_NPOINT)

      IF (TWO_D) THEN
      
         SL%ND_STD  = 1
         SL%NLZ_STD = 2
         SL%NLX_STD = 3
         SL%NUX_STD = 4
         SL%NUZ_STD = 5
      
         ID  => SL%ND_STD
         ILZ => SL%NLZ_STD
         ILX => SL%NLX_STD
         IUX => SL%NUX_STD
         IUZ => SL%NUZ_STD
      
         SL%N_COUPLED_STD = 5
      
         ALLOCATE (SL%STENCIL_STD(SL%N_COUPLED_STD), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_STENCIL', 'STENCIL', IERR)
      
         SL%STENCIL_STD(ID ) =  0
         SL%STENCIL_STD(ILZ) = -SL%IBAR
         SL%STENCIL_STD(ILX) = -1         
         SL%STENCIL_STD(IUX) =  1
         SL%STENCIL_STD(IUZ) =  SL%IBAR
      
      
      ELSE
      
         SL%ND_STD  = 1
         SL%NLZ_STD = 2
         SL%NLY_STD = 3
         SL%NLX_STD = 4
         SL%NUX_STD = 5
         SL%NUY_STD = 6
         SL%NUZ_STD = 7
      
         ID  => SL%ND_STD
         ILZ => SL%NLZ_STD
         ILY => SL%NLY_STD
         ILX => SL%NLX_STD
         IUX => SL%NUX_STD
         IUY => SL%NUY_STD
         IUZ => SL%NUZ_STD
      
         SL%N_COUPLED_STD = 7
      
         ALLOCATE (SL%STENCIL_STD(SL%N_COUPLED_STD), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_STENCIL', 'STENCIL', IERR)
      
         SL%STENCIL_STD(ID ) =  0
         SL%STENCIL_STD(ILZ) = -SL%IBAR * SL%JBAR
         SL%STENCIL_STD(ILY) = -SL%IBAR
         SL%STENCIL_STD(ILX) = -1
         SL%STENCIL_STD(IUX) =  1
         SL%STENCIL_STD(IUY) =  SL%IBAR
         SL%STENCIL_STD(IUZ) =  SL%IBAR * SL%JBAR
      
      ENDIF
      
      
   !!!-------------------------------------------------------------------------------------------------
   !!! Define structure of 'aggressive' matrix stencil on finest level
   !!! N_COUPLED_AGG is the number of couplings in the stencil
   !!! ID, ILZ, ILY, ILX, IUX, IUY, IUZ are pointers to offsets in diagonal,  lower x-, y-, z- 
   !!! and upper x-, y-, z-directions
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_STENCIL_AGGRESSIVE)

      IF (TWO_D) THEN
      
         SL%ND_AGG  = 1
         SL%NLZ_AGG = 2
         SL%NLX_AGG = 3
         SL%NUX_AGG = 4
         SL%NUZ_AGG = 5
      
         ID  => SL%ND_AGG
         ILZ => SL%NLZ_AGG
         ILX => SL%NLX_AGG
         IUX => SL%NUX_AGG
         IUZ => SL%NUZ_AGG
      
         SL%N_COUPLED_AGG = 5
      
         ALLOCATE (SL%STENCIL_AGG(SL%N_COUPLED_AGG), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_STENCIL', 'STENCIL', IERR)
      
         SL%STENCIL_AGG(ID ) =  0
         SL%STENCIL_AGG(ILZ) = -2* SL%IBAR
         SL%STENCIL_AGG(ILX) = -2         
         SL%STENCIL_AGG(IUX) =  2
         SL%STENCIL_AGG(IUZ) =  2* SL%IBAR
      
      
      ELSE
      
         SL%ND_AGG  = 1
         SL%NLZ_AGG = 2
         SL%NLY_AGG = 3
         SL%NLX_AGG = 4
         SL%NUX_AGG = 5
         SL%NUY_AGG = 6
         SL%NUZ_AGG = 7
      
         ID  => SL%ND_AGG
         ILZ => SL%NLZ_AGG
         ILY => SL%NLY_AGG
         ILX => SL%NLX_AGG
         IUX => SL%NUX_AGG
         IUY => SL%NUY_AGG
         IUZ => SL%NUZ_AGG
      
         SL%N_COUPLED_AGG = 7
      
         ALLOCATE (SL%STENCIL_AGG(SL%N_COUPLED_AGG), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_STENCIL', 'STENCIL', IERR)
      
         SL%STENCIL_AGG(ILZ) = -2* SL%IBAR * SL%JBAR
         SL%STENCIL_AGG(ILY) = -2* SL%IBAR
         SL%STENCIL_AGG(ILX) = -2
         SL%STENCIL_AGG(IUX) =  2
         SL%STENCIL_AGG(IUY) =  2* SL%IBAR
         SL%STENCIL_AGG(IUZ) =  2* SL%IBAR * SL%JBAR
         SL%STENCIL_AGG(ID ) =  0
      
      ENDIF
      
END SELECT

END SUBROUTINE SCARC_SETUP_STENCIL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Return matrix value for given connection in standard npoint-matrix stencil
!!! depending on the used storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION MATRIX_VALUE(IC, ICPL, NM, NL)
INTEGER, INTENT(IN) :: IC, ICPL, NM, NL
INTEGER :: ICOL, ICO

SL  => SCARC(NM)%LEVEL(NL)

IF (TYPE_STORAGE == NSCARC_STORAGE_COMPACT) THEN
   DO ICOL = SL%ROW(IC), SL%ROW(IC+1)-1
      ICO = IC + SL%STENCIL_STD(ICPL)
      IF (SL%COL(ICOL) == ICO) EXIT
   ENDDO
   MATRIX_VALUE = SL%CA(ICOL)
ELSE
   MATRIX_VALUE = SL%BA(IC, ICPL)
ENDIF

RETURN
END FUNCTION MATRIX_VALUE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set initial measure of each cell, i.e. the number of connections in the
!!! matrix graph and determine adjacency graph 
!!! Note: This is all related to the 5-point stencil in 2D and the 7-point
!!! stencil in 3D and won't work for general matrices !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_AMG_MEASURE_STD(MEASURE_STD, COUPLED_STD, NM, NL)
INTEGER, INTENT(IN)  :: NM, NL
INTEGER, DIMENSION(:)  , INTENT(OUT)  :: MEASURE_STD
LOGICAL, DIMENSION(:,2:), INTENT(OUT) :: COUPLED_STD
INTEGER  :: IC, ICPL, JCPL, IROW 
REAL(EB) :: AROW(5), EPS, AMG_TOL

EPS     =  1.0E-12_EB
AMG_TOL =  0.25_EB

SL  => SCARC(NM)%LEVEL(NL)

!!! for each subdiagonal entry check size in comparion to remaining subdiagonal entries
CELL_LOOP: DO IC = 1, SL%N_CELLS

   COUPLED_LOOP: DO ICPL = 2, SL%N_COUPLED_STD    ! omit diagonal entry
      
      IF (MATRIX_VALUE(IC, ICPL, NM, NL) > EPS) THEN
         IROW = 1
         DO JCPL = 2, SL%N_COUPLED_STD
            IF (JCPL /= ICPL) THEN
               AROW(IROW) = MATRIX_VALUE(IC,JCPL,NM,NL)
               IROW = IROW + 1
            ENDIF
         ENDDO
         IF (MATRIX_VALUE(IC,ICPL,NM,NL) >= AMG_TOL * MAXVAL(AROW)) THEN
            MEASURE_STD(IC) = MEASURE_STD(IC)+1
            COUPLED_STD(IC, ICPL) = .TRUE.
         ENDIF
      ENDIF
   
   ENDDO COUPLED_LOOP

ENDDO CELL_LOOP

END SUBROUTINE SCARC_AMG_MEASURE_STD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform standard coarsening:
!!! Determine the set of strongly coupled C- and F-points (S1-coupling)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_AMG_COARSENING_STD(COUPLED_STD, CPOINTS_STD, FPOINTS_STD, MEASURE_STD, NMEASURE, NM, NL)
LOGICAL, DIMENSION(:,2:), INTENT(IN)   :: COUPLED_STD
INTEGER, DIMENSION(:)  , INTENT(OUT)   :: CPOINTS_STD, FPOINTS_STD
INTEGER, DIMENSION(:)  , INTENT(INOUT) :: MEASURE_STD
INTEGER, INTENT(IN)   :: NM, NL
INTEGER, INTENT(OUT)  :: NMEASURE
INTEGER :: IC, JC, JCA, ICPL, JCPL, N_POINTS

SL  => SCARC(NM)%LEVEL(NL)

COARSENING_LOOP1: DO

   NMEASURE = MAXVAL(MEASURE_STD)
   IF (NMEASURE == 0) EXIT COARSENING_LOOP1

   CELL_LOOP: DO IC = 1, SL%N_CELLS

      !!! Take first cell with maximum measure as next CPOINT
      IF (MEASURE_STD(IC) /= NMEASURE) THEN
         CYCLE
      ELSE

         CPOINTS_STD(IC) = 1
         MEASURE_STD(IC) = 0

         !!! Determine set of F-points 
         C_COUPLED_LOOP: DO ICPL = 2, SL%N_COUPLED_STD
   
            IF (COUPLED_STD(IC, ICPL)) THEN

               !!! JC is set to be a FPOINT which is no longer measured
               JC = IC + SL%STENCIL_STD(ICPL)

               FPOINTS_STD(JC) = 1
               MEASURE_STD(JC) = 0

               !!!  increase measures of cells JCA adjacent to FPOINTS JC based on strong S1-couplings
               F_COUPLED_LOOP: DO JCPL = 2, SL%N_COUPLED_STD
                  IF (COUPLED_STD(JC, JCPL)) THEN
                     JCA = JC + SL%STENCIL_STD(JCPL)
                     IF (JCA /= IC .AND. (FPOINTS_STD(JCA)==0) .AND. (CPOINTS_STD(JCA)==0)) THEN
                        MEASURE_STD(JCA) = MEASURE_STD(JCA)+1
                        NMEASURE = MAX(NMEASURE, MEASURE_STD(JCA))
                     ENDIF
                  ENDIF
               ENDDO F_COUPLED_LOOP

            ENDIF

         ENDDO C_COUPLED_LOOP

         EXIT CELL_LOOP
      ENDIF
   ENDDO CELL_LOOP

ENDDO COARSENING_LOOP1

!!! determine number of C- and F-Points and check correctness of computation
SL%N_CPOINTS = SCARC_NUMBER_OF_POINTS(CPOINTS_STD)
SL%N_FPOINTS = SCARC_NUMBER_OF_POINTS(FPOINTS_STD)

N_POINTS = SL%N_CPOINTS + SL%N_FPOINTS
IF (N_POINTS /= SL%N_CELLS) THEN
   WRITE(*,*) 'Error in AMG standard coarsening, N_POINTS = ', N_POINTS,' differs from N_CELLS = ', SL%N_CELLS
   STOP
ENDIF

END SUBROUTINE SCARC_AMG_COARSENING_STD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set measure of each cell for aggressive coarsening:
!!! Get connections only for CPOINTS_STD via neighboring FPOINTS_STD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_AMG_COARSENING_AGG(COUPLED_STD, FPOINTS_STD, COUPLED_AGG, CPOINTS_AGG, FPOINTS_AGG, NM, NL)
LOGICAL, DIMENSION(:,2:), INTENT(IN) :: COUPLED_STD
LOGICAL, DIMENSION(:,2:), INTENT(OUT):: COUPLED_AGG
INTEGER, DIMENSION(:)  , INTENT(IN)  :: FPOINTS_STD
INTEGER, DIMENSION(:)  , INTENT(OUT) :: CPOINTS_AGG, FPOINTS_AGG
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, JC, JCA, ICPL, JCPL, N_POINTS

SL => SCARC(NM)%LEVEL(NL)

!!! Take first CPOINT1 as CPOINT2
CELL_LOOP: DO IC = 1, SL%N_CELLS

   IF (CPOINTS_AGG(IC) == 0) CYCLE CELL_LOOP

   !!! For a given CPOINT analyze all surrounding FPOINTS 
   CPOINT_COUPLED_LOOP: DO ICPL = 2, SL%N_COUPLED_AGG

      IF (.NOT.COUPLED_STD(IC, ICPL)) CYCLE CPOINT_COUPLED_LOOP

      JC = IC + SL%STENCIL_STD(ICPL)
      IF (FPOINTS_STD(JC) == 0) THEN
         WRITE(*,*) 'Error in SCARC_AMG_COARSENING_AGG, ',JC,' must be a FPOINT !'
         STOP
      ENDIF

      !!! Analyze all couplings of given FPOINT
      FPOINT_COUPLED_LOOP: DO JCPL = 2, SL%N_COUPLED_STD

         IF (.NOT.COUPLED_STD(JC, JCPL)) CYCLE FPOINT_COUPLED_LOOP

         !!! If F-point has same coupling direction as C-point take adjacent cell as next C-point
         !!! else set adjacent cell to be an weakly coupled FPOINT2 (strongly couplings not changed)
         JCA = JC + SL%STENCIL_STD(JCPL)
         IF (JCPL == ICPL) THEN
            IF (JCA /= IC + SL%STENCIL_AGG(ICPL)) THEN
               WRITE(*,*) 'Error in SCARC_AMG_COARSENING_AGG, ',JCA,' must have offset ',&
                          SL%STENCIL_AGG(ICPL),' !'
            ENDIF
            COUPLED_AGG(IC, ICPL) = .TRUE.
         ELSE IF (JCA /= IC) THEN
            CPOINTS_AGG(JCA) = 0
            FPOINTS_AGG(JCA) = 1
         ENDIF
   
      ENDDO FPOINT_COUPLED_LOOP


   ENDDO CPOINT_COUPLED_LOOP
ENDDO CELL_LOOP

!!! determine number of C- and F-Points and check correctness of computation
SL%N_CPOINTS = SCARC_NUMBER_OF_POINTS(CPOINTS_AGG)
SL%N_FPOINTS = SCARC_NUMBER_OF_POINTS(FPOINTS_AGG)

N_POINTS = SL%N_CPOINTS + SL%N_FPOINTS
IF (N_POINTS /= SL%N_CELLS) THEN
   WRITE(*,*) 'Error in AMG aggressive coarsening, N_POINTS = ', N_POINTS,' differs from N_CELLS = ', SL%N_CELLS
   STOP
ENDIF

END SUBROUTINE SCARC_AMG_COARSENING_AGG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Define interpolation weights for standard coarsening
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_AMG_INTERPOLATION_STD(COUPLED_STD, CPOINTS_STD, FPOINTS_STD, NM, NL)
LOGICAL, DIMENSION(:,2:),INTENT(INOUT) :: COUPLED_STD
INTEGER, DIMENSION(:)  , INTENT(INOUT) :: CPOINTS_STD, FPOINTS_STD
INTEGER, INTENT(IN) :: NM, NL
INTEGER  :: ICP, IFP, IC, ICPL, IP
REAL(EB) :: DSCAL

SL => SCARC(NM)%LEVEL(NL)


!!! store numbers of cpoints and fpoints 
!!! update arrays CPOiNTS_STD and FPOINTS_STD (needed for following loop)
ICP   = 0
IFP   = 0
DO IC = 1, SL%N_CELLS
   IF (CPOINTS_STD(IC) == 1) THEN
      ICP = ICP +1
      SL%CPOINTS(ICP) = IC           
      CPOINTS_STD(IC) = ICP
   ELSE IF (FPOINTS_STD(IC) == 1) THEN
      IFP = IFP + 1
      SL%FPOINTS(IFP) = IC       
      FPOINTS_STD(IC) = IFP
   ELSE
      WRITE(*,*) 'Error in SCARC_AMG_INTERPOLATION_STD: ',IC,' neither CPOINT nor FPOINT !!!' 
      STOP
   ENDIF
ENDDO

IF (ICP /= SL%N_CPOINTS .OR. IFP /= SL%N_FPOINTS) THEN
   WRITE(*,*) 'Error while transition of C- and F-points in level ', NL
   STOP
ENDIF

!!! take full weight for cpoints, compute corresponding interpolation weights for fpoints via coupled C-points
!!! generate corresponding interpolation matrix which is stored using comact storage technique
IP = 1
DO IC = 1, SL%N_CELLS

   IF (CPOINTS_STD(IC) /= 0) THEN

      SL%WEIGHTS2_ROW(IC)  = IP
      SL%WEIGHTS2_COL(IP) = CPOINTS_STD(IC)
      SL%WEIGHTS2(IP)      = 1.0_EB
      IP = IP +1

   ELSE IF (FPOINTS_STD(IC) /= 0) THEN

      SL%WEIGHTS2_ROW(IC) = IP

      DSCAL = 1.0_EB/MATRIX_VALUE(IC, 1, NM, NL)

      DO ICPL = 2, SL%N_COUPLED_STD
         IF (COUPLED_STD(IC, ICPL)) THEN
            SL%WEIGHTS2_COL(IP) = CPOINTS_STD(IC+SL%STENCIL_STD(ICPL))
            SL%WEIGHTS2(IP)      = -MATRIX_VALUE(IC, ICPL, NM, NL) * DSCAL
            IP = IP +1
         ENDIF
      ENDDO

   ENDIF
ENDDO
SL%WEIGHTS2_ROW(SL%N_CELLS+1) = IP

END SUBROUTINE SCARC_AMG_INTERPOLATION_STD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute coarser matrix A_coarse by use of interpolation matrices:
!!!
!!!      A_coarse = I_fine^coarse * A_fine * I_coarse_fine
!!! 
!!! Note the different storage techniques (compact storage technique for 
!!! transfer matrices and coarse matrix, bandwise for finest matrix)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_AMG_TRANSFER (COUPLED_STD, NM, NL)
LOGICAL, DIMENSION(:,2:), INTENT(IN) :: COUPLED_STD
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, ICO, ICP, IERR, ICPL, ICP1, ICP2, IP
REAL (EB), ALLOCATABLE, DIMENSION(:) :: AROW
REAL:: AUX1, AUX2

IERR = 0
SLHI => SCARC(NM)%LEVEL(NL)
SLLO => SCARC(NM)%LEVEL(NL-1)

ALLOCATE (AROW(SLHI%N_CPOINTS), STAT=IERR)
CALL CHKMEMERR ('SCARC_AMG_TRANSFER', 'AROW', IERR)
AROW = 0.0_EB

ALLOCATE (SLLO%CA(SLHI%N_CPOINTS*9), STAT=IERR)
CALL CHKMEMERR ('SCARC_AMG_TRANSFER', 'SLLO%CA', IERR)
SLLO%CA = 0.0_EB

ALLOCATE (SLLO%ROW(SLHI%N_CPOINTS+1), STAT=IERR)
CALL CHKMEMERR ('SCARC_AMG_TRANSFER', 'SLLO%ROW', IERR)
SLLO%ROW = 0

ALLOCATE (SLLO%COL(SLHI%N_CPOINTS*9), STAT=IERR)
CALL CHKMEMERR ('SCARC_AMG_TRANSFER', 'SLLO%COL', IERR)
SLLO%COL = 0


IP  = 1
ICP = 1
DO ICP1 = 1, SLHI%N_CPOINTS
   DO ICP2 = 1, SLHI%N_CPOINTS
      AUX2 = 0.0_EB
      ICP = (ICP1-1)*SLHI%N_CPOINTS + ICP2
      DO IC = 1, SLHI%N_CELLS

         AUX1 = MATRIX_VALUE(IC,1,NM,NL) * CPOINT_WEIGHT(IC, ICP2, NM, NL)

         COUPLED_LOOP: DO ICPL = 2, SLHI%N_COUPLED_STD                
            IF (COUPLED_STD(IC,ICPL)) THEN
               ICO = IC + SLHI%STENCIL_STD(ICPL)
               AUX1 = AUX1 + MATRIX_VALUE(IC,ICPL,NM,NL) * CPOINT_WEIGHT(ICO, ICP2, NM, NL)
            ENDIF
         ENDDO COUPLED_LOOP

         AUX2 = AUX2 + CPOINT_WEIGHT(IC, ICP1, NM, NL) * AUX1

      ENDDO
      AROW(ICP2) = AUX2
      ICP = ICP + 1
   ENDDO


   !!! analyze new matrix line and store it corresponding to compact storage technique:
   !!! (diagonal entry first)
   SLLO%CA(IP) = AROW(ICP1)
   SLLO%ROW(ICP1) = IP
   SLLO%COL(IP)  = ICP1
   IP  = IP + 1
   DO ICP2 = 1, SLHI%N_CPOINTS
      IF (ICP2 /= ICP1 .AND. ABS(AROW(ICP2)) >= 1.0E-12_EB) THEN
         SLLO%CA(IP) = AROW(ICP2)
         SLLO%COL(IP) = ICP2
         IP  = IP + 1
      ENDIF
   ENDDO

ENDDO
SLLO%ROW(SL%N_CPOINTS+1) = IP

DEALLOCATE(AROW)

END SUBROUTINE SCARC_AMG_TRANSFER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if corresponding CPOINT contributes a non-zero interpoation 
!!! weight or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION CPOINT_WEIGHT(IC, ICP, NM, NL)
INTEGER, INTENT(IN):: IC, ICP, NM, NL
INTEGER :: ICOL
LOGICAL :: FOUND

SL => SCARC(NM)%LEVEL(NL)

FOUND = .FALSE.
CPOINT_COL_LOOP: DO ICOL = SL%WEIGHTS2_ROW(IC), SL%WEIGHTS2_ROW(IC+1)-1
   IF (SL%WEIGHTS2_COL(ICOL) == ICP) THEN
      FOUND = .TRUE.
      EXIT CPOINT_COL_LOOP
   ENDIF
ENDDO CPOINT_COL_LOOP

IF (FOUND)  THEN
   CPOINT_WEIGHT =  SL%WEIGHTS2(ICOL)
ELSE
   CPOINT_WEIGHT =  0.0_EB
ENDIF

RETURN
END FUNCTION CPOINT_WEIGHT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine number of activated elements in logical array POINTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION SCARC_NUMBER_OF_POINTS(POINTS)
INTEGER, DIMENSION(:), INTENT(IN) :: POINTS
INTEGER:: INUM, I

INUM=0
DO I = 1, SIZE(POINTS)
  IF (POINTS(I) /= 0 ) INUM=INUM+1
ENDDO
SCARC_NUMBER_OF_POINTS = INUM
RETURN

END FUNCTION SCARC_NUMBER_OF_POINTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Interface for the call of ScaRC-solver with requested storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SOLVER
REAL (EB) :: TNOW_SOLVER

TNOW_SOLVER = SECOND()

SELECT CASE (TYPE_STORAGE)

   !!!-------------------------------------------------------------------------------------------------
   !!! Bandwise storage technique
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_STORAGE_BANDWISE)

      SELECT CASE (TYPE_METHOD)

         !!! Krylov-type method
         CASE (NSCARC_METHOD_KRYLOV)

            SELECT CASE (TYPE_KRYLOV)
               CASE (NSCARC_KRYLOV_CG)                            !!! Krylov-CG
                  CALL SCARC_CG_BANDWISE(NSCARC_SCOPE_MAIN)
               CASE (NSCARC_KRYLOV_BICG)                          !!! Krylov-BICG
                  CALL SCARC_BICG_BANDWISE(NSCARC_SCOPE_MAIN)
            END SELECT

         !!! Multigrid-type method
         CASE (NSCARC_METHOD_MULTIGRID)

            CALL SCARC_MULTIGRID_BANDWISE(NSCARC_SCOPE_MAIN, NSCARC_DUMMY)

      END SELECT

   !!!-------------------------------------------------------------------------------------------------
   !!! Compact storage technique
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_STORAGE_COMPACT)

      SELECT CASE (TYPE_METHOD)

         !!! Krylov-type method
         CASE (NSCARC_METHOD_KRYLOV)

            SELECT CASE (TYPE_KRYLOV)
               CASE (NSCARC_KRYLOV_CG)                            !!! Krylov-CG
                  CALL SCARC_CG_COMPACT(NSCARC_SCOPE_MAIN)
               CASE (NSCARC_KRYLOV_BICG)                          !!! Krylov-BICG
                  CALL SCARC_BICG_COMPACT(NSCARC_SCOPE_MAIN)
            END SELECT

         !!! Multigrid-type method
         CASE (NSCARC_METHOD_MULTIGRID)

            CALL SCARC_MULTIGRID_COMPACT(NSCARC_SCOPE_MAIN, NSCARC_DUMMY)

      END SELECT

END SELECT

TUSED_SCARC(NSCARC_TIME_SOLVER  ,:)=TUSED_SCARC(NSCARC_TIME_SOLVER  ,:)+SECOND()-TNOW_SOLVER
TUSED_SCARC(NSCARC_TIME_COMPLETE,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE,:)+SECOND()-TNOW_SOLVER
END SUBROUTINE SCARC_SOLVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform global CG-method based on global possion-matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CG_BANDWISE(NTYPE)
INTEGER, INTENT(IN) :: NTYPE
INTEGER   :: NM, NL, ITE, NIT, ISTATE, IERR
REAL (EB) :: SIGMA0, SIGMA1, ALPHA0, GAMMA0
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_KRYLOV
CHARACTER(30) :: CROUTINE
LOGICAL, SAVE :: BFIRST = .TRUE.
TYPE (BANDWISE_SYSTEM_POINTER_TYPE), SAVE, POINTER, DIMENSION(:,:) :: P

TNOW_KRYLOV = SECOND()
TYPE_SCOPE  = NTYPE

!!!----------------------------------------------------------------------------------------------------
!!! Initialize data structures
!!!----------------------------------------------------------------------------------------------------
!!! Define name of routine corresponding to current scope and set working level
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CROUTINE = 'SCARC_GLOBAL_CG_BANDWISE'
   NL = NLEVEL_MAX
ELSE
   CROUTINE = 'SCARC_COARSE_CG_BANDWISE'
   NL = NLEVEL_MIN
ENDIF

!!! Point to all needed vectors of corresponding level NL 
IF (BFIRST) THEN

   ALLOCATE(P(NMESHES_MIN:NMESHES_MAX, NL:NL), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'CG', IERR)

   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%NC => SCARC(NM)%LEVEL(NL)%N_CELLS
      P(NM,NL)%NA => SCARC(NM)%LEVEL(NL)%N_MATRIX_ENTRIES
      P(NM,NL)%NW => SCARC(NM)%LEVEL(NL)%N_EXTERNAL_WALL_CELLS
      P(NM,NL)%NX => SCARC(NM)%LEVEL(NL)%IBAR
      P(NM,NL)%NY => SCARC(NM)%LEVEL(NL)%JBAR
      P(NM,NL)%NZ => SCARC(NM)%LEVEL(NL)%KBAR
      P(NM,NL)%A  => SCARC(NM)%LEVEL(NL)%BA
      P(NM,NL)%X  => SCARC(NM)%LEVEL(NL)%BX
      P(NM,NL)%F  => SCARC(NM)%LEVEL(NL)%BF
      P(NM,NL)%D  => SCARC(NM)%LEVEL(NL)%BD
      P(NM,NL)%G  => SCARC(NM)%LEVEL(NL)%BG
      P(NM,NL)%Y  => SCARC(NM)%LEVEL(NL)%BY
      P(NM,NL)%R  => SCARC(NM)%LEVEL(NL)%BR
   ENDDO

   BFIRST = .FALSE.

ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! (Re-)Initiaize working vectors and right hand side
!!!  - In case of CG as main solver for pressure solution (only finest level): 
!!!      ---> initialize rhs with data from calling routine including boundary conditions
!!!  - In case of CG as caorse grid solver (only coarsest level):
!!!      ---> take rhs as given from restriction routine
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%D = 0.0_EB
      P(NM,NL)%G = 0.0_EB                   
      P(NM,NL)%Y = 0.0_EB
      P(NM,NL)%R = 0.0_EB
      CALL SCARC_INITIALIZE_SOLUTION(P(NM,NL)%X, P(NM,NL)%F, NLEVEL_MAX)
   ENDDO
ENDIF

!!! define iterations parameters
EPS = SCARC_KRYLOV_ACCURACY
NIT = SCARC_KRYLOV_ITERATIONS
ITE = 0

!!!----------------------------------------------------------------------------------------------------
!!! Compute initial defect and perform initial preconditioning
!!!----------------------------------------------------------------------------------------------------
!!! Calculate initial matrix vector product R: = A*X 
DO NM = NMESHES_MIN, NMESHES_MAX
   CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%X, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
ENDDO
CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_R, NL)

!!! Calculate initial residual R := B - D and get its global L2-norm 
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%R  = -P(NM,NL)%F + P(NM,NL)%R
   CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
ENDDO
RESIN = SCARC_GLOBAL_L2NORM(NL)


!!! Print information about norm of initial residual
CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)
   
!!! Perform initial preconditioning
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%G = P(NM,NL)%R
ENDDO

IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
   CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_G)
ELSE
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_PRECONDITIONER(P(NM,NL)%A, P(NM,NL)%R, P(NM,NL)%G, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
ENDIF

!!! get initial scaling factor and direction
DO NM = NMESHES_MIN, NMESHES_MAX
   CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%G, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   P(NM,NL)%D = -P(NM,NL)%G
ENDDO
SIGMA0 = SCARC_GLOBAL_SCALPROD(NL)


!!!----------------------------------------------------------------------------------------------------
!!! Conjugate gradient looping
!!!----------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, NIT
 
   !!! get new matrix vector product Y = A*D
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%D, P(NM,NL)%Y, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_D, NSCARC_VECTOR_Y, NL)

   !!! get scaling factor for descent direction
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%Y, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
   ALPHA0 = SCARC_GLOBAL_SCALPROD(NL)
   ALPHA0 = SIGMA0 / ALPHA0

   !!! get new descent direction
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%X = ALPHA0 * P(NM,NL)%D + P(NM,NL)%X
      P(NM,NL)%R = ALPHA0 * P(NM,NL)%Y + P(NM,NL)%R
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
   RES = SCARC_GLOBAL_L2NORM(NL)
 
   !!! 
   !!! Convergence or divergence ?
   !!!
   ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
   IF (ISTATE /= NSCARC_LOOP_PROCEED) EXIT CG_LOOP
 
   !!! Perform preconditioning
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM, NL)%G = P(NM, NL)%R
   ENDDO

   IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
      CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_G)
   ELSE
      DO NM = NMESHES_MIN, NMESHES_MAX
         CALL SCARC_PRECONDITIONER(P(NM,NL)%A, P(NM,NL)%R, P(NM,NL)%G, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
      ENDDO
   ENDIF

   !!! Compute new scaling factors
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%G, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
   SIGMA1 = SCARC_GLOBAL_SCALPROD(NL)
   GAMMA0 = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1

   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%D = -P(NM,NL)%G + GAMMA0 * P(NM,NL)%D
   ENDDO

ENDDO CG_LOOP
 
!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

!!!----------------------------------------------------------------------------------------------------
!!! In case of CG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_TRANSFER_DATA(P(NM,NLEVEL_MAX)%X, NM)
      CALL SCARC_UPDATE_GHOSTCELLS(NM)
   ENDDO
   CALL SCARC_EXCHANGE_BDRY(NLEVEL_MAX)
ENDIF

TUSED_SCARC(NSCARC_TIME_KRYLOV  ,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV  ,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_COMPLETE,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_CG_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform global CG-method based on global possion-matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CG_COMPACT(NTYPE)
INTEGER, INTENT(IN) :: NTYPE
INTEGER   :: NM, NL, ITE, NIT, ISTATE, IERR
REAL (EB) :: SIGMA0, SIGMA1, ALPHA0, GAMMA0
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_KRYLOV
CHARACTER(30) :: CROUTINE
LOGICAL, SAVE :: BFIRST = .TRUE.
TYPE (COMPACT_SYSTEM_POINTER_TYPE), SAVE, POINTER, DIMENSION(:,:) :: P

TNOW_KRYLOV = SECOND()
TYPE_SCOPE  = NTYPE

!!!----------------------------------------------------------------------------------------------------
!!! Initialize data structures
!!!----------------------------------------------------------------------------------------------------
!!! Define name of routine corresponding to current scope and set working level
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CROUTINE = 'SCARC_GLOBAL_CG_COMPACT'
   NL = NLEVEL_MAX
ELSE
   CROUTINE = 'SCARC_COARSE_CG_COMPACT'
   NL = NLEVEL_MIN
ENDIF

!!! Point to all needed vectors of corresponding level NL and
IF (BFIRST) THEN

   ALLOCATE(P(NMESHES_MIN:NMESHES_MAX, NL:NL), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'CG', IERR)

   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%NC  => SCARC(NM)%LEVEL(NL)%N_CELLS
      P(NM,NL)%NA  => SCARC(NM)%LEVEL(NL)%N_MATRIX_ENTRIES
      P(NM,NL)%NW  => SCARC(NM)%LEVEL(NL)%N_EXTERNAL_WALL_CELLS
      P(NM,NL)%NX  => SCARC(NM)%LEVEL(NL)%IBAR
      P(NM,NL)%NY  => SCARC(NM)%LEVEL(NL)%JBAR
      P(NM,NL)%NZ  => SCARC(NM)%LEVEL(NL)%KBAR
      P(NM,NL)%A   => SCARC(NM)%LEVEL(NL)%CA
      P(NM,NL)%COL => SCARC(NM)%LEVEL(NL)%COL
      P(NM,NL)%ROW => SCARC(NM)%LEVEL(NL)%ROW
      P(NM,NL)%X   => SCARC(NM)%LEVEL(NL)%CX
      P(NM,NL)%F   => SCARC(NM)%LEVEL(NL)%CF
      P(NM,NL)%D   => SCARC(NM)%LEVEL(NL)%CD
      P(NM,NL)%G   => SCARC(NM)%LEVEL(NL)%CG
      P(NM,NL)%Y   => SCARC(NM)%LEVEL(NL)%CY
      P(NM,NL)%R   => SCARC(NM)%LEVEL(NL)%CR
   ENDDO

   BFIRST = .FALSE.

ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! (Re-)Initiaize working vectors and right hand side
!!!  - In case of CG as main solver for pressure solution (only finest level): 
!!!      ---> initialize rhs with data from calling routine including boundary conditions
!!!  - In case of CG as caorse grid solver (only coarsest level):
!!!      ---> take rhs as given from restriction routine
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%D = 0.0_EB
      P(NM,NL)%G = 0.0_EB                   
      P(NM,NL)%Y = 0.0_EB
      P(NM,NL)%R = 0.0_EB
      CALL SCARC_INITIALIZE_SOLUTION(P(NM,NL)%X, P(NM,NL)%F, NLEVEL_MAX)
   ENDDO
ENDIF

!!! define iterations parameters
EPS = SCARC_KRYLOV_ACCURACY
NIT = SCARC_KRYLOV_ITERATIONS
ITE = 0

!!!----------------------------------------------------------------------------------------------------
!!! Compute initial defect and perform initial preconditioning
!!!----------------------------------------------------------------------------------------------------
!!! Calculate initial matrix vector D: = A*X product 
DO NM = NMESHES_MIN, NMESHES_MAX
   CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%X, P(NM,NL)%R, &
                            P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
ENDDO
CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_R, NL)

!!! Calculate initial residual R := B - D and get its global L2-norm 
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%R  = -P(NM,NL)%F + P(NM,NL)%R
   CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
ENDDO
RESIN = SCARC_GLOBAL_L2NORM(NL)

!!! Print information about norm of initial residual
CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)
   
!!! Perform initial preconditioning
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%G = P(NM,NL)%R
ENDDO

IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
   CALL SCARC_MULTIGRID_COMPACT (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_G)
ELSE
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_PRECONDITIONER(P(NM,NL)%A, P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%R, P(NM,NL)%G, &
                                P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
ENDIF


!!! get initial scaling factor and direction
DO NM = NMESHES_MIN, NMESHES_MAX
   CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%G, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   P(NM,NL)%D = -P(NM,NL)%G
ENDDO
SIGMA0 = SCARC_GLOBAL_SCALPROD(NL)


!!!----------------------------------------------------------------------------------------------------
!!! Defect correction loop
!!!----------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, NIT
 
   !!! get new matrix vector product Y = A*D
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%D, P(NM,NL)%Y, &
                               P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_D, NSCARC_VECTOR_Y, NL)

   !!! get scaling factor for descent direction
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%Y, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
   ALPHA0 = SCARC_GLOBAL_SCALPROD(NL)
   ALPHA0 = SIGMA0 / ALPHA0

   !!! get new descent direction
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%X = ALPHA0 * P(NM,NL)%D + P(NM,NL)%X
      P(NM,NL)%R = ALPHA0 * P(NM,NL)%Y + P(NM,NL)%R
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
   RES = SCARC_GLOBAL_L2NORM(NL)
 
   !!! 
   !!! Convergence or divergence ?
   !!! 
   ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
   IF (ISTATE /= NSCARC_LOOP_PROCEED) EXIT CG_LOOP
 
   !!! Perform preconditioning
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM, NL)%G = P(NM, NL)%R
   ENDDO

   IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
      CALL SCARC_MULTIGRID_COMPACT (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_G)
   ELSE
      DO NM = NMESHES_MIN, NMESHES_MAX
         CALL SCARC_PRECONDITIONER(P(NM,NL)%A , P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%R, P(NM,NL)%G, &
                                   P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ , NM) 
      ENDDO
   ENDIF

   !!! Compute new scaling factors
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%G, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
   SIGMA1 = SCARC_GLOBAL_SCALPROD(NL)
   GAMMA0 = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1

   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%D = -P(NM,NL)%G + GAMMA0 * P(NM,NL)%D
   ENDDO

ENDDO CG_LOOP
 
!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

!!!----------------------------------------------------------------------------------------------------
!!! In case of CG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_TRANSFER_DATA(P(NM,NLEVEL_MAX)%X, NM)
      CALL SCARC_UPDATE_GHOSTCELLS(NM)
   ENDDO
   CALL SCARC_EXCHANGE_BDRY(NLEVEL_MAX)
ENDIF

TUSED_SCARC(NSCARC_TIME_KRYLOV  ,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV  ,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_COMPLETE,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_CG_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform global BICGstab-method based on global possion-matrix - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_BICG_BANDWISE(NTYPE)
INTEGER, INTENT(IN) :: NTYPE
INTEGER   :: NM, NL, ITE, NIT, ISTATE, IERR
REAL (EB) :: ALPHA0, ALPHA1, ALPHA2, RHO0, RHO1, DTHETA, DBETA
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_KRYLOV
CHARACTER(30) :: CROUTINE
LOGICAL, SAVE :: BFIRST = .TRUE.
TYPE (BANDWISE_SYSTEM_POINTER_TYPE), SAVE, POINTER, DIMENSION(:,:) :: P
 
TNOW_KRYLOV = SECOND()
TYPE_SCOPE  = NTYPE

!!!----------------------------------------------------------------------------------------------------
!!! Initialize data structures
!!!----------------------------------------------------------------------------------------------------
!!! Define name of routine corresponding to current scope and set working level
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CROUTINE = 'SCARC_GLOBAL_BICG_BANDWISE'
   NL = NLEVEL_MAX
ELSE
   CROUTINE = 'SCARC_COARSE_BICG_BANDWISE'
   NL = NLEVEL_MIN
ENDIF

!!! Point to all needed vectors of corresponding level NL and
IF (BFIRST) THEN

   ALLOCATE(P(NMESHES_MIN:NMESHES_MAX, NL:NL), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'CG', IERR)

   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%NC => SCARC(NM)%LEVEL(NL)%N_CELLS
      P(NM,NL)%NA => SCARC(NM)%LEVEL(NL)%N_MATRIX_ENTRIES
      P(NM,NL)%NW => SCARC(NM)%LEVEL(NL)%N_EXTERNAL_WALL_CELLS
      P(NM,NL)%NX => SCARC(NM)%LEVEL(NL)%IBAR
      P(NM,NL)%NY => SCARC(NM)%LEVEL(NL)%JBAR
      P(NM,NL)%NZ => SCARC(NM)%LEVEL(NL)%KBAR
      P(NM,NL)%A  => SCARC(NM)%LEVEL(NL)%BA
      P(NM,NL)%X  => SCARC(NM)%LEVEL(NL)%BX
      P(NM,NL)%F  => SCARC(NM)%LEVEL(NL)%BF
      P(NM,NL)%D  => SCARC(NM)%LEVEL(NL)%BD
      P(NM,NL)%G  => SCARC(NM)%LEVEL(NL)%BG
      P(NM,NL)%Y  => SCARC(NM)%LEVEL(NL)%BY
      P(NM,NL)%R  => SCARC(NM)%LEVEL(NL)%BR
      P(NM,NL)%Z  => SCARC(NM)%LEVEL(NL)%BZ
   ENDDO

   BFIRST = .FALSE.

ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! (Re-)Initiaize working vectors and right hand side
!!!  - In case of CG as main solver for pressure solution (only finest level): 
!!!      ---> initialize rhs with data from calling routine including boundary conditions
!!!  - In case of CG as caorse grid solver (only coarsest level):
!!!      ---> take rhs as given from restriction routine
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%D =  0.0_EB
      P(NM,NL)%G =  0.0_EB                   
      P(NM,NL)%Y =  0.0_EB
      P(NM,NL)%R =  0.0_EB
      P(NM,NL)%Z =  0.0_EB
      CALL SCARC_INITIALIZE_SOLUTION(P(NM,NL)%X, P(NM,NL)%F, NLEVEL_MAX)
   ENDDO
ENDIF

!!! define iterations parameters
EPS    = SCARC_KRYLOV_ACCURACY
NIT    = SCARC_KRYLOV_ITERATIONS
ITE    = 0
ALPHA0 = 1.0_EB
RHO0   = 1.0_EB
RHO1   = 0.0_EB
DBETA  = 0.0_EB
DTHETA = 1.0_EB
   
!!!----------------------------------------------------------------------------------------------------
!!! Compute initial defect and perform initial preconditioning
!!!----------------------------------------------------------------------------------------------------
!!! Perform first initial preconditioning R := PRECON(F)
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%R = P(NM,NL)%F
ENDDO

IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
   CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_R)
ELSE
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_PRECONDITIONER(P(NM,NL)%A, P(NM,NL)%R, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
ENDIF

!!! Compute initial matrix-vector product R := A*R
DO NM = NMESHES_MIN, NMESHES_MAX
   CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%X, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
ENDDO
CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_R, NL)

!!! Perform second initial preconditioning R := PRECON(F)
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%R = P(NM,NL)%F - P(NM,NL)%R
ENDDO

IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
   CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_R)
ELSE
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_PRECONDITIONER(P(NM,NL)%A , P(NM,NL)%R , P(NM,NL)%R , &
                                P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
ENDIF

!!! compute norm of initial residual
DO NM = NMESHES_MIN, NMESHES_MAX
   CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
ENDDO
RESIN = SCARC_GLOBAL_L2NORM(NL)

!!! print initial convergence info
CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)
   
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%G=P(NM,NL)%R
ENDDO

!!!----------------------------------------------------------------------------------------------------
!!! Start BICG-loop 
!!!----------------------------------------------------------------------------------------------------
BICG_LOOP: DO ITE = 1, NIT

   !!! Compute global scalar product
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%G, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   RHO1  = SCARC_GLOBAL_SCALPROD(NL)
   DBETA = (RHO1*DTHETA)/(RHO0*ALPHA0)
   RHO0  = RHO1

   !!! Compute new directions and global matrix-vector product
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%Z = P(NM,NL)%R + DBETA * P(NM,NL)%Z
      P(NM,NL)%Z = -DBETA * ALPHA0 * P(NM,NL)%Y + P(NM,NL)%Z
      CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%Z, P(NM,NL)%Y, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_Z, NSCARC_VECTOR_Y, NL)

   !!! Perform first preconditioning and compute global scalarproduct
   IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
      CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_Y)
   ELSE
      DO NM = NMESHES_MIN, NMESHES_MAX
         CALL SCARC_PRECONDITIONER(P(NM,NL)%A , P(NM,NL)%Y , P(NM,NL)%Y , &
                                   P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
      ENDDO
   ENDIF

   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%G, P(NM,NL)%Y, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   DTHETA = SCARC_GLOBAL_SCALPROD(NL)
   DTHETA = RHO1/DTHETA

   !!! Set new direction and perform second matrix-vector multiplication 
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%R = -DTHETA * P(NM,NL)%Y + P(NM,NL)%R
      CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%R, P(NM,NL)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_R, NSCARC_VECTOR_D, NL)

   !!! perform second preconditioning and compute global scalarproducts
   IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
      CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_D)
   ELSE
      DO NM = NMESHES_MIN, NMESHES_MAX
         CALL SCARC_PRECONDITIONER(P(NM,NL)%A , P(NM,NL)%D , P(NM,NL)%D , &
                                   P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
      ENDDO
   ENDIF

   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   ALPHA1 = SCARC_GLOBAL_SCALPROD(NL)

   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   ALPHA2 = SCARC_GLOBAL_SCALPROD(NL)
   ALPHA0 = ALPHA1/ALPHA2

   !!! Determine new descent directions and get norm of new residual
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%X =  DTHETA * P(NM,NL)%Z + P(NM,NL)%X
      P(NM,NL)%X =  ALPHA0 * P(NM,NL)%R + P(NM,NL)%X
      P(NM,NL)%R = -ALPHA0 * P(NM,NL)%D + P(NM,NL)%R
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   RES = SCARC_GLOBAL_L2NORM(NL)

   !!! Leave BICG-loop in case of convergence or divergence 
   ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
   IF (ISTATE /= NSCARC_LOOP_PROCEED) EXIT BICG_LOOP
 
ENDDO BICG_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

!!!----------------------------------------------------------------------------------------------------
!!! In case of BICG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_TRANSFER_DATA(P(NM,NLEVEL_MAX)%X, NM)
      CALL SCARC_UPDATE_GHOSTCELLS(NM)
   ENDDO
   CALL SCARC_EXCHANGE_BDRY(NLEVEL_MAX)
ENDIF

TUSED_SCARC(NSCARC_TIME_KRYLOV  ,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV  ,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_COMPLETE,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_BICG_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform global BICGstab-method based on global possion-matrix - compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_BICG_COMPACT(NTYPE)
INTEGER, INTENT(IN) :: NTYPE
INTEGER   :: NM, NL, ITE, NIT, ISTATE, IERR
REAL (EB) :: ALPHA0, ALPHA1, ALPHA2, RHO0, RHO1, DTHETA, DBETA
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_KRYLOV
CHARACTER(30) :: CROUTINE
LOGICAL, SAVE :: BFIRST = .TRUE.
TYPE (COMPACT_SYSTEM_POINTER_TYPE), SAVE, POINTER, DIMENSION(:,:) :: P
 
TNOW_KRYLOV = SECOND()
TYPE_SCOPE  = NTYPE

!!!----------------------------------------------------------------------------------------------------
!!! Initialize data structures
!!!----------------------------------------------------------------------------------------------------
!!! Define name of routine corresponding to current scope and set working level
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CROUTINE = 'SCARC_GLOBAL_BICG_BANDWISE'
   NL = NLEVEL_MAX
ELSE
   CROUTINE = 'SCARC_COARSE_BICG_BANDWISE'
   NL = NLEVEL_MIN
ENDIF

!!! Point to all needed vectors of corresponding level NL and
IF (BFIRST) THEN

   ALLOCATE(P(NMESHES_MIN:NMESHES_MAX, NL:NL), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'CG', IERR)

   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%NC => SCARC(NM)%LEVEL(NL)%N_CELLS
      P(NM,NL)%NA => SCARC(NM)%LEVEL(NL)%N_MATRIX_ENTRIES
      P(NM,NL)%NW => SCARC(NM)%LEVEL(NL)%N_EXTERNAL_WALL_CELLS
      P(NM,NL)%NX => SCARC(NM)%LEVEL(NL)%IBAR
      P(NM,NL)%NY => SCARC(NM)%LEVEL(NL)%JBAR
      P(NM,NL)%NZ => SCARC(NM)%LEVEL(NL)%KBAR
      P(NM,NL)%A  => SCARC(NM)%LEVEL(NL)%CA
      P(NM,NL)%COL=> SCARC(NM)%LEVEL(NL)%COL
      P(NM,NL)%ROW=> SCARC(NM)%LEVEL(NL)%ROW
      P(NM,NL)%X  => SCARC(NM)%LEVEL(NL)%CX
      P(NM,NL)%F  => SCARC(NM)%LEVEL(NL)%CF
      P(NM,NL)%D  => SCARC(NM)%LEVEL(NL)%CD
      P(NM,NL)%G  => SCARC(NM)%LEVEL(NL)%CG
      P(NM,NL)%Y  => SCARC(NM)%LEVEL(NL)%CY
      P(NM,NL)%R  => SCARC(NM)%LEVEL(NL)%CR
      P(NM,NL)%Z  => SCARC(NM)%LEVEL(NL)%CZ
   ENDDO

   BFIRST = .FALSE.

ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! (Re-)Initiaize working vectors and right hand side
!!!  - In case of CG as main solver for pressure solution (only finest level): 
!!!      ---> initialize rhs with data from calling routine including boundary conditions
!!!  - In case of CG as caorse grid solver (only coarsest level):
!!!      ---> take rhs as given from restriction routine
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%D =  0.0_EB
      P(NM,NL)%G =  0.0_EB                   
      P(NM,NL)%Y =  0.0_EB
      P(NM,NL)%R =  0.0_EB
      P(NM,NL)%Z =  0.0_EB
      CALL SCARC_INITIALIZE_SOLUTION(P(NM,NL)%X, P(NM,NL)%F, NLEVEL_MAX)
   ENDDO
ENDIF

!!! Define iteration parameters
EPS    = SCARC_KRYLOV_ACCURACY
NIT    = SCARC_KRYLOV_ITERATIONS
ITE    = 0
ALPHA0 = 1.0_EB
RHO0   = 1.0_EB
RHO1   = 0.0_EB
DBETA  = 0.0_EB
DTHETA = 1.0_EB
   
!!!----------------------------------------------------------------------------------------------------
!!! Compute initial defect and perform initial preconditioning
!!!----------------------------------------------------------------------------------------------------
!!! Perform first initial preconditioning R := PRECON(F)
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%R = P(NM,NL)%F
ENDDO

IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
   CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_R)
ELSE
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_PRECONDITIONER(P(NM,NL)%A , P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%R, P(NM,NL)%R, &
                                P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ , NM)
   ENDDO
ENDIF

!!! Compute initial matrix-vector product R := A*R
DO NM = NMESHES_MIN, NMESHES_MAX
   CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%X, P(NM,NL)%R, &
                            P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
ENDDO
CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_R, NL)


!!! Perform second initial preconditioning R := PRECON(F)
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%R = P(NM,NL)%F - P(NM,NL)%R
ENDDO

IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
   CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_R)
ELSE
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_PRECONDITIONER(P(NM,NL)%A, P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%R, P(NM,NL)%R, &
                                P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
ENDIF

!!! compute norm of initial residual
DO NM = NMESHES_MIN, NMESHES_MAX
   CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
ENDDO
RESIN = SCARC_GLOBAL_L2NORM(NL)

!!! print initial convergence info
CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)
   
DO NM = NMESHES_MIN, NMESHES_MAX
   P(NM,NL)%G=P(NM,NL)%R
ENDDO

!!!----------------------------------------------------------------------------------------------------
!!! Start BICG-loop 
!!!----------------------------------------------------------------------------------------------------
BICG_LOOP: DO ITE = 1, NIT

   !!! Compute global scalar product
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%G, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   RHO1  = SCARC_GLOBAL_SCALPROD(NL)
   DBETA = (RHO1*DTHETA)/(RHO0*ALPHA0)
   RHO0  = RHO1

   !!! Compute new directions and global matrix-vector product
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%Z = P(NM,NL)%R + DBETA * P(NM,NL)%Z
      P(NM,NL)%Z = -DBETA * ALPHA0 * P(NM,NL)%Y + P(NM,NL)%Z
      CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%Z, P(NM,NL)%Y, &
                               P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_Z, NSCARC_VECTOR_Y, NL)

   !!! Perform first preconditioning and compute global scalarproduct
   IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
      CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_Y)
   ELSE
      DO NM = NMESHES_MIN, NMESHES_MAX
         CALL SCARC_PRECONDITIONER(P(NM,NL)%A , P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%Y, P(NM,NL)%Y, &
                                   P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ , NM)
      ENDDO
   ENDIF

   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%G, P(NM,NL)%Y, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   DTHETA = SCARC_GLOBAL_SCALPROD(NL)
   DTHETA = RHO1/DTHETA

   !!! Set new direction and perform second matrix-vector multiplication 
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%R = -DTHETA * P(NM,NL)%Y + P(NM,NL)%R
      CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A , P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%R, P(NM,NL)%D, &
                               P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_R, NSCARC_VECTOR_D, NL)

   !!! perform second preconditioning and compute global scalarproducts
   IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
      CALL SCARC_MULTIGRID_BANDWISE (NSCARC_SCOPE_PRECON, NSCARC_VECTOR_D)
   ELSE
      DO NM = NMESHES_MIN, NMESHES_MAX
         CALL SCARC_PRECONDITIONER(P(NM,NL)%A, P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%D, P(NM,NL)%D, &
                                   P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
      ENDDO
   ENDIF

   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   ALPHA1 = SCARC_GLOBAL_SCALPROD(NL)

   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   ALPHA2 = SCARC_GLOBAL_SCALPROD(NL)
   ALPHA0 = ALPHA1/ALPHA2

   !!! Determine new descent directions and get norm of new residual
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%X =  DTHETA * P(NM,NL)%Z + P(NM,NL)%X
      P(NM,NL)%X =  ALPHA0 * P(NM,NL)%R + P(NM,NL)%X
      P(NM,NL)%R = -ALPHA0 * P(NM,NL)%D + P(NM,NL)%R
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%R, P(NM,NL)%R, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   RES = SCARC_GLOBAL_L2NORM(NL)

   !!! Leave BICG-loop in case of convergence or divergence 
   ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
   IF (ISTATE /= NSCARC_LOOP_PROCEED) EXIT BICG_LOOP
 
ENDDO BICG_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

!!!----------------------------------------------------------------------------------------------------
!!! In case of BICG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_TRANSFER_DATA(P(NM,NLEVEL_MAX)%X, NM)
      CALL SCARC_UPDATE_GHOSTCELLS(NM)
   ENDDO
   CALL SCARC_EXCHANGE_BDRY(NLEVEL_MAX)
ENDIF

TUSED_SCARC(NSCARC_TIME_KRYLOV  ,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV  ,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_COMPLETE,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_BICG_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform geometric multigrid method based on global possion-matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MULTIGRID_BANDWISE(NTYPE, NRHS)
INTEGER, INTENT(IN) :: NTYPE, NRHS
INTEGER   :: NM, NL, ITE, NIT, ISTATE, IERR, ICYCLE
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_MULTIGRID
CHARACTER(30):: CROUTINE
LOGICAL, SAVE :: BFIRST = .TRUE.
TYPE (BANDWISE_SYSTEM_POINTER_TYPE), SAVE, POINTER, DIMENSION(:,:) :: P

TNOW_MULTIGRID = SECOND()
TYPE_SCOPE     = NTYPE

!!!----------------------------------------------------------------------------------------------------
!!! Initialize pointers to system variables corresponding to the type of the method
!!! (main solver / preconditioner)
!!! Note: This is only done at the very first call of the routine
!!!----------------------------------------------------------------------------------------------------
IF (BFIRST) THEN

   ALLOCATE(P(NMESHES_MIN:NMESHES_MAX, NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
   CALL CHKMEMERR ('SCARC_MULTIGRID_BANDWISE', 'P', IERR)

   SELECT CASE (TYPE_SCOPE)

      CASE(NSCARC_SCOPE_MAIN)
         CROUTINE = 'SCARC_GLOBAL_MG_BANDWISE'
         DO NM = NMESHES_MIN, NMESHES_MAX
            DO NL=NLEVEL_MIN,NLEVEL_MAX
               P(NM,NL)%NC => SCARC(NM)%LEVEL(NL)%N_CELLS
               P(NM,NL)%NA => SCARC(NM)%LEVEL(NL)%N_MATRIX_ENTRIES
               P(NM,NL)%NW => SCARC(NM)%LEVEL(NL)%N_EXTERNAL_WALL_CELLS
               P(NM,NL)%NX => SCARC(NM)%LEVEL(NL)%IBAR
               P(NM,NL)%NY => SCARC(NM)%LEVEL(NL)%JBAR
               P(NM,NL)%NZ => SCARC(NM)%LEVEL(NL)%KBAR
               P(NM,NL)%A  => SCARC(NM)%LEVEL(NL)%BA
               P(NM,NL)%X  => SCARC(NM)%LEVEL(NL)%BX
               P(NM,NL)%F  => SCARC(NM)%LEVEL(NL)%BF
               P(NM,NL)%D  => SCARC(NM)%LEVEL(NL)%BD
            ENDDO
         ENDDO

      CASE(NSCARC_SCOPE_PRECON)
         CROUTINE = 'SCARC_PRECON_MG_BANDWISE'
         DO NM = NMESHES_MIN, NMESHES_MAX
            DO NL=NLEVEL_MIN,NLEVEL_MAX
               P(NM,NL)%NC => SCARC(NM)%LEVEL(NL)%N_CELLS
               P(NM,NL)%NA => SCARC(NM)%LEVEL(NL)%N_MATRIX_ENTRIES
               P(NM,NL)%NW => SCARC(NM)%LEVEL(NL)%N_EXTERNAL_WALL_CELLS
               P(NM,NL)%NX => SCARC(NM)%LEVEL(NL)%IBAR
               P(NM,NL)%NY => SCARC(NM)%LEVEL(NL)%JBAR
               P(NM,NL)%NZ => SCARC(NM)%LEVEL(NL)%KBAR
               P(NM,NL)%A  => SCARC(NM)%LEVEL(NL)%BA
               P(NM,NL)%X  => SCARC(NM)%LEVEL(NL)%BX2
               P(NM,NL)%D  => SCARC(NM)%LEVEL(NL)%BD2
               P(NM,NL)%X  = 0.0_EB
               SELECT CASE(NRHS)
                  CASE(NSCARC_VECTOR_D)
                     P(NM,NL)%F => SCARC(NM)%LEVEL(NL)%BD
                  CASE(NSCARC_VECTOR_G)
                     P(NM,NL)%F => SCARC(NM)%LEVEL(NL)%BG
                  CASE(NSCARC_VECTOR_R)
                     P(NM,NL)%F => SCARC(NM)%LEVEL(NL)%BR
                  CASE(NSCARC_VECTOR_Y)
                     P(NM,NL)%F => SCARC(NM)%LEVEL(NL)%BY
               END SELECT
            ENDDO
         ENDDO
   END SELECT
            
   BFIRST = .FALSE.

ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! In case of global multigrid:
!!!     ---> initialize rhs with data from calling routine including boundary conditions
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_INITIALIZE_SOLUTION(P(NM,NLEVEL_MAX)%X, P(NM,NLEVEL_MAX)%F, NLEVEL_MAX)
   ENDDO
ENDIF

!!! define iteration parameters
EPS = SCARC_MULTIGRID_ACCURACY
NIT = SCARC_MULTIGRID_ITERATIONS
ITE = 0
NL  = NLEVEL_MAX

!!!----------------------------------------------------------------------------------------------------
!!! Only one level: solve problem exactly (either by CG or Gaussian elimination)
!!!----------------------------------------------------------------------------------------------------
ONLY_ONE_LEVEL_IF: IF (SCARC_MULTIGRID_LEVEL==1) THEN   

      IF (NLEVEL_MIN /= NLEVEL_MAX) THEN
         WRITE(*,*) 'Error in 1-level-MG, NLEVEL_MIN=',NLEVEL_MIN,' not equal to NLEVEL_MAX=',NLEVEL_MAX
         STOP
      ENDIF

      IF (SCARC_COARSE == 'CG') THEN
         CALL SCARC_CG_BANDWISE(NSCARC_SCOPE_COARSE)
      ELSE
         CALL SCARC_GE_BANDWISE(NSCARC_SCOPE_COARSE)
      ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! More than one level: start MG-cycling
!!!----------------------------------------------------------------------------------------------------
ELSE

   !!! initialize cycle counts for MG-iteration
   ICYCLE = SCARC_CYCLE_STATE(NSCARC_CYCLE_INIT, NLEVEL_MAX)

   !!! perform initial matrix-vector product on finest level
   NL = NLEVEL_MAX
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%X, P(NM,NL)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_D, NL)

   !!! calculate norm of initial residual on finest level
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%D = P(NM,NL)%F - P(NM,NL)%D
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
   RESIN = SCARC_GLOBAL_L2NORM(NL)

   !!! print initial residual information 
   CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)

    
   !!!-------------------------------------------------------------------------------------------------
   !!! start MG-iteration
   !!!-------------------------------------------------------------------------------------------------
   MULTIGRID_LOOP: DO ITE = 1, NIT
    
      !!! reset cycling-information at beginning of each single mg-loop, start cycling in finest level
      NL = NLEVEL_MAX
      ICYCLE = SCARC_CYCLE_STATE(NSCARC_CYCLE_RESET, NL)

      CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLE_EXIT)
   
         !!!-------------------------------------------------------------------------------------------
         !!!  Presmoothing
         !!!-------------------------------------------------------------------------------------------
         PRESMOOTH_LOOP: DO WHILE (NL > NLEVEL_MIN)
   
            !!! perform presmoothing on level NL
            DO NM = NMESHES_MIN, NMESHES_MAX
               CALL SCARC_SMOOTHER(P(NM,NL)%A , P(NM,NL)%X , P(NM,NL)%D,  P(NM,NL)%F, &
                                   P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, &
                                   NSCARC_CYCLE_PRESMOOTH, NL)
            ENDDO

            !!! perform restriction to coarser level NL-1,   F:=restriction(D)
            NL=NL-1
            DO NM = NMESHES_MIN, NMESHES_MAX
               CALL SCARC_RESTRICTION(P(NM,NL)%F, P(NM,NL+1)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ)
            ENDDO

         ENDDO PRESMOOTH_LOOP
   
         !!!-------------------------------------------------------------------------------------------
         !!! Coarse grid solver (either by CG or Gaussian elimination)
         !!!-------------------------------------------------------------------------------------------
         IF (SCARC_COARSE == 'CG') THEN
            CALL SCARC_CG_BANDWISE(NSCARC_SCOPE_COARSE)
         ELSE
            CALL SCARC_GE_BANDWISE(NSCARC_SCOPE_COARSE)
         ENDIF
         TYPE_SCOPE = NTYPE
   
         !!!-------------------------------------------------------------------------------------------
         !!! Postsmoothing
         !!!-------------------------------------------------------------------------------------------
         POSTSMOOTH_LOOP: DO WHILE (NL < NLEVEL_MAX) 
       
            !!! perform prolongation to finer grid: D:=prolongation(X)
            DO NM = NMESHES_MIN, NMESHES_MAX
               CALL SCARC_PROLONGATION(P(NM,NL)%X, P(NM,NL+1)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ)
            ENDDO
   
            !!! compute new iterate and perform postsmoothing on finer level 
            NL=NL+1
            DO NM = NMESHES_MIN, NMESHES_MAX
               P(NM,NL)%X = P(NM,NL)%D + P(NM,NL)%X
               CALL SCARC_SMOOTHER(P(NM,NL)%A , P(NM,NL)%X , P(NM,NL)%D, P(NM,NL)%F, &
                                   P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ,&
                                   NSCARC_CYCLE_POSTSMOOTH, NL)
            ENDDO
   
            !!! determine how to proceed with cycling depending on chosen cycle-type (F/V/W)
            ICYCLE = SCARC_CYCLE_STATE(NSCARC_CYCLE_PROCEED, NL)
            IF (ICYCLE /= NSCARC_CYCLE_POSTSMOOTH) CYCLE CYCLE_LOOP

         ENDDO POSTSMOOTH_LOOP

      ENDDO CYCLE_LOOP

      IF (NL /= NLEVEL_MAX) THEN
         WRITE(*,*) 'ERROR in SCARC_GMG, wrong level ', NL
         STOP
      ENDIF

      !!!----------------------------------------------------------------------------------------------
      !!! Compute norm of new residual on finest level and
      !!! leave loop in case of convergence or divergence
      !!!----------------------------------------------------------------------------------------------
      DO NM = NMESHES_MIN, NMESHES_MAX
         CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A, P(NM,NL)%X, P(NM,NL)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM)
      ENDDO
      CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_D, NL)

      DO NM = NMESHES_MIN, NMESHES_MAX
         P(NM,NL)%D = P(NM,NL)%F - P(NM,NL)%D
         CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
      ENDDO
      RES = SCARC_GLOBAL_L2NORM(NL)

      ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
      IF (ISTATE /= NSCARC_LOOP_PROCEED) EXIT MULTIGRID_LOOP
 
   ENDDO MULTIGRID_LOOP

   !!!-------------------------------------------------------------------------------------------------
   !!! Determine convergence rate and print corresponding information
   !!!-------------------------------------------------------------------------------------------------
   CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

ENDIF ONLY_ONE_LEVEL_IF
 
!!!----------------------------------------------------------------------------------------------------
!!! In case of MG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_TRANSFER_DATA(P(NM,NLEVEL_MAX)%X, NM)
      CALL SCARC_UPDATE_GHOSTCELLS(NM)
   ENDDO
   CALL SCARC_EXCHANGE_BDRY(NLEVEL_MAX)
ENDIF
 
TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)=TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)+SECOND()-TNOW_MULTIGRID
TUSED_SCARC(NSCARC_TIME_COMPLETE ,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE ,:)+SECOND()-TNOW_MULTIGRID
END SUBROUTINE SCARC_MULTIGRID_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform geometric multigrid method based on global possion-matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MULTIGRID_COMPACT(NTYPE, NRHS)
INTEGER, INTENT(IN) :: NTYPE, NRHS
INTEGER   :: NM, NL, ITE, NIT, ISTATE, IERR, ICYCLE
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_MULTIGRID
CHARACTER(30):: CROUTINE
LOGICAL, SAVE :: BFIRST = .TRUE.
TYPE (COMPACT_SYSTEM_POINTER_TYPE), SAVE, POINTER, DIMENSION(:,:) :: P

TNOW_MULTIGRID = SECOND()
TYPE_SCOPE     = NTYPE

!!!----------------------------------------------------------------------------------------------------
!!! Initialize pointers to system variables corresponding to the type of the method
!!! (main solver / preconditioner)
!!! Note: This is only done at the very first call of the routine
!!!----------------------------------------------------------------------------------------------------
IF (BFIRST) THEN

   ALLOCATE(P(NMESHES_MIN:NMESHES_MAX, NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
   CALL CHKMEMERR ('SCARC_MULTIGRID_BANDWISE', 'P', IERR)

   SELECT CASE (TYPE_SCOPE)

      CASE(NSCARC_SCOPE_MAIN)
         CROUTINE = 'SCARC_GLOBAL_MG_COMPACT'
         DO NM = NMESHES_MIN, NMESHES_MAX
            DO NL=NLEVEL_MIN,NLEVEL_MAX
               P(NM,NL)%NC  => SCARC(NM)%LEVEL(NL)%N_CELLS
               P(NM,NL)%NA  => SCARC(NM)%LEVEL(NL)%N_MATRIX_ENTRIES
               P(NM,NL)%NW  => SCARC(NM)%LEVEL(NL)%N_EXTERNAL_WALL_CELLS
               P(NM,NL)%NX  => SCARC(NM)%LEVEL(NL)%IBAR
               P(NM,NL)%NY  => SCARC(NM)%LEVEL(NL)%JBAR
               P(NM,NL)%NZ  => SCARC(NM)%LEVEL(NL)%KBAR
               P(NM,NL)%A   => SCARC(NM)%LEVEL(NL)%CA
               P(NM,NL)%COL => SCARC(NM)%LEVEL(NL)%COL
               P(NM,NL)%ROW => SCARC(NM)%LEVEL(NL)%ROW
               P(NM,NL)%X   => SCARC(NM)%LEVEL(NL)%CX
               P(NM,NL)%F   => SCARC(NM)%LEVEL(NL)%CF
               P(NM,NL)%D   => SCARC(NM)%LEVEL(NL)%CD
               !P(NM,NL)%X   = 0.0_EB
            ENDDO
         ENDDO

      CASE(NSCARC_SCOPE_PRECON)
         CROUTINE = 'SCARC_PRECON_MG_COMPACT'
         DO NM = NMESHES_MIN, NMESHES_MAX
            DO NL=NLEVEL_MIN,NLEVEL_MAX
               P(NM,NL)%NC  => SCARC(NM)%LEVEL(NL)%N_CELLS
               P(NM,NL)%NA  => SCARC(NM)%LEVEL(NL)%N_MATRIX_ENTRIES
               P(NM,NL)%NW  => SCARC(NM)%LEVEL(NL)%N_EXTERNAL_WALL_CELLS
               P(NM,NL)%NX  => SCARC(NM)%LEVEL(NL)%IBAR
               P(NM,NL)%NY  => SCARC(NM)%LEVEL(NL)%JBAR
               P(NM,NL)%NZ  => SCARC(NM)%LEVEL(NL)%KBAR
               P(NM,NL)%A   => SCARC(NM)%LEVEL(NL)%CA
               P(NM,NL)%COL => SCARC(NM)%LEVEL(NL)%COL
               P(NM,NL)%ROW => SCARC(NM)%LEVEL(NL)%ROW
               P(NM,NL)%X   => SCARC(NM)%LEVEL(NL)%CX2
               P(NM,NL)%D   => SCARC(NM)%LEVEL(NL)%CD2
               P(NM,NL)%X   = 0.0_EB
               SELECT CASE(NRHS)
                  CASE(NSCARC_VECTOR_D)
                     P(NM,NL)%F => SCARC(NM)%LEVEL(NL)%CD
                  CASE(NSCARC_VECTOR_G)
                     P(NM,NL)%F => SCARC(NM)%LEVEL(NL)%CG
                  CASE(NSCARC_VECTOR_R)
                     P(NM,NL)%F => SCARC(NM)%LEVEL(NL)%CR
                  CASE(NSCARC_VECTOR_Y)
                     P(NM,NL)%F => SCARC(NM)%LEVEL(NL)%CY
               END SELECT
            ENDDO
         ENDDO
   END SELECT
            
   BFIRST = .FALSE.

ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! In case of global multigrid:
!!!     ---> initialize rhs with data from calling routine including boundary conditions
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_INITIALIZE_SOLUTION(P(NM,NLEVEL_MAX)%X, P(NM,NLEVEL_MAX)%F, NLEVEL_MAX)
   ENDDO
ENDIF

!!! define iteration parameters
EPS = SCARC_MULTIGRID_ACCURACY
NIT = SCARC_MULTIGRID_ITERATIONS
ITE = 0
NL  = NLEVEL_MAX

!!!----------------------------------------------------------------------------------------------------
!!! Only one level: solve problem exactly (either by CG or Gaussian elimination)
!!!----------------------------------------------------------------------------------------------------
ONLY_ONE_LEVEL_IF: IF (SCARC_MULTIGRID_LEVEL==1) THEN   

      IF (NLEVEL_MIN /= NLEVEL_MAX) THEN
         WRITE(*,*) 'Error in 1-level-MG, NLEVEL_MIN=',NLEVEL_MIN,' not equal to NLEVEL_MAX=',NLEVEL_MAX
         STOP
      ENDIF

      IF (SCARC_COARSE == 'CG') THEN
         CALL SCARC_CG_COMPACT(NSCARC_SCOPE_COARSE)
      ELSE
         CALL SCARC_GE_COMPACT(NSCARC_SCOPE_COARSE)
      ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! More than one level: start MG-cycling
!!!----------------------------------------------------------------------------------------------------
ELSE

   !!! initialize cycle counts for MG-iteration
   ICYCLE = SCARC_CYCLE_STATE(NSCARC_CYCLE_INIT, NLEVEL_MAX)

   !!! perform initial matrix-vector product on finest level
   NL = NLEVEL_MAX
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A , P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%X, P(NM,NL)%D, &
                               P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_D, NL)

   !!! calculate norm of initial residual on finest level
   DO NM = NMESHES_MIN, NMESHES_MAX
      P(NM,NL)%D = P(NM,NL)%F - P(NM,NL)%D
      CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
   ENDDO
   RESIN = SCARC_GLOBAL_L2NORM(NL)

   !!! print initial residual information 
   CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)

    
   !!!-------------------------------------------------------------------------------------------------
   !!! start MG-iteration
   !!!-------------------------------------------------------------------------------------------------
   MULTIGRID_LOOP: DO ITE = 1, NIT
    
      !!! reset cycling-information at beginning of each single mg-loop, start cycling in finest level
      NL = NLEVEL_MAX
      ICYCLE = SCARC_CYCLE_STATE(NSCARC_CYCLE_RESET, NL)

      CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLE_EXIT)
   
         !!!-------------------------------------------------------------------------------------------
         !!!  Presmoothing
         !!!-------------------------------------------------------------------------------------------
         PRESMOOTH_LOOP: DO WHILE (NL > NLEVEL_MIN)
   
            !!! perform presmoothing on level NL
            DO NM = NMESHES_MIN, NMESHES_MAX
               CALL SCARC_SMOOTHER(P(NM,NL)%A , P(NM,NL)%ROW, P(NM,NL)%COL, &
                                   P(NM,NL)%X , P(NM,NL)%D  , P(NM,NL)%F  , &
                                   P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ , &
                                   NSCARC_CYCLE_PRESMOOTH, NL)
            ENDDO

            !!! perform restriction to coarser level NL-1,   F:=restriction(D)
            NL=NL-1
            DO NM = NMESHES_MIN, NMESHES_MAX
               CALL SCARC_RESTRICTION(P(NM,NL)%F , P(NM,NL+1)%D, P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ)
            ENDDO

         ENDDO PRESMOOTH_LOOP
   
         !!!-------------------------------------------------------------------------------------------
         !!! Coarse grid solver (either by CG or Gaussian elimination)
         !!!-------------------------------------------------------------------------------------------
         IF (SCARC_COARSE == 'CG') THEN
            CALL SCARC_CG_COMPACT(NSCARC_SCOPE_COARSE)
         ELSE
            CALL SCARC_GE_COMPACT(NSCARC_SCOPE_COARSE)
         ENDIF
         TYPE_SCOPE = NTYPE

         !!!-------------------------------------------------------------------------------------------
         !!! Postsmoothing
         !!!-------------------------------------------------------------------------------------------
         POSTSMOOTH_LOOP: DO WHILE (NL < NLEVEL_MAX) 
       
            !!! perform prolongation to finer grid: D:=prolongation(X)
            DO NM = NMESHES_MIN, NMESHES_MAX
               CALL SCARC_PROLONGATION(P(NM,NL)%X , P(NM,NL+1)%D, P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ)
            ENDDO
   
            !!! compute new iterate and perform postsmoothing on finer level 
            NL=NL+1
            DO NM = NMESHES_MIN, NMESHES_MAX
               P(NM,NL)%X = P(NM,NL)%D + P(NM,NL)%X
               CALL SCARC_SMOOTHER(P(NM,NL)%A , P(NM,NL)%ROW, P(NM,NL)%COL, &
                                   P(NM,NL)%X , P(NM,NL)%D  , P(NM,NL)%F  , &
                                   P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ ,&
                                   NSCARC_CYCLE_POSTSMOOTH, NL)
            ENDDO
   
            !!! determine how to proceed with cycling depending on chosen cycle-type (F/V/W)
            ICYCLE = SCARC_CYCLE_STATE(NSCARC_CYCLE_PROCEED, NL)
            SELECT CASE(ICYCLE)
               CASE (NSCARC_CYCLE_POSTSMOOTH)
                  CYCLE POSTSMOOTH_LOOP
               CASE DEFAULT
                  CYCLE CYCLE_LOOP
            END SELECT
   
         ENDDO POSTSMOOTH_LOOP

      ENDDO CYCLE_LOOP

      IF (NL /= NLEVEL_MAX) THEN
         WRITE(*,*) 'ERROR in SCARC_GMG, wrong level ', NL
         STOP
      ENDIF


      !!!----------------------------------------------------------------------------------------------
      !!! Compute norm of new residual on finest level and
      !!! leave loop in case of convergence or divergence
      !!!----------------------------------------------------------------------------------------------
      DO NM = NMESHES_MIN, NMESHES_MAX
         CALL SCARC_LOCAL_MATVEC (P(NM,NL)%A , P(NM,NL)%ROW, P(NM,NL)%COL, P(NM,NL)%X, P(NM,NL)%D, &
                                  P(NM,NL)%NX, P(NM,NL)%NY , P(NM,NL)%NZ, NM)
      ENDDO
      CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_D, NL)

      DO NM = NMESHES_MIN, NMESHES_MAX
         P(NM,NL)%D = P(NM,NL)%F - P(NM,NL)%D
         CALL SCARC_LOCAL_SCALPROD (P(NM,NL)%D, P(NM,NL)%D, P(NM,NL)%NX, P(NM,NL)%NY, P(NM,NL)%NZ, NM) 
      ENDDO
      RES = SCARC_GLOBAL_L2NORM(NL)

      ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
      IF (ISTATE /= NSCARC_LOOP_PROCEED) EXIT MULTIGRID_LOOP
 
   ENDDO MULTIGRID_LOOP

   !!!-------------------------------------------------------------------------------------------------
   !!! Determine convergence rate and print corresponding information
   !!!-------------------------------------------------------------------------------------------------
   CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

ENDIF ONLY_ONE_LEVEL_IF
 
!!!----------------------------------------------------------------------------------------------------
!!! In case of MG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_TRANSFER_DATA(P(NM,NLEVEL_MAX)%X, NM)
      CALL SCARC_UPDATE_GHOSTCELLS(NM)
   ENDDO
   CALL SCARC_EXCHANGE_BDRY(NLEVEL_MAX)
ENDIF
 
TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)=TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)+SECOND()-TNOW_MULTIGRID
TUSED_SCARC(NSCARC_TIME_COMPLETE ,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE ,:)+SECOND()-TNOW_MULTIGRID
END SUBROUTINE SCARC_MULTIGRID_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Control multigrid cycling (F/V/W)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION SCARC_CYCLE_STATE(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NL0, ICYCLE

SELECT CASE(NTYPE)

   !!!-------------------------------------------------------------------------------------------------
   !!! initialize cycle counts at beginning of multigrid method
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_INIT)

      DO NM = NMESHES_MIN, NMESHES_MAX
         SCARC(NM)%CYCLE_COUNT(2,NL)=1
         DO NL0=NLEVEL_MIN+1,NLEVEL_MAX-1
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
         DO NL0 = NLEVEL_MIN,NLEVEL_MAX
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
            IF (NL == NLEVEL_MAX) THEN
               ICYCLE = NSCARC_CYCLE_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLE_POSTSMOOTH
            ENDIF
         ELSE
            IF (NL == NLEVEL_MAX) THEN
               ICYCLE = NSCARC_CYCLE_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLE_NEXT
            ENDIF
         ENDIF
      ENDDO

END SELECT

SCARC_CYCLE_STATE = ICYCLE
RETURN

END FUNCTION SCARC_CYCLE_STATE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform smoothing - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SMOOTHER_BANDWISE(A, X, D, F, NX, NY, NZ, NTYPE, NL)
REAL(EB), POINTER, DIMENSION(:, :),    INTENT(INOUT) :: A
REAL(EB), POINTER, DIMENSION(:, :, :), INTENT(INOUT) :: X, D, F
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER, POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER :: NM, ITE, NIT, ISTATE
REAL(EB):: RES, RESIN, EPS, OMEGA
REAL(EB):: TNOW_SMOOTH
LOGICAL :: BMATVEC, BL2NORM
CHARACTER(30) :: CROUTINE

TNOW_SMOOTH = SECOND()

!!!----------------------------------------------------------------------------------------------------
!!! Initialization
!!!----------------------------------------------------------------------------------------------------
NIT      = SCARC_SMOOTH_ITERATIONS
EPS      = SCARC_SMOOTH_ACCURACY
OMEGA    = SCARC_SMOOTH_OMEGA
ITE      = 0
RESIN    = 1.0_EB
CROUTINE = 'SCARC_SMOOTHER_BANDWISE'
BL2NORM  = .TRUE.
BMATVEC  = .TRUE.
IF (NTYPE == NSCARC_CYCLE_PRESMOOTH .AND. NL == NLEVEL_MAX) BMATVEC = .FALSE.

!!! Calculate initial defect on level NL (only if BMATVEC = .TRUE.)
!!! Because initial vector is set to zero, this defect corresponds to F
IF (BMATVEC) THEN

   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_MATVEC (A, X, D, NX, NY, NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_D, NL)

   DO NM = NMESHES_MIN, NMESHES_MAX
      D = F - D
      CALL SCARC_LOCAL_SCALPROD (D, D, NX, NY, NZ, NM)
   ENDDO
   IF (BL2NORM) THEN
      RESIN = SCARC_GLOBAL_L2NORM(NL)
      CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)
   ENDIF

ENDIF


!!!----------------------------------------------------------------------------------------------------
!!! Smoothing loop
!!!----------------------------------------------------------------------------------------------------
SMOOTH_LOOP: DO ITE=1, NIT
 
   !!! Perform preconditioning
   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_PRECONDITIONER(A, D, D, NX, NY, NZ, NM)
   ENDDO

   !!! get new iterate and compute matrix-vector product
   DO NM = NMESHES_MIN, NMESHES_MAX
      X = X + OMEGA * D 
      CALL SCARC_LOCAL_MATVEC (A, X, D, NX, NY, NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_D, NL)

   !!! Compute norm of new residual if NIT isn't predefined
   DO NM = NMESHES_MIN, NMESHES_MAX
      D = F - D
      CALL SCARC_LOCAL_SCALPROD (D, D, NX, NY, NZ, NM)
   ENDDO
   IF (BL2NORM) THEN
      RES = SCARC_GLOBAL_L2NORM(NL)
      ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
      IF (ISTATE /= NSCARC_LOOP_PROCEED) EXIT SMOOTH_LOOP
   ENDIF

ENDDO SMOOTH_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!!----------------------------------------------------------------------------------------------------
IF (BL2NORM) CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

TUSED_SCARC(NSCARC_TIME_SMOOTH  ,:)=TUSED_SCARC(NSCARC_TIME_SMOOTH  ,:)+SECOND()-TNOW_SMOOTH
TUSED_SCARC(NSCARC_TIME_COMPLETE,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE,:)+SECOND()-TNOW_SMOOTH
END SUBROUTINE SCARC_SMOOTHER_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform smoothing - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SMOOTHER_COMPACT(A, ROW, COL, X, D, F, NX, NY, NZ, NTYPE, NL)
INTEGER,  POINTER, DIMENSION(:), INTENT(IN)    :: ROW, COL
REAL(EB), POINTER, DIMENSION(:), INTENT(IN)    :: A
REAL(EB), POINTER, DIMENSION(:), INTENT(INOUT) :: X, D, F
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER, POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER :: NM, ITE, NIT, ISTATE
REAL(EB):: RES, RESIN, EPS, OMEGA
REAL(EB):: TNOW_SMOOTH
LOGICAL :: BMATVEC, BL2NORM
CHARACTER(30) :: CROUTINE

TNOW_SMOOTH = SECOND()

!!!----------------------------------------------------------------------------------------------------
!!! Initialization
!!!----------------------------------------------------------------------------------------------------
NIT      = SCARC_SMOOTH_ITERATIONS
EPS      = SCARC_SMOOTH_ACCURACY
OMEGA    = SCARC_SMOOTH_OMEGA
ITE      = 0
RESIN    = 1.0_EB
CROUTINE = 'SCARC_SMOOTHER_COMPACT'
BL2NORM  = .TRUE.
BMATVEC  = .TRUE.
IF (NTYPE == NSCARC_CYCLE_PRESMOOTH .AND. NL == NLEVEL_MAX) BMATVEC = .FALSE.

!!! Calculate initial defect on level NL (only if BMATVEC = .TRUE.)
!!! Because initial vector is set to zero, this defect corresponds to F
IF (BMATVEC) THEN

   DO NM = NMESHES_MIN, NMESHES_MAX
      CALL SCARC_LOCAL_MATVEC (A, ROW, COL, X, D, NX, NY, NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_D, NL)

   DO NM = NMESHES_MIN, NMESHES_MAX
      D = F - D
      CALL SCARC_LOCAL_SCALPROD (D, D, NX, NY, NZ, NM)
   ENDDO
   IF (BL2NORM) THEN
      RESIN = SCARC_GLOBAL_L2NORM(NL)
      CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)
   ENDIF

ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! Smoothing loop
!!!----------------------------------------------------------------------------------------------------
SMOOTH_LOOP: DO ITE=1, NIT
 
   !!! Perform preconditioning
   DO NM = NMESHES_MIN, NMESHES_MAX
       CALL SCARC_PRECONDITIONER(A, ROW, COL, D, D, NX, NY, NZ, NM)
   ENDDO

   !!! get new iterate and perform matrix-vector product
   DO NM = NMESHES_MIN, NMESHES_MAX
      X = X + OMEGA * D 
      CALL SCARC_LOCAL_MATVEC (A, ROW, COL, X, D, NX, NY, NZ, NM)
   ENDDO
   CALL SCARC_GLOBAL_MATVEC(NSCARC_VECTOR_X, NSCARC_VECTOR_D, NL)

   !!! Compute norm of new residual if NIT isn't predefined
   DO NM = NMESHES_MIN, NMESHES_MAX
      D = F - D
      CALL SCARC_LOCAL_SCALPROD (D, D, NX, NY, NZ, NM)
   ENDDO
   IF (BL2NORM) THEN
      RES = SCARC_GLOBAL_L2NORM(NL)
      ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
      IF (ISTATE /= NSCARC_LOOP_PROCEED) EXIT SMOOTH_LOOP
   ENDIF

ENDDO SMOOTH_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!!----------------------------------------------------------------------------------------------------
IF (BL2NORM) CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

TUSED_SCARC(NSCARC_TIME_SMOOTH  ,:)=TUSED_SCARC(NSCARC_TIME_SMOOTH  ,:)+SECOND()-TNOW_SMOOTH
TUSED_SCARC(NSCARC_TIME_COMPLETE,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE,:)+SECOND()-TNOW_SMOOTH
END SUBROUTINE SCARC_SMOOTHER_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Gaussian elimination for coarse grid solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_GE_BANDWISE(NTYPE)
INTEGER, INTENT(IN) :: NTYPE
IF (NTYPE == NSCARC_SCOPE_COARSE) THEN
   WRITE(*,*) 'global GE not yet implemented'
   STOP
ENDIF
END SUBROUTINE SCARC_GE_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Gaussian elimination for coarse grid solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_GE_COMPACT(NTYPE)
INTEGER, INTENT(IN) :: NTYPE
IF (NTYPE == NSCARC_SCOPE_COARSE) THEN
   WRITE(*,*) 'global GE not yet implemented'
   STOP
ENDIF
END SUBROUTINE SCARC_GE_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Select preconditioning method - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECONDITIONER_BANDWISE(A, R, G, NX, NY, NZ, NM)
INTEGER  , INTENT(IN) :: NM
INTEGER  , POINTER, INTENT(IN) :: NX, NY, NZ
REAL (EB), POINTER, DIMENSION (:, :),    INTENT(IN)    :: A
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(IN)    :: R
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(INOUT) :: G
REAL(EB):: TNOW_PRECON

TNOW_PRECON = SECOND()

SELECT CASE(TYPE_PRECON)
   CASE(NSCARC_PRECON_JACOBI)
      CALL SCARC_JACOBI_BANDWISE(A, G, NX, NY, NZ)
   CASE(NSCARC_PRECON_SSOR)
      CALL SCARC_SSOR_BANDWISE(A, G, NX, NY, NZ)
   CASE(NSCARC_PRECON_GSTRIX)
      CALL SCARC_GSTRIX_BANDWISE(G, NX, NY, NZ, NM)
   CASE(NSCARC_PRECON_FFT)
      CALL SCARC_FFT_BANDWISE(R, G, NX, NY, NZ, NM)
END SELECT

TUSED_SCARC(NSCARC_TIME_PRECON  ,NM)=TUSED_SCARC(NSCARC_TIME_PRECON  ,NM)+SECOND()-TNOW_PRECON
TUSED_SCARC(NSCARC_TIME_COMPLETE,NM)=TUSED_SCARC(NSCARC_TIME_COMPLETE,NM)+SECOND()-TNOW_PRECON
END SUBROUTINE SCARC_PRECONDITIONER_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Select preconditioning method - compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECONDITIONER_COMPACT(A, ROW, COL, R, G, NX, NY, NZ, NM)
INTEGER  , INTENT(IN) :: NM
INTEGER  , POINTER, INTENT(IN) :: NX, NY, NZ
REAL (EB), POINTER, DIMENSION (:), INTENT(IN)    :: A
INTEGER  , POINTER, DIMENSION (:), INTENT(IN)    :: ROW, COL
REAL (EB), POINTER, DIMENSION (:), INTENT(IN)    :: R
REAL (EB), POINTER, DIMENSION (:), INTENT(INOUT) :: G
REAL(EB):: TNOW_PRECON

TNOW_PRECON = SECOND()

SELECT CASE(TYPE_PRECON)
   CASE(NSCARC_PRECON_JACOBI)
      CALL SCARC_JACOBI_COMPACT(A, ROW, G, NX, NY, NZ)
   CASE(NSCARC_PRECON_SSOR)
      CALL SCARC_SSOR_COMPACT(A, ROW, COL, G, NX, NY, NZ)
   CASE(NSCARC_PRECON_GSTRIX)
      CALL SCARC_GSTRIX_COMPACT(G, NX, NY, NZ, NM)
   CASE(NSCARC_PRECON_FFT)
      CALL SCARC_FFT_COMPACT(R, G, NX, NY, NZ, NM)
END SELECT

TUSED_SCARC(NSCARC_TIME_PRECON  ,NM)=TUSED_SCARC(NSCARC_TIME_PRECON  ,NM)+SECOND()-TNOW_PRECON
TUSED_SCARC(NSCARC_TIME_COMPLETE,NM)=TUSED_SCARC(NSCARC_TIME_COMPLETE,NM)+SECOND()-TNOW_PRECON
END SUBROUTINE SCARC_PRECONDITIONER_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Jacobi preconditioner/smoother - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_JACOBI_BANDWISE(A, G, NX, NY, NZ)
REAL (EB), POINTER, DIMENSION (:, :),    INTENT(IN)    :: A
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(INOUT) :: G
INTEGER  , POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER :: I, J, K, IC
DO K = 1, NZ
   DO J = 1, NY
      DO I = 1, NX
         IC = (K-1) * NX * NY + (J-1) * NX + I
         G (I, J, K) = G (I, J, K) / A(IC, ID)
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE SCARC_JACOBI_BANDWISE
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Jacobi preconditioner/smoother - compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_JACOBI_COMPACT(A, ROW, G, NX, NY, NZ)
REAL (EB), POINTER, DIMENSION (:), INTENT(IN)    :: A
INTEGER  , POINTER, DIMENSION (:), INTENT(IN)    :: ROW
REAL (EB), POINTER, DIMENSION (:), INTENT(INOUT) :: G
INTEGER  , POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER :: I, J, K, IC
DO K = 1, NZ
   DO J = 1, NY
      DO I = 1, NX
         IC = (K-1) * NX * NY + (J-1) * NX + I
         G (IC) = G (IC) / A(ROW(IC))
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE SCARC_JACOBI_COMPACT
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SSOR preconditioner/smoother - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SSOR_BANDWISE(A, G, NX, NY, NZ)
REAL (EB), POINTER, DIMENSION (:, :),    INTENT(IN)    :: A
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(INOUT) :: G
INTEGER  , POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER :: I, J, K, IC
REAL (EB) :: AUX, OMEGA=1.5_EB

DIMENSION_IF: IF (TWO_D) THEN
   DO K = 1, NZ
      DO I = 1, NX
         IC = (K-1) * NX + I
         AUX =    A(IC,ILZ) * G (I  , 1, K-1)  &
                + A(IC,ILX) * G (I-1, 1, K  )
         G (I, 1, K) = (G(I, 1, K) - AUX*OMEGA) / A(IC,ID)
      ENDDO
   ENDDO
   DO K = NZ, 1, - 1
      DO I = NX, 1, - 1
         IC = (K-1) * NX + I
         AUX =    A(IC,IUZ) * G (I  , 1, K+1) &
                + A(IC,IUX) * G (I+1, 1, K  )
         G (I, 1, K) = G (I, 1, K) - AUX * OMEGA / A(IC,ID)
      ENDDO
   ENDDO
ELSE
   DO K = 1, NZ
      DO J = 1, NY
         DO I = 1, NX
            IC = (K-1) * NX * NY + (J-1) * NX + I
            AUX =    A(IC,ILZ) * G (I  , J  , K-1) &
                   + A(IC,ILY) * G (I  , J-1, K  ) &
                   + A(IC,ILX) * G (I-1, J  , K  )
            G (I, J, K) = (G(I, J, K) - AUX * OMEGA) / A(IC,ID)
         ENDDO
      ENDDO
   ENDDO
   DO K = NZ, 1, - 1
      DO J = NY, 1, - 1
         DO I = NX, 1, - 1
            IC = (K-1) * NX * NY + (J-1) * NX + I
            AUX =    A(IC,IUZ) * G (I  , J  , K+1) &
                   + A(IC,IUY) * G (I  , J+1, K  ) &
                   + A(IC,IUZ) * G (I+1, J  , K  )
            G (I, J, K) = G (I, J, K) - AUX * OMEGA / A(IC,ID)
         ENDDO
      ENDDO
   ENDDO
ENDIF DIMENSION_IF
END SUBROUTINE SCARC_SSOR_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SSOR preconditioner/smoother - compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SSOR_COMPACT(A, ROW, COL, G, NX, NY, NZ)
REAL (EB), POINTER, DIMENSION (:), INTENT(IN)    :: A
INTEGER  , POINTER, DIMENSION (:), INTENT(IN)    :: ROW, COL
REAL (EB), POINTER, DIMENSION (:), INTENT(INOUT) :: G
INTEGER  , POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER :: IC, ICOL, NC
REAL (EB) :: AUX, OMEGA=1.5_EB

NC = NX * NY * NZ

!!! use only matrix superdiagonals
FORWARD_CELL_LOOP: DO IC = 1, NC

   AUX = 0.0_EB
   LOWER_DIAG_LOOP: DO ICOL = ROW(IC)+1, ROW(IC+1)-1
      IF (COL(ICOL) >= IC) EXIT LOWER_DIAG_LOOP
      AUX = AUX + A(ICOL) * G(COL(ICOL))
   ENDDO LOWER_DIAG_LOOP
 
   G(IC) = (G(IC) - AUX * OMEGA) / A(ROW(IC))

ENDDO FORWARD_CELL_LOOP

!!! use only matrix subdiagonals
BACKWARD_CELL_LOOP: DO IC = NC-1, 1, -1
   
   AUX = 0.0_EB
   UPPER_DIAG_LOOP: DO ICOL = ROW(IC)+1, ROW(IC+1)-1
      IF (COL(ICOL) <= IC) CYCLE
      AUX = AUX + A(ICOL) * G(COL(ICOL))
   ENDDO UPPER_DIAG_LOOP
   
   G(IC) = G(IC) - AUX * OMEGA / A(ROW(IC))

ENDDO BACKWARD_CELL_LOOP

END SUBROUTINE SCARC_SSOR_COMPACT
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FFT preconditioner/smoother - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_FFT_BANDWISE (R, G, NX, NY, NZ, NM)
USE POIS, ONLY: H3CZSS, H2CZSS
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(IN)    :: R
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(INOUT) :: G
INTEGER  , POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER  , INTENT(IN) :: NM

M  => MESHES(NM)
SL => SCARC(NM)%LEVEL(NLEVEL_MAX)

SL%FFT(1:NX, 1:NY, 1:NZ) = R(1:NX, 1:NY, 1:NZ)
DIMENSION_IF: IF (TWO_D) THEN
   CALL H2CZSS (M%BXS, M%BXF, M%BZS, M%BZF,&
                NX+1, SL%FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
ELSE
   CALL H3CZSS (M%BXS, M%BXF, M%BYS, M%BYF, M%BZS, M%BZF, &
                NX+1, NY+1, SL%FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
ENDIF DIMENSION_IF
G(1:NX, 1:NY, 1:NZ) = SL%FFT(1:NX, 1:NY, 1:NZ)

END SUBROUTINE SCARC_FFT_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FFT preconditioner/smoother - compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_FFT_COMPACT (R, G, NX, NY, NZ, NM)
USE POIS, ONLY: H3CZSS, H2CZSS
REAL (EB), POINTER, DIMENSION (:), INTENT(IN)    :: R
REAL (EB), POINTER, DIMENSION (:), INTENT(INOUT) :: G
INTEGER  , POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER  , INTENT(IN) :: NM
INTEGER :: I, J, K, IC

M  => MESHES(NM)
SL => SCARC(NM)%LEVEL(NLEVEL_MAX)

DO K = 1, NZ
   DO J = 1, NY
      DO I = 1, NX
         IC = (K-1) * NX * NY + (J-1) * NX + I
WRITE(*,*) 'RICHTIG ??????????????'
         SL%FFT(I, J, K) = R(IC)
      ENDDO
   ENDDO
ENDDO
DIMENSION_IF: IF (TWO_D) THEN
   CALL H2CZSS (M%BXS, M%BXF, M%BZS, M%BZF,&
                NX+1, SL%FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
ELSE
   CALL H3CZSS (M%BXS, M%BXF, M%BYS, M%BYF, M%BZS, M%BZF, &
                NX+1, NY+1, SL%FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
ENDIF DIMENSION_IF
DO K = 1, NZ
   DO J = 1, NY
      DO I = 1, NX
         IC = (K-1) * NX * NY + (J-1) * NX + I
WRITE(*,*) 'RICHTIG ??????????'
WRITE(SCARC_LU,*) 'RICHTIG ??????????'
         G(IC) = SL%FFT(I, J, K)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_FFT_COMPACT

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize GSTRIX preconditioner/smoother
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECONDITIONER_INIT_GSTRIX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, NC, IERR

IERR=0

SL => SCARC(NM)%LEVEL(NL)
NC =  SL%N_CELLS

!!!----------------------------------------------------------------------------------------------------
!!! 2D-version
!!!----------------------------------------------------------------------------------------------------
DIMENSION_IF: IF (TWO_D) THEN

   ALLOCATE (SL%MDX(1:NC), STAT=IERR)
   CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SL%MDX', IERR)
   
   ALLOCATE (SL%UDX(1:NC), STAT=IERR)
   CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SL%UDX', IERR)
   
   ALLOCATE (SL%LDX(1:NC), STAT=IERR)
   CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SL%LDX', IERR)
   
   ALLOCATE (SL%MDZ(1:NC), STAT=IERR)
   CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SL%MDZ', IERR)
   
   ALLOCATE (SL%DWORK(1:NC), STAT=IERR)
   CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SL%DWORK', IERR)
   
   !!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
   IF (TYPE_STORAGE == NSCARC_STORAGE_BANDWISE) THEN
      DO IC = 1, NC
         SL%MDX(IC) = SL%BA(IC,ID)         ! main  diagonal in x-direction
         SL%MDZ(IC) = SL%BA(IC,ILZ)        ! main  diagonal in z-direction
         SL%LDX(IC) = SL%BA(IC,ILX)        ! lower diagonal in x-direction
         SL%UDX(IC) = SL%BA(IC,IUX)        ! upper diagonal in x-direction
      ENDDO
   ELSE
      WRITE(*,*) 'HIER NOCH ANPASSEN !!'
      WRITE(SCARC_LU,*) 'HIER NOCH ANPASSEN !!'
      DO IC = 1, NC
!          SL%MDX(IC) = SL%CA(IC,ID)         ! main  diagonal in x-direction
!         SL%MDZ(IC) = SL%CA(IC,ILZ)        ! main  diagonal in z-direction
!         SL%LDX(IC) = SL%CA(IC,ILX)        ! lower diagonal in x-direction
!         SL%UDX(IC) = SL%CA(IC,IUX)        ! upper diagonal in x-direction
      ENDDO
   ENDIF
    
   !!! perform LU-factorization of matrix AG according to bandwise storage technique
   DO IC = 2, NC
      SL%LDX (IC) = SL%LDX(IC) / SL%MDX(IC-1)
      SL%MDX (IC) = SL%MDX(IC) - SL%LDX(IC) * SL%UDX(IC-1)
   ENDDO
    
   !!! replace diagonal values diag(i) by 1/diag(i) and multiply UDX with 1/diag(i)
   DO IC = 1, NC
      SL%MDX (IC) = 1.0_EB / SL%MDX(IC)
      SL%UDX (IC) = SL%MDX(IC) * SL%UDX(IC)
   ENDDO
   
!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE

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
   
   ALLOCATE (SL%DWORK(1:NC), STAT=IERR)
   CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DWORK', IERR)
   
   !!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
   IF (TYPE_STORAGE == NSCARC_STORAGE_BANDWISE) THEN
      DO IC = 1, NC
         SL%MDX(IC) = SL%BA(IC,ID)        ! main  diagonal in x-direction
         SL%MDZ(IC) = SL%BA(IC,ILZ)       ! main  diagonal in z-direction
         SL%MDY(IC) = SL%BA(IC,ILY)       ! main  diagonal in y-direction
         SL%LDX(IC) = SL%BA(IC,ILX)       ! lower diagonal in x-direction
         SL%UDX(IC) = SL%BA(IC,IUX)       ! upper diagonal in x-direction
      ENDDO
   ELSE
      WRITE(*,*) 'HIER NOCH ANPASSEN !!'
      WRITE(SCARC_LU,*) 'HIER NOCH ANPASSEN !!'
      DO IC = 1, NC
!         SL%MDX(IC) = SL%CA(IC,ID)        ! main  diagonal in x-direction
!         SL%MDZ(IC) = SL%CA(IC,ILZ)       ! main  diagonal in z-direction
!         SL%MDY(IC) = SL%CA(IC,ILY)       ! main  diagonal in y-direction
!         SL%LDX(IC) = SL%CA(IC,ILX)       ! lower diagonal in x-direction
!         SL%UDX(IC) = SL%CA(IC,IUX)       ! upper diagonal in x-direction
      ENDDO
   ENDIF
    
   !!! perform LU-factorization of matrix AG according to bandwise storage technique
   DO IC = 2, NC
      SL%LDX (IC) = SL%LDX(IC) / SL%MDX(IC-1)
      SL%MDX (IC) = SL%MDX(IC) - SL%LDX(IC) * SL%UDX(IC-1)
   ENDDO
    
   !!! replace diagonal values diag(i) by 1/diag(i) and multiply UDX with 1/diag(i)
   DO IC = 1, NC
      SL%MDX (IC) = 1.0_EB / SL%MDX(IC)
      SL%UDX (IC) = SL%MDX(IC) * SL%UDX(IC)
   ENDDO
   
ENDIF DIMENSION_IF

END SUBROUTINE SCARC_PRECONDITIONER_INIT_GSTRIX
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GSTRIX preconditioner/smoother
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_GSTRIX_BANDWISE (G, NX, NY, NZ, NM)
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(INOUT) :: G
INTEGER, POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER, INTENT(IN) :: NM
INTEGER :: I, J, K, IC

WRITE(SCARC_LU,*) 'ACHTUNG, NOCH ANPASSEN!'
WRITE(*,*) 'ACHTUNG, NOCH ANPASSEN!'
SL => SCARC(NM)%LEVEL(NLEVEL_MAX)

!!!----------------------------------------------------------------------------------------------------
!!! 2D-version
!!!----------------------------------------------------------------------------------------------------
DIMENSION_IF: IF (TWO_D) THEN

   ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
   DO I=2,NX
      G(I,1,1) = G(I,1,1)-SL%LDX(I-1)*G(I-1,1,1)
   ENDDO
   DO I=1,NX
      G(I,1,1) = G(I,1,1)*SL%MDX(I)
   ENDDO
   DO I=NX-1,1,-1
      G(I,1,1) = G(I,1,1)-SL%UDX(I)*G(I+1,1,1)
   ENDDO
     
   
   ! backward elimination of following unknowns (here the subdiagonal in z-direction must be taken into account)
   DO K=2,NZ
   
      IC = (K-1)*NX + I
   
      ! eliminate subdiagonal in z-direction into right hand side (related solution component already known)
      DO I=1,NX
         G(I,1,K) = G(I,1,K) - SL%MDZ(IC)*G(I,1,K-1) 
      ENDDO
   
      ! perform elimination of matrix lines corresponding to K
      DO I=2,NX
         G(I,1,K) = G(I,1,K)-SL%LDX(IC-1)*G(I-1,1,K)
      ENDDO
      DO I=1,NX
         G(I,1,K) = G(I,1,K)*SL%MDX(IC)
      ENDDO
      DO I=NX-1,1,-1
         G(I,1,K) = G(I,1,K)-SL%UDX(IC)*G(I+1,1,K)
      ENDDO
     
   ENDDO
   
!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE

   !!! NOT WORKING YET, has still to be adopted to 3D !!!!!
   DO J=1,NY
   
      ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
      DO I=2,NX
         G(I,J,1) = G(I,J,1)-SL%LDX(I-1)*G(I-1,J,1)
      ENDDO
      DO I=1,NX
         G(I,J,1) = G(I,J,1)*SL%MDX(I)
      ENDDO
      DO I=NX-1,1,-1
         G(I,J,1) = G(I,J,1)-SL%UDX(I)*G(I+1,J,1)
      ENDDO
        
      ! backward elimination of following unknowns 
      ! (here the subdiagonals in y- and z-direction must be taken into account)
      DO K=2,NZ
      
         IC = (K-1)*NX + I
      
         ! bring subdiagonal in z-direction to right hand side (corresponding solution component already known)
         DO I=1,NX
            G(I,J,K) = G(I,J,K) - SL%MDZ(IC)*G(I,J,K-1) 
         ENDDO
      
         ! perform elimination of matrix lines corresponding to K
         DO I=2,NX
            G(I,J,K) = G(I,J,K)-SL%LDX(IC-1)*G(I-1,J,K)
         ENDDO
         DO I=1,NX
            G(I,J,K) = G(I,J,K)*SL%MDX(IC)
         ENDDO
         DO I=NX-1,1,-1
            G(I,J,K) = G(I,J,K)-SL%UDX(IC)*G(I+1,J,K)
         ENDDO
        
      ENDDO
   
   ENDDO
   
ENDIF DIMENSION_IF
 
END SUBROUTINE SCARC_GSTRIX_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GSTRIX preconditioner/smoother
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_GSTRIX_COMPACT (G, NX, NY, NZ, NM)
REAL(EB), POINTER, DIMENSION (:), INTENT(INOUT) :: G
INTEGER , POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER , INTENT(IN) :: NM
INTEGER :: I, J, K, IC, JC

WRITE(SCARC_LU,*) 'ACHTUNG, NOCH ANPASSEN!'
WRITE(*,*) 'ACHTUNG, NOCH ANPASSEN!'
SL => SCARC(NM)%LEVEL(NLEVEL_MAX)

!!!----------------------------------------------------------------------------------------------------
!!! 2D-version
!!!----------------------------------------------------------------------------------------------------
DIMENSION_IF: IF (TWO_D) THEN

   ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
   DO IC=2,NX
      G(IC) = G(IC)-SL%LDX(I-1)*G(IC-1)
   ENDDO
   DO I=1,NX
      G(IC) = G(IC)*SL%MDX(I)
   ENDDO
   DO I=NX-1,1,-1
      G(IC) = G(IC)-SL%UDX(I)*G(IC+1)
   ENDDO
     
   
   ! backward elimination of following unknowns (here the subdiagonal in z-direction must be taken into account)
   DO K=2,NZ
   
      IC = (K-1)*NX + I
   
      ! eliminate subdiagonal in z-direction into right hand side (related solution component already known)
      DO I=1,NX
         G(IC) = G(IC) - SL%MDZ(IC)*G(IC-1) 
      ENDDO
   
      ! perform elimination of matrix lines corresponding to K
      DO I=2,NX
         G(IC) = G(IC)-SL%LDX(IC-1)*G(IC-1)
      ENDDO
      DO I=1,NX
         G(IC) = G(IC)*SL%MDX(IC)
      ENDDO
      DO I=NX-1,1,-1
         G(IC) = G(IC)-SL%UDX(IC)*G(IC+1)
      ENDDO
     
   ENDDO
   
!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE

   !!! NOT WORKING YET, has still to be adopted to 3D !!!!!
   DO J=1,NY
   
      JC = J*NX

      ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
      DO I=2,NX
         G(JC+I) = G(JC+I)-SL%LDX(I-1)*G(JC+I-1)
      ENDDO
      DO I=1,NX
         G(JC+I) = G(JC+I)*SL%MDX(I)
      ENDDO
      DO I=NX-1,1,-1
         G(JC+I) = G(JC+I)-SL%UDX(I)*G(JC+I+1)
      ENDDO
        
      ! backward elimination of following unknowns 
      ! (here the subdiagonals in y- and z-direction must be taken into account)
      DO K=2,NZ
      
         IC = (K-1)*NX*NY + (J-1)*NX + I
      
         ! bring subdiagonal in z-direction to right hand side (corresponding solution component already known)
         DO I=1,NX
            G(IC) = G(IC) - SL%MDZ(IC)*G(IC-1) 
         ENDDO
      
         ! perform elimination of matrix lines corresponding to K
         DO I=2,NX
            G(IC) = G(IC)-SL%LDX(IC-1)*G(IC-1)
         ENDDO
         DO I=1,NX
            G(IC) = G(IC)*SL%MDX(IC)
         ENDDO
         DO I=NX-1,1,-1
            G(IC) = G(IC)-SL%UDX(IC)*G(IC+1)
         ENDDO
        
      ENDDO
   
   ENDDO
   
ENDIF DIMENSION_IF
 
END SUBROUTINE SCARC_GSTRIX_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set initial solution corresponding to H and HS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_SOLUTION_BANDWISE(X, F, NL)
INTEGER, INTENT(IN) :: NL
REAL(EB), POINTER, DIMENSION(:,:,:) :: X, F
INTEGER :: NM, IW, IOR0, I, J, K

!!!----------------------------------------------------------------------------------------------------
!!!  2D-version
!!!----------------------------------------------------------------------------------------------------
DIMENSION_BANDWISE_IF: IF (TWO_D) THEN

   MESH_LOOP_BANDWISE2D: DO NM = NMESHES_MIN, NMESHES_MAX
   
      M  => MESHES(NM)
      SL => SCARC(NM)%LEVEL(NL)
   
      ! get right hand side and initial vector
      F (1:SL%IBAR, 1, 1:SL%KBAR) = M%PRHS (1:M%IBAR, 1, 1:M%KBAR)
      IF (PREDICTOR) THEN
         X (1:SL%IBAR, 1, 1:SL%KBAR) = M%H (1:M%IBAR, 1, 1:M%KBAR)
      ELSE
         X (1:SL%IBAR, 1, 1:SL%KBAR) = M%HS(1:M%IBAR, 1, 1:M%KBAR)
      ENDIF

      ! set boundary conditions at exterior boundaries (corresponding to pois.f90)
      WALL_CELL_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
   
         I = M%IJKW(6,IW)
         J = M%IJKW(7,IW)
         K = M%IJKW(8,IW)

         IF (J /= 1) THEN
            WRITE(*,*) 'Wrong index for J =',J,' in SCARC_INITIALIZE_SOLUTION !!!'
            STOP
         ENDIF
      
         IOR0 = M%IJKW(4,IW)
      
         SELECT CASE(IOR0)
            CASE(1)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DXI2 * M%BXS(1,K)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) + SL%DXI * M%BXS(1,K)                   ! Neumann
               ENDIF
            CASE(-1)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(1,K)          ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) - SL%DXI *M%BXF(1,K)                    ! Neumann
               ENDIF
            CASE(3)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,1)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) + SL%DZI * M%BZS(I,1)                   ! Neumann
               ENDIF
            CASE(-3)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,1)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) - SL%DZI  * M%BZF(I,1)                  ! Neumann
               ENDIF
         END SELECT
      
      ENDDO WALL_CELL_LOOP2D

   ENDDO MESH_LOOP_BANDWISE2D


!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE

   MESH_LOOP_BANDWISE3D: DO NM = NMESHES_MIN, NMESHES_MAX
   
      M  => MESHES(NM)
      SL => SCARC(NM)%LEVEL(NL)

      ! get right hand side and initial vector
      F (1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR) = M%PRHS (1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR)
      IF (PREDICTOR) THEN
         SL%BX (1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR) = M%H (1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR)
      ELSE
         SL%BX (1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR) = M%HS(1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR)
      ENDIF

      WALL_CELL_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
      
         I = M%IJKW(6,IW)
         J = M%IJKW(7,IW)
         K = M%IJKW(8,IW)
      
         IOR0 = M%IJKW(4,IW)
      
         SELECT CASE(IOR0)
            CASE(1)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DXI2 * M%BXS(J,K)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) + SL%DXI * M%BXS(J,K)                   ! Neumann
               ENDIF
            CASE(-1)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(J,K)          ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) - SL%DXI *M%BXF(J,K)                    ! Neumann
               ENDIF
            CASE(2)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DYI2 * M%BYS(I,K)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) + SL%DYI * M%BYS(I,K)                   ! Neumann
               ENDIF
            CASE(-2)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DYI2 *M%BYF(I,K)          ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) - SL%DYI *M%BYF(I,K)                    ! Neumann
               ENDIF
            CASE(3)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,J)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) + SL%DZI * M%BZS(I,J)                   ! Neumann
               ENDIF
            CASE(-3)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(I,J,K) = F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,J)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(I,J,K) = F(I,J,K) - SL%DZI  * M%BZF(I,J)                  ! Neumann
               ENDIF
         END SELECT
      
      ENDDO WALL_CELL_LOOP3D
   ENDDO MESH_LOOP_BANDWISE3D
ENDIF DIMENSION_BANDWISE_IF

END SUBROUTINE SCARC_INITIALIZE_SOLUTION_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set initial solution corresponding to H and HS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_SOLUTION_COMPACT(X, F, NL)
INTEGER, INTENT(IN) :: NL
REAL(EB), POINTER, DIMENSION(:) :: X, F
INTEGER :: NM, IW, IOR0, I, J, K, IC

!!!----------------------------------------------------------------------------------------------------
!!!  2D-version
!!!----------------------------------------------------------------------------------------------------
DIMENSION_COMPACT_IF: IF (TWO_D) THEN

   MESH_LOOP_COMPACT2D: DO NM = NMESHES_MIN, NMESHES_MAX
   
      M  => MESHES(NM)
   
      ! get right hand side and initial vector
      DO K = 1, M%KBAR
         DO I = 1, M%IBAR
            IC = (K-1)*SL%IBAR + I
            F(IC) = M%PRHS (I, 1, K)
            IF (PREDICTOR) THEN
               X(IC) = M%H (I, 1, K)
            ELSE
               X(IC) = M%HS(I, 1, K)
            ENDIF
         ENDDO
      ENDDO

      ! set boundary conditions at exterior boundaries (corresponding to pois.f90)
      WALL_CELL_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
   
         I = M%IJKW(6,IW)
         J = M%IJKW(7,IW)
         K = M%IJKW(8,IW)

         IF (J /= 1) THEN
            WRITE(*,*) 'Wrong index for J =',J,' in SCARC_INITIALIZE_SOLUTION !!!'
            STOP
         ENDIF
      
         IOR0 = M%IJKW(4,IW)
         IC = (K-1)*M%IBAR + I
         
         SELECT CASE(IOR0)
            CASE(1)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DXI2 * M%BXS(1,K)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) + SL%DXI * M%BXS(1,K)                   ! Neumann
               ENDIF
            CASE(-1)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DXI2 *M%BXF(1,K)          ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) - SL%DXI *M%BXF(1,K)                    ! Neumann
               ENDIF
            CASE(3)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DZI2 * M%BZS(I,1)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) + SL%DZI * M%BZS(I,1)                   ! Neumann
               ENDIF
            CASE(-3)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DZI2 * M%BZF(I,1)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) - SL%DZI  * M%BZF(I,1)                  ! Neumann
               ENDIF
         END SELECT
      
      ENDDO WALL_CELL_LOOP2D

   ENDDO MESH_LOOP_COMPACT2D

!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE

   MESH_LOOP_COMPACT3D: DO NM = NMESHES_MIN, NMESHES_MAX
   
      M  => MESHES(NM)
      SL => SCARC(NM)%LEVEL(NL)

      ! get right hand side and initial vector
      DO K = 1, M%KBAR
         DO J = 1, M%JBAR
            DO I = 1, M%IBAR
               IC = (K-1) * M%IBAR * JBAR + (J-1) * M%IBAR + I
               F(IC) = M%PRHS (I, J, K)
               IF (PREDICTOR) THEN
                  X(IC) = M%H (I, J, K)
               ELSE
                  X(IC) = M%HS(I, J, K)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      WALL_CELL_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
      
         I = M%IJKW(6,IW)
         J = M%IJKW(7,IW)
         K = M%IJKW(8,IW)
      
         IOR0 = M%IJKW(4,IW)
         IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
      
         SELECT CASE(IOR0)
            CASE(1)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DXI2 * M%BXS(J,K)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) + SL%DXI * M%BXS(J,K)                   ! Neumann
               ENDIF
            CASE(-1)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DXI2 *M%BXF(J,K)          ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) - SL%DXI *M%BXF(J,K)                    ! Neumann
               ENDIF
            CASE(2)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DYI2 * M%BYS(I,K)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) + SL%DYI * M%BYS(I,K)                   ! Neumann
               ENDIF
            CASE(-2)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DYI2 *M%BYF(I,K)          ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) - SL%DYI *M%BYF(I,K)                    ! Neumann
               ENDIF
            CASE(3)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DZI2 * M%BZS(I,J)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) + SL%DZI * M%BZS(I,J)                   ! Neumann
               ENDIF
            CASE(-3)
               IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                  F(IC) = F(IC) - 2.0_EB * SL%DZI2 * M%BZF(I,J)         ! Dirichlet
               ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                  F(IC) = F(IC) - SL%DZI  * M%BZF(I,J)                  ! Neumann
               ENDIF
         END SELECT
      
      ENDDO WALL_CELL_LOOP3D
   ENDDO MESH_LOOP_COMPACT3D
ENDIF DIMENSION_COMPACT_IF

END SUBROUTINE SCARC_INITIALIZE_SOLUTION_COMPACT


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

1000 FORMAT (5X,A30,': MG-level=',i4,': #ite= ',i4,': res =',e14.6)
END SUBROUTINE SCARC_CONVERGENCE_INFO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Check if solver converges or diverges and print out residual information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)

INTEGER, INTENT(IN) :: ITE, NL
INTEGER :: NM, ISTATE
REAL(EB), INTENT(IN) :: RESIN, RES, EPS
CHARACTER(*), INTENT(IN) :: CROUTINE

ISTATE = NSCARC_LOOP_PROCEED

DO NM = NMESHES_MIN, NMESHES_MAX
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS)          WRITE(SCARC_LU,1000) TRIM(CROUTINE), NL, ITE, RES
ENDDO
IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS.AND.MYID==0) WRITE(*       ,1000) TRIM(CROUTINE), NL, ITE, RES

IF (TYPE_ACCURACY == NSCARC_ACCURACY_RELATIVE) THEN
   IF (RES <= RESIN*EPS)  ISTATE = NSCARC_LOOP_CONV
ELSE
   IF (RES <= EPS .AND. RES <= RESIN*SCARC_ACCURACY_RELATIVE) ISTATE = NSCARC_LOOP_CONV
ENDIF
IF (RES > SCARC_ACCURACY_DIVERGENCE) ISTATE = NSCARC_LOOP_DIVG

SCARC_CONVERGENCE_STATE = ISTATE
RETURN

1000 FORMAT (5X,A30,': MG-level=',i4,': #ite= ',i4,': res =',e14.6)
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
IF (ISTATE == NSCARC_LOOP_DIVG) THEN                       
   SCARC_ITERATIONS   = - 1
   SCARC_CAPPA = 1.0_EB
ELSE 
   IF (ISTATE == NSCARC_LOOP_CONV) THEN
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
!!! Perform matrix-vector multiplication Y = A*X  - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_LOCAL_MATVEC_BANDWISE (A, X, Y, NX, NY, NZ, NM)
REAL (EB), POINTER, DIMENSION (:, :),    INTENT(IN)    :: A
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(IN)    :: X
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(INOUT) :: Y
INTEGER,   POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER,   INTENT(IN) :: NM
INTEGER :: I, J, K, IC
REAL (EB):: TNOW_LOCAL_MATVEC
TNOW_LOCAL_MATVEC = SECOND()
 
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
CALL SCARC_SHOW_LEVEL (X, NX, NY, NZ, 'MATVEC', 'X matvec ')
CALL SCARC_SHOW_LEVEL (Y, NX, NY, NZ, 'MATVEC', 'Y matvec ')
ENDIF

IF (TWO_D) THEN
   DO K = 1, NZ
      DO I = 1, NX
         IC = (K-1) * NX + I
         Y (I, 1, K) =   (  A(IC, ID ) * X(I  , 1, K  )   &          ! diagonal component
                          + A(IC, ILZ) * X(I  , 1, K-1)   &          ! lower z-component
                          + A(IC, ILX) * X(I-1, 1, K  )   &          ! lower x-component
                          + A(IC, IUX) * X(I+1, 1, K  )   &          ! upper x-component
                          + A(IC, IUZ) * X(I  , 1, K+1) )            ! upper z-component
      ENDDO
   ENDDO
ELSE
   DO K = 1, NZ
      DO J = 1, NY
         DO I = 1, NX
            IC = (K-1) * NX * NY + (J-1) * NX + I
            Y(I, J, K) =   ( A(IC, ID ) * X (I  , J  , K  )   &      ! diagonal component
                         +   A(IC, ILZ) * X (I  , J  , K-1)   &      ! lower z-component
                         +   A(IC, ILY) * X (I  , J-1, K  )   &      ! lower y-component
                         +   A(IC, ILX) * X (I-1, J  , K  )   &      ! lower x-component
                         +   A(IC, IUX) * X (I+1, J  , K  )   &      ! upper x-component
                         +   A(IC, IUY) * X (I  , J+1, K  )   &      ! upper y-component
                         +   A(IC, IUZ) * X (I  , J  , K+1) )        ! upper z-component
         ENDDO
      ENDDO
   ENDDO
ENDIF 

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
CALL SCARC_SHOW_LEVEL (Y, NX, NY, NZ, 'MATVEC', 'Y matvec2 ')
ENDIF

TUSED_SCARC(NSCARC_TIME_LOCAL_MATVEC,NM)=TUSED_SCARC(NSCARC_TIME_LOCAL_MATVEC,NM)+SECOND()-TNOW_LOCAL_MATVEC
TUSED_SCARC(NSCARC_TIME_COMPLETE    ,NM)=TUSED_SCARC(NSCARC_TIME_COMPLETE    ,NM)+SECOND()-TNOW_LOCAL_MATVEC
END SUBROUTINE SCARC_LOCAL_MATVEC_BANDWISE
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform matrix-vector multiplication Y = A*X  - compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_LOCAL_MATVEC_COMPACT (A, ROW, COL, X, Y, NX, NY, NZ, NM)

REAL (EB), POINTER, DIMENSION (:), INTENT(IN)    :: A
INTEGER,   POINTER, DIMENSION (:), INTENT(IN)    :: ROW, COL
REAL (EB), POINTER, DIMENSION (:), INTENT(IN)    :: X
REAL (EB), POINTER, DIMENSION (:), INTENT(INOUT) :: Y
INTEGER,   POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER,   INTENT(IN) :: NM
INTEGER :: I, J, K, IC, ICOL
REAL (EB):: TNOW_LOCAL_MATVEC
TNOW_LOCAL_MATVEC = SECOND()

DO K = 1, NZ
   DO J = 1, NY
      DO I = 1, NX

         IC = (K-1) * NX * NY + (J-1) * NX + I
        
         ICOL = ROW(IC)                                    ! diagonal entry
         Y (IC) =  A(ICOL)* X(COL(ICOL))

         DO ICOL = ROW(IC)+1, ROW(IC+1)-1                  ! subdiagonal entries
            Y (IC) =  Y(IC) + A(ICOL)* X(COL(ICOL))
         ENDDO

      ENDDO
   ENDDO
ENDDO
      
TUSED_SCARC(NSCARC_TIME_LOCAL_MATVEC,NM)=TUSED_SCARC(NSCARC_TIME_LOCAL_MATVEC,NM)+SECOND()-TNOW_LOCAL_MATVEC
TUSED_SCARC(NSCARC_TIME_COMPLETE    ,NM)=TUSED_SCARC(NSCARC_TIME_COMPLETE    ,NM)+SECOND()-TNOW_LOCAL_MATVEC
END SUBROUTINE SCARC_LOCAL_MATVEC_COMPACT
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute global matrix-vector products by exchanging corresponding matrix-vector
!!! information along internal boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_GLOBAL_MATVEC(TYPE1, TYPE2, NL)
INTEGER, INTENT(IN) :: TYPE1, TYPE2, NL
REAL (EB):: TNOW_GLOBAL_MATVEC
 
TNOW_GLOBAL_MATVEC = SECOND()

IF (NMESHES > 1) THEN
   NREQ_SCARC = 0
   CALL SCARC_RECEIVE  (NSCARC_EXCHANGE_MATVEC, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATVEC, NL, TYPE1, TYPE2)
ENDIF

TUSED_SCARC(NSCARC_TIME_GLOBAL_MATVEC,:)=TUSED_SCARC(NSCARC_TIME_GLOBAL_MATVEC,:)+SECOND()-TNOW_GLOBAL_MATVEC
TUSED_SCARC(NSCARC_TIME_COMPLETE     ,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE     ,:)+SECOND()-TNOW_GLOBAL_MATVEC
END SUBROUTINE SCARC_GLOBAL_MATVEC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform update of internal boundaries to get neighbouring data 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_EXCHANGE_BDRY(NL)
INTEGER, INTENT(IN) :: NL
REAL (EB):: TNOW_EXCHANGE_BDRY
 
TNOW_EXCHANGE_BDRY = SECOND()

IF (NMESHES > 1) THEN
   NREQ_SCARC = 0
   CALL SCARC_RECEIVE  (NSCARC_EXCHANGE_BDRY, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_BDRY, NL, NSCARC_DUMMY, NSCARC_DUMMY)
ENDIF

TUSED_SCARC(NSCARC_TIME_EXCHANGE_BDRY,:)=TUSED_SCARC(NSCARC_TIME_EXCHANGE_BDRY,:)+SECOND()-TNOW_EXCHANGE_BDRY
TUSED_SCARC(NSCARC_TIME_COMPLETE     ,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE     ,:)+SECOND()-TNOW_EXCHANGE_BDRY
END SUBROUTINE SCARC_EXCHANGE_BDRY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute local scalarproduct of vector X and vector Y in 3D: SP=(X,Y)
!!! - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_LOCAL_SCALPROD_BANDWISE (X, Y, NX, NY, NZ, NM)
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(IN) :: X, Y
INTEGER, POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER, INTENT(IN) :: NM
INTEGER :: I, J, K
REAL(EB):: TNOW_LOCAL_SCALPROD
 
TNOW_LOCAL_SCALPROD = SECOND()
 
SP_LOCAL(NM) = 0.0_EB
DO K = 1, NZ
   DO J = 1, NY
      DO I = 1, NX
         SP_LOCAL(NM) = SP_LOCAL(NM) + X (I, J, K) * Y (I, J, K)
      ENDDO
   ENDDO
ENDDO
 
TUSED_SCARC(NSCARC_TIME_LOCAL_SCALPROD,NM)=TUSED_SCARC(NSCARC_TIME_LOCAL_SCALPROD,NM)+SECOND()-TNOW_LOCAL_SCALPROD
TUSED_SCARC(NSCARC_TIME_COMPLETE      ,NM)=TUSED_SCARC(NSCARC_TIME_COMPLETE      ,NM)+SECOND()-TNOW_LOCAL_SCALPROD
END SUBROUTINE SCARC_LOCAL_SCALPROD_BANDWISE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute local scalarproduct of vector X and vector Y in 3D: SP=(X,Y)
!!! - compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_LOCAL_SCALPROD_COMPACT (X, Y, NX, NY, NZ, NM)
REAL (EB), POINTER, DIMENSION (:), INTENT(IN) :: X, Y
INTEGER, POINTER, INTENT(IN) :: NX, NY, NZ
INTEGER, INTENT(IN) :: NM
INTEGER :: I, J, K, IC
REAL(EB):: TNOW_LOCAL_SCALPROD
 
TNOW_LOCAL_SCALPROD = SECOND()
 
SP_LOCAL(NM) = 0.0_EB
DO K = 1, NZ
   DO J = 1, NY
      DO I = 1, NX
         IC = (K-1) * NX * NY + (J-1) * NX + I
         SP_LOCAL(NM) = SP_LOCAL(NM) + X(IC) * Y(IC)
      ENDDO
   ENDDO
ENDDO
 
TUSED_SCARC(NSCARC_TIME_LOCAL_SCALPROD,NM)=TUSED_SCARC(NSCARC_TIME_LOCAL_SCALPROD,NM)+SECOND()-TNOW_LOCAL_SCALPROD
TUSED_SCARC(NSCARC_TIME_COMPLETE      ,NM)=TUSED_SCARC(NSCARC_TIME_COMPLETE      ,NM)+SECOND()-TNOW_LOCAL_SCALPROD
END SUBROUTINE SCARC_LOCAL_SCALPROD_COMPACT
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get global sum of local scalar products by global data exchange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_GLOBAL_SCALPROD(NL)
INTEGER, INTENT(IN):: NL
REAL(EB):: SP
INTEGER :: IERR, NM, NL0
REAL(EB):: TNOW_GLOBAL_SCALPROD
 
TNOW_GLOBAL_SCALPROD = SECOND()

IERR = 0
NL0  = NL
SP   = 0.0_EB
IF (NMESHES>1 .AND. USE_MPI) THEN
   CALL MPI_ALLREDUCE (SP_LOCAL, SP, NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   DO NM=1,NMESHES
      SP = SP + SP_LOCAL(NM)
   ENDDO
ENDIF

SCARC_GLOBAL_SCALPROD=SP
RETURN

TUSED_SCARC(NSCARC_TIME_GLOBAL_SCALPROD,:)=TUSED_SCARC(NSCARC_TIME_GLOBAL_SCALPROD,:)+SECOND()-TNOW_GLOBAL_SCALPROD
TUSED_SCARC(NSCARC_TIME_COMPLETE       ,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE       ,:)+SECOND()-TNOW_GLOBAL_SCALPROD
END FUNCTION SCARC_GLOBAL_SCALPROD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get global L2-norm by global data exchange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_GLOBAL_L2NORM(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: IERR, NM
REAL(EB):: SP
REAL(EB):: TNOW_GLOBAL_L2NORM
TNOW_GLOBAL_L2NORM = SECOND()

IERR = 0
SP   = 0.0_EB
IF (NMESHES>1 .AND. USE_MPI) THEN
   CALL MPI_ALLREDUCE (SP_LOCAL, SP, NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   DO NM=1,NMESHES
      SP = SP + SP_LOCAL(NM)
   ENDDO
ENDIF
SCARC_GLOBAL_L2NORM = SQRT (SP/REAL(NC_GLOBAL(NL), EB))   
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'GLOBAL SP=',SP
RETURN

TUSED_SCARC(NSCARC_TIME_GLOBAL_L2NORM,:)=TUSED_SCARC(NSCARC_TIME_GLOBAL_L2NORM,:)+SECOND()-TNOW_GLOBAL_L2NORM
TUSED_SCARC(NSCARC_TIME_COMPLETE     ,:)=TUSED_SCARC(NSCARC_TIME_COMPLETE     ,:)+SECOND()-TNOW_GLOBAL_L2NORM
END FUNCTION SCARC_GLOBAL_L2NORM

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Restrict vector X from level NL to vector Y on level NL-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RESTRICTION_BANDWISE (F_LO, D_HI, NX_LO, NY_LO, NZ_LO)
REAL(EB), POINTER, DIMENSION(:, :, :), INTENT(IN)  :: D_HI
REAL(EB), POINTER, DIMENSION(:, :, :), INTENT(OUT) :: F_LO
INTEGER, INTENT(IN) :: NX_LO, NY_LO, NZ_LO
INTEGER :: NM
INTEGER :: IX_HI, IY_HI, IZ_HI
INTEGER :: IX_LO, IY_LO, IZ_LO

!!!----------------------------------------------------------------------------------------------------
!!! 2D-version
!!!----------------------------------------------------------------------------------------------------
DIMENSION_IF: IF (TWO_D) THEN

   MESH_LOOP2D: DO NM = NMESHES_MIN, NMESHES_MAX

      DO IZ_LO = 1, NZ_LO
         DO IX_LO = 1, NX_LO

            IX_HI = 2*IX_LO
            IY_HI = 1
            IZ_HI = 2*IZ_LO

            F_LO(IX_LO, 1, IZ_LO) = 0.25_EB * (  D_HI(IX_HI  , IY_HI, IZ_HI-1)  &
                                               + D_HI(IX_HI-1, IY_HI, IZ_HI-1)  &
                                               + D_HI(IX_HI  , IY_HI, IZ_HI)  &
                                               + D_HI(IX_HI-1, IY_HI, IZ_HI) )
         ENDDO
      ENDDO

   ENDDO MESH_LOOP2D

!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE
   
   MESH_LOOP3D: DO NM = NMESHES_MIN, NMESHES_MAX

      DO IZ_LO = 1, NZ_LO
         DO IY_LO = 1, NY_LO
            DO IX_LO = 1, NX_LO

               IX_HI = 2*IX_LO
               IY_HI = 2*IY_HI
               IZ_HI = 2*IZ_LO

               F_LO(IX_LO, IY_LO , IZ_LO) = 0.125_EB * (  D_HI(IX_HI-1, IY_HI-1, IZ_HI-1)  &
                                                        + D_HI(IX_HI  , IY_HI-1, IZ_HI-1)  &
                                                        + D_HI(IX_HI-1, IY_HI  , IZ_HI-1)  &
                                                        + D_HI(IX_HI  , IY_HI  , IZ_HI-1)  &
                                                        + D_HI(IX_HI-1, IY_HI-1, IZ_HI  )  &
                                                        + D_HI(IX_HI  , IY_HI-1, IZ_HI  )  &
                                                        + D_HI(IX_HI-1, IY_HI  , IZ_HI  )  &
                                                        + D_HI(IX_HI  , IY_HI  , IZ_HI  ) )
            ENDDO
         ENDDO
      ENDDO

   ENDDO MESH_LOOP3D
ENDIF DIMENSION_IF

END SUBROUTINE SCARC_RESTRICTION_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Restrict vector X from level NL to vector Y on level NL-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RESTRICTION_COMPACT (F_LO, D_HI, NX_LO, NY_LO, NZ_LO)
REAL(EB), POINTER, DIMENSION(:), INTENT(IN)  :: D_HI
REAL(EB), POINTER, DIMENSION(:), INTENT(OUT) :: F_LO
INTEGER, INTENT(IN) :: NX_LO, NY_LO, NZ_LO
INTEGER :: NM, IC_LO, IC_HI(8)
INTEGER :: NX_HI, NY_HI, NZ_HI
INTEGER :: IX_LO, IY_LO, IZ_LO
INTEGER :: IX_HI, IY_HI, IZ_HI

NX_HI = 2*NX_LO
NY_HI = 2*NY_LO
NZ_HI = 2*NZ_LO

!!!----------------------------------------------------------------------------------------------------
!!! 2D-version
!!!----------------------------------------------------------------------------------------------------
DIMENSION_IF: IF (TWO_D) THEN

   MESH_LOOP2D: DO NM = NMESHES_MIN, NMESHES_MAX

      DO IZ_LO = 1, NZ_LO
         DO IX_LO = 1, NX_LO
 
            IX_HI = 2*IX_LO
            IZ_HI = 2*IZ_LO

            IC_LO     = (IZ_LO-1)*NX_LO + IX_LO

            IC_HI(1) = (IZ_HI-1)*NX_HI + IX_HI - 1
            IC_HI(2) = (IZ_HI-1)*NX_HI + IX_HI   
            IC_HI(3) =  IZ_HI   *NX_HI + IX_HI - 1
            IC_HI(4) =  IZ_HI   *NX_HI + IX_HI    

            F_LO(IC_LO) = 0.25_EB * (  D_HI(IC_HI(1)) &
                                     + D_HI(IC_HI(2)) &
                                     + D_HI(IC_HI(3)) &
                                     + D_HI(IC_HI(4)) )

         ENDDO
      ENDDO

   ENDDO MESH_LOOP2D

!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE
   
   MESH_LOOP3D: DO NM = NMESHES_MIN, NMESHES_MAX

      DO IZ_LO = 1, NZ_LO
         DO IY_LO = 1, NY_LO
            DO IX_LO = 1, NX_LO
 
               IX_HI = 2*IX_LO
               IY_HI = 2*IY_LO
               IZ_HI = 2*IZ_LO

               IC_LO    = (IZ_LO-1)*NX_LO*NY_LO + (IY_LO-1)*NX_LO + IX_LO

               IC_HI(1) = (IZ_HI-1)*NX_HI*NY_HI + (IY_HI-1)*NX_HI + IX_HI - 1
               IC_HI(2) = (IZ_HI-1)*NX_HI*NY_HI + (IY_HI-1)*NX_HI + IX_HI    
               IC_HI(3) = (IZ_HI-1)*NX_HI*NY_HI +  IY_HI   *NX_HI + IX_HI - 1
               IC_HI(4) = (IZ_HI-1)*NX_HI*NY_HI +  IY_HI   *NX_HI + IX_HI    
               IC_HI(5) =  IZ_HI   *NX_HI*NY_HI + (IY_HI-1)*NX_HI + IX_HI - 1
               IC_HI(6) =  IZ_HI   *NX_HI*NY_HI + (IY_HI-1)*NX_HI + IX_HI    
               IC_HI(7) =  IZ_HI   *NX_HI*NY_HI +  IY_HI   *NX_HI + IX_HI - 1
               IC_HI(8) =  IZ_HI   *NX_HI*NY_HI +  IY_HI   *NX_HI + IX_HI    

               F_LO(IC_LO) = 0.125_EB * (  D_HI(IC_HI(1)) &
                                         + D_HI(IC_HI(2)) &
                                         + D_HI(IC_HI(3)) &
                                         + D_HI(IC_HI(4)) &
                                         + D_HI(IC_HI(5)) &
                                         + D_HI(IC_HI(6)) &
                                         + D_HI(IC_HI(7)) &
                                         + D_HI(IC_HI(8)) )
            ENDDO
         ENDDO
      ENDDO

   ENDDO MESH_LOOP3D
ENDIF DIMENSION_IF

END SUBROUTINE SCARC_RESTRICTION_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Restrict vector X from level NL to vector Y on level NL-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PROLONGATION_BANDWISE (X_LO, D_HI, NX_LO, NY_LO, NZ_LO)
REAL(EB), POINTER, DIMENSION(:, :, :), INTENT(IN) :: X_LO
REAL(EB), POINTER, DIMENSION(:, :, :), INTENT(OUT):: D_HI
INTEGER, INTENT(IN) :: NX_LO, NY_LO, NZ_LO
INTEGER :: NM
INTEGER :: IX_LO, IY_LO, IZ_LO
INTEGER :: IX_HI, IY_HI, IZ_HI

!!!----------------------------------------------------------------------------------------------------
!!! 2D-version
!!!----------------------------------------------------------------------------------------------------
DIMENSION_IF: IF (TWO_D) THEN

   MESH_LOOP2D: DO NM = NMESHES_MIN, NMESHES_MAX

      DO IZ_LO = 1, NZ_LO
         DO IX_LO = 1, NX_LO

            IX_HI = 2*IX_LO
            IY_HI = 1
            IZ_HI = 2*IZ_LO

            D_HI(IX_HI-1, 1, IZ_HI-1) = X_LO(IX_LO, 1, IZ_LO)
            D_HI(IX_HI  , 1, IZ_HI-1) = X_LO(IX_LO, 1, IZ_LO)
            D_HI(IX_HI-1, 1, IZ_HI  ) = X_LO(IX_LO, 1, IZ_LO)
            D_HI(IX_HI  , 1, IZ_HI  ) = X_LO(IX_LO, 1, IZ_LO)

         ENDDO
      ENDDO
   ENDDO MESH_LOOP2D

!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE
   
   MESH_LOOP3D: DO NM = NMESHES_MIN, NMESHES_MAX

      DO IZ_LO = 1, NZ_LO
         DO IY_LO = 1, NY_LO
            DO IX_LO = 1, NX_LO

               IX_HI = 2*IX_LO
               IY_HI = 2*IZ_LO
               IZ_HI = 2*IZ_LO

               D_HI(IX_HI-1, IY_HI-1, IZ_HI-1) = X_LO(IX_LO, IY_LO, IZ_LO)
               D_HI(IX_HI  , IY_HI-1, IZ_HI-1) = X_LO(IX_LO, IY_LO, IZ_LO)
               D_HI(IX_HI-1, IY_HI  , IZ_HI-1) = X_LO(IX_LO, IY_LO, IZ_LO)
               D_HI(IX_HI  , IY_HI  , IZ_HI-1) = X_LO(IX_LO, IY_LO, IZ_LO)
               D_HI(IX_HI-1, IY_HI-1, IZ_HI  ) = X_LO(IX_LO, IY_LO, IZ_LO)
               D_HI(IX_HI  , IY_HI-1, IZ_HI  ) = X_LO(IX_LO, IY_LO, IZ_LO)
               D_HI(IX_HI-1, IY_HI  , IZ_HI  ) = X_LO(IX_LO, IY_LO, IZ_LO)
               D_HI(IX_HI  , IY_HI  , IZ_HI  ) = X_LO(IX_LO, IY_LO, IZ_LO)

            ENDDO
         ENDDO
      ENDDO
   ENDDO MESH_LOOP3D
ENDIF DIMENSION_IF

END SUBROUTINE SCARC_PROLONGATION_BANDWISE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Restrict vector X from level NL to vector Y on level NL-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PROLONGATION_COMPACT (X_LO, D_HI, NX_LO, NY_LO, NZ_LO)
REAL(EB), POINTER, DIMENSION(:), INTENT(IN) :: X_LO
REAL(EB), POINTER, DIMENSION(:), INTENT(OUT):: D_HI
INTEGER,  INTENT(IN) :: NX_LO, NY_LO, NZ_LO
INTEGER :: NM, IC_LO, IC_HI(8)
INTEGER :: NX_HI, NY_HI, NZ_HI
INTEGER :: IX_HI, IY_HI, IZ_HI
INTEGER :: IX_LO, IY_LO, IZ_LO

NX_HI = 2*NX_LO
NY_HI = 2*NY_LO
NZ_HI = 2*NZ_LO

!!!----------------------------------------------------------------------------------------------------
!!! 2D-version
!!!----------------------------------------------------------------------------------------------------
DIMENSION_IF: IF (TWO_D) THEN

   MESH_LOOP2D: DO NM = NMESHES_MIN, NMESHES_MAX

      DO IZ_LO = 1, NZ_LO
         DO IX_LO = 1, NX_LO

            IX_HI = 2*IX_LO
            IY_HI = 1
            IZ_HI = 2*IZ_LO

            IC_LO     = (IZ_LO-1)*NX_LO + IX_LO

            IC_HI(1) = (IZ_HI-1)*NX_HI + IX_HI - 1
            IC_HI(2) = (IZ_HI-1)*NX_HI + IX_HI   
            IC_HI(3) =  IZ_HI   *NX_HI + IX_HI - 1
            IC_HI(4) =  IZ_HI   *NX_HI + IX_HI    

            D_HI(IC_HI(1)) = X_LO(IC_LO)
            D_HI(IC_HI(2)) = X_LO(IC_LO)
            D_HI(IC_HI(3)) = X_LO(IC_LO)
            D_HI(IC_HI(4)) = X_LO(IC_LO)

         ENDDO
      ENDDO

   ENDDO MESH_LOOP2D

!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
ELSE
   
   MESH_LOOP3D: DO NM = NMESHES_MIN, NMESHES_MAX

      DO IZ_LO = 1, NZ_LO
         DO IY_LO = 1, NY_LO
            DO IX_LO = 1, NX_LO

               IX_HI = 2*IX_LO
               IY_HI = 2*IY_LO
               IZ_HI = 2*IZ_LO

               IC_LO    = (IZ_LO-1)*NX_LO*NY_LO + (IY_LO-1)*NX_LO + IX_LO

               IC_HI(1) = (IZ_HI-1)*NX_HI*NY_HI + (IY_HI-1)*NX_HI + IX_HI - 1
               IC_HI(2) = (IZ_HI-1)*NX_HI*NY_HI + (IY_HI-1)*NX_HI + IX_HI    
               IC_HI(3) = (IZ_HI-1)*NX_HI*NY_HI +  IY_HI   *NX_HI + IX_HI - 1
               IC_HI(4) = (IZ_HI-1)*NX_HI*NY_HI +  IY_HI   *NX_HI + IX_HI    
               IC_HI(5) =  IZ_HI   *NX_HI*NY_HI + (IY_HI-1)*NX_HI + IX_HI - 1
               IC_HI(6) =  IZ_HI   *NX_HI*NY_HI + (IY_HI-1)*NX_HI + IX_HI    
               IC_HI(7) =  IZ_HI   *NX_HI*NY_HI +  IY_HI   *NX_HI + IX_HI - 1
               IC_HI(8) =  IZ_HI   *NX_HI*NY_HI +  IY_HI   *NX_HI + IX_HI    

               D_HI(IC_HI(1)) = X_LO(IC_LO)
               D_HI(IC_HI(2)) = X_LO(IC_LO)
               D_HI(IC_HI(3)) = X_LO(IC_LO)
               D_HI(IC_HI(4)) = X_LO(IC_LO)
               D_HI(IC_HI(5)) = X_LO(IC_LO)
               D_HI(IC_HI(6)) = X_LO(IC_LO)
               D_HI(IC_HI(7)) = X_LO(IC_LO)
               D_HI(IC_HI(8)) = X_LO(IC_LO)

            ENDDO
         ENDDO
      ENDDO

   ENDDO MESH_LOOP3D
ENDIF DIMENSION_IF

END SUBROUTINE SCARC_PROLONGATION_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Finalize data - bandwise storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_TRANSFER_SOLUTION_BANDWISE(X, NM)
REAL(EB), POINTER, DIMENSION(:, :, :), INTENT(IN):: X
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()

M => MESHES(NM)
IF (PREDICTOR) THEN
WRITE(SCARC_LU,*) 'PREDICTOR'
   HP => M%H
ELSE
WRITE(SCARC_LU,*) 'CORRECTOR'
   HP => M%HS
ENDIF

!!! Overwrite internal values of H or HS by corresponding data of X
HP(1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = X(1:M%IBAR, 1:M%JBAR, 1:M%KBAR)

END SUBROUTINE SCARC_TRANSFER_SOLUTION_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Finalize data - compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_TRANSFER_SOLUTION_COMPACT(X, NM)
REAL(EB), POINTER, DIMENSION(:), INTENT(IN):: X
INTEGER, INTENT(IN) :: NM
INTEGER :: I, J, K, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()

M => MESHES(NM)
IF (PREDICTOR) THEN
   HP => M%H
ELSE
   HP => M%HS
ENDIF

!!! Overwrite internal values of H or HS by corresponding data of X
DO K = 1, M%KBAR
   DO J = 1, M%JBAR
      DO I = 1, M%IBAR
         IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
         HP(I, J, K) = X(IC)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_TRANSFER_SOLUTION_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set correct boundary values at external and internal boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE_GHOSTCELLS(NM)
INTEGER, INTENT(IN) :: NM
INTEGER :: IW, IOR0, I0, J0, K0, I1, J1, K1
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()

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

END SUBROUTINE SCARC_UPDATE_GHOSTCELLS

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Receive data from neighbors (corresponds to POST_RECEIVES)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RECEIVE (CODE, NL)

INTEGER, INTENT(IN) :: CODE, NL
INTEGER :: NM, NOM, IERR, ILEN

IERR=0
RECEIVE_MESH_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   RECEIVE_OMESH_LOOP: DO NOM = 1, NMESHES
    
      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_OMESH_LOOP

      SNM   => SCARC(NM)                         ! corresponds to M
      OSNM  => SCARC(NM)%OSCARC(NOM)             ! corresponds to M3

      IF (OSNM%NICMAX_S==0 .AND. OSNM%NICMAX_R==0) CYCLE RECEIVE_OMESH_LOOP
    
      SNML  => SCARC(NM)%LEVEL(NL)               ! corresponds to M  for the level 'NL'
      OSNML => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)   ! corresponds to M3 for the level 'NL'

      !!!----------------------------------------------------------------------------------------------
      !!! Initialize communication structures for the receiving of data
      !!!---------------------------------------------------------------------------------------------------
      RECEIVE_INIT_IF: IF (CODE==NSCARC_EXCHANGE_INIT) THEN

         IF (NL/=NLEVEL_MAX) THEN  ! for lower levels get neighboring IJKW's
            IF (USE_MPI) THEN
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_IRECV(OSNML%IJKW(1,1),15*OSNML%N_EXTERNAL_WALL_CELLS,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF
         ELSE                                ! on maximum level neighboring IJKW already exists
            OSNML%IJKW => MESHES(NM)%OMESH(NOM)%IJKW     
         ENDIF

         IF (NL==NLEVEL_MAX .AND. OSNML%NIC_R > 0) THEN
            ILEN=(MAX(OSNML%NIC_R, OSNML%NIC_S)+2)*2+1
            ALLOCATE (OSNM%RECV_BUF(ILEN))
            OSNM%RECV_BUF = 0.0_EB
         ENDIF

      ENDIF RECEIVE_INIT_IF

      !!!----------------------------------------------------------------------------------------------
      !!! Perform exchange for matrix-vector multiplication 
      !!!----------------------------------------------------------------------------------------------
      RECEIVE_MATVEC_IF: IF (CODE==NSCARC_EXCHANGE_MATVEC) THEN

         NREQ_SCARC = NREQ_SCARC+1
         IF (USE_MPI) CALL MPI_IRECV(OSNM%RECV_BUF(1),SIZE(OSNM%RECV_BUF),MPI_DOUBLE_PRECISION,&
                                     SNODE,TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

      ENDIF RECEIVE_MATVEC_IF

      !!!----------------------------------------------------------------------------------------------
      !!! Perform full exchange including edge 
      !!!----------------------------------------------------------------------------------------------
      RECEIVE_UPDATE_IF: IF (CODE==NSCARC_EXCHANGE_BDRY) THEN

         NREQ_SCARC = NREQ_SCARC+1
         IF (USE_MPI) CALL MPI_IRECV(OSNM%RECV_BUF(1),SIZE(OSNM%RECV_BUF),MPI_DOUBLE_PRECISION,&
                                     SNODE,TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

      ENDIF RECEIVE_UPDATE_IF

   ENDDO RECEIVE_OMESH_LOOP
ENDDO RECEIVE_MESH_LOOP

END SUBROUTINE SCARC_RECEIVE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pack send buffer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PACK_SEND_BUF (X0, SEND_BUF0, IJKW0, NEWC0, NLEN0, NM)
REAL (EB), POINTER, DIMENSION (:), INTENT(OUT) :: SEND_BUF0
REAL (EB), POINTER, DIMENSION (:,:,:) :: X0
INTEGER  , POINTER, DIMENSION (:,:)   :: IJKW0
INTEGER, POINTER     :: NEWC0
INTEGER, INTENT(OUT) :: NLEN0
INTEGER, INTENT(IN)  :: NM
INTEGER ::  IW, IWW, LL, II, JJ, KK

LL = 0
IWW = 0
PACK_SEND: DO IW=1,NEWC0
   IF (IJKW0(9,IW)/=NM) CYCLE PACK_SEND
   DO KK=IJKW0(12,IW),IJKW0(15,IW)
      DO JJ=IJKW0(11,IW),IJKW0(14,IW)
         DO II=IJKW0(10,IW),IJKW0(13,IW)
            IWW = IWW + 1
            SEND_BUF0(LL+1) = REAL(IW,EB)
            SEND_BUF0(LL+2) = X0(II,JJ,KK)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO
ENDDO PACK_SEND
NLEN0=2*IWW+1

SEND_BUF0(NLEN0) = -999.0_EB

END SUBROUTINE SCARC_PACK_SEND_BUF
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for matrix-vector multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_RECV_BUF1 (X0,RECV_BUF0,IJKW0,ASUB)
REAL (EB), POINTER, DIMENSION (:,:,:), INTENT(OUT) :: X0
REAL (EB), POINTER, DIMENSION (:)    , INTENT(IN)  :: RECV_BUF0
INTEGER  , POINTER, DIMENSION (:,:)  , INTENT(IN)  :: IJKW0
REAL (EB), INTENT(IN) :: ASUB(3)
REAL (EB):: ZSUM
INTEGER IW, LL, ISUM, I, J, K, II, JJ, KK

LL = 0
UNPACK_RECV1: DO

   IW = NINT(RECV_BUF0(LL+1))
   IF (IW==-999) EXIT UNPACK_RECV1
   ZSUM=0.0_EB
   DO KK=IJKW0(12,IW),IJKW0(15,IW)
      DO JJ=IJKW0(11,IW),IJKW0(14,IW)
         DO II=IJKW0(10,IW),IJKW0(13,IW)
            ZSUM=ZSUM+RECV_BUF0(LL+2)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO

   ISUM = (IJKW0(13,IW)-IJKW0(10,IW)+1) * &
          (IJKW0(14,IW)-IJKW0(11,IW)+1) * &
          (IJKW0(15,IW)-IJKW0(12,IW)+1)

   I=IJKW0(6,IW)
   J=IJKW0(7,IW)
   K=IJKW0(8,IW)
   X0(I, J, K) = X0(I, J, K) + ASUB(ABS(IJKW0(4,IW))) * ZSUM/REAL(ISUM,EB)

ENDDO UNPACK_RECV1

END SUBROUTINE SCARC_UNPACK_RECV_BUF1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for internal update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_RECV_BUF2 (X0,RECV_BUF0,IJKW0)
REAL (EB), POINTER, DIMENSION (:,:,:), INTENT(OUT) :: X0
REAL (EB), POINTER, DIMENSION (:)    , INTENT(IN)  :: RECV_BUF0
INTEGER  , POINTER, DIMENSION (:,:)  , INTENT(IN)  :: IJKW0
REAL (EB):: ZSUM
INTEGER IW, LL, ISUM, I, J, K, II, JJ, KK

LL = 0
UNPACK_RECV2: DO

   IW = NINT(RECV_BUF0(LL+1))
   IF (IW==-999) EXIT UNPACK_RECV2
   ZSUM=0.0_EB
   DO KK=IJKW0(12,IW),IJKW0(15,IW)
      DO JJ=IJKW0(11,IW),IJKW0(14,IW)
         DO II=IJKW0(10,IW),IJKW0(13,IW)
            ZSUM=ZSUM+RECV_BUF0(LL+2)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO

   ISUM = (IJKW0(13,IW)-IJKW0(10,IW)+1) * &
          (IJKW0(14,IW)-IJKW0(11,IW)+1) * &
          (IJKW0(15,IW)-IJKW0(12,IW)+1)

   I=IJKW0(1,IW)
   J=IJKW0(2,IW)
   K=IJKW0(3,IW)
   X0(I, J, K) = ZSUM/REAL(ISUM,EB)

ENDDO UNPACK_RECV2

END SUBROUTINE SCARC_UNPACK_RECV_BUF2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Send data to neighbors (corresponds to MESH_EXCHANGE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_EXCHANGE (CODE, NL, TYPE1, TYPE2)

INTEGER, INTENT(IN) :: CODE, NL, TYPE1, TYPE2
INTEGER :: NM, NOM
INTEGER :: IERR
INTEGER :: ILEN, NLEN0
INTEGER , POINTER :: NEWC0
INTEGER , POINTER, DIMENSION(:,:)  :: IJKW0
REAL(EB), POINTER, DIMENSION(:)    :: BUFFER0
REAL(EB), POINTER, DIMENSION(:,:,:):: DATA0

IERR = 0

EXCHANGE_SEND_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX

   IF (PROCESS(NM)/=MYID)  CYCLE EXCHANGE_SEND_LOOP1

   EXCHANGE_RECV_LOOP1: DO NOM = 1, NMESHES
    
      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      SNM   => SCARC(NM)                            ! corresponds to M
      OSNM  => SCARC(NM)%OSCARC(NOM)                ! corresponds to M3

      IF (OSNM%NICMAX_S == 0 .AND. OSNM%NICMAX_R == 0) CYCLE EXCHANGE_RECV_LOOP1

      SNML  => SCARC(NM)%LEVEL(NL)                  ! corresponds to M  for the level 'NL'
      OSNML => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)      ! corresponds to M3 for the level 'NL'

      !!!----------------------------------------------------------------------------------------------
      !!! Initialize the communication structures for sending data
      !!!----------------------------------------------------------------------------------------------
      EXCHANGE_INIT_IF: IF (CODE == NSCARC_EXCHANGE_INIT) THEN
 
         ! in case of multigrid send IJKW's from lower levels (parallel) or taken neighboring ones (serial)
         ! (on max level existing M%IJKW is used)
         IF (RNODE /= SNODE) THEN
            IF (NL /= NLEVEL_MAX) THEN
               NREQ_SCARC = NREQ_SCARC+1
               IF (USE_MPI) CALL MPI_ISEND(SNML%IJKW(1,1),15*SNML%N_EXTERNAL_WALL_CELLS,MPI_INTEGER,SNODE, &
                                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF
         ELSE
            SNOML => SCARC(NOM)%LEVEL(NL)
            OSNML%IJKW => SNOML%IJKW(:,1:SNOML%N_EXTERNAL_WALL_CELLS)
         ENDIF

         ! send buffer from maximum level is used for all levels
         IF (NL==NLEVEL_MAX) THEN
            ILEN=(MAX(OSNML%NIC_R, OSNML%NIC_S)+2)*2+1   
            ALLOCATE (OSNM%SEND_BUF(ILEN))
            OSNM%SEND_BUF = 0.0_EB
         ENDIF

      ENDIF EXCHANGE_INIT_IF

      !!!----------------------------------------------------------------------------------------------
      !!! Perform exchange for matrix-vector multiplication 
      !!!----------------------------------------------------------------------------------------------
      EXCHANGE_MATVEC_IF: IF (CODE==NSCARC_EXCHANGE_MATVEC) THEN

         SELECT CASE(TYPE1)
            CASE(NSCARC_VECTOR_X)
               DATA0 => SNML%BX
            CASE(NSCARC_VECTOR_Y)
               DATA0 => SNML%BY
            CASE(NSCARC_VECTOR_G)
               DATA0 => SNML%BG
            CASE(NSCARC_VECTOR_R)
               DATA0 => SNML%BR
            CASE(NSCARC_VECTOR_D)
               DATA0 => SNML%BD
            CASE(NSCARC_VECTOR_Z)
               DATA0 => SNML%BZ
            CASE(NSCARC_VECTOR_X2)
               DATA0 => SNML%BX2
            CASE(NSCARC_VECTOR_D2)
               DATA0 => SNML%BD2
            CASE(NSCARC_VECTOR_R2)
               DATA0 => SNML%BR2
            CASE(NSCARC_VECTOR_Y2)
               DATA0 => SNML%BY2
         END SELECT

         BUFFER0 => OSNM%SEND_BUF
         IJKW0   => OSNML%IJKW
         NEWC0   => OSNML%N_EXTERNAL_WALL_CELLS

         CALL SCARC_PACK_SEND_BUF(DATA0, BUFFER0, IJKW0, NEWC0, NLEN0, NM)

         IF (RNODE/=SNODE) THEN
            NREQ_SCARC=NREQ_SCARC+1
            IF (USE_MPI) CALL MPI_ISEND(BUFFER0,NLEN0,MPI_DOUBLE_PRECISION,SNODE, &
                                        TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
         ENDIF

      ENDIF  EXCHANGE_MATVEC_IF


   !!! Perform update of internal boundaries
      EXCHANGE_UPDATE_IF: IF (CODE==NSCARC_EXCHANGE_BDRY) THEN

         IF (PREDICTOR) THEN
            DATA0   => MESHES(NM)%H
         ELSE
            DATA0   => MESHES(NM)%HS
         ENDIF
         BUFFER0 => OSNM%SEND_BUF
         IJKW0   => OSNML%IJKW
         NEWC0   => OSNML%N_EXTERNAL_WALL_CELLS

         CALL SCARC_PACK_SEND_BUF(DATA0, BUFFER0, IJKW0, NEWC0, NLEN0, NM)

         IF (RNODE/=SNODE) THEN
            NREQ_SCARC=NREQ_SCARC+1
            IF (USE_MPI) CALL MPI_ISEND(OSNM%SEND_BUF(1),NLEN0,MPI_DOUBLE_PRECISION,SNODE, &
                                        TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
         ENDIF

      ENDIF EXCHANGE_UPDATE_IF

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
    
      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      SNOML => SCARC(NOM)%LEVEL(NL)  
      OSNOM => SCARC(NOM)%OSCARC(NM)
      OSNM  => SCARC(NM)%OSCARC(NOM)

      EXCHANGE_RECV_IF: IF (OSNOM%NICMAX_S/=0 .AND. OSNOM%NICMAX_R/=0) THEN

         IF (RNODE/=SNODE) THEN
            BUFFER0 => OSNOM%RECV_BUF
         ELSE
            BUFFER0 => OSNM%SEND_BUF
         ENDIF

         RECV_BUF_MATV_IF: IF (CODE==NSCARC_EXCHANGE_MATVEC) THEN
     
            SELECT CASE(TYPE2)
               CASE(NSCARC_VECTOR_X)
                  DATA0 => SNOML%BX
               CASE(NSCARC_VECTOR_Y)
                  DATA0 => SNOML%BY
               CASE(NSCARC_VECTOR_G)
                  DATA0 => SNOML%BG
               CASE(NSCARC_VECTOR_R)
                  DATA0 => SNOML%BR
               CASE(NSCARC_VECTOR_D)
                  DATA0 => SNOML%BD
               CASE(NSCARC_VECTOR_X2)
                  DATA0 => SNOML%BX2
               CASE(NSCARC_VECTOR_D2)
                  DATA0 => SNOML%BD2
               CASE(NSCARC_VECTOR_R2)
                  DATA0 => SNOML%BR2
               CASE(NSCARC_VECTOR_Y2)
                  DATA0 => SNOML%BY2
            END SELECT

            CALL SCARC_UNPACK_RECV_BUF1(DATA0, BUFFER0, SNOML%IJKW, SNOML%ASUB)

         ENDIF RECV_BUF_MATV_IF

         !!! Extract data for full communication including edge (3D!) and vertex data
         RECV_BUF_FULL_IF: IF (CODE==NSCARC_EXCHANGE_BDRY) THEN

            IF (PREDICTOR) THEN
               DATA0 => MESHES(NOM)%H
            ELSE
               DATA0 => MESHES(NOM)%HS
            ENDIF
            CALL SCARC_UNPACK_RECV_BUF2(DATA0, BUFFER0, SNOML%IJKW)

         ENDIF RECV_BUF_FULL_IF
      ENDIF EXCHANGE_RECV_IF
   ENDDO EXCHANGE_RECV_LOOP2
ENDDO EXCHANGE_SEND_LOOP2

END SUBROUTINE SCARC_EXCHANGE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! check if two values differ only by a given tolerance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB), INTENT(IN) :: VAL1, VAL2
REAL (EB) :: TOL
TOL = 1.0E-8_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
RETURN
END FUNCTION MATCH
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only for debugging reasons: print out vector information on level NL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW_BANDWISE (X, CROUTINE, CNAME, NM)
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(IN) :: X
INTEGER, INTENT(IN):: NM
REAL (EB):: VALUES(10)
INTEGER :: II, JJ, KK, IBAR8,JBAR8,KBAR8
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
 
IBAR8=MIN(8,MESHES(NM)%IBAR)
JBAR8=MIN(8,MESHES(NM)%JBAR)
KBAR8=MIN(8,MESHES(NM)%KBAR)

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '===================================================================================='
   WRITE(SCARC_LU,2000) CROUTINE, CNAME
   WRITE(SCARC_LU,*) '===================================================================================='
   IF (TWO_D) THEN
      DO KK = KBAR8, 1, - 1
         DO II=1,IBAR8
            IF (ABS(X(II,1,KK))<1.0E-14_EB) THEN
               VALUES(II)=0.0_EB
            ELSE
               VALUES(II)=X(II,1,KK)
            ENDIF
         ENDDO
         WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, IBAR8)
      ENDDO
      WRITE(SCARC_LU,*)  '------------------------------------------------',&
                         '---------------------------------------------------'
      WRITE(SCARC_LU,1000) ('I = ',II,II=1,IBAR8)
   ELSE
      DO KK = KBAR8, 1, - 1
         WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                     '----------------------------------------'
         DO JJ = JBAR8, 1, - 1
            DO II=1,IBAR8
               IF (ABS(X(II,JJ,KK))<1.0E-14_EB) THEN
                  VALUES(II)=0.0_EB
               ELSE
                  VALUES(II)=X(II,JJ,KK)
               ENDIF
            ENDDO
            WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, IBAR8)
         ENDDO
      ENDDO
      WRITE(SCARC_LU,*)  '------------------------------------------------',&
                         '---------------------------------------------------'
      WRITE(SCARC_LU,1000) ('I = ',II,II=1,IBAR8)
      WRITE(SCARC_LU,*)  '--------------------------------------------------',&
                         '-------------------------------------------------'
      WRITE(SCARC_LU,*)
   ENDIF
ENDIF

1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4) 
END SUBROUTINE SCARC_SHOW_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only for debugging reasons: print out vector information on level NL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW_COMPACT (X, CROUTINE, CNAME, NM)
REAL (EB), POINTER, DIMENSION (:), INTENT(IN) :: X
INTEGER, INTENT(IN):: NM
REAL (EB):: VALUES(8)
INTEGER :: II, JJ, KK, IBAR8,JBAR8,KBAR8, IC
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
 
M => MESHES(NM)

IBAR8=MIN(8,M%IBAR)
JBAR8=MIN(8,M%JBAR)
KBAR8=MIN(8,M%KBAR)

WRITE(SCARC_LU,*) '===================================================================================='
WRITE(SCARC_LU,2000) CROUTINE, CNAME
WRITE(SCARC_LU,*) '===================================================================================='
IF (TWO_D) THEN
   DO KK = KBAR8, 1, - 1
      DO II=1,IBAR8
         IC = (KK-1)*M%IBAR + II
         IF (ABS(X(IC))<1.0E-14_EB) THEN
            VALUES(II)=0.0_EB
         ELSE
            VALUES(II)=X(IC)
         ENDIF
      ENDDO
      WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',1,' : ',(VALUES(II), II=1, IBAR8)
   ENDDO
   WRITE(SCARC_LU,*)  '------------------------------------------------',&
                      '---------------------------------------------------'
   WRITE(SCARC_LU,1000) ('I = ',II,II=1,IBAR8)
ELSE
   DO KK = KBAR8, 1, - 1
      WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                  '----------------------------------------'
      DO JJ = JBAR8, 1, - 1
         DO II=1,IBAR8
            IC = (KK-1)*M%IBAR*M%JBAR + (JJ-1)*M%IBAR + II
            IF (ABS(X(IC))<1.0E-14_EB) THEN
               VALUES(II)=0.0_EB
            ELSE
               VALUES(II)=X(IC)
            ENDIF
         ENDDO
         WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, IBAR8)
      ENDDO
   ENDDO
   WRITE(SCARC_LU,*)  '------------------------------------------------',&
                      '---------------------------------------------------'
   WRITE(SCARC_LU,1000) ('I = ',II,II=1,IBAR8)
   WRITE(SCARC_LU,*)  '--------------------------------------------------',&
                      '-------------------------------------------------'
   WRITE(SCARC_LU,*)
ENDIF

1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4) 
END SUBROUTINE SCARC_SHOW_COMPACT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only for debugging reasons: print out vector information on level NL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW_LEVEL_BANDWISE (X, NX, NY, NZ, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :), INTENT(IN) :: X
INTEGER, INTENT(IN):: NX, NY, NZ
REAL (EB):: VALUES(10)
INTEGER :: II, JJ, KK, IBAR8,JBAR8,KBAR8
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
 
IBAR8=MIN(8,NX)
JBAR8=MIN(8,NY)
KBAR8=MIN(8,NZ)

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '===================================================================================='
   WRITE(SCARC_LU,2000) CROUTINE, CNAME
   WRITE(SCARC_LU,*) '===================================================================================='
   IF (TWO_D) THEN
      DO KK = KBAR8, 1, - 1
         DO II=1,IBAR8
            IF (ABS(X(II,1,KK))<1.0E-14_EB) THEN
               VALUES(II)=0.0_EB
            ELSE
               VALUES(II)=X(II,1,KK)
            ENDIF
         ENDDO
         WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',1,' : ',(VALUES(II), II=1, IBAR8)
      ENDDO
      WRITE(SCARC_LU,*)  '------------------------------------------------',&
                         '---------------------------------------------------'
      !WRITE(SCARC_LU,1000) ('I = ',II,II=1,IBAR8)
   ELSE
      DO KK = KBAR8, 1, - 1
         WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                     '----------------------------------------'
         DO JJ = JBAR8, 1, - 1
            DO II=1,IBAR8
               IF (ABS(X(II,JJ,KK))<1.0E-14_EB) THEN
                  VALUES(II)=0.0_EB
               ELSE
                  VALUES(II)=X(II,JJ,KK)
               ENDIF
            ENDDO
            WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, IBAR8)
         ENDDO
      ENDDO
      WRITE(SCARC_LU,*)  '------------------------------------------------',&
                         '---------------------------------------------------'
      WRITE(SCARC_LU,1000) ('I = ',II,II=1,IBAR8)
      WRITE(SCARC_LU,*)  '--------------------------------------------------',&
                         '-------------------------------------------------'
      WRITE(SCARC_LU,*)
   ENDIF
ENDIF

1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4) 
END SUBROUTINE SCARC_SHOW_LEVEL_BANDWISE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only for debugging reasons: print out vector information on level NL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW_LEVEL_COMPACT (X, NX, NY, NZ, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:), INTENT(IN) :: X
INTEGER, INTENT(IN):: NX, NY, NZ
REAL (EB):: VALUES(8)
INTEGER :: II, JJ, KK, IBAR8,JBAR8,KBAR8, IC
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
 

IBAR8=MIN(8,NX)
JBAR8=MIN(8,NY)
KBAR8=MIN(8,NZ)

WRITE(SCARC_LU,*) '===================================================================================='
WRITE(SCARC_LU,2000) CROUTINE, CNAME
WRITE(SCARC_LU,*) '===================================================================================='
IF (TWO_D) THEN
   DO KK = KBAR8, 1, - 1
      DO II=1,IBAR8
         IC = (KK-1)*NX + II
         IF (ABS(X(IC))<1.0E-14_EB) THEN
            VALUES(II)=0.0_EB
         ELSE
            VALUES(II)=X(IC)
         ENDIF
      ENDDO
      WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',1,' : ',(VALUES(II), II=1, IBAR8)
   ENDDO
   WRITE(SCARC_LU,*)  '------------------------------------------------',&
                      '---------------------------------------------------'
   !WRITE(SCARC_LU,1000) ('I = ',II,II=1,IBAR8)
ELSE
   DO KK = KBAR8, 1, - 1
      WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                  '----------------------------------------'
      DO JJ = JBAR8, 1, - 1
         DO II=1,IBAR8
            IC = (KK-1)*NX*NY + (JJ-1)*NX + II
            IF (ABS(X(IC))<1.0E-14_EB) THEN
               VALUES(II)=0.0_EB
            ELSE
               VALUES(II)=X(IC)
            ENDIF
         ENDDO
         WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, IBAR8)
      ENDDO
   ENDDO
   WRITE(SCARC_LU,*)  '------------------------------------------------',&
                      '---------------------------------------------------'
   WRITE(SCARC_LU,1000) ('I = ',II,II=1,IBAR8)
   WRITE(SCARC_LU,*)  '--------------------------------------------------',&
                      '-------------------------------------------------'
   WRITE(SCARC_LU,*)
ENDIF

1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4) 
END SUBROUTINE SCARC_SHOW_LEVEL_COMPACT

 
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
   NAME_SCARC                              = 'null'
   NAME_SCARC(NSCARC_TIME_COMPLETE)        = 'ScaRC complete'
   NAME_SCARC(NSCARC_TIME_SETUP)           = 'ScaRC setup'
   NAME_SCARC(NSCARC_TIME_SOLVER)          = 'ScaRC solver'
   NAME_SCARC(NSCARC_TIME_KRYLOV)          = 'ScaRC krylov method'
   NAME_SCARC(NSCARC_TIME_MULTIGRID)       = 'ScaRC multigrid method'
   NAME_SCARC(NSCARC_TIME_PRECON)          = 'ScaRC preconditioner'
   NAME_SCARC(NSCARC_TIME_SMOOTH)          = 'ScaRC smoother'
   NAME_SCARC(NSCARC_TIME_COARSE)          = 'ScaRC coarse grid solver'
   NAME_SCARC(NSCARC_TIME_LOCAL_MATVEC)    = 'ScaRC local  matrix-vector product'
   NAME_SCARC(NSCARC_TIME_GLOBAL_MATVEC)   = 'ScaRC global matrix-vector product'
   NAME_SCARC(NSCARC_TIME_LOCAL_SCALPROD)  = 'ScaRC local  scalar product'
   NAME_SCARC(NSCARC_TIME_GLOBAL_SCALPROD) = 'ScaRC global scalar product'
   NAME_SCARC(NSCARC_TIME_GLOBAL_L2NORM)   = 'ScaRC global L2-norm'
   NAME_SCARC(NSCARC_TIME_EXCHANGE_BDRY)   = 'ScaRC exchange of internal boundaries'
   
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Stops the code gracefully after writing a message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
 


