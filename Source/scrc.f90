#undef WITH_MKL_FB
#ifdef WITH_MKL
#define __MKL_PARDISO_F90
#define __MKL_CLUSTER_SPARSE_SOLVER_F90
!include "include/mkl_pardiso.f90"
!include "include/mkl_cluster_sparse_solver.f90"
#endif

MODULE SCRC

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME, GET_FILE_NUMBER, SHUTDOWN
USE MPI
USE TYPES, ONLY: MULTIPLIER_TYPE

#ifdef WITH_MKL
  USE MKL_PARDISO
  USE MKL_CLUSTER_SPARSE_SOLVER
#endif 

IMPLICIT NONE

!> ------------------------------------------------------------------------------------------------
!> Public subroutines (initialization, solver, time measurement and revisioning)
!> ------------------------------------------------------------------------------------------------
PUBLIC SCARC_SETUP                             !> Setup routine for ScaRC, needed in main.f90
PUBLIC SCARC_SOLVER                            !> Call of basic ScaRC solver
PUBLIC SCARC_TIMINGS                           !> Call of time measurements for ScaRC

!> ------------------------------------------------------------------------------------------------
!> Public variables   (explanations in declaration part below)
!> Note: For input parameters in character format corresponding INTEGER type-parameters will
!> Be introduced later to simplify inquiries
!> ------------------------------------------------------------------------------------------------
PUBLIC SCARC_METHOD                            !> ScaRC method 
PUBLIC SCARC_DISCRETIZATION                    !> ScaRC discretization type 
PUBLIC SCARC_TWOLEVEL                          !> Type of Twolevel-Krylov-method 

PUBLIC SCARC_DEBUG                             !> Level for plotting debug information 
PUBLIC SCARC_VERBOSE                           !> Level for plotting verbose information 
PUBLIC SCARC_CSV                               !> Level for plotting csv information 

PUBLIC SCARC_RESIDUAL                          !> Residual of iterative solver
PUBLIC SCARC_ITERATIONS                        !> Number of iterations
PUBLIC SCARC_CAPPA                             !> Convergence rate
PUBLIC SCARC_ACCURACY                          !> Chosen accuracy type (relative/absolute)
PUBLIC SCARC_PRECISION                         !> Single/double precision for preconditioning or LU-decomposition

PUBLIC SCARC_KRYLOV                            !> Type of Krylov method
PUBLIC SCARC_KRYLOV_ITERATIONS                 !> Maximum number of iterations for Krylov method
PUBLIC SCARC_KRYLOV_INTERPOL                   !> Interpolation type for twolevel methods
PUBLIC SCARC_KRYLOV_ACCURACY                   !> Requested accuracy for Krylov method

PUBLIC SCARC_MULTIGRID                         !> Type of multigrid method
PUBLIC SCARC_MULTIGRID_LEVEL                   !> Multigrid level
PUBLIC SCARC_MULTIGRID_ITERATIONS              !> Maximum number of iterations for multigrid method
PUBLIC SCARC_MULTIGRID_INTERPOL                !> Interpolation method for multigrid (AMG only)
PUBLIC SCARC_MULTIGRID_ACCURACY                !> Requested accuracy for multigrid method
PUBLIC SCARC_MULTIGRID_CYCLE                   !> Type of multigrid cycle 
PUBLIC SCARC_MULTIGRID_COARSENING              !> Coarsening method for multigrid (AMG only)

PUBLIC SCARC_SMOOTH                            !> Smoother for multigrid method
PUBLIC SCARC_SMOOTH_ITERATIONS                 !> Maximum number of iterations for smoothing method
PUBLIC SCARC_SMOOTH_ACCURACY                   !> Requested accuracy for smoothing method
PUBLIC SCARC_SMOOTH_OMEGA                      !> Damping parameter for smoothing method

PUBLIC SCARC_PRECON                            !> Preconditioner of defect correction method
PUBLIC SCARC_PRECON_ITERATIONS                 !> Maximum number of iterations for perconditioning method
PUBLIC SCARC_PRECON_ACCURACY                   !> Requested accuracy for preconditioning method
PUBLIC SCARC_PRECON_OMEGA                      !> Relaxation parameter for preconditioning method

PUBLIC SCARC_COARSE                            !> Coarse grid solver for multigrid method
PUBLIC SCARC_COARSE_ITERATIONS                 !> Maximum number of iterations for coarse grid solver
PUBLIC SCARC_COARSE_ACCURACY                   !> Requested accuracy for coarse grid solver
PUBLIC SCARC_COARSE_OMEGA                      !> Relaxation parameter for coarse grid solver 
PUBLIC SCARC_COARSE_LEVEL                      !> Coarse grid level


!> ------------------------------------------------------------------------------------------------
!> Miscellaneous declarations 
!> ------------------------------------------------------------------------------------------------
!> General definitions
CHARACTER(40) :: SCARC_METHOD   = 'NONE'                    !> Requested solver method (KRYLOV/MULTIGRID)
CHARACTER(40) :: SCARC_TWOLEVEL = 'NONE'                    !> Type of two-level method (NONE/ADDITIVE/MULTIPLICATIVE)
CHARACTER(40) :: SCARC_DISCRETIZATION   = 'STRUCTURED'      !> Type of discretization (STRUCTURED/UNSTRUCTURED)

!> General iteration parameters
INTEGER       :: SCARC_ITERATIONS           =  0            !> Number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_RESIDUAL             =  0.0_EB       !> Residual of global selected solver
REAL (EB)     :: SCARC_CAPPA                =  0.0_EB       !> Convergence rate of selected ScarC solver
CHARACTER(40) :: SCARC_ACCURACY             = 'ABSOLUTE'    !> Accuracy type (ABSOLUTE/RELATIVE)
CHARACTER(6)  :: SCARC_PRECISION            = 'DOUBLE'      !> Single/double precision for preconditioning or LU-decomposition

!> Parameters for multigrid-type methods
CHARACTER(40) :: SCARC_MULTIGRID            = 'GEOMETRIC'   !> Type of MG-method (GEOMETRIC/ALGEBRAIC)
CHARACTER(40) :: SCARC_MULTIGRID_COARSENING = 'FALGOUT'     !> Coarsening strategy  (Falgout/RS3/A1/A2/...)
CHARACTER(1)  :: SCARC_MULTIGRID_CYCLE      = 'V'           !> Cycling type  (F/V/W)
INTEGER       :: SCARC_MULTIGRID_LEVEL      = -1            !> User defined number of MG-levels (optionally, otherwise maximum)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS = 100           !> Max number of iterations
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY   = 1.E-8_EB      !> Requested accuracy for convergence
CHARACTER(40) :: SCARC_MULTIGRID_INTERPOL   = 'CONSTANT'    !> Interpolation strategy (CONSTANT/BILINEAR)

!> Parameters for Krylov-type methods
CHARACTER(40) :: SCARC_KRYLOV            = 'CG'             !> Type of Krylov-method (CG/BICG)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 1000             !> Max number of iterations
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-8_EB         !> Requested accuracy for convergence
CHARACTER(40) :: SCARC_KRYLOV_INTERPOL   = 'CONSTANT'       !> twolevel-interpolation (CONSTANT/BILINEAR)

!> Parameters for smoothing method (used in multigrids-methods)
CHARACTER(40) :: SCARC_SMOOTH            = 'SSOR'           !> Smoother for MG (JACOBI/SSOR)
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 5                !> Max number of iterations
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-8_EB         !> Requested accuracy for convergence
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.80E+0_EB       !> Relaxation parameter

!> Parameters for preconditioning method (used in Krylov-methods)
CHARACTER(40) :: SCARC_PRECON            = 'NONE'           !> Preconditioner for CG/BICG (JACOBI/SSOR/FFT/PARDISO/MG)
INTEGER       :: SCARC_PRECON_ITERATIONS = 100              !> Max number of iterations
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-10_EB        !> Requested accuracy for convergence
REAL (EB)     :: SCARC_PRECON_OMEGA      = 0.80E+0_EB       !> Relaxation parameter

!> Parameters for coarse grid method
CHARACTER(40) :: SCARC_COARSE            = 'ITERATIVE'      !> Type of coarse grid solver (ITERATIVE/DIRECT)
INTEGER       :: SCARC_COARSE_ITERATIONS = 100              !> Max number of iterations for iterative solver
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-14_EB        !> Requested accuracy for iterative solver
REAL (EB)     :: SCARC_COARSE_OMEGA      = 0.80E+0_EB       !> Relaxation parameter
INTEGER       :: SCARC_COARSE_LEVEL      =  1               !> Coarse grid level for twolevel-Krylov method (default minimum level)

#ifdef WITH_MKL
!> Parameter for MKL solver
CHARACTER(40) :: SCARC_MKL       = 'GLOBAL'                 !> Type of MKL solver (LOCAL->Pardiso/GLOBAL->Cluster_Sparse_solver)
CHARACTER(40) :: SCARC_MKL_MTYPE = 'SYMMETRIC'              !> Type of MKL matrix (SYMMETRIC/UNSYMMETRIC)
#endif

!> Debugging parameters
CHARACTER(40) :: SCARC_DEBUG   = 'NONE'                     !> Debugging level (NONE/LESS/MUCH)
CHARACTER(40) :: SCARC_VERBOSE = 'NONE'                     !> Verbose level (NONE/LESS/MUCH)
CHARACTER(40) :: SCARC_CSV     = 'NONE'                     !> CSV writing level (NONE/MAIN/MEDIUM/FULL)

PRIVATE

!> ------------------------------------------------------------------------------------------------
!> Global constants
!> ------------------------------------------------------------------------------------------------
INTEGER, PARAMETER :: NSCARC_DISCRET_STRUCTURED      =  1, &    !> structured discretization
                      NSCARC_DISCRET_UNSTRUCTURED    =  2, &    !> unstructured discretization
                      NSCARC_DISCRET_GASPHASE        =  3, &    !> gasphase cell
                      NSCARC_DISCRET_SOLID           =  4       !> solid cell

INTEGER, PARAMETER :: NSCARC_SCOPE_ONE               =  1, &    !> primary scope for solution vectors
                      NSCARC_SCOPE_TWO               =  2       !> secondary scope for solution vectors

INTEGER, PARAMETER :: NSCARC_METHOD_KRYLOV           =  1, &    !> Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID        =  2, &    !> multigrid-method as global solver
                      NSCARC_METHOD_LUDECOMP         =  3       !> MKL method as global solver

INTEGER, PARAMETER :: NSCARC_KRYLOV_CG               =  1, &    !> CG   as Krylov solver
                      NSCARC_KRYLOV_BICG             =  2, &    !> BICG as Krylov solver
                      NSCARC_KRYLOV_MAIN             =  1, &    !> Krylov solver as main solver
                      NSCARC_KRYLOV_COARSE           =  2       !> Krylov solver as coarse grid solver

INTEGER, PARAMETER :: NSCARC_MULTIGRID_GEOMETRIC     =  1, &    !> geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC     =  2, &    !> algebraic multigrid
                      NSCARC_MULTIGRID_MAIN          =  1, &    !> multigrid solver as main solver
                      NSCARC_MULTIGRID_PRECON        =  2       !> multigrid solver as preconditioner

INTEGER, PARAMETER :: NSCARC_MKL_LOCAL               =  1, &    !> local use of MKL solver (Pardiso)
                      NSCARC_MKL_GLOBAL              =  2, &    !> global use of MKL solver (Cluster_Sparse_Solver)
                      NSCARC_MKL_COARSE              =  3       !> MKL solver on coarse grid level

INTEGER, PARAMETER :: NSCARC_EXCHANGE_BASIC          =  1, &    !> initialize wall information
                      NSCARC_EXCHANGE_WALLINFO       =  2, &    !> initialize wall information
                      NSCARC_EXCHANGE_VECTOR         =  3, &    !> matrix-vector communication
                      NSCARC_EXCHANGE_PRESSURE       =  4, &    !> vector values along internal boundaries
                      NSCARC_EXCHANGE_MEASURE        =  5, &    !> measure values along internal boundaries
                      NSCARC_EXCHANGE_MEASURE_ADD    =  6, &    !> measure values along internal boundaries
                      NSCARC_EXCHANGE_CELL_TYPE      =  7, &    !> cell types along internal boundaries
                      NSCARC_EXCHANGE_CELL_INDEX     =  8, &    !> internal transfer weights
                      NSCARC_EXCHANGE_PROLONGATION   =  9, &    !> internal transfer weights
                      NSCARC_EXCHANGE_RESTRICTION    = 10, &    !> internal transfer weights
                      NSCARC_EXCHANGE_STENCIL        = 11, &    !> internal subdiagonal matrix values
                      NSCARC_EXCHANGE_MATRIX_SIZE    = 12, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_SUBDIAG = 13, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_VALUE   = 14, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_STENCIL = 15, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_SYSTEM  = 16, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_PROL    = 17, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_REST    = 18, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_TRANSFER_SIZE  = 19, &    !> neighboring transfer matrix size
                      NSCARC_EXCHANGE_WIDTHINFO      = 20, &    !> neighboring grid resolution
                      NSCARC_EXCHANGE_MESHINFO       = 21, &    !> neighboring mesh information
                      NSCARC_EXCHANGE_GRAPH          = 22, &    !> graph along internal boundaries
                      NSCARC_EXCHANGE_DISCRETIZATION = 23       !> exchange discretization information

INTEGER, PARAMETER :: NSCARC_BROADCAST_SUM          =  1, &    !> broadcast local value and deliver sum of all
                      NSCARC_BROADCAST_PRODUCT      =  2, &    !> broadcast local value and deliver product of all
                      NSCARC_BROADCAST_MEAN         =  3, &    !> broadcast local value and deliver mean of all
                      NSCARC_BROADCAST_FIRST        =  4, &    !> broadcast local value and deliver first 
                      NSCARC_BROADCAST_LAST         =  5       !> broadcast local value and deliver last

INTEGER, PARAMETER :: NSCARC_RELAX_JACOBI           =  1, &    !> preconditioning by JACOBI-method
                      NSCARC_RELAX_SSOR             =  2, &    !> preconditioning by SSOR-method
                      NSCARC_RELAX_FFT              =  3, &    !> preconditioning by FFT-method
                      NSCARC_RELAX_FFT_OVERLAP      =  4, &    !> preconditioning by FFT-method
                      NSCARC_RELAX_MULTIGRID        =  5, &    !> preconditioning by MG-method
                      NSCARC_RELAX_PARDISO          =  6, &    !> preconditioning by PARDISO-method
                      NSCARC_RELAX_CLUSTER          =  7       !> preconditioning by PARDISO-method

INTEGER, PARAMETER :: NSCARC_TWOLEVEL_NONE           =  0, &    !> no two levels, only one level
                      NSCARC_TWOLEVEL_ADD            =  1, &    !> additive 2-level method
                      NSCARC_TWOLEVEL_MUL            =  2, &    !> multiplicative 2-level method
                      NSCARC_TWOLEVEL_MUL2           =  3, &    !> multiplicative 2-level method
                      NSCARC_TWOLEVEL_COARSE         =  4       !> only coarse grid

INTEGER, PARAMETER :: NSCARC_CYCLING_F               =  0, &    !> F-cycle for mg-method
                      NSCARC_CYCLING_V               =  1, &    !> V-cycle for mg-method
                      NSCARC_CYCLING_W               =  2, &    !> W-cycle for mg-method
                      NSCARC_CYCLING_SETUP           =  3, &    !> initialize cycle counts
                      NSCARC_CYCLING_RESET           =  4, &    !> reset cycle counts
                      NSCARC_CYCLING_PROCEED         =  5, &    !> proceed cycle counts
                      NSCARC_CYCLING_NEXT            =  6, &    !> perform next cycling loop
                      NSCARC_CYCLING_EXIT            =  7, &    !> exit cycling loop
                      NSCARC_CYCLING_PRESMOOTH       = -1, &    !> presmoothing cycle
                      NSCARC_CYCLING_POSTSMOOTH      =  1       !> postsmoothing cycle

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED           =  0, &    !> proceed loop
                      NSCARC_STATE_CONV0             =  1, &    !> check convergence already for initial residual
                      NSCARC_STATE_CONV              =  2, &    !> check convergence for residual
                      NSCARC_STATE_DIVG              =  3       !> check divergence for residual

INTEGER, PARAMETER :: NSCARC_CSV_NONE                =  0, &    !> don't print csv-information about iteration parameters
                      NSCARC_CSV_MAIN                =  1, &    !> print information only related to main solvers
                      NSCARC_CSV_MEDIUM              =  2, &    !> print information also related to relaxation solvers
                      NSCARC_CSV_FULL                =  3       !> print all iteration information

INTEGER, PARAMETER :: NSCARC_VERBOSE_NONE            =  0, &    !> info0  level of debugging requested
                      NSCARC_VERBOSE_LESS            =  1, &    !> info0  level of debugging requested
                      NSCARC_VERBOSE_MUCH            =  2       !> info1  level of debugging requested

INTEGER, PARAMETER :: NSCARC_DEBUG_NONE              =  0, &    !> info0  level of debugging requested
                      NSCARC_DEBUG_LESS              =  1, &    !> low    level of debugging requested
                      NSCARC_DEBUG_MUCH              =  2, &    !> strong level of debugging requested
                      NSCARC_DEBUG_MATRIX            = 11, &    !> show matrix
                      NSCARC_DEBUG_MATRIXS           = 12, &    !> show matrix
                      NSCARC_DEBUG_MATRIXE           = 13, &    !> show matrix
                      NSCARC_DEBUG_IJKW              = 14, &    !> show IJKW
                      NSCARC_DEBUG_WALLINFO          = 15, &    !> show WALLINFO
                      NSCARC_DEBUG_FACEINFO          = 16, &    !> show FACEINFO
                      NSCARC_DEBUG_BCINDEX           = 17, &    !> show PRESSURE_BC_INDEX
                      NSCARC_DEBUG_ACELL             = 18, &    !> show WALL_CELL
                      NSCARC_DEBUG_GCELL             = 19, &    !> show GHOST_CELL
                      NSCARC_DEBUG_NCELL             = 20, &    !> show NOM_CELL
                      NSCARC_DEBUG_SUBDIVISION       = 21, &    !> show SUBDIVISION
                      NSCARC_DEBUG_MEASURE           = 22, &    !> show MEASURE
                      NSCARC_DEBUG_CELL_TYPE         = 23, &    !> show CELL_TYPE
                      NSCARC_DEBUG_GRAPH             = 24, &    !> show CELL_TYPE
                      NSCARC_DEBUG_COARSE            = 25, &    !> show coarse grid
                      NSCARC_DEBUG_PROLONGATION      = 26, &    !> show prolongation matrix
                      NSCARC_DEBUG_RESTRICTION       = 27, &    !> show restriction matrix
                      NSCARC_DEBUG_STACK             = 28       !> show restriction matrix

INTEGER, PARAMETER :: NSCARC_COARSENING_BASIC        =  1, &    !> basic coarsening
                      NSCARC_COARSENING_FALGOUT      =  2, &    !> parallel Falgout
                      NSCARC_COARSENING_RS3          =  3, &    !> parallel RS3
                      NSCARC_COARSENING_A1           =  4, &    !> aggressive 1 (path =1, length =2)
                      NSCARC_COARSENING_A2           =  5, &    !> aggressive 2 (path =2, length =2)
                      NSCARC_COARSENING_BDRY         =  6       !> FDSA2  : FDS variant similar to A2

INTEGER, PARAMETER :: NSCARC_COARSE_ITERATIVE        =  1, &    !> iterative solution of coarse grid problem
                      NSCARC_COARSE_DIRECT           =  2, &    !> direct solution of coarse grid problem
                      NSCARC_COARSE_KRYLOV           =  3, &    !> direct solution of coarse grid problem
                      NSCARC_COARSE_PARDISO          =  4, &    !> direct solution of coarse grid problem
                      NSCARC_COARSE_CLUSTER          =  5       !> direct solution of coarse grid problem

INTEGER, PARAMETER :: NSCARC_SIZE_MATRIX             =  1, &    !> size of system matrix for compact system
                      NSCARC_SIZE_TRANSFER           =  2       !> size of transfer matrices for compact system

INTEGER, PARAMETER :: NSCARC_SOLVER_MAIN             =  1, &    !> Krylov solver as main solver
                      NSCARC_SOLVER_PRECON           =  2, &    !> Multigrid solver as main solver
                      NSCARC_SOLVER_SMOOTH           =  3, &    !> Cluster sparse solver as main solver
                      NSCARC_SOLVER_COARSE           =  4       !> Cluster sparse solver as main solver

INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_X            =  1, &    !> selection parameter for 1-scope vector X
                      NSCARC_VECTOR_ONE_F            =  2, &    !> selection parameter for 1-scope vector F
                      NSCARC_VECTOR_ONE_D            =  3, &    !> selection parameter for 1-scope vector D
                      NSCARC_VECTOR_ONE_G            =  4, &    !> selection parameter for 1-scope vector G
                      NSCARC_VECTOR_ONE_W            =  5, &    !> selection parameter for 1-scope vector R
                      NSCARC_VECTOR_ONE_Y            =  6, &    !> selection parameter for 1-scope vector Y
                      NSCARC_VECTOR_ONE_Z            =  7, &    !> selection parameter for 1-scope vector Z
                      NSCARC_VECTOR_ONE_E            = 18, &    !> selection parameter for 1-scope vector Z
                      NSCARC_VECTOR_TWO_X            =  8, &    !> selection parameter for 2-scope vector X
                      NSCARC_VECTOR_TWO_F            =  9, &    !> selection parameter for 2-scope vector F
                      NSCARC_VECTOR_TWO_D            = 10, &    !> selection parameter for 2-scope vector D
                      NSCARC_VECTOR_TWO_G            = 11, &    !> selection parameter for 2-scope vector G
                      NSCARC_VECTOR_TWO_W            = 12, &    !> selection parameter for 2-scope vector R
                      NSCARC_VECTOR_TWO_Y            = 13, &    !> selection parameter for 2-scope vector Y
                      NSCARC_VECTOR_TWO_Z            = 14, &    !> selection parameter for 2-scope vector Z
                      NSCARC_VECTOR_TWO_E            = 19, &    !> selection parameter for 2-scope vector Z
                      NSCARC_VECTOR_H                = 15, &    !> selection parameter for vector H
                      NSCARC_VECTOR_HS               = 16, &    !> selection parameter for vector HS
                      NSCARC_VECTOR_MEASURE          = 17       !> selection parameter for vector Z

INTEGER, PARAMETER :: NSCARC_MATRIX_SYSTEM           =  1, &    !> exchange subdiagonal  matrix entries
                      NSCARC_MATRIX_SUBDIAG          =  2, &    !> exchange subdiagonal  matrix entries
                      NSCARC_MATRIX_SUBDIAG_LOWER    =  3, &    !> exchange subdiagonal  matrix entries
                      NSCARC_MATRIX_SUBDIAG_UPPER    =  4, &    !> exchange subdiagonal  matrix entries
                      NSCARC_MATRIX_TRANSFER         =  5, &    !> exchange prolongation matrix
                      NSCARC_MATRIX_STENCIL          =  6, &    !> exchange prolongation matrix
                      NSCARC_MATRIX_PROLONGATION     =  7, &    !> exchange prolongation matrix
                      NSCARC_MATRIX_RESTRICTION      =  8       !> exchange restriction  matrix

INTEGER, PARAMETER :: NSCARC_ACCURACY_ABSOLUTE       =  1, &    !> absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE       =  2       !> relative accuracy must be reached

INTEGER, PARAMETER :: NSCARC_PRECISION_FB            =  1, &    !> single precision for preconditioning or LU-decomposition
                      NSCARC_PRECISION_EB            =  2       !> double precision for preconditioning or LU-decomposition

REAL(EB), PARAMETER:: NSCARC_MEASURE_NONE            =  0.0_EB, &
                      NSCARC_MEASURE_ONE             =  1.0_EB, &  !> coarse-grid cell
                      NSCARC_MEASURE_COARSE          =  6.0_EB, &  !> coarse-grid cell
                      NSCARC_MEASURE_FINE            =  5.0_EB, &  !> fine-grid cell
                      NSCARC_MEASURE_SFINE           =  5.0_EB, &  !> strongly coupled fine-grid cell
                      NSCARC_MEASURE_WFINE           =  4.0_EB, &  !> weakly   coupled fine-grid cell
                      NSCARC_MEASURE_BDRY            =  1.0_EB     !> boundry weight

INTEGER, PARAMETER :: NSCARC_CELL_TYPE_COARSE        =  1, &    !> coarse-grid cell
                      NSCARC_CELL_TYPE_COARSE0       =  3, &    !> special coarse-grid cell
                      NSCARC_CELL_TYPE_COMMON        =  3, &    !> common cell
                      NSCARC_CELL_TYPE_FINE          = -1, &    !> fine-grid cell
                      NSCARC_CELL_TYPE_FINE0         = -3, &    !> special fine-grid cell
                      NSCARC_CELL_TYPE_SFINE         = -1, &    !> strongly coupled fine-grid cell
                      NSCARC_CELL_TYPE_WFINE         = -2, &    !> weakly   coupled fine-grid cell
                      NSCARC_CELL_TYPE_FPNT          = -1, &    !> special f-point
                      NSCARC_CELL_TYPE_ZPNT          = -2, &    !> special z-point
                      NSCARC_CELL_TYPE_SFPNT         = -3, &    !> special sf-point
                      NSCARC_CELL_TYPE_CPNT          =  2       !> special c-point

INTEGER, PARAMETER :: NSCARC_INTERPOL_STANDARD       =  1, &    !> standard interpolation
                      NSCARC_INTERPOL_CONSTANT       =  2, &    !> standard interpolation
                      NSCARC_INTERPOL_BILINEAR       =  3, &    !> standard interpolation
                      NSCARC_INTERPOL_CLASSICAL      =  4, &    !> classical interpolation
                      NSCARC_INTERPOL_CLASSICAL2     =  5, &    !> classical interpolation
                      NSCARC_INTERPOL_DIRECT         =  6, &    !> direct interpolation
                      NSCARC_INTERPOL_DIRECT_BDRY    =  7       !> direct interpolation with special boundary

INTEGER, PARAMETER :: NSCARC_LEVEL_MIN               =  0, &    !> minimum multigrid level
                      NSCARC_LEVEL_MAX               = 15, &    !> maximum multigrid level
                      NSCARC_LEVEL_SINGLE            =  1, &    !> only one grid level needed
                      NSCARC_LEVEL_MULTI             =  2, &    !> multi grid levels needed
                      NSCARC_LEVEL_AMG               =  3       !> maximum number of grid levels 

INTEGER, PARAMETER :: NSCARC_DUMP_RHS                =  1, &    !> dump rhs
                      NSCARC_DUMP_PRES               =  2       !> dump pressure

INTEGER, PARAMETER :: NSCARC_STENCIL_CENTRAL         =  1, &    !> standard 5- or 7-point stencil
                      NSCARC_STENCIL_AMG             =  2       !> arbitrary AMG-stencil

INTEGER, PARAMETER :: NSCARC_DIAG_MAIN               =  1, &    !> standard 5- or 7-point stencil
                      NSCARC_DIAG_LOWER              =  2, &    !> standard 5- or 7-point stencil
                      NSCARC_DIAG_UPPER              =  3       !> standard 5- or 7-point stencil

INTEGER, PARAMETER :: NSCARC_COUPLING_MAX            = 10       !> maximum of possible couplings in stencil

INTEGER, PARAMETER :: NSCARC_N_FACES                 =  6       !> number of faces per mesh
INTEGER, PARAMETER :: NSCARC_MAX_FACE_NEIGHBORS      = 10       !> max number neighbors per mesh face
INTEGER, PARAMETER :: NSCARC_MAX_STENCIL             =  7

INTEGER, PARAMETER :: NSCARC_UNDEFINED_INT           = -1, &         !> undefined integer value
                      NSCARC_ZERO_INT                =  0, &         !> zero integer value
                      NSCARC_ONE_INT                 =  1            !> one integer value

REAL(EB), PARAMETER:: NSCARC_UNDEFINED_REAL_EB       = -1.0_EB, &    !> undefined real value
                      NSCARC_ZERO_REAL_EB            =  0.0_EB, &    !> zero real value
                      NSCARC_ONE_REAL_EB             =  1.0_EB       !> one real value

REAL(EB), PARAMETER:: NSCARC_UNDEFINED_REAL_FB       = -1.0_FB, &    !> undefined real value
                      NSCARC_ZERO_REAL_FB            =  0.0_FB, &    !> zero real value
                      NSCARC_ONE_REAL_FB             =  1.0_FB       !> one real value

INTEGER, PARAMETER :: NSCARC_INIT_UNDEFINED          = -99, &        !> initialize allocated array as undefined
                      NSCARC_INIT_NONE               =  -1, &        !> do not initialize allocated arrays
                      NSCARC_INIT_ZERO               =   0, &        !> initialize allocated array with zero
                      NSCARC_INIT_ONE                =   1           !> initialize allocated array with one

INTEGER, PARAMETER :: NSCARC_HUGE_INT                = -999999999    !> undefined integer value
REAL(EB), PARAMETER:: NSCARC_HUGE_REAL               = -999999999_EB !> undefined integer value

INTEGER, PARAMETER :: NSCARC_STACK_ZERO              =  0, &         !> root stage of stack
                      NSCARC_STACK_ROOT              =  1, &         !> maximum number of consecutive solvers
                      NSCARC_STACK_MAX               = 10            !> maximum number of consecutive solvers

REAL(EB), PARAMETER:: NSCARC_THRESHOLD_CONVERGENCE   = 1.0E-15_EB, & !> threshold for convergence
                      NSCARC_THRESHOLD_DIVGERGENCE   = 1.0E+15_EB    !> threshold for divergence

INTEGER :: IERROR = 0

!> --------------------------------------------------------------------------------------------
!> Some auxiliary parameters used for transfer operators
!> --------------------------------------------------------------------------------------------
REAL(EB), PARAMETER :: SCALR  = 0.015625_EB
REAL(EB), PARAMETER :: SCALP  = 0.0625_EB
REAL(EB), PARAMETER :: W1     =  1.0_EB
REAL(EB), PARAMETER :: W3     =  3.0_EB
REAL(EB), PARAMETER :: W4     =  4.0_EB
REAL(EB), PARAMETER :: W9     =  9.0_EB
REAL(EB), PARAMETER :: W12    = 12.0_EB
REAL(EB), PARAMETER :: W16    = 16.0_EB

!> --------------------------------------------------------------------------------------------
!> Logical indicators for different methods and mechanisms
!> --------------------------------------------------------------------------------------------
LOGICAL :: BCG         = .FALSE.         ! Krylov-method (different preconditioners possible) used?
LOGICAL :: BCGGMG      = .FALSE.         ! Krylov-method with GMG-preconditioning used?
LOGICAL :: BCGADD      = .FALSE.         ! additive Twolevel-Krylov-method  used?
LOGICAL :: BCGMUL      = .FALSE.         ! multiplicative Twolevel-Krylov-method  used?
LOGICAL :: BCGCOARSE   = .FALSE.         ! only coarse grid preconditiner used ?
LOGICAL :: BMG         = .FALSE.         ! Multigrid-method (different smoothers possible) used?
LOGICAL :: BGMG        = .FALSE.         ! Geometric Multigrid-method used?
LOGICAL :: BAMG        = .FALSE.         ! Algebraic Multigrid-method used?
LOGICAL :: BFFT        = .FALSE.         ! FFT-method used?
LOGICAL :: BMKL        = .FALSE.         ! MKL-method used?
LOGICAL :: BTWOLEVEL   = .FALSE.         ! Method with two grid levels used?
LOGICAL :: BMULTILEVEL = .FALSE.         ! Method with multiple grid levels used?
LOGICAL :: BSTRUCTURED = .TRUE.          ! Structured/Unstructured discretization used?

LOGICAL :: BDEBUG_NONE = .TRUE.          ! no debug messages
LOGICAL :: BDEBUG_LESS = .FALSE.         ! less debug messages
LOGICAL :: BDEBUG_MUCH = .FALSE.         ! much debug messages

LOGICAL :: BVERBOSE_NONE = .TRUE.        ! no verbose messages
LOGICAL :: BVERBOSE_LESS = .FALSE.       ! less verbose messages
LOGICAL :: BVERBOSE_MUCH = .FALSE.       ! much verbose messages

LOGICAL :: BCSV_NONE   = .TRUE.          ! no csv file 
LOGICAL :: BCSV_MAIN   = .FALSE.         ! only main solver related csv information
LOGICAL :: BCSV_MEDIUM = .FALSE.         ! also relaxation solver related csv information
LOGICAL :: BCSV_FULL   = .FALSE.         ! full csv information

!> --------------------------------------------------------------------------------------------
!> Globally used types for description of different solvers
!> --------------------------------------------------------------------------------------------
INTEGER :: TYPE_DISCRET     = NSCARC_DISCRET_STRUCTURED    !> Type of discretization (structured/unstructured)
INTEGER :: TYPE_METHOD      = NSCARC_METHOD_KRYLOV         !> Type of ScaRC method
INTEGER :: TYPE_SOLVER      = NSCARC_SOLVER_MAIN           !> Type of surrounding solver scope
INTEGER :: TYPE_SCOPE       = NSCARC_SCOPE_ONE             !> Type of surrounding solver scope
INTEGER :: TYPE_TWOLEVEL    = NSCARC_TWOLEVEL_NONE         !> Type of two-level method
INTEGER :: TYPE_INTERPOL    = NSCARC_INTERPOL_CONSTANT     !> Type of interpolation method
INTEGER :: TYPE_KRYLOV      = NSCARC_KRYLOV_CG             !> Type of Krylov method (CG/BICG)
INTEGER :: TYPE_MULTIGRID   = NSCARC_MULTIGRID_GEOMETRIC   !> Type of multigrid method (GMG/AMG)
INTEGER :: TYPE_ACCURACY    = NSCARC_ACCURACY_ABSOLUTE     !> Type of requested accuracy
INTEGER :: TYPE_CYCLING     = NSCARC_CYCLING_V             !> Type of cycling for multigrid method
INTEGER :: TYPE_COARSE      = NSCARC_COARSE_ITERATIVE      !> Type of coarse grid solver for multigrid method
INTEGER :: TYPE_LUDECOMP    = NSCARC_UNDEFINED_INT         !> Type of MKL method (PARDISO/CLUSTER_SPARSE_SOLVER)
INTEGER :: TYPE_COARSENING  = NSCARC_UNDEFINED_INT         !> Type of coarsening algorithm for AMG
INTEGER :: TYPE_RELAX       = NSCARC_UNDEFINED_INT         !> Type of preconditioner for iterative solver
INTEGER :: TYPE_PRECON      = NSCARC_UNDEFINED_INT         !> Type of preconditioner for iterative solver
INTEGER :: TYPE_PRECISION   = NSCARC_UNDEFINED_INT         !> Type of preconditioner for iterative solver
INTEGER :: TYPE_SMOOTH      = NSCARC_UNDEFINED_INT         !> Type of smoother for multigrid method
INTEGER :: TYPE_DEBUG       = NSCARC_UNDEFINED_INT         !> Type of debugging level
INTEGER :: TYPE_VERBOSE     = NSCARC_UNDEFINED_INT         !> Type of verbose level
INTEGER :: TYPE_EXCHANGE    = NSCARC_UNDEFINED_INT         !> Type of data exchange
INTEGER :: TYPE_VECTOR      = NSCARC_UNDEFINED_INT         !> Type of vector to point to
INTEGER :: TYPE_PARENT      = NSCARC_UNDEFINED_INT         !> Type of parent (calling) solver
INTEGER :: TYPE_LU_LEVEL(NSCARC_LEVEL_MAX) = NSCARC_UNDEFINED_INT

!> ------------------------------------------------------------------------------------------------
!> Globally used parameters
!> ------------------------------------------------------------------------------------------------
!> discretization information
INTEGER :: NLEVEL_MAX, NLEVEL_MIN                         !> Total, minimum and maximum number of multigrid levels
INTEGER :: N_CELLS_GLOBAL(NSCARC_LEVEL_MAX)     = 0       !> number of global cells 
INTEGER :: N_CELLS_GLOBAL_DOF(NSCARC_LEVEL_MAX) = 0       !> number of degrees of freedom (unstructured case)
INTEGER :: N_DIRIC_GLOBAL(NSCARC_LEVEL_MAX) = 0           !> global number of Dirichlet BCs

!> stack information
INTEGER :: N_STACK_TOTAL                                  !> maximum number of used solvers in stack

!> communication parameters
INTEGER :: N_REQ, N_EXCHANGE, TAG                         !> Variables for data exchange
INTEGER :: SNODE, RNODE                                   !> Process identifier for data exchange

INTEGER,  ALLOCATABLE, DIMENSION (:)  :: REQ              !> Request array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: COUNTS           !> Counter array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: DISPLS           !> Displacement array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: LOCAL_INT        !> Local integer data array for data exchange
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: LOCAL_REAL       !> Local real data array for data exchange

INTEGER :: GLOBAL_INT
REAL(EB):: GLOBAL_REAL

INTEGER  :: FACE_ORIENTATION(6) = (/1,-1,2,-2,3,-3/)           

!> -----------------------------------------------------------------------------------------------
!> Iteration parameters vectors and iteration parameters
!> -----------------------------------------------------------------------------------------------
INTEGER   :: X, F, D, G, W, Y, Z, E
INTEGER   :: NIT, ITE
REAL (EB) :: EPS, RES, RESIN, OMEGA, CAPPA, ERR
INTEGER   :: ITE_PRES=0, ITE_TOTAL=0, ITE_CG=0, ITE_MG=0, ITE_LU=0, ITE_SMOOTH=0, ITE_COARSE=0
CHARACTER(60) :: CNAME

!> -----------------------------------------------------------------------------------------------
!> Time measurement
!> ------------------------------------------------------------------------------------------------
TYPE SCARC_TIME_TYPE
REAL(EB) :: OVERALL   = 0.0_EB               !> complete time
REAL(EB) :: SOLVER    = 0.0_EB               !> time for solver (general version)
REAL(EB) :: CLUSTER   = 0.0_EB               !> time for cluster solver
REAL(EB) :: PARDISO   = 0.0_EB               !> time for pardiso solver
REAL(EB) :: KRYLOV    = 0.0_EB               !> time for krylov solver
REAL(EB) :: MULTIGRID = 0.0_EB               !> time for multigrid solver
REAL(EB) :: MATVEC    = 0.0_EB               !> time for matrix vector multiplication
REAL(EB) :: SCALPROD  = 0.0_EB               !> time for scalar product 
REAL(EB) :: L2NORM    = 0.0_EB               !> time for l2-norm
REAL(EB) :: PRECON    = 0.0_EB               !> time for preconditioner
REAL(EB) :: SMOOTH    = 0.0_EB               !> time for smoother
REAL(EB) :: COARSE    = 0.0_EB               !> time for coarse grid solver
REAL(EB) :: EXCHANGE  = 0.0_EB               !> time for data exchange
END TYPE SCARC_TIME_TYPE


!> -----------------------------------------------------------------------------------------------
!> Messaging and debugging mechanisms
!> ------------------------------------------------------------------------------------------------
TYPE SCARC_MESSAGE_TYPE
CHARACTER(60), ALLOCATABLE, DIMENSION(:) :: HISTORY
CHARACTER(100):: TEXT
CHARACTER(60) :: FILE_DEBUG, FILE_CSV, FILE_DUMP
INTEGER :: LU_DEBUG = 0, LU_CSV = 0, LU_DUMP = 0
INTEGER :: NCURRENT, NHISTORY_MAX = 20
END TYPE SCARC_MESSAGE_TYPE

!> ------------------------------------------------------------------------------------------------
!> Data exchange type
!> ------------------------------------------------------------------------------------------------
TYPE SCARC_EXCHANGE_TYPE
REAL (EB), ALLOCATABLE, DIMENSION (:) :: SEND_REAL, RECV_REAL    !> main real send and receive buffers
INTEGER  , ALLOCATABLE, DIMENSION (:) :: SEND_INT, RECV_INT      !> main integer send and receive buffers
REAL(EB) :: SEND_REAL_BASIC(50), RECV_REAL_BASIC(50)             !> initial real send and receive buffers
INTEGER  :: SEND_INT_BASIC(50) , RECV_INT_BASIC(50)              !> initial integer send and receive buffers
INTEGER  :: NICMAX_R=0, NICMAX_S=0, NIC_R=0, NIC_S=0
END TYPE SCARC_EXCHANGE_TYPE

!> --------------------------------------------------------------------------------------------
!> Face information related to wall cells and neighbors
!> --------------------------------------------------------------------------------------------
TYPE SCARC_FACE_TYPE
INTEGER :: NFC, NFW                                !> number of cells and wall cells along face
INTEGER :: NFX, NFY, NFZ                           !> local face dimensions
INTEGER :: NCPL                                    !> number of adjacent couplings
INTEGER :: N_NEIGHBORS = 0                         !> number of adjacent neighbors
INTEGER :: IWG_PTR                                 !> first (global) IW number to that face
INTEGER :: IOFFSET_WALL   = 0                      !> counter for wall cells over all faces
INTEGER , ALLOCATABLE, DIMENSION(:) :: NEIGHBORS   !> adjacent neighbors
REAL(EB), POINTER, DIMENSION(:) :: DH              !> adjacent grid sizes
END TYPE SCARC_FACE_TYPE

!> --------------------------------------------------------------------------------------------
!> Wall information related to neighbors and BC's
!> --------------------------------------------------------------------------------------------
TYPE SCARC_WALL_TYPE
INTEGER :: DOF, STATE                             !> Degree of freedom and state of related cell (gasphase/solid)
INTEGER :: BTYPE                                  !> type of wall cell (Dirichlet/Neumann/Internal)
INTEGER :: BOUNDARY_TYPE = 0                      !> state of wall cell (Solid/Interpolated/Open))
INTEGER :: IOR = 0                                !> orientation of wall cell
INTEGER :: IXG, IYG, IZG                          !> x-, y- and z-indices of ghost cells
INTEGER :: IXW, IYW, IZW                          !> x-, y- and z-indices of (internal) wall cells
INTEGER :: IXN(2), IYN(2), IZN(2)                 !> x-, y- and z-indices of neighboring cells
INTEGER :: NCPL = 1                               !> number of couplings at wall cell (depending on resolution of neighbor)
INTEGER :: NOM = 0                                !> adjacent neighbor at wall cell
INTEGER :: ICW = NSCARC_UNDEFINED_INT             !> internal wall cell for IW
INTEGER :: ICO = NSCARC_UNDEFINED_INT             !> overlapping cell for IW
INTEGER :: IWL = NSCARC_UNDEFINED_INT             !> corresponding local wall cell number for neighbor NOM
INTEGER, ALLOCATABLE, DIMENSION(:) :: ICE         !> extended cell for IW
INTEGER, ALLOCATABLE, DIMENSION(:) :: ICG         !> ghost cell for IW
END TYPE SCARC_WALL_TYPE

!> --------------------------------------------------------------------------------------------
!> Mappings between different discretization description arrays
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MAPPING_TYPE
INTEGER :: ICG_PTR = 0                                     !> ghost cell pointer
INTEGER :: ICO_PTR = 0                                     !> overlapping cell pointer
INTEGER :: ICE_PTR = 0                                     !> extended cell pointer
INTEGER :: IWL_PTR = 0                                     !> local wall cell pointer
!INTEGER , ALLOCATABLE, DIMENSION (:)  :: ICN_TO_ICE        !> mapping from ICN to ICE    ! ACHTUNG: RIESIG
INTEGER , ALLOCATABLE, DIMENSION (:)  :: ICE_TO_IWG        !> mapping from ICE to IWG
INTEGER , ALLOCATABLE, DIMENSION (:)  :: ICE_TO_IWL        !> mapping from ICE to IWL
INTEGER , ALLOCATABLE, DIMENSION (:)  :: ICE_TO_ICG        !> mapping from ICE to ICG
INTEGER , ALLOCATABLE, DIMENSION (:)  :: ICE_TO_ICN        !> mapping from ICE to ICN
INTEGER , ALLOCATABLE, DIMENSION (:)  :: ICG_TO_IWG        !> mapping from ICG to IWG
INTEGER , ALLOCATABLE, DIMENSION (:)  :: ICG_TO_ICE        !> mapping from ICG to ICE
INTEGER , ALLOCATABLE, DIMENSION (:)  :: ICG_TO_ICO        !> mapping from ICG to ICE
INTEGER , ALLOCATABLE, DIMENSION (:)  :: IWL_TO_IWG        !> mapping from IWL to IWG
INTEGER , ALLOCATABLE, DIMENSION (:)  :: IWL_TO_ICW        !> mapping from IWL to ICW
INTEGER , ALLOCATABLE, DIMENSION (:)  :: IWL_TO_ICO        !> mapping from IWL to ICO
INTEGER , ALLOCATABLE, DIMENSION (:,:):: IWL_TO_ICG        !> mapping from IWL to ICG
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: ICE_TO_VAL        !> mapping from ICE to VAL
END TYPE SCARC_MAPPING_TYPE

!> --------------------------------------------------------------------------------------------
!> Obstruction information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_OBST_TYPE
INTEGER :: I1, I2, J1, J2, K1, K2                              !> cell indices of obstructions
END TYPE SCARC_OBST_TYPE

!> --------------------------------------------------------------------------------------------
!> Preconditioners information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_AMG_TYPE
INTEGER :: NP, NPE, NPS                                        !> number of elements for prolongation matrix
INTEGER :: NR, NRE, NRS                                        !> number of elements for restriction matrix
INTEGER :: NCF, NCFE, NCFS                                     !> number of internal and extended fine cells 
INTEGER :: NCC, NCCE, NCCS                                     !> number of internal and extended coarse cells 
INTEGER :: NCE, NCCI                                           !> number of internal coarse cells
INTEGER :: ICG0 = 0                                            !> auxiliary counter
REAL(EB):: MAX_ROWSUM=0.9_EB                                   !> row sum threshold for AMG
INTEGER,  ALLOCATABLE, DIMENSION (:)   :: CTYPE                !> type for coarsening
INTEGER,  ALLOCATABLE, DIMENSION (:)   :: GRAPH                !> graph vector 
REAL(EB), ALLOCATABLE, DIMENSION (:)   :: MEASURE              !> measure of different cells
INTEGER,  ALLOCATABLE, DIMENSION (:)   :: MARKER               !> cell markers for coarsening
INTEGER,  ALLOCATABLE, DIMENSION (:,:) :: MAPPING              !> mapping for coarsening
INTEGER,  ALLOCATABLE, DIMENSION (:)   :: INTERN               !> index of internal boundary cells
INTEGER,  ALLOCATABLE, DIMENSION (:)   :: EXTERN               !> index of external cells
END TYPE SCARC_AMG_TYPE

!> --------------------------------------------------------------------------------------------
!> Coordinates information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_COORD_TYPE
REAL(EB) :: DX , DY , DZ                                       !> step sizes in x-, y- and z-direction
REAL(EB) :: DXI, DYI, DZI                                      !> inversed of step sizes in x-, y- and z-direction
REAL(EB) :: DXI2, DYI2, DZI2                                   !> squared and inversed step sizes in x-, y- and z-direction
REAL(EB) :: DH = 0.0_EB                                        !> local step sizes
REAL(EB), ALLOCATABLE, DIMENSION (:) :: XCOR, YCOR, ZCOR       !> coordinate vectors in x-, y- and z-direction
REAL(EB), ALLOCATABLE, DIMENSION (:) :: XMID, YMID, ZMID       !> midpoint vectors in x-, y- and z-direction
REAL(EB), ALLOCATABLE, DIMENSION (:) :: DXL, DYL, DZL          !> step size vectors in x-, y- and z-direction
END TYPE SCARC_COORD_TYPE

!> --------------------------------------------------------------------------------------------
!> Wall information related to neighbors and BC's
!> --------------------------------------------------------------------------------------------
TYPE SCARC_CELL_TYPE
INTEGER :: NC_GLOBAL = 0                                        !> number of global cells 
INTEGER :: NC_GLOBAL_DOF = 0                                    !> number of global cells 
INTEGER,  ALLOCATABLE, DIMENSION (:)     :: NC_LOCAL            !> number of local cells, 
INTEGER,  ALLOCATABLE, DIMENSION (:)     :: NC_LOCAL_DOF        !> number of local degrees of freedom (unstructured)
INTEGER,  ALLOCATABLE, DIMENSION (:)     :: NC_OFFSET           !> offset in cell numbering 
INTEGER,  ALLOCATABLE, DIMENSION (:)     :: NC_OFFSET_DOF       !> offset for local degrees of freedom (unstructured)
INTEGER,  ALLOCATABLE, DIMENSION (:,:,:) :: DOF                 !> indicates if cell is degree of freedom
INTEGER,  ALLOCATABLE, DIMENSION (:,:,:) :: STATE               !> state of single cells (gasphase/solid)
INTEGER,  ALLOCATABLE, DIMENSION (:,:,:) :: CINDEX              !> cell index list
INTEGER,  ALLOCATABLE, DIMENSION (:,:)   :: WINDEX              !> wall index list 
END TYPE SCARC_CELL_TYPE

!> --------------------------------------------------------------------------------------------
!> Collection of grid level related information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_LEVEL_TYPE
INTEGER :: NOM                                              !> number of adjacent neighbor
INTEGER :: IOR = 0                                          !> local orientations
INTEGER :: SUBDIVISION(3,-3:3)=0                            !> basic information related to single faces
INTEGER :: N_OBST                                           !> number of obstructions
INTEGER :: N_CELLS                                          !> number of cells
INTEGER :: N_WALL_CELLS                                     !> number of wall cells
INTEGER :: N_WALL_CELLS_EXT                                 !> number of external cells
INTEGER :: N_WALL_CELLS_INT                                 !> number of internal cells
INTEGER :: N_DIRIC = 0                                      !> number of Dirichlet BCs 
INTEGER :: N_NEUMANN = 0                                    !> number of Neumann BCs 
INTEGER :: NX, NY, NZ                                       !> number of grid cells in x-, y- and z-direction
INTEGER :: NC, NCS, NCW                                     !> number of cells
INTEGER :: NW, NWL                                          !> number of global and local wall cells
INTEGER :: NCG, NCE, NCO                                    !> number of ghost, extended, overallping cells
INTEGER :: NCPL=1, NCPL_MAX=-NSCARC_UNDEFINED_INT           !> number of couplings
INTEGER :: NCPLS, NCPLR                                     !> number of couplings to send and read
TYPE (SCARC_CELL_TYPE)    :: CELL                           !> cell information
TYPE (SCARC_COORD_TYPE)   :: COORD                          !> coordinates information
TYPE (SCARC_AMG_TYPE)     :: AMG                            !> coordinates information
TYPE (SCARC_MAPPING_TYPE) :: MAP                            !> mappings between different meshes 
TYPE (SCARC_WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: WALL   !> wall information
TYPE (SCARC_FACE_TYPE), ALLOCATABLE, DIMENSION(:) :: FACE   !> face information
TYPE (SCARC_OBST_TYPE), ALLOCATABLE, DIMENSION(:) :: OBST   !> obstruction information
END TYPE SCARC_LEVEL_TYPE

!> --------------------------------------------------------------------------------------------
!> Matrix entries which will be stored and exchanged during generation of condensed system
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MATRIX_CONDENSING_TYPE                      !> only needed for condensed system
REAL(EB) :: VAL_ORIG(NSCARC_MAX_STENCIL) = 0.0_EB      !> save original values
REAL(EB) :: VAL_COND(NSCARC_MAX_STENCIL) = 0.0_EB      !> save condensed values
INTEGER  :: COL(NSCARC_MAX_STENCIL) = 0                !> save column pointers
INTEGER  :: PTR(NSCARC_MAX_STENCIL) = 0                !> save storage pointer
INTEGER  :: NROW, NCOL
END TYPE SCARC_MATRIX_CONDENSING_TYPE

!> --------------------------------------------------------------------------------------------
!> Matrices including storage pointers and lengths - double precision version
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MATRIX_TYPE
INTEGER :: NAV, NAC, NAR, NAS, NAE                     !> number of elements for different parts of structure
INTEGER :: NSTENCIL                                    !> number of points in matrix stencil
REAL(EB), ALLOCATABLE, DIMENSION (:) :: VAL            !> value of matrix (double precision)
INTEGER,  ALLOCATABLE, DIMENSION (:) :: ROW            !> row pointer 
INTEGER,  ALLOCATABLE, DIMENSION (:) :: COL            !> column pointer
INTEGER,  ALLOCATABLE, DIMENSION (:) :: COLG           !> column pointer for global numbering
INTEGER,  ALLOCATABLE, DIMENSION (:) :: STENCIL        !> matrix stencil information
INTEGER,  ALLOCATABLE, DIMENSION (:) :: POS            !> position in stencil
INTEGER,  ALLOCATABLE, DIMENSION (:) :: TAG            !> marking tag
INTEGER :: NSTORE = 0
TYPE (SCARC_MATRIX_CONDENSING_TYPE) :: STORE(NSCARC_MAX_STENCIL)
END TYPE SCARC_MATRIX_TYPE

#ifdef WITH_MKL_FB
!> --------------------------------------------------------------------------------------------
!> Matrices including storage pointers and lengths - single precision version
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MATRIX_FB_TYPE
INTEGER :: NAV, NAC, NAR, NAS, NAE                        !> number of elements for different parts of structure
INTEGER :: NSTENCIL                                       !> number of points in matrix stencil
REAL (FB), ALLOCATABLE, DIMENSION (:) :: VAL              !> value of matrix (single precision)
INTEGER,   ALLOCATABLE, DIMENSION (:) :: ROW              !> row pointer 
INTEGER,   ALLOCATABLE, DIMENSION (:) :: COL              !> column pointer
INTEGER,   ALLOCATABLE, DIMENSION (:) :: COLG             !> column pointer for global numbering
INTEGER,   ALLOCATABLE, DIMENSION (:) :: STENCIL          !> matrix stencil information
INTEGER,   ALLOCATABLE, DIMENSION (:) :: POS              !> position in stencil
INTEGER,   ALLOCATABLE, DIMENSION (:) :: TAG              !> marking tag
END TYPE SCARC_MATRIX_FB_TYPE
#endif

!> --------------------------------------------------------------------------------------------
!> Collection of different matrices used for system of equations
!> --------------------------------------------------------------------------------------------
TYPE SCARC_SYSTEM_TYPE
TYPE (SCARC_MATRIX_TYPE) :: A                       !> Poisson matrix, Poisson matrix only symmetric part
TYPE (SCARC_MATRIX_TYPE) :: P, R                    !> Prolongation, Restriction matrix
TYPE (SCARC_MATRIX_TYPE) :: S, ST                   !> Strength matrix, Strength matrix transpose
#ifdef WITH_MKL_FB
TYPE (SCARC_MATRIX_FB_TYPE) :: ASYM_FB              !> Poisson matrix, only symmetric part, single precision
#else
TYPE (SCARC_MATRIX_TYPE) :: ASYM                    !> Poisson matrix, only symmetric part, double precision
#endif
END TYPE SCARC_SYSTEM_TYPE

!> --------------------------------------------------------------------------------------------
!> Different scopes for solution, rhs and auxiliary vectors of different solvers
!> --------------------------------------------------------------------------------------------
TYPE SCARC_SCOPE_TYPE
REAL (EB), ALLOCATABLE, DIMENSION (:) :: X          !> solution vector
REAL (EB), ALLOCATABLE, DIMENSION (:) :: F          !> right hand side vector
REAL (EB), ALLOCATABLE, DIMENSION (:) :: D          !> auxiliary vector
REAL (EB), ALLOCATABLE, DIMENSION (:) :: G          !> auxiliary vector
REAL (EB), ALLOCATABLE, DIMENSION (:) :: W          !> auxiliary vector
REAL (EB), ALLOCATABLE, DIMENSION (:) :: Y          !> auxiliary vector
REAL (EB), ALLOCATABLE, DIMENSION (:) :: Z          !> auxiliary vector
REAL (EB), ALLOCATABLE, DIMENSION (:) :: E          !> auxiliary vector
#ifdef WITH_MKL_FB
REAL (FB), ALLOCATABLE, DIMENSION (:) :: X_FB       !> auxiliary vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: F_FB       !> auxiliary vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: D_FB       !> auxiliary vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: G_FB       !> auxiliary vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: W_FB       !> auxiliary vector single precision
#endif
END TYPE SCARC_SCOPE_TYPE

!> --------------------------------------------------------------------------------------------
!> Pointers to iteration vectors
!> --------------------------------------------------------------------------------------------
TYPE SCARC_POINTERS_TYPE
INTEGER  :: X = NSCARC_UNDEFINED_INT                !> reference to local X-vector
INTEGER  :: F = NSCARC_UNDEFINED_INT                !> reference to local F-vector
INTEGER  :: D = NSCARC_UNDEFINED_INT                !> reference to local D-vector
INTEGER  :: G = NSCARC_UNDEFINED_INT                !> reference to local G-vector
INTEGER  :: W = NSCARC_UNDEFINED_INT                !> reference to local W-vector
INTEGER  :: Y = NSCARC_UNDEFINED_INT                !> reference to local Y-vector
INTEGER  :: Z = NSCARC_UNDEFINED_INT                !> reference to local Z-vector
INTEGER  :: E = NSCARC_UNDEFINED_INT                !> reference to local E-vector
#ifdef WITH_MKL_FB
INTEGER  :: X_FB = NSCARC_UNDEFINED_INT             !> reference to local X-vector, single precision
INTEGER  :: F_FB = NSCARC_UNDEFINED_INT             !> reference to local F-vector, single precision
INTEGER  :: D_FB = NSCARC_UNDEFINED_INT             !> reference to local D-vector, single precision
INTEGER  :: G_FB = NSCARC_UNDEFINED_INT             !> reference to local G-vector, single precision
INTEGER  :: W_FB = NSCARC_UNDEFINED_INT             !> reference to local W-vector, single precision
#endif
END TYPE SCARC_POINTERS_TYPE

!> --------------------------------------------------------------------------------------------
!> Store parameter types of different solvers
!> --------------------------------------------------------------------------------------------
TYPE SCARC_TYPES_TYPE
INTEGER :: TYPE_METHOD    = NSCARC_UNDEFINED_INT    !> type of current solver
INTEGER :: TYPE_SOLVER    = NSCARC_UNDEFINED_INT    !> type of current solver
INTEGER :: TYPE_PARENT    = NSCARC_UNDEFINED_INT    !> parent (calling) solver
INTEGER :: TYPE_SCOPE     = NSCARC_UNDEFINED_INT    !> scope for working vectors
INTEGER :: TYPE_NLMIN     = NSCARC_UNDEFINED_INT    !> minimum level for that solver
INTEGER :: TYPE_NLMAX     = NSCARC_UNDEFINED_INT    !> maximum level for that solver
INTEGER :: TYPE_RELAX     = NSCARC_UNDEFINED_INT    !> relaxation method
INTEGER :: TYPE_TWOLEVEL  = NSCARC_UNDEFINED_INT    !> schwarz method?
INTEGER :: TYPE_INTERPOL  = NSCARC_UNDEFINED_INT    !> interpolation type
INTEGER :: TYPE_ACCURACY  = NSCARC_UNDEFINED_INT    !> accuracy requirements 
INTEGER :: TYPE_PRECISION = NSCARC_UNDEFINED_INT    !> precision type for preconditioning or LU-decomposition
INTEGER :: TYPE_CYCLING   = NSCARC_UNDEFINED_INT    !> multigrid cycle
END TYPE SCARC_TYPES_TYPE

!> --------------------------------------------------------------------------------------------
!> Settings of iterations parameters for different solvers
!> --------------------------------------------------------------------------------------------
TYPE SCARC_CONVERGENCE_TYPE
INTEGER  :: NIT   = NSCARC_UNDEFINED_INT            !> maximum iteration number 
INTEGER  :: ITE   = NSCARC_UNDEFINED_INT            !> current iteration number 
REAL(EB) :: EPS   = NSCARC_UNDEFINED_REAL_EB        !> required accuracy
REAL(EB) :: RES   = NSCARC_UNDEFINED_REAL_EB        !> current residual
REAL(EB) :: RESIN = NSCARC_UNDEFINED_REAL_EB        !> initial residual
REAL(EB) :: ERR   = NSCARC_UNDEFINED_REAL_EB        !> initial residual
REAL(EB) :: OMEGA = NSCARC_UNDEFINED_REAL_EB        !> relaxation parameter
REAL(EB) :: CAPPA = NSCARC_UNDEFINED_REAL_EB        !> convergence rate
END TYPE SCARC_CONVERGENCE_TYPE

!> --------------------------------------------------------------------------------------------
!> Cycling information for multigrid methods
!> --------------------------------------------------------------------------------------------
TYPE SCARC_CYCLING_TYPE
INTEGER :: COUNTER(2) = 0                           !> Counter for multigrid cycling
END TYPE SCARC_CYCLING_TYPE

!> --------------------------------------------------------------------------------------------
!> Workspace for FFT preconditioners
!> --------------------------------------------------------------------------------------------
TYPE SCARC_FFT_TYPE
INTEGER  :: LSAVE, LWORK
INTEGER  :: LBC, MBC, NBC
INTEGER  :: ITRN, JTRN, KTRN
INTEGER  :: IBAR0, JBAR0, KBAR0
INTEGER  :: ITRN0, JTRN0, KTRN0
INTEGER  :: IS0=0, IF0=0, JS0=0, JF0=0, KS0=0, KF0=0
REAL(EB) :: XS0, XF0, YS0, YF0, ZS0, ZF0
REAL(EB) :: POIS_PTB = 0.0_EB, XLM = 0.0_EB
REAL(EB), ALLOCATABLE, DIMENSION (:)       :: SAVE1, WORK, HX
REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BXS, BXF, BYS, BYF, BZS, BZF
REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: PRHS
END TYPE SCARC_FFT_TYPE

#ifdef WITH_MKL
!> --------------------------------------------------------------------------------------------
!> MKL information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MKL_TYPE
INTEGER, ALLOCATABLE :: IPARM(:)                      
INTEGER :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR, MSGLVL
INTEGER :: PERM(1)
TYPE(MKL_PARDISO_HANDLE),               ALLOCATABLE :: PT_H(:), PT(:)
TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), ALLOCATABLE :: CT_H(:), CT(:)
END TYPE SCARC_MKL_TYPE
#endif

!> --------------------------------------------------------------------------------------------
!> Basic Information about current solver
!> --------------------------------------------------------------------------------------------
TYPE SCARC_SOLVER_TYPE
CHARACTER(60):: CNAME = 'NONE'                      !> name of current solver
TYPE (SCARC_CONVERGENCE_TYPE):: CONVREQS            !> convergence requirements
TYPE (SCARC_POINTERS_TYPE)   :: POINTERS            !> references to vectors in current scope
TYPE (SCARC_TYPES_TYPE)      :: TYPES               !> current settings for solver
END TYPE SCARC_SOLVER_TYPE

!> --------------------------------------------------------------------------------------------
!> Sample sequence of used solvers in stack
!> --------------------------------------------------------------------------------------------
TYPE SCARC_STACK_TYPE
TYPE (SCARC_SOLVER_TYPE), POINTER :: SOLVER
END TYPE SCARC_STACK_TYPE

!> --------------------------------------------------------------------------------------------
!> Administration other mesh data needed for the coupling of adjacent neighbors
!> --------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE
TYPE (SCARC_EXCHANGE_TYPE) :: EXCHANGE                                !> Variables/arrays needed for data exchange
TYPE (SCARC_LEVEL_TYPE) ,  ALLOCATABLE, DIMENSION(:) :: LEVEL         !> Description of level related information
TYPE (SCARC_SYSTEM_TYPE),  ALLOCATABLE, DIMENSION(:) :: SYSTEM        !> Descriptions of systems of equation
END TYPE OSCARC_TYPE

!> --------------------------------------------------------------------------------------------
!> Basic administration type for ScaRC-method
!> --------------------------------------------------------------------------------------------
TYPE SCARC_TYPE
INTEGER  :: N_NEIGHBORS = 0                                           !> number of adjacent neighbors of whole mesh
INTEGER  :: N_CELLS = 0                                               !> number of cells on that mesh
REAL(EB) :: RHS_END = 0.0_EB
INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS                       !> List of adjacent neighbors of whole mesh
INTEGER  :: IBAR, JBAR, KBAR                                          !> number of cells (corresponding to main prg)
REAL(EB) :: XS, XF, YS, YF, ZS, ZF                                    !> x-, y- and z-bounds (corresponding to main prg)
TYPE (OSCARC_TYPE)       , ALLOCATABLE, DIMENSION(:)   :: OSCARC      !> ScaRC type on other mesh
TYPE (SCARC_LEVEL_TYPE)  , ALLOCATABLE, DIMENSION(:)   :: LEVEL       !> of level related information
TYPE (SCARC_SYSTEM_TYPE) , ALLOCATABLE, DIMENSION(:)   :: SYSTEM      !> system matrices (Poisson/Transfer/AMG)
TYPE (SCARC_CYCLING_TYPE), ALLOCATABLE, DIMENSION(:)   :: CYCLING     !> cycling information for  multigrid (V/W/F)
TYPE (SCARC_FFT_TYPE)    , ALLOCATABLE, DIMENSION(:)   :: FFT         !> FFT preconditioner
TYPE (SCARC_SCOPE_TYPE)  , ALLOCATABLE, DIMENSION(:,:) :: SCOPE       !> different scopes for solution vector
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE)    , ALLOCATABLE, DIMENSION(:)   :: MKL         !> MKL solver
#endif
END TYPE SCARC_TYPE


!> --------------------------------------------------------------------------------------------
!> Basic ScaRC type, different solver types and stack of used solvers
!> --------------------------------------------------------------------------------------------
TYPE (SCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: SCARC
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_KRYLOV, MAIN_MULTIGRID
#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_LUDECOMP
#endif
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: COARSE_KRYLOV
#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: COARSE_CLUSTER, COARSE_PARDISO
#endif
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_JACOBI, PRECON_SSOR, PRECON_FFT, PRECON_MULTIGRID
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: SMOOTH_JACOBI, SMOOTH_SSOR, SMOOTH_FFT
#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_PARDISO, PRECON_CLUSTER
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: SMOOTH_PARDISO, SMOOTH_CLUSTER
#endif
TYPE (SCARC_STACK_TYPE)  , SAVE, DIMENSION(:), ALLOCATABLE :: STACK
TYPE (SCARC_TIME_TYPE)   , SAVE, DIMENSION(:), ALLOCATABLE :: TSETUP, TSUM, TSTEP
TYPE (SCARC_MESSAGE_TYPE), SAVE :: MSG

CONTAINS

!> ------------------------------------------------------------------------------------------------
!> Initialize ScaRC structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP
REAL(EB):: TNOW

TNOW = CURRENT_TIME()

CALL SCARC_SETUP_MESSAGING                                                 !> setup messaging/debugging mechanisms
CALL SCARC_SETUP_TIMING                                                    !> setup time measurment for different parts of code

CALL SCARC_PARSE_INPUT          ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> parse input parameters for ScaRC

CALL SCARC_SETUP_LEVELS         ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup different grid levels
CALL SCARC_SETUP_TYPES          ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup basic ScaRC-types for all used grid levels
CALL SCARC_SETUP_MESHES         ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup mesh information
CALL SCARC_SETUP_DISCRETIZATION ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup discretization information
CALL SCARC_SETUP_INTERFACES     ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup structures related to mesh interfaces
CALL SCARC_SETUP_GLOBALS        ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup global variables on fine level
CALL SCARC_SETUP_WALLS          ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup information along neighboring walls
CALL SCARC_SETUP_EXCHANGE       ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup information for data exchange
CALL SCARC_SETUP_SYSTEM         ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup linear system of equations
CALL SCARC_SETUP_AMG            ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup algebraic coarsening if requested (AMG only)
CALL SCARC_SETUP_METHODS        ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup types and parameters for all needed solvers
CALL SCARC_SETUP_VECTORS        ; IF (STOP_STATUS==SETUP_STOP) RETURN     !> setup vectors for all needed solvers

TSETUP(MYID+1)%OVERALL = TSETUP(MYID+1)%OVERALL + CURRENT_TIME() - TNOW
END SUBROUTINE SCARC_SETUP

!> ------------------------------------------------------------------------------------------------
!> Setup debug file if requested
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TIMING
CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_TIMING')

ALLOCATE (TSETUP(LOWER_MESH_INDEX:UPPER_MESH_INDEX), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_TIMING', 'TSETUP', IERROR)

ALLOCATE (TSTEP(LOWER_MESH_INDEX:UPPER_MESH_INDEX), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_TIMING', 'TSTEP', IERROR)

ALLOCATE (TSUM(LOWER_MESH_INDEX:UPPER_MESH_INDEX), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_TIMING', 'TSUM', IERROR)

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_TIMING

!> ------------------------------------------------------------------------------------------------
!> Setup debug file if requested
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MESSAGING
INTEGER:: NM, LASTID

BDEBUG_NONE = TRIM(SCARC_DEBUG) == 'NONE'
BDEBUG_LESS = TRIM(SCARC_DEBUG) == 'LESS'
BDEBUG_MUCH = TRIM(SCARC_DEBUG) == 'MUCH'
IF (BDEBUG_MUCH) BDEBUG_LESS = .TRUE.

BVERBOSE_NONE = TRIM(SCARC_VERBOSE) == 'NONE'
BVERBOSE_LESS = TRIM(SCARC_VERBOSE) == 'LESS'
BVERBOSE_MUCH = TRIM(SCARC_VERBOSE) == 'MUCH'
IF (BVERBOSE_MUCH) BVERBOSE_LESS = .TRUE.

BCSV_NONE   = TRIM(SCARC_CSV) == 'NONE'
BCSV_MAIN   = TRIM(SCARC_CSV) == 'MAIN'
BCSV_MEDIUM = TRIM(SCARC_CSV) == 'MEDIUM'
BCSV_FULL   = TRIM(SCARC_CSV) == 'FULL'
IF (BCSV_MEDIUM) BCSV_MAIN = .TRUE.
IF (BCSV_FULL) THEN
   BCSV_MAIN   = .TRUE.
   BCSV_MEDIUM = .TRUE.
ENDIF

!> Allocate history stack for called routines
ALLOCATE (MSG%HISTORY(MSG%NHISTORY_MAX), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_MESSAGING', 'MSG%HISTORY', IERROR)

MSG%HISTORY  = 'NONE'
MSG%NCURRENT = 0

!> If requested, open file for debug messages
IF (.NOT.BDEBUG_NONE) THEN
   LASTID = -99999
   DO NM=LOWER_MESH_INDEX, UPPER_MESH_INDEX
      IF (MYID == LASTID) CYCLE
      WRITE (MSG%FILE_DEBUG, '(A,A,i3.3)') TRIM(CHID),'.scarc',MYID+1
      WRITE(*,*) 'Opening ', MSG%LU_DEBUG
      MSG%LU_DEBUG = GET_FILE_NUMBER()
      OPEN (MSG%LU_DEBUG, FILE=MSG%FILE_DEBUG, ACTION = 'readwrite')
      !WRITE(*,*) 'MYID=',MYID,': MSG%LU_DEBUG=', MSG%LU_DEBUG
      LASTID = MYID
   ENDDO
ENDIF

!> If requested, open file for CSV-information about convergence of different solvers
IF (.NOT.BCSV_NONE) THEN
   IF (MYID == 0) THEN
      WRITE (MSG%FILE_CSV, '(A,A)') TRIM(CHID),'_scarc.csv'
      MSG%LU_CSV = GET_FILE_NUMBER()       
      OPEN (MSG%LU_CSV, FILE=MSG%FILE_CSV)
      WRITE(MSG%LU_CSV,*) 'ITE_PRES, ITE_TOTAL, ITE_CG, ITE_MG, ITE_LU, ITE_COARSE, ITE_SMOOTH, TSMOOTH, LEVEL, SOLVER, RES'
   ENDIF
ENDIF

END SUBROUTINE SCARC_SETUP_MESSAGING


!> ------------------------------------------------------------------------------------------------
!> Shutdown ScaRC with error message
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SHUTDOWN(CTEXT, CPARAM, NPARAM)
CHARACTER(*), INTENT(IN) :: CTEXT, CPARAM
INTEGER, INTENT(IN) :: NPARAM

IF (TRIM(CPARAM) /= 'NONE') THEN
   WRITE(MSG%TEXT,'(5A)') TRIM(MSG%HISTORY(MSG%NCURRENT)),' : ', TRIM(CTEXT), ' : ',TRIM(CPARAM)
   IF (MYID == 0) WRITE(LU_ERR,'(/A,A,A,A,A,A)') TRIM(CTEXT),' : ',TRIM(CPARAM),' (CHID: ',TRIM(CHID),')'
ELSE IF (NPARAM /= -999) THEN
   WRITE(MSG%TEXT,'(4A, I8)') TRIM(MSG%HISTORY(MSG%NCURRENT)), ' : ',TRIM(CTEXT), ' : ',NPARAM
   IF (MYID == 0) WRITE(LU_ERR,'(/A,A,I8,A,A,A)') TRIM(CTEXT),' : ',NPARAM,' (CHID: ',TRIM(CHID),')'
ELSE 
   WRITE(MSG%TEXT,'(3A)') TRIM(MSG%HISTORY(MSG%NCURRENT)), ' : ',TRIM(CTEXT)
   IF (MYID == 0) WRITE(LU_ERR,'(/A,A,A,A)') TRIM(CTEXT),' (CHID: ',TRIM(CHID),')'
ENDIF

STOP_STATUS = SETUP_STOP
RETURN

END SUBROUTINE SCARC_SHUTDOWN

!> ----------------------------------------------------------------------------------------------------
!> Determine types of input parameters
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PARSE_INPUT
CHARACTER(80) :: MKL_ERROR_PARDISO, MKL_ERROR_CLUSTER, FFT_ERROR_UNSTRUCTURED

CALL SCARC_ENTER_ROUTINE('SCARC_PARSE_INPUT')

ITERATE_PRESSURE = .TRUE.  ! Although there is no need to do pressure iterations to drive down velocity error
                           ! leave it .TRUE. to write out velocity error diagnostics.

MKL_ERROR_PARDISO      = 'Error: MKL Library compile flag not defined, Pardiso solver not available'
MKL_ERROR_CLUSTER      = 'Error: MKL Library compile flag not defined, Cluster_Sparse_Solver not available'
FFT_ERROR_UNSTRUCTURED = 'Error: FFT-preconditioning only possible for structured discretization '

!> 
!> ------------- set type of discretization
!> 
SELECT CASE (TRIM(SCARC_DISCRETIZATION))
   CASE ('STRUCTURED')
      TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED
   CASE ('UNSTRUCTURED')
      TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED
   CASE DEFAULT
      CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_DISCRETIZATION, -999)
END SELECT
PRES_ON_WHOLE_DOMAIN = (TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED)

!> 
!> ------------ set type of global solver
!> 
SELECT CASE (TRIM(SCARC_METHOD))

   CASE ('KRYLOV')

      TYPE_METHOD = NSCARC_METHOD_KRYLOV

      !> set type of Krylov-method (CG/BICG)
      SELECT CASE (TRIM(SCARC_KRYLOV))
         CASE ('CG')
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
         CASE ('BICG')
            TYPE_KRYLOV = NSCARC_KRYLOV_BICG
         CASE DEFAULT
            CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_KRYLOV, -999)
      END SELECT

      !> set type of two-level method
      SELECT CASE (TRIM(SCARC_TWOLEVEL))
         CASE ('NONE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE
         CASE ('ADDITIVE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_ADD
         CASE ('MULTIPLICATIVE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_MUL
         CASE ('MULTIPLICATIVE2')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_MUL2
         CASE ('COARSE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_COARSE
         CASE DEFAULT
            CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_TWOLEVEL, -999)
      END SELECT

      !> set type of interpolation for two-level Krylov method
      SELECT CASE (TRIM(SCARC_KRYLOV_INTERPOL))
         CASE ('NONE')
            TYPE_INTERPOL = NSCARC_UNDEFINED_INT
         CASE ('CONSTANT')
            TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
         CASE ('BILINEAR')
            TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
         CASE DEFAULT
            CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_KRYLOV_INTERPOL, -999)
      END SELECT

      !> set type of preconditioner (JACOBI/SSOR/MG)
      SELECT CASE (TRIM(SCARC_PRECON))
         CASE ('JACOBI')
            TYPE_PRECON = NSCARC_RELAX_JACOBI
         CASE ('SSOR')
            TYPE_PRECON = NSCARC_RELAX_SSOR
         CASE ('MULTIGRID')
            TYPE_PRECON = NSCARC_RELAX_MULTIGRID
            SELECT CASE (TRIM(SCARC_SMOOTH))
               CASE ('JACOBI')
                  TYPE_SMOOTH = NSCARC_RELAX_JACOBI
               CASE ('SSOR')
                  TYPE_SMOOTH = NSCARC_RELAX_SSOR
               CASE ('FFT')
                  IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) &
                     CALL SCARC_SHUTDOWN(FFT_ERROR_UNSTRUCTURED, 'NONE', -999)
                  TYPE_PRECON = NSCARC_RELAX_FFT
               CASE ('FFT_OVERLAP')
                  IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) &
                     CALL SCARC_SHUTDOWN(FFT_ERROR_UNSTRUCTURED, 'NONE', -999)
                  TYPE_PRECON = NSCARC_RELAX_FFT_OVERLAP
               CASE ('PARDISO')
#ifdef WITH_MKL
                  TYPE_SMOOTH = NSCARC_RELAX_PARDISO
#else
                  CALL SCARC_SHUTDOWN(MKL_ERROR_PARDISO, 'NONE', -999)
#endif

               CASE ('CLUSTER')
#ifdef WITH_MKL
                  TYPE_SMOOTH = NSCARC_RELAX_CLUSTER
#else
                  CALL SCARC_SHUTDOWN(MKL_ERROR_CLUSTER, 'NONE', -999)
#endif
            END SELECT
         CASE ('FFT')
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) &
               CALL SCARC_SHUTDOWN(FFT_ERROR_UNSTRUCTURED, 'NONE', -999)
            TYPE_PRECON = NSCARC_RELAX_FFT
         CASE ('FFT_OVERLAP')
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) &
               CALL SCARC_SHUTDOWN(FFT_ERROR_UNSTRUCTURED, 'NONE', -999)
            TYPE_PRECON = NSCARC_RELAX_FFT_OVERLAP
         CASE ('PARDISO')
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_PARDISO
            TYPE_LUDECOMP = NSCARC_MKL_LOCAL
#else
            CALL SCARC_SHUTDOWN(MKL_ERROR_PARDISO, 'NONE', -999)
#endif
         CASE ('CLUSTER')
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_CLUSTER
            TYPE_LUDECOMP = NSCARC_MKL_GLOBAL
#else
            CALL SCARC_SHUTDOWN(MKL_ERROR_CLUSTER, 'NONE', -999)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_PRECON, -999)
      END SELECT

   CASE ('MULTIGRID')

      TYPE_METHOD = NSCARC_METHOD_MULTIGRID

      !> set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_MULTIGRID))
         CASE ('GEOMETRIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
         CASE ('ALGEBRAIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
         CASE DEFAULT
            CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_MULTIGRID, -999)
      END SELECT

      !> set type of smoother (JACOBI/SSOR)
      SELECT CASE (TRIM(SCARC_SMOOTH))                        !> use same parameters as for preconditioner
         CASE ('JACOBI')
            TYPE_SMOOTH = NSCARC_RELAX_JACOBI
         CASE ('SSOR')
            TYPE_SMOOTH = NSCARC_RELAX_SSOR
         CASE ('FFT')
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) &
               CALL SCARC_SHUTDOWN(FFT_ERROR_UNSTRUCTURED, 'NONE', -999)
            TYPE_SMOOTH = NSCARC_RELAX_FFT
         CASE ('FFT_OVERLAP')
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) &
               CALL SCARC_SHUTDOWN(FFT_ERROR_UNSTRUCTURED, 'NONE', -999)
            TYPE_SMOOTH = NSCARC_RELAX_FFT_OVERLAP
         CASE ('PARDISO')
#ifdef WITH_MKL
            TYPE_SMOOTH = NSCARC_RELAX_PARDISO
#else
            CALL SCARC_SHUTDOWN(MKL_ERROR_PARDISO, 'NONE', -999)
#endif
         CASE ('CLUSTER')
#ifdef WITH_MKL
            TYPE_SMOOTH = NSCARC_RELAX_CLUSTER
#else
            CALL SCARC_SHUTDOWN(MKL_ERROR_CLUSTER, 'NONE', -999)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_SMOOTH, -999)
      END SELECT

#ifdef WITH_MKL
   CASE ('MKL')

      TYPE_METHOD  = NSCARC_METHOD_LUDECOMP

      !> set type of MKL method (global/local)
      SELECT CASE (TRIM(SCARC_MKL))                      !Achtung, hier noch nacharbeiten!
         CASE ('GLOBAL')
#ifdef WITH_MKL
            TYPE_LUDECOMP = NSCARC_MKL_GLOBAL
#else
            CALL SCARC_SHUTDOWN(MKL_ERROR_CLUSTER, 'NONE', -999)
#endif
         CASE ('LOCAL')
#ifdef WITH_MKL
            TYPE_LUDECOMP = NSCARC_MKL_LOCAL
#else
            CALL SCARC_SHUTDOWN(MKL_ERROR_PARDISO, 'NONE', -999)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_MKL, -999)
      END SELECT
#endif

   CASE DEFAULT
      CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_METHOD, -999)

END SELECT

!> 
!> if a multigrid solver is used (either as main solver or as preconditioner)
!> set types for multigrid, coarse grid solver and cycling pattern
!> 
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_RELAX_MULTIGRID) THEN

   !> set type of multigrid (GEOMETRIC/ALGEBRAIC with corresponding coarsening strategy)
   SELECT CASE (TRIM(SCARC_MULTIGRID))

      CASE ('GEOMETRIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC

      CASE ('ALGEBRAIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC

         !> set type of coarsening strategy (STANDARD/AGGRESSIVE)
         SELECT CASE (TRIM(SCARC_MULTIGRID_COARSENING))
            CASE ('BASIC')
               TYPE_COARSENING = NSCARC_COARSENING_BASIC
            CASE ('FALGOUT')
               TYPE_COARSENING = NSCARC_COARSENING_FALGOUT
            CASE ('RS3')
               TYPE_COARSENING = NSCARC_COARSENING_RS3
            CASE ('A1')
               TYPE_COARSENING = NSCARC_COARSENING_A1
            CASE ('A2')
               TYPE_COARSENING = NSCARC_COARSENING_A2
            CASE DEFAULT
            CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_MULTIGRID_COARSENING, -999)
         END SELECT

      CASE DEFAULT
         CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_MULTIGRID, -999)
   END SELECT

   !> set type of cycling pattern (F/V/W)
   SELECT CASE (TRIM(SCARC_MULTIGRID_CYCLE))
      CASE ('F')
         TYPE_CYCLING = NSCARC_CYCLING_F
      CASE ('V')
         TYPE_CYCLING = NSCARC_CYCLING_V
      CASE ('W')
         TYPE_CYCLING = NSCARC_CYCLING_W
      CASE DEFAULT
         CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_MULTIGRID_CYCLE, -999)
   END SELECT

   !> set type of interpolation (STANDARD/DIRECT/MULTIPASS)
   SELECT CASE (TRIM(SCARC_MULTIGRID_INTERPOL))
      CASE ('STANDARD')
         TYPE_INTERPOL = NSCARC_INTERPOL_STANDARD
      CASE ('CONSTANT')
         TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
      CASE ('BILINEAR')
         TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
      CASE ('CLASSICAL')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL
      CASE ('CLASSICAL2')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL2
      CASE ('DIRECT')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT
      CASE ('DIRECT_BDRY')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT_BDRY
      CASE DEFAULT
         CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_MULTIGRID_INTERPOL, -999)
   END SELECT

ENDIF

!> set type of coarse grid solver
SELECT CASE (TRIM(SCARC_COARSE))
   CASE ('ITERATIVE')
      TYPE_COARSE = NSCARC_COARSE_ITERATIVE
      TYPE_KRYLOV = NSCARC_KRYLOV_CG
   CASE ('DIRECT')
#ifdef WITH_MKL
      TYPE_COARSE   = NSCARC_COARSE_DIRECT
      TYPE_LUDECOMP = NSCARC_MKL_COARSE
#else
      CALL SCARC_SHUTDOWN(MKL_ERROR_CLUSTER, 'NONE', -999)
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_COARSE, -999)
END SELECT

!> 
!> set type of accuracy (ABSOLUTE/RELATIVE)
!> 
SELECT CASE (TRIM(SCARC_ACCURACY))
   CASE ('ABSOLUTE')
      TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
   CASE ('RELATIVE')
      TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_ACCURACY, -999)
END SELECT

!> 
!> set type of precision for preconditioner (SINGLE/DOUBLE)
!> 
SELECT CASE (TRIM(SCARC_PRECISION))
   CASE ('SINGLE')
      TYPE_PRECISION = NSCARC_PRECISION_FB
   CASE ('DOUBLE')
      TYPE_PRECISION = NSCARC_PRECISION_EB
   CASE DEFAULT
      CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_PRECISION, -999)
END SELECT

!> 
!> set level of debugging (NONE/LESS/MUCH)
!> 
SELECT CASE (TRIM(SCARC_DEBUG))
   CASE ('NONE')
      TYPE_DEBUG = NSCARC_UNDEFINED_INT
   CASE ('LESS')
      TYPE_DEBUG = NSCARC_DEBUG_LESS
   CASE ('MUCH')
      TYPE_DEBUG = NSCARC_DEBUG_MUCH
   CASE DEFAULT
      CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_DEBUG, -999)
END SELECT

!> 
!> set level for verbose messages (NONE/LESS/MUCH)
!> 
SELECT CASE (TRIM(SCARC_VERBOSE))
   CASE ('NONE')
      TYPE_VERBOSE = NSCARC_UNDEFINED_INT
   CASE ('LESS')
      TYPE_VERBOSE = NSCARC_VERBOSE_LESS
   CASE ('MUCH')
      TYPE_VERBOSE = NSCARC_VERBOSE_MUCH
   CASE DEFAULT
      CALL SCARC_SHUTDOWN('Error with input parameter ', SCARC_VERBOSE, -999)
END SELECT

!>
!> Define some logical variables, simply for notational convenience
!>
BCG = TYPE_METHOD == NSCARC_METHOD_KRYLOV
BCGGMG  = BCG .AND. TYPE_PRECON == NSCARC_RELAX_MULTIGRID .AND. &
          TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC

BMG = TYPE_METHOD == NSCARC_METHOD_MULTIGRID
BGMG = BMG .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC
BAMG = BMG .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC

BTWOLEVEL   = BCG .AND. &
              TYPE_PRECON /= NSCARC_RELAX_MULTIGRID .AND. &
              TYPE_TWOLEVEL > NSCARC_TWOLEVEL_NONE
BMULTILEVEL = BGMG .OR. BCGGMG .OR. BTWOLEVEL

BCGADD    = BTWOLEVEL .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_ADD
BCGMUL    = BTWOLEVEL .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_MUL
BCGCOARSE = BTWOLEVEL .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_COARSE

BFFT =  TYPE_PRECON == NSCARC_RELAX_FFT         .OR. &
        TYPE_SMOOTH == NSCARC_RELAX_FFT         .OR. &
        TYPE_PRECON == NSCARC_RELAX_FFT_OVERLAP .OR. &
        TYPE_SMOOTH == NSCARC_RELAX_FFT_OVERLAP

BMKL = (TYPE_PRECON >= NSCARC_RELAX_PARDISO) .OR. &
       (TYPE_SMOOTH >= NSCARC_RELAX_PARDISO) .OR. &
       (BMULTILEVEL .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT)

BSTRUCTURED = TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED 

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_PARSE_INPUT


!> ------------------------------------------------------------------------------------------------
!> Determine number of grid levels  (1 for CG/BICG-method, NLEVEL for MG-method)
!> Note: NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LEVELS
#ifdef WITH_MKL
INTEGER :: NL
#endif

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_LEVELS')

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !> ----------------------- Krylov method --------------------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV: SELECT CASE (TYPE_PRECON)

#ifdef WITH_MKL
         CASE (NSCARC_RELAX_PARDISO)
            IF (BTWOLEVEL) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_LU_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF
            TYPE_LU_LEVEL(NLEVEL_MIN) = NSCARC_MKL_LOCAL

         CASE (NSCARC_RELAX_CLUSTER)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            TYPE_LU_LEVEL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
#endif

         CASE (NSCARC_RELAX_MULTIGRID)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
            IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_LU_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif

         CASE DEFAULT
            IF (BTWOLEVEL) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_LU_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF
      END SELECT SELECT_KRYLOV

   !> ----------------------- Multigrid method -----------------------------------------
   CASE (NSCARC_METHOD_MULTIGRID)

      SELECT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)
      
         !> predefined hierarchy of levels in case of geometric multigrid-method
         CASE (NSCARC_MULTIGRID_GEOMETRIC)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
      
#ifdef WITH_MKL
            SELECT_SMOOTHER: SELECT CASE (TYPE_SMOOTH)
               CASE (NSCARC_RELAX_PARDISO)
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     TYPE_LU_LEVEL(NL) = NSCARC_MKL_LOCAL
                  ENDDO
               CASE (NSCARC_RELAX_CLUSTER)
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     TYPE_LU_LEVEL(NL) = NSCARC_MKL_GLOBAL
                  ENDDO
            END SELECT SELECT_SMOOTHER
      
            IF (TYPE_LUDECOMP == NSCARC_MKL_COARSE) TYPE_LU_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif
      
         !> first, only finest level is set, further levels are defined during coarsening process
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_AMG)
      
      END SELECT SELECT_MULTIGRID
      
   !> ----------------------- MKL method -----------------------------------------
   CASE (NSCARC_METHOD_LUDECOMP)

      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)

      SELECT_MKL: SELECT CASE (TYPE_LUDECOMP)
         CASE (NSCARC_MKL_LOCAL)
            TYPE_LU_LEVEL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
         CASE (NSCARC_MKL_GLOBAL)
            TYPE_LU_LEVEL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
      END SELECT SELECT_MKL

END SELECT SELECT_METHOD

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_LEVELS


!> ------------------------------------------------------------------------------------------------
!> Setup single level in case of default Krylov method
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS(NTYPE)
INTEGER, INTENT(IN) :: NTYPE
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM, NLEVEL

SELECT_LEVEL_TYPE: SELECT CASE (NTYPE)

   !> only use the finest grid level
   CASE(NSCARC_LEVEL_SINGLE)

      NLEVEL     = 1
      NLEVEL_MIN = 1
      NLEVEL_MAX = 1

   !> determine maximum number of possible levels based on number of grid cells
   CASE(NSCARC_LEVEL_MULTI)

      NLEVEL = NSCARC_LEVEL_MAX
      KLEVEL = NSCARC_LEVEL_MAX       
      
      DO NM=1,NMESHES
         KLEVEL(1)=SCARC_GET_MAX_LEVEL(MESHES(NM)%IBAR,1)
         IF (.NOT.TWO_D) KLEVEL(2)=SCARC_GET_MAX_LEVEL(MESHES(NM)%JBAR,2)
         KLEVEL(3)=SCARC_GET_MAX_LEVEL(MESHES(NM)%KBAR,3)
         KLEVEL_MIN = MINVAL(KLEVEL)
         IF (KLEVEL_MIN<NLEVEL) NLEVEL=KLEVEL_MIN
      ENDDO
      
      NLEVEL_MIN  = 1
      IF (BGMG.OR.BCGGMG) THEN
         IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
            NLEVEL_MAX  = NLEVEL_MIN + SCARC_MULTIGRID_LEVEL - 1
         ELSE
            NLEVEL_MAX  = NLEVEL
         ENDIF
      ELSE IF (BTWOLEVEL) THEN
         IF (SCARC_COARSE_LEVEL /= -1) THEN
            NLEVEL_MAX  = NLEVEL_MIN + SCARC_COARSE_LEVEL - 1
         ELSE
            NLEVEL_MAX  = NLEVEL
         ENDIF
      ENDIF

   !> use user specified number of grid levels
   CASE(NSCARC_LEVEL_AMG)

      NLEVEL_MIN = 1
      IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
         NLEVEL_MAX  = SCARC_MULTIGRID_LEVEL
      ELSE
         NLEVEL_MAX  = NSCARC_LEVEL_MAX
      ENDIF
      NLEVEL = NLEVEL_MAX

   END SELECT SELECT_LEVEL_TYPE

END SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS


!> ------------------------------------------------------------------------------------------------
!> Determine maximum number of possible levels on direction IOR0 of mesh NM
!> In case of GMG- or 2-Level-method, NC must be divisable by 2 at least one time
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_GET_MAX_LEVEL(NC, IOR0)
INTEGER, INTENT(IN) :: NC, IOR0
INTEGER :: NC0, NL

CALL SCARC_ENTER_ROUTINE('SCARC_GET_MAX_LEVEL')

IF (BMULTILEVEL .AND.  MOD(NC,2)/=0) THEN
   SELECT CASE (IOR0)
      CASE (1)
         CALL SCARC_SHUTDOWN('Step size in x-direction not divisable by 2 ', 'NONE', NC)
      CASE (2)
         CALL SCARC_SHUTDOWN('Step size in y-direction not divisable by 2 ', 'NONE', NC)
      CASE (3)
         CALL SCARC_SHUTDOWN('Step size in z-direction not divisable by 2 ', 'NONE', NC)
   END SELECT
ENDIF

!> Divide by 2 as often as possible or until user defined max-level is reached
NC0=NC
DO NL=1,NSCARC_LEVEL_MAX
   NC0=NC0/2
   IF (MOD(NC0,2)/=0) EXIT                !> NC no longer divisable by two
   IF (NL==SCARC_MULTIGRID_LEVEL) EXIT    !> max number of levels defined by user
   IF (NC0==1) EXIT                       !> NC is power of two, minimum has been reached
ENDDO

SCARC_GET_MAX_LEVEL=NL
CALL SCARC_LEAVE_ROUTINE()

RETURN
END FUNCTION SCARC_GET_MAX_LEVEL

!> ------------------------------------------------------------------------------------------------
!> Allocate basic ScaRC-structures for all needed levels
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TYPES
INTEGER :: NM

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_TYPES')

!> Basic information for all requested grid levels
ALLOCATE (SCARC(NMESHES), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERROR)

!> Basic solver stack
ALLOCATE (STACK(NSCARC_STACK_MAX), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'STACK', IERROR)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   !> Needed information about other meshes
   ALLOCATE (SCARC(NM)%OSCARC(NMESHES), STAT=IERROR) 
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSCARC', IERROR)

   !> Information for single grid levels
   ALLOCATE (SCARC(NM)%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'LEVEL', IERROR)

   !> Vector scopes (1: for main solver, 2: for preconditioner)
   ALLOCATE (SCARC(NM)%SCOPE(2, NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'ENVIRON', IERROR)

   !> Matrices on different grid levels
   ALLOCATE (SCARC(NM)%SYSTEM(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'SYSTEM', IERROR)

   !> Multigrid type
   IF (BMG .OR. BCGGMG) THEN
      ALLOCATE (SCARC(NM)%CYCLING(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'CYCLING', IERROR)
   ENDIF

   !> Discretizations for different grid leves
   IF (BFFT) THEN
      ALLOCATE (SCARC(NM)%FFT(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'MULTIGRID', IERROR)
   ENDIF

#ifdef WITH_MKL
   !> Information for Intel MKL routines on different grid levels
   IF (BMKL) THEN
      ALLOCATE (SCARC(NM)%MKL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'MKL', IERROR)
   ENDIF
#endif

ENDDO MESHES_LOOP
CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_TYPES


!> ----------------------------------------------------------------------------------------------------
!> Setup geometry information for mesh NM
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MESHES
INTEGER  :: NL, NM, IX, IY, IZ, IO
INTEGER  :: NX, NY, NZ
TYPE (MESH_TYPE)         , POINTER :: M
TYPE (SCARC_TYPE)        , POINTER :: S
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_COORD_TYPE)  , POINTER :: C
REAL(EB), DIMENSION(:), POINTER :: XCOR, YCOR, ZCOR, XMID, YMID, ZMID

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_MESHES')
LEVEL_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   S => SCARC(NM)

   !> store bounds of mesh in SCARC-structure
   S%XS = M%XS                 
   S%XF = M%XF                 
   S%YS = M%YS                 
   S%YF = M%YF                 
   S%ZS = M%ZS                 
   S%ZF = M%ZF                 

   S%IBAR = M%IBAR
   S%JBAR = M%JBAR
   S%KBAR = M%KBAR

   NX  = M%IBAR
   NY  = M%JBAR
   NZ  = M%KBAR
 
   LEVEL_LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

      IF (NL > NLEVEL_MIN .AND. BAMG) EXIT LEVEL_LEVEL_LOOP
      L => S%LEVEL(NL)

      L%NX = NX
      L%NY = NY
      L%NZ = NZ

      L%NC = L%NX * L%NY * L%NZ

      L%N_WALL_CELLS_EXT = M%N_EXTERNAL_WALL_CELLS
      L%N_WALL_CELLS_INT = M%N_INTERNAL_WALL_CELLS

      L%NW = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT 

      !> Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(L%CELL%NC_LOCAL     , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(L%CELL%NC_LOCAL_DOF , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL_DOF')
      CALL SCARC_ALLOCATE_INT1(L%CELL%NC_OFFSET    , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET')
      CALL SCARC_ALLOCATE_INT1(L%CELL%NC_OFFSET_DOF, 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET_DOF')

      IF (NL == NLEVEL_MIN) THEN
         L%N_OBST = M%N_OBST
         ALLOCATE(L%OBST(L%N_OBST), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','OBST',IERROR)
         DO IO = 1, L%N_OBST
            L%OBST(IO)%I1  = M%OBSTRUCTION(IO)%I1
            L%OBST(IO)%I2  = M%OBSTRUCTION(IO)%I2
            L%OBST(IO)%J1  = M%OBSTRUCTION(IO)%J1
            L%OBST(IO)%J2  = M%OBSTRUCTION(IO)%J2
            L%OBST(IO)%K1  = M%OBSTRUCTION(IO)%K1
            L%OBST(IO)%K2  = M%OBSTRUCTION(IO)%K2
         ENDDO
      ENDIF

      !> get coordination information
      C => L%COORD

      C%DX = (S%XF-S%XS)/REAL(L%NX,EB)
      C%DY = (S%YF-S%YS)/REAL(L%NY,EB)
      C%DZ = (S%ZF-S%ZS)/REAL(L%NZ,EB)

      C%DXI = 1.0_EB/C%DX
      C%DYI = 1.0_EB/C%DY
      C%DZI = 1.0_EB/C%DZ

      C%DXI2 = C%DXI**2
      C%DYI2 = C%DYI**2
      C%DZI2 = C%DZI**2

      !> needed in case of GMG with multiple grid levels
      NX=NX/2
      IF (.NOT.TWO_D) NY=NY/2
      NZ=NZ/2

      !> Allocate vectors for coordinate information
      IF (NL == NLEVEL_MIN) THEN
   
         XCOR => M%X
         YCOR => M%Y
         ZCOR => M%Z

         XMID => M%XC
         YMID => M%YC
         ZMID => M%ZC

      ELSE

         CALL SCARC_ALLOCATE_REAL1(C%XCOR, 0, L%NX, NSCARC_INIT_NONE, 'XCOR')
         CALL SCARC_ALLOCATE_REAL1(C%YCOR, 0, L%NY, NSCARC_INIT_NONE, 'YCOR')
         CALL SCARC_ALLOCATE_REAL1(C%ZCOR, 0, L%NZ, NSCARC_INIT_NONE, 'ZCOR')

         !> compute coordinates in x-, y- and z-direction
         DO IX = 0, L%NX
            C%XCOR(IX) = S%XS + IX*C%DX
         ENDDO
         DO IY = 0, L%NY
            C%YCOR(IY) = S%YS + IY*C%DY
         ENDDO
         DO IZ = 0, L%NZ
            C%ZCOR(IZ) = S%ZS + IZ*C%DZ
         ENDDO

         !> compute midpoints in x-, y- and z-direction
         CALL SCARC_ALLOCATE_REAL1(C%XMID, 0, L%NX+1, NSCARC_INIT_NONE, 'XMID')
         CALL SCARC_ALLOCATE_REAL1(C%YMID, 0, L%NY+1, NSCARC_INIT_NONE, 'YMID')
         CALL SCARC_ALLOCATE_REAL1(C%ZMID, 0, L%NZ+1, NSCARC_INIT_NONE, 'ZMID')

         C%XMID(0) = S%XS - 0.5_EB*C%DX
         DO IX = 1, L%NX
            C%XMID(IX) = 0.5_EB*(C%XCOR(IX-1) + C%XCOR(IX))
         ENDDO
         C%XMID(L%NX+1) = S%XF + 0.5_EB*C%DX

         C%YMID(0) = S%YS - 0.5_EB*C%DY
         DO IY = 1, L%NY
            C%YMID(IY) = 0.5_EB*(C%YCOR(IY-1) + C%YCOR(IY))
         ENDDO
         C%YMID(L%NY+1) = S%YF + 0.5_EB*C%DY

         C%ZMID(0) = S%ZS - 0.5_EB*C%DZ
         DO IZ = 1, L%NZ
            C%ZMID(IZ) = 0.5_EB*(C%ZCOR(IZ-1) + C%ZCOR(IZ))
         ENDDO
         C%ZMID(L%NZ+1) = S%ZF + 0.5_EB*C%DZ

         XCOR => C%XCOR
         YCOR => C%YCOR
         ZCOR => C%ZCOR

         XMID => C%XMID
         YMID => C%YMID
         ZMID => C%ZMID
         
      ENDIF

      !> Allocate vectors for step sizes in different directions
      CALL SCARC_ALLOCATE_REAL1(C%DXL, 0, L%NX, NSCARC_INIT_ZERO, 'DXL')
      CALL SCARC_ALLOCATE_REAL1(C%DYL, 0, L%NY, NSCARC_INIT_ZERO, 'DYL')
      CALL SCARC_ALLOCATE_REAL1(C%DZL, 0, L%NZ, NSCARC_INIT_ZERO, 'DZL')

      !> set step sizes between cell midpoints, use interior step sizes for ghost cells as initial values
      !> correct sizes for ghost cells are exchanged later
      DO IX = 1, L%NX-1
         C%DXL(IX) = XMID(IX+1) - XMID(IX)
      ENDDO
      C%DXL(0)    = C%DXL(1)
      C%DXL(L%NX) = C%DXL(L%NX-1)

      DO IY = 1, L%NY-1
         C%DYL(IY) = YMID(IY+1) - YMID(IY)
      ENDDO
      C%DYL(0)    = C%DYL(1)
      C%DYL(L%NY) = C%DYL(L%NY-1)

      DO IZ = 1, L%NZ-1
         C%DZL(IZ) = ZMID(IZ+1) - ZMID(IZ)
      ENDDO
      C%DZL(0)    = C%DZL(1)
      C%DZL(L%NZ) = C%DZL(L%NZ-1)

   ENDDO LEVEL_LEVEL_LOOP
ENDDO LEVEL_MESHES_LOOP

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_MESHES

!> ------------------------------------------------------------------------------------------------
!> Setup communication structure for data exchange along mesh interfaces
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERFACES
INTEGER :: NM, NOM, NL
TYPE (MESH_TYPE)          , POINTER :: M
TYPE (SCARC_TYPE)         , POINTER :: S
TYPE (SCARC_LEVEL_TYPE)   , POINTER :: LF, LC, OLF, OLC
TYPE (SCARC_MATRIX_TYPE)  , POINTER :: OFA, OCA
TYPE (SCARC_EXCHANGE_TYPE), POINTER :: OX

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_INTERFACES')

!> Initialize communication counter for ScaRC, use same TAG for all communications
TAG   = 99
N_REQ =  0
N_EXCHANGE =  0

!> Allocate basic WALL and FACE types on mesh NM for all requested grid levels
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M  => MESHES(NM)
   S  => SCARC(NM)
   LF => S%LEVEL(NLEVEL_MIN)                         !> level fine 

   ALLOCATE(LF%WALL(LF%NW), STAT=IERROR)
   CALL ChkMemErr('SCARC_SETUP_INTERFACES','WALL',IERROR)

   ALLOCATE(LF%FACE(-3:3), STAT=IERROR)
   CALL ChkMemErr('SCARC_SETUP_INTERFACES','FACE',IERROR)

   IF (NLEVEL_MAX > NLEVEL_MIN .AND. .NOT.BAMG) THEN
      DO NL=NLEVEL_MIN+1,NLEVEL_MAX
         LC => S%LEVEL(NL)                           !> level coarse
         ALLOCATE(LC%FACE(-3:3), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_INTERFACES','WALL',IERROR)
      ENDDO
   ENDIF

   !> Get communication lengths for other meshes from main program
   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES
      !IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
      OX => S%OSCARC(NOM)%EXCHANGE                   !> exchange structures for other mesh
      OX%NIC_R    = M%OMESH(NOM)%NIC_R
      OX%NIC_S    = M%OMESH(NOM)%NIC_S
      OX%NICMAX_R = M%OMESH(NOM)%NIC_R
      OX%NICMAX_S = M%OMESH(NOM)%NIC_S
      IF (OX%NICMAX_R==0 .AND. OX%NICMAX_S==0) CYCLE OTHER_MESHES_LOOP
      N_EXCHANGE  = N_EXCHANGE+1
   ENDDO OTHER_MESHES_LOOP
ENDDO MESHES_LOOP

!> Initialize level structures on neighboring meshes
LEVEL_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S => SCARC(NM)
   LEVEL_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      !IF (NOM == NM) CYCLE LEVEL_OTHER_MESHES_LOOP

      OX => S%OSCARC(NOM)%EXCHANGE
      IF (OX%NICMAX_S==0 .AND. OX%NICMAX_R==0)  CYCLE LEVEL_OTHER_MESHES_LOOP

      ALLOCATE (S%OSCARC(NOM)%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OS%LEVEL', IERROR)

      ALLOCATE (S%OSCARC(NOM)%SYSTEM(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OS%SYSTEM', IERROR)

      OLF => S%OSCARC(NOM)%LEVEL(NLEVEL_MIN)                       ! pointer to fine level on other mesh
      OFA => S%OSCARC(NOM)%SYSTEM(NLEVEL_MIN)%A                    ! pointer to matrix on other mesh

      OLF%NX = MESHES(NOM)%IBAR                                    !> number of cells in x-direction on other mesh
      OLF%NY = MESHES(NOM)%JBAR                                    !> number of cells in y-direction on other mesh
      OLF%NZ = MESHES(NOM)%KBAR                                    !> number of cells in z-direction on other mesh

      OLF%N_WALL_CELLS_EXT = MESHES(NOM)%N_EXTERNAL_WALL_CELLS     !> number of external wall cells on other mesh
      OLF%N_WALL_CELLS_INT = MESHES(NOM)%N_INTERNAL_WALL_CELLS     !> number of external wall cells on other mesh

      OLF%NC  = OLF%NX*OLF%NY*OLF%NZ                               !> number of cells on other mesh
      OLF%NW  = OLF%N_WALL_CELLS_EXT + OLF%N_WALL_CELLS_INT        !> number of walls cell on other mesh
      OLF%NCG = 0                                                  !> number of ghost cells on other mesh

      OFA%NAV = 0                                                  !> parameters for size of matrix on other mesh
      OFA%NAC = 0
      OFA%NAR = 0

      IF (OX%NICMAX_S == 0 .AND. OX%NICMAX_R == 0) CYCLE LEVEL_OTHER_MESHES_LOOP

      !> In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN .AND. .NOT.BAMG) THEN

         DO NL=NLEVEL_MIN+1,NLEVEL_MAX

            OLF => S%OSCARC(NOM)%LEVEL(NL-1)                       !> pointer to fine level on other mesh
            OLC => S%OSCARC(NOM)%LEVEL(NL)                         !> pointer to coarse level on other mesh
            OCA => S%OSCARC(NOM)%SYSTEM(NL)%A                      !> pointer to coarse matrix of other meshes

            OLC%NX = OLF%NX/2                          
            IF (TWO_D) THEN                            
               OLC%NY = 1
            ELSE
               OLC%NY = OLF%NY/2
            ENDIF
            OLC%NZ = OLF%NZ/2

            OLC%N_WALL_CELLS_EXT = 4 * (OLC%NX * OLC%NZ + OLC%NX * OLC%NY + OLC%NY * OLC%NZ)   !> ACHTUNG: warum 4?

            OLC%NC  = OLC%NX * OLC%NY * OLC%NZ                     !> see above
            OLC%NW  = OLC%N_WALL_CELLS_EXT 
            OLC%NCG = 0

            OCA%NAV = 0                                            !> see above
            OCA%NAC = 0
            OCA%NAR = 0

         ENDDO
      ENDIF

   ENDDO LEVEL_OTHER_MESHES_LOOP
ENDDO LEVEL_MESHES_LOOP

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_INTERFACES

!> ------------------------------------------------------------------------------------------------
!> Initialize arrays for data exchange
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE
INTEGER :: NL, NM, NOM, NLEN_SEND, NLEN_RECV
INTEGER :: NLMIN, NLMAX
INTEGER :: INBR
TYPE (SCARC_TYPE)         , POINTER :: S
TYPE (SCARC_LEVEL_TYPE)   , POINTER :: L, OL
TYPE (SCARC_EXCHANGE_TYPE), POINTER :: OX

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_EXCHANGE')
!>
!> Allocate request array for data exchanges
!> Exchange basic information about wall sizes (needed for the dimensioning of the exchange buffers)
!>
IF (N_MPI_PROCESSES>1) THEN
   ALLOCATE (REQ(N_EXCHANGE*40), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'REQ', IERROR)
   REQ = MPI_REQUEST_NULL
ENDIF

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_BASIC, NLEVEL_MIN)

!>
!> Allocate send and receive buffers (real and integer) in correct lengths
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S => SCARC(NM)
   DO INBR = 1, S%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)

      L  => SCARC(NM)%LEVEL(NLEVEL_MIN)
      OX => SCARC(NM)%OSCARC(NOM)%EXCHANGE
      OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)

      !> allocate real buffers in maximum needed length
      NLEN_SEND = NSCARC_MAX_STENCIL*MAX(OL%NWL, OL%NCG)+10
      NLEN_RECV = NLEN_SEND

      CALL SCARC_ALLOCATE_REAL1(OX%SEND_REAL, 1, NLEN_SEND, NSCARC_INIT_ZERO, 'SEND_REAL')
      CALL SCARC_ALLOCATE_REAL1(OX%RECV_REAL, 1, NLEN_RECV, NSCARC_INIT_ZERO, 'RECV_REAL')

      !> allocate integer buffers in maximum needed length
      OL%NCPLS = OL%NCPL
      IF (OL%NCG == OL%NWL) THEN
         NLEN_SEND = 15*OL%NWL
         NLEN_RECV = 15*OL%NCG
         OL%NCPLR = OL%NCPL
      ELSE IF (OL%NCG == 2*OL%NWL) THEN
         NLEN_SEND = 13*OL%NWL + 2*OL%NCG
         NLEN_RECV = 15*OL%NCG
         OL%NCPLR =  1
      ELSE IF (OL%NWL == 2*OL%NCG) THEN
         NLEN_SEND = 15*OL%NWL
         NLEN_RECV = 13*OL%NCG + 2*OL%NWL
         OL%NCPLR =  2
      ENDIF

      CALL SCARC_ALLOCATE_INT1(OX%SEND_INT, 1, NLEN_SEND, NSCARC_INIT_ZERO, 'SEND_INT')
      CALL SCARC_ALLOCATE_INT1(OX%RECV_INT, 1, NLEN_RECV, NSCARC_INIT_ZERO, 'RECV_INT')

      !> neighboring wall (and face?) structures for common wall cells
      ALLOCATE (OL%WALL(OL%NCG), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OL%WALL', IERROR)

      !ALLOCATE (OL%FACE(-3:3), STAT=IERROR)                                     !> needed ?
      !CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OL%FACE', IERROR)

      !> In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN .AND. .NOT.BAMG) THEN
         DO NL=NLEVEL_MIN+1,NLEVEL_MAX
            OL => S%OSCARC(NOM)%LEVEL(NL)                          
            ALLOCATE (OL%WALL(OL%NCG), STAT=IERROR)
            CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OL%WALL', IERROR)
            !ALLOCATE (OL%FACE(-3:3), STAT=IERROR)
            !CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OL%FACE', IERROR)         !> needed ?
         ENDDO
      ENDIF

   ENDDO
ENDDO

!>
!> Initialize communication structures on finest level (if there is more than 1 mesh)
!>
IF (N_MPI_PROCESSES > 1) THEN
   NLMIN = NLEVEL_MIN
   NLMAX = NLEVEL_MIN
   IF (.NOT.BAMG) NLMAX = NLEVEL_MAX
   DO NL = NLMIN, NLMAX
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_DISCRETIZATION, NL)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WIDTHINFO, NL)
!      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WALLINFO, NL)
   ENDDO
ENDIF

!>
!> Correct boundary types for cells adjacent to obstructions on ghost cells
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (.NOT.PRES_ON_WHOLE_DOMAIN) CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NLEVEL_MIN)
   IF (.NOT.BAMG) THEN
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX
         CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
      ENDDO
   ENDIF
ENDDO

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_EXCHANGE

!> ----------------------------------------------------------------------------------------------------
!> Setup neighborship structures and boundary conditions
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLS
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES
INTEGER :: NL, NM, NM2
INTEGER :: IREFINE, IFACE
INTEGER :: INBR
INTEGER :: IWG, IWL, IWC
INTEGER :: NOM_LAST , NOM
INTEGER :: NCPL_LAST, NCPL
INTEGER :: IOR_LAST, IOR0, JOR0
INTEGER :: FACE_NEIGHBORS(-3:3, NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: MESH_NEIGHBORS(6*NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: NUM_FACE_NEIGHBORS(-3:3)
INTEGER :: NUM_MESH_NEIGHBORS
INTEGER :: FACE_ORDER_XYZ(6) = (/1,-1,2,-2,3,-3/)           !> Coordinate direction related order of mesh faces
LOGICAL :: BKNOWN(-3:3), IS_BC_DIRICHLET, IS_OPEN_BOUNDARY
TYPE (SCARC_TYPE)        , POINTER :: S
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L, OL, LF, OLF, LC, OLC
TYPE (SCARC_COORD_TYPE)  , POINTER :: C
TYPE (SCARC_MAPPING_TYPE), POINTER :: M, MC
TYPE (WALL_TYPE)         , POINTER :: WC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_WALLS')

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S => SCARC(NM)

   FACE_NEIGHBORS = -1
   MESH_NEIGHBORS = -1

   NUM_FACE_NEIGHBORS = 0
   NUM_MESH_NEIGHBORS = 0

   !> For all solvers: Determine array WALL and PRESSURE_BC_INDEX on finest level 
   L => S%LEVEL(NLEVEL_MIN)
   C => L%COORD

   L%NCE = L%NCS
   L%NCO = L%NCS

   L%N_WALL_CELLS_EXT = MESHES(NM)%N_EXTERNAL_WALL_CELLS
   L%N_WALL_CELLS_INT = MESHES(NM)%N_INTERNAL_WALL_CELLS

   !> Determine basic data for single faces (orientation, dimensions)
   FACES_OF_MESH_INDEX1: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE FACES_OF_MESH_INDEX1

      !> information about face orientation and local dimensions
      SELECT CASE (ABS(IOR0))
         CASE (1)
            L%FACE(IOR0)%NFC  =  L%NX                   !> number of cells between opposite mesh faces
            L%FACE(IOR0)%NFX  =  1                      !> number of cells in x-direction
            L%FACE(IOR0)%NFY  =  L%NY                   !> number of cells in y-direction
            L%FACE(IOR0)%NFZ  =  L%NZ                   !> number of cells in z-direction
            L%FACE(IOR0)%NFW  =  L%NY*L%NZ              !> number of wall cells at that face
            L%FACE(IOR0)%DH   => C%DXL                  !> step size vector between opposite mesh faces
         CASE (2)
            L%FACE(IOR0)%NFC  =  L%NY                   !> see above
            L%FACE(IOR0)%NFX  =  L%NX
            L%FACE(IOR0)%NFY  =  1
            L%FACE(IOR0)%NFZ  =  L%NZ
            L%FACE(IOR0)%NFW  =  L%NX*L%NZ
            L%FACE(IOR0)%DH   => C%DYL
         CASE (3)
            L%FACE(IOR0)%NFC  =  L%NZ                   !> see above
            L%FACE(IOR0)%NFX  =  L%NX
            L%FACE(IOR0)%NFY  =  L%NY
            L%FACE(IOR0)%NFZ  =  1
            L%FACE(IOR0)%NFW  =  L%NX*L%NY
            L%FACE(IOR0)%DH   => C%DZL
      END SELECT

   ENDDO FACES_OF_MESH_INDEX1

   !> Store local IWG-number for each face
   IWG = 1
   FACE_ORDER_LOOP: DO IFACE = 1, 6
      IOR0 = FACE_ORDER_XYZ(IFACE)
      L%FACE(IOR0)%IWG_PTR = IWG
      IWG = IWG + L%FACE(IOR0)%NFW
   ENDDO FACE_ORDER_LOOP

   !> loop over global IW's:
   !> store basic data and determine number of adajacent neighbors to each
   !> face with corresponding number of IW's
   IOR_LAST  =  0
   NOM_LAST  = -1
   NCPL_LAST = -1
   IWL = 0

   !> process external wall cells
   EXTERNAL_WALL_CELLS_LOOP1: DO IWG = 1, L%N_WALL_CELLS_EXT

      WC  => MESHES(NM)%WALL(IWG)
      EWC => MESHES(NM)%EXTERNAL_WALL(IWG)

      !> Determine and store neighbors, orientation and number of couplings for a single wall cell
      NOM  =  EWC%NOM
      IOR0 =  WC%ONE_D%IOR
      NCPL = (EWC%IIO_MAX - EWC%IIO_MIN + 1) * &
             (EWC%JJO_MAX - EWC%JJO_MIN + 1) * &
             (EWC%KKO_MAX - EWC%KKO_MIN + 1)

      L%WALL(IWG)%NOM  = NOM                            !> store number of neighbor in wall cell
      L%WALL(IWG)%IOR  = IOR0                           !> store orientation of that cell

      IWL = IWL + 1                                     !> count local wall cells for that face

      IF (NOM /= 0) THEN
         L%WALL(IWG)%NCPL = NCPL                        !> store number of couplings for that cell
         BKNOWN = .FALSE.
         DO JOR0 = -3, 3
            IF (JOR0 == 0) CYCLE
            DO INBR = 1, NUM_FACE_NEIGHBORS(JOR0)
               IF (FACE_NEIGHBORS(JOR0, INBR) == NOM) THEN
                  BKNOWN(JOR0) = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDDO
      ELSE
         L%WALL(IWG)%NCPL = 0                           !> no couplings
      ENDIF

      !> for wall cells with a neighbor compute extended and overlapping structures
      IF (NOM /= 0) THEN

         L%NCE = L%NCE + NCPL                           !> increase number of extended grid cells
         L%NCO = L%NCO + 1                              !> increase number of overlapping grid cells

         OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)
         OL%NCPL = NCPL                                 !> initialize own counter for local wall cells
         OL%IOR  = -IOR0                                !> initialize own orientation variable

         IF (ANY(BKNOWN)) THEN
            OL%NWL = OL%NWL + 1                         !> increase own counter for local wall cells
            OL%NCG = OL%NCG + NCPL                      !> increase counter for local ghost cells
         ELSE
            OL%NWL = 1                                  !> initialize own counter for local wall cells
            OL%NCG = NCPL                               !> initialize counter for local ghost cells
            L%NCPL_MAX  = MAX(L%NCPL_MAX, NCPL)         !> get max NCPL ever used on this mesh
         ENDIF
      ENDIF

      IF (NOM /= 0) THEN
         IF (.NOT.BKNOWN(IOR0)) THEN
            NUM_FACE_NEIGHBORS(IOR0) = NUM_FACE_NEIGHBORS(IOR0) + 1     !> increase neighbor counter for face
            FACE_NEIGHBORS(IOR0, NUM_FACE_NEIGHBORS(IOR0)) = NOM        !> store number of neighbor for face
         ENDIF
         IF (.NOT.ANY(BKNOWN)) THEN
            NUM_MESH_NEIGHBORS = NUM_MESH_NEIGHBORS + 1                 !> increase neighbor counter for mesh
            MESH_NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0)) = NOM              !> store number of neighbor for mesh
         ENDIF
      ENDIF

      IOR_LAST  = IOR0                                                  !> save former values
      NOM_LAST  = NOM
      NCPL_LAST = NCPL

   ENDDO EXTERNAL_WALL_CELLS_LOOP1
 
   !> Then process internal wall cells
   INTERNAL_WALL_CELLS_LOOP1: DO IWG = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT

      WC => MESHES(NM)%WALL(IWG)

      L%WALL(IWG)%IOR  = MESHES(NM)%WALL(IWG)%ONE_D%IOR
      L%WALL(IWG)%NOM  = 0

      L%WALL(IWG)%BTYPE  = NEUMANN
      L%WALL(IWG)%BOUNDARY_TYPE  = MESHES(NM)%WALL(IWG)%BOUNDARY_TYPE

      L%WALL(IWG)%IXG =  WC%ONE_D%II                        !> ghost cell indices
      L%WALL(IWG)%IYG =  WC%ONE_D%JJ
      L%WALL(IWG)%IZG =  WC%ONE_D%KK

      L%WALL(IWG)%IXW =  WC%ONE_D%IIG                       !> (internal) wall cell indices
      L%WALL(IWG)%IYW =  WC%ONE_D%JJG
      L%WALL(IWG)%IZW =  WC%ONE_D%KKG

   ENDDO INTERNAL_WALL_CELLS_LOOP1

   !> Allocate array which stores numbers of all neighboring meshes
   IF (NUM_MESH_NEIGHBORS /= 0) &
      CALL SCARC_ALLOCATE_INT1(S%NEIGHBORS, 1, NUM_MESH_NEIGHBORS, NSCARC_INIT_UNDEFINED, 'NEIGHBORS')
   S%N_NEIGHBORS = NUM_MESH_NEIGHBORS

   !> Store information about adjacent neighbors on different faces
   !> Allocate corresponding index arrays in OSCARC-structures
   !> First allocate administrative mapping arrays for own mesh
   CALL SCARC_SETUP_MAPPINGS(NM, NLEVEL_MIN)

   NEIGHBORS_OF_FACE_LOOP: DO IOR0 = -3, 3
      IF (IOR0 == 0) CYCLE NEIGHBORS_OF_FACE_LOOP

      !> if there are neighbors at face IOR0 store information about them
      IF (NUM_FACE_NEIGHBORS(IOR0) /= 0) THEN

         L%FACE(IOR0)%N_NEIGHBORS = NUM_FACE_NEIGHBORS(IOR0)        !> store number of neighbors on face

         !> allocate array for storing the numbers of the single neighbors
         CALL SCARC_ALLOCATE_INT1(L%FACE(IOR0)%NEIGHBORS, 1, NUM_FACE_NEIGHBORS(IOR0), NSCARC_INIT_NONE, 'NEIGHBORS')

         !> store every neighbor and allocate corresponding administration arrays
         DO INBR = 1, NUM_FACE_NEIGHBORS(IOR0)

            NOM = FACE_NEIGHBORS(IOR0, INBR)

            L%FACE(IOR0)%NEIGHBORS(INBR) = NOM             !> store NOM as a neighbor of that face and if
            CALL SCARC_UPDATE_NEIGHBORS(NM, NOM)           !> not already done also as mesh neighbor itself

            !> allocate administrative arrays for neighboring meshes
            CALL SCARC_SETUP_OMAPPINGS(NM, NOM, NLEVEL_MIN)

         ENDDO
      ENDIF
   ENDDO NEIGHBORS_OF_FACE_LOOP

   !> Second loop over external wall cells:
   !> Store detailed coordinate and cell data and get type of boundary condition
   M => L%MAP
   M%ICE_PTR = L%NCS
   M%ICO_PTR = L%NCS

   IOR_LAST = 0
   WALL_CELLS_LOOP2: DO IWG = 1, L%N_WALL_CELLS_EXT

      !> Determine neighbors and orientation again
      NOM  = L%WALL(IWG)%NOM
      IOR0 = L%WALL(IWG)%IOR
      NCPL = L%WALL(IWG)%NCPL

      WC  => MESHES(NM)%WALL(IWG)
      EWC => MESHES(NM)%EXTERNAL_WALL(IWG)

      !>
      !> Preset ScaRC's boundary type indicator BTYPE
      !> INTERNAL  : the global Poisson problem is solved, no need to impose BC's along mesh interfaces
      !> DIRICHLET : in the structured case face-wise BC-settings are used ccording to FFT-settings
      !>             (this also allows to use FFT as local preconditioner)
      !>             in the unstructured case Dirichlet BCs are only used for open boundary cells 
      !> NEUMANN   : is used for the rest
      !>
      IS_BC_DIRICHLET  = WC%PRESSURE_BC_INDEX == DIRICHLET
      IS_OPEN_BOUNDARY = WC%BOUNDARY_TYPE     == OPEN_BOUNDARY

      IF (EWC%NOM /= 0) THEN
         L%WALL(IWG)%BTYPE = INTERNAL
      ELSE IF ((    PRES_ON_WHOLE_DOMAIN .AND. IS_BC_DIRICHLET) .OR. &
              (.NOT.PRES_ON_WHOLE_DOMAIN .AND. IS_OPEN_BOUNDARY)) THEN
         L%WALL(IWG)%BTYPE = DIRICHLET
         L%N_DIRIC = L%N_DIRIC + 1
      ELSE
         L%WALL(IWG)%BTYPE = NEUMANN
         L%N_NEUMANN = L%N_NEUMANN + 1
      ENDIF

      L%WALL(IWG)%BOUNDARY_TYPE  = WC%BOUNDARY_TYPE

      L%WALL(IWG)%IXG =  WC%ONE_D%II                                 !> ghost cell indices
      L%WALL(IWG)%IYG =  WC%ONE_D%JJ
      L%WALL(IWG)%IZG =  WC%ONE_D%KK

      L%WALL(IWG)%IXW =  WC%ONE_D%IIG                                !> (internal) wall cell indices
      L%WALL(IWG)%IYW =  WC%ONE_D%JJG
      L%WALL(IWG)%IZW =  WC%ONE_D%KKG

      !> If there exists a neighbor for that wall cell, setup corresponding neighborship information
      IF (NOM /= 0) THEN
         CALL SCARC_SETUP_WALLCELL_NEIGHBOR(EWC%IIO_MIN, EWC%IIO_MAX, &
                                            EWC%JJO_MIN, EWC%JJO_MAX, &
                                            EWC%KKO_MIN, EWC%KKO_MAX, &
                                            IWG, IOR0, NM, NOM, NLEVEL_MIN)
      ENDIF

      NOM_LAST  = NOM
      IOR_LAST = IOR0

   ENDDO WALL_CELLS_LOOP2

ENDDO MESHES_LOOP1

CALL SCARC_SETUP_SUBDIVISION(NLEVEL_MIN)

!> 
!> Set DISCRET information on finest level and if requested also on coarser levels
!> 
CALL SCARC_SETUP_GLOBALS_UNSTRUCTURED(NLEVEL_MIN)
IF (.NOT.BAMG) THEN
   DO NL = NLEVEL_MIN+1, NLEVEL_MAX
      CALL SCARC_SETUP_DISCRETIZATION_LEVEL(NL)
      CALL SCARC_SETUP_GLOBALS_STRUCTURED(NL)
      CALL SCARC_SETUP_GLOBALS_UNSTRUCTURED(NL)
   ENDDO
ENDIF

!> Count number of Dirichlet BCs on finest level and, in case that there are none,
!> allocate mapping information for overlapped cells for later definition of condensed matrix system
LOCAL_INT = 0
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L => SCARC(NM)%LEVEL(NLEVEL_MIN)
   DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      LOCAL_INT(NM2) = SCARC(NM2)%LEVEL(NLEVEL_MIN)%N_DIRIC
   ENDDO
ENDDO
N_DIRIC_GLOBAL(NLEVEL_MIN) = SCARC_BROADCAST_INT(NSCARC_BROADCAST_SUM)

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM)%LEVEL(NLEVEL_MIN)
      M => L%MAP
      CALL SCARC_ALLOCATE_INT1 (L%MAP%ICE_TO_ICN, L%NCS+1, L%NCE, NSCARC_INIT_ZERO, 'ICE_TO_ICN')
      CALL SCARC_ALLOCATE_REAL1(L%MAP%ICE_TO_VAL, L%NCS+1, L%NCE, NSCARC_INIT_ZERO, 'ICE_TO_VAL')
   ENDDO
ENDIF

!>
!> Only in case of Twolevel-CG- or GMG-method (as main solver or preconditioner):
!> Determine WALL, FACE and OSCARC types for coarser levels
!> 
MULTI_LEVEL_IF: IF (BMULTILEVEL) THEN

   MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      S => SCARC(NM)

      IREFINE=1
      LEVEL_GMG_LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX

         LF => S%LEVEL(NL-1)                   !> level fine
         LC => S%LEVEL(NL)                     !> level coarse 
         MC => LC%MAP                          !> mappings on coarse level

         IREFINE=IREFINE*2

         CALL SCARC_CHECK_DIVISIBILITY(LF%NCE-LF%NCS, 'LF%NCE')
         CALL SCARC_CHECK_DIVISIBILITY(LF%NCO-LF%NCS, 'LF%NCO')

         LC%NCE = LC%NCS + (LF%NCE-LF%NCS)/2
         LC%NCO = LC%NCS + (LF%NCO-LF%NCS)/2

         LC%MAP%ICE_PTR = LC%NCS
         LC%MAP%ICO_PTR = LC%NCS

         LC%N_WALL_CELLS_EXT = SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
         LC%N_WALL_CELLS_INT = SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)

         LC%NW = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT
         ALLOCATE(LC%WALL(LC%NW), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_INTERFACES','WALL',IERROR)

         CALL SCARC_SETUP_MAPPINGS(NM, NL)

         !> compute FACE and WALL information for all faces of coarser level
         IWC = 1
         IWG = 1
         DO IFACE = 1, 6

            IOR0 = FACE_ORDER_XYZ(IFACE)

            !> compute mesh dimensions of coarser mesh level
            CALL SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)

            !> for every neighbor do:
            IF (LF%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
               DO INBR = 1, LF%FACE(IOR0)%N_NEIGHBORS

                  NOM = LF%FACE(IOR0)%NEIGHBORS(INBR)

                  OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL-1)

                  OLC%NCPL = OLF%NCPL

                  CALL SCARC_CHECK_DIVISIBILITY(OLF%NWL, 'OLF%NWL')
                  CALL SCARC_CHECK_DIVISIBILITY(OLF%NCG, 'OLF%NCG')

                  IF (.NOT.TWO_D) THEN
                     OLC%NWL = OLF%NWL/4
                     OLC%NCG = OLF%NCG/4
                  ELSE
                     OLC%NWL = OLF%NWL/2
                     OLC%NCG = OLF%NCG/2
                  ENDIF

                  CALL SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)
                  CALL SCARC_SETUP_OMAPPINGS(NM, NOM, NL)

               ENDDO
            ENDIF

            !> setup complete face information for coarser mesh
            CALL SCARC_SETUP_FACE(IOR0, IWC, IREFINE, NM, NL)

         ENDDO
         CALL SCARC_SETUP_CELL_INDEX(NM, NL)
         CALL SCARC_SETUP_INTERNAL_WALL_COORDS(NM, NL)
         CALL SCARC_SETUP_WALL_INDEX(NM, NL)

      ENDDO LEVEL_GMG_LEVEL_LOOP
   ENDDO MESHES_LOOP2
ENDIF MULTI_LEVEL_IF

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_WALLS

!> -----------------------------------------------------------------------------------------
!> --- Store neighbors of mesh
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_NEIGHBORS(NM, NOM)
INTEGER, INTENT(IN) :: NM, NOM
INTEGER:: INBR
DO INBR = 1, SCARC(NM)%N_NEIGHBORS
   IF (SCARC(NM)%NEIGHBORS(INBR) == NSCARC_UNDEFINED_INT) EXIT      !> not found, to be stored
   IF (SCARC(NM)%NEIGHBORS(INBR) == NOM) RETURN                     !> nothing to do, already stored
ENDDO
SCARC(NM)%NEIGHBORS(INBR) = NOM 
RETURN
END SUBROUTINE SCARC_UPDATE_NEIGHBORS


!> -----------------------------------------------------------------------------------------
!> --- Setup CELL_INDEX array on coarser grid levels in case of MG-method
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CELL_INDEX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: I, J, K, NO
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_OBST_TYPE)   , POINTER :: OB

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_CELL_INDEX')

L => SCARC(NM)%LEVEL(NL)

CALL SCARC_ALLOCATE_INT3(L%CELL%CINDEX, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'CELL_INDEX')
L%N_CELLS = 0

!>
!> Preset it for all grid cells
!>
DO K=0,L%NZ+1
   DO J=0,L%NY+1
      DO I=0,1
         IF (L%CELL%CINDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL%CINDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
      DO I=L%NX,L%NX+1
         IF (L%CELL%CINDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL%CINDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO K=0,L%NZ+1
   DO I=0,L%NX+1
      DO J=0,1
         IF (L%CELL%CINDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL%CINDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
      DO J=L%NY,L%NY+1
         IF (L%CELL%CINDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL%CINDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO J=0,L%NY+1
   DO I=0,L%NX+1
      DO K=0,1
         IF (L%CELL%CINDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL%CINDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
      DO K=L%NZ,L%NZ+1
         IF (L%CELL%CINDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL%CINDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
   ENDDO
ENDDO

!>
!> Consider cells in obstructions
!>
DO NO=1,L%N_OBST
   OB=>SCARC(NM)%LEVEL(NL)%OBST(NO)
   DO K=OB%K1,OB%K2+1
      DO J=OB%J1,OB%J2+1
         DO I=OB%I1,OB%I2+1
            IF (L%CELL%CINDEX(I,J,K)==0) THEN
               L%N_CELLS = L%N_CELLS + 1
               L%CELL%CINDEX(I,J,K) = L%N_CELLS
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_CELL_INDEX


!> -----------------------------------------------------------------------------------------
!> --- Setup WALL_INDEX array on coarser grid levels in case of MG-method
!> --- corresponding to cell index information
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_INDEX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: I, J, K, ICG, IW, IOR0
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_WALL_INDEX')

L => SCARC(NM)%LEVEL(NL)

CALL SCARC_ALLOCATE_INT2(L%CELL%WINDEX, 1, L%N_CELLS, -3, 3, NSCARC_INIT_ZERO, 'WINDEX')

DO IW = 1, L%NW

   I = L%WALL(IW)%IXW
   J = L%WALL(IW)%IYW
   K = L%WALL(IW)%IZW

   IOR0 = L%WALL(IW)%IOR
   ICG  = L%CELL%CINDEX(I,J,K)
  
   L%CELL%WINDEX(ICG,-IOR0) = IW

ENDDO

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_WALL_INDEX


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor in case of MG-method
!> Number of obstructions on coarse level is the same as on fine level
!> ACHTUNG: FUNKTIONIERT NUR FUER SPEZIALFAELLE, DIE AUCH FUER GMG LAUFEN !>!
!> ACHTUNG: MUSS DRINGEND NOCH ERWEITERT WERDEN
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERNAL_WALL_COORDS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, IO, IWC
INTEGER :: I, J, K
TYPE (SCARC_LEVEL_TYPE), POINTER :: LC, LF
TYPE (SCARC_CELL_TYPE) , POINTER :: CC, CF
TYPE (SCARC_OBST_TYPE) , POINTER :: OB

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_INTERNAL_WALL_COORDS')

LF => SCARC(NM)%LEVEL(NL-1)
LC => SCARC(NM)%LEVEL(NL)

CF => SCARC(NM)%LEVEL(NL-1)%CELL
CC => SCARC(NM)%LEVEL(NL)%CELL

IWC = LC%N_WALL_CELLS_EXT + 1

DO IO = 1, LC%N_OBST

   OB => LC%OBST(IO)

   !> Analyze IOR = 1
   I = OB%I1
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = CC%DOF(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         LC%WALL(IWC)%IXW = I+1
         LC%WALL(IWC)%IYW = J
         LC%WALL(IWC)%IZW = K
         LC%WALL(IWC)%IXG = I
         LC%WALL(IWC)%IYG = J
         LC%WALL(IWC)%IZG = K
         LC%WALL(IWC)%IOR = 1
         LC%WALL(IWC)%BTYPE = NEUMANN
         LC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = CC%DOF(I+1, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         LC%WALL(IWC)%IXW = I
         LC%WALL(IWC)%IYW = J
         LC%WALL(IWC)%IZW = K
         LC%WALL(IWC)%IXG = I+1
         LC%WALL(IWC)%IYG = J
         LC%WALL(IWC)%IZG = K
         LC%WALL(IWC)%IOR =-1
         LC%WALL(IWC)%BTYPE = NEUMANN
         LC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 2
   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = CC%DOF(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         LC%WALL(IWC)%IXW = I
         LC%WALL(IWC)%IYW = J+1
         LC%WALL(IWC)%IZW = K
         LC%WALL(IWC)%IXG = I
         LC%WALL(IWC)%IYG = J
         LC%WALL(IWC)%IZG = K
         LC%WALL(IWC)%IOR = 2
         LC%WALL(IWC)%BTYPE = NEUMANN
         LC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -2
   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = CC%DOF(I, J+1, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         LC%WALL(IWC)%IXW = I
         LC%WALL(IWC)%IYW = J
         LC%WALL(IWC)%IZW = K
         LC%WALL(IWC)%IXG = I
         LC%WALL(IWC)%IYG = J+1
         LC%WALL(IWC)%IZG = K
         LC%WALL(IWC)%IOR =-2
         LC%WALL(IWC)%BTYPE = NEUMANN
         LC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 3
   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = CC%DOF(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         LC%WALL(IWC)%IXW = I
         LC%WALL(IWC)%IYW = J
         LC%WALL(IWC)%IZW = K+1
         LC%WALL(IWC)%IXG = I
         LC%WALL(IWC)%IYG = J
         LC%WALL(IWC)%IZG = K
         LC%WALL(IWC)%IOR = 3
         LC%WALL(IWC)%BTYPE = NEUMANN
         LC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -3
   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = CC%DOF(I, J, K+1)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         LC%WALL(IWC)%IXW = I
         LC%WALL(IWC)%IYW = J
         LC%WALL(IWC)%IZW = K
         LC%WALL(IWC)%IXG = I
         LC%WALL(IWC)%IYG = J
         LC%WALL(IWC)%IZG = K+1
         LC%WALL(IWC)%IOR =-3
         LC%WALL(IWC)%BTYPE = NEUMANN
         LC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO
         
ENDDO

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_INTERNAL_WALL_COORDS


!> -------------------------------------------------------------------------------------------------
!> Correct BTYPE related to internal obstructions on ghost cells
!> ACHTUNG: NOCHMAL UEBERARBEITEN !
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWG
!INTEGER, POINTER, DIMENSION(:,:,:):: CELL_INDEX
!INTEGER, POINTER, DIMENSION(:,:)  :: WALL_INDEX
INTEGER :: IX, IY, IZ, IOR0, BTYPE0
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_INTERNAL_NEUMANNS')
L => SCARC(NM)%LEVEL(NL)

!IF (NL == NLEVEL_MIN) THEN
!   CELL_INDEX => MESHES(NM)%CELL_INDEX
!   WALL_INDEX => MESHES(NM)%WALL_INDEX
!ELSE
!   CELL_INDEX => L%CELL%CINDEX
!   WALL_INDEX => L%CELL%WINDEX
!ENDIF

DO IWG = 1, L%N_WALL_CELLS_EXT

   ! ACHTUNG ACHTUNG ACHTUNG: HIER NOCHMAL SCHAUEN GLEICH ODER UNGLEICH NULL????????
   IF (L%WALL(IWG)%NOM == 0) CYCLE

   IX = L%WALL(IWG)%IXW
   IY = L%WALL(IWG)%IYW
   IZ = L%WALL(IWG)%IZW

   IOR0   = L%WALL(IWG)%IOR
   BTYPE0 = L%WALL(IWG)%BTYPE

  ! ICG = CELL_INDEX(IX, IY, IZ)
  ! IWG = WALL_INDEX(ICG, IOR0)
  ! IF (L%WALL(IWG)%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY) L%WALL(IWG)%BTYPE=NEUMANN
  ! IF (L%CELL%STATE(IX, IY, IZ) == NSCARC_DISCRET_SOLID) L%WALL(IWG)%BTYPE=NEUMANN

   IF (L%WALL(IWG)%BOUNDARY_TYPE == SOLID_BOUNDARY) L%WALL(IWG)%BTYPE=NEUMANN

ENDDO

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IXC, IYC, IZC
INTEGER :: IXF, IYF, IZF
INTEGER :: IWC, ICF(4)=0, IWF(4)=0, IOR0
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
INTEGER, POINTER, DIMENSION(:,:)   :: WALL_INDEX
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: LF , LC

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_EXTERNAL_WALL_CELLS')

LC => SCARC(NM)%LEVEL(NL)
LF => SCARC(NM)%LEVEL(NL-1)

IF (NL == NLEVEL_MIN+1) THEN
   CELL_INDEX => MESHES(NM)%CELL_INDEX
   WALL_INDEX => MESHES(NM)%WALL_INDEX
ELSE
   CELL_INDEX => LF%CELL%CINDEX
   WALL_INDEX => LF%CELL%WINDEX
ENDIF

ICF = 0
IWC = 0
IWF = 0

IF (TWO_D) THEN
   IYC = 1
   IYF = 1

   !> IOR = 1
   IOR0 = 1
   DO IZC = 1, LC%NZ 
      IZF = 2*IZC - 1
      ICF(1) = CELL_INDEX(1  , IYF  , IZF  )
      ICF(2) = CELL_INDEX(1  , IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
   ENDDO

   !> IOR = -1
   IOR0 = -1
   DO IZC = 1, LC%NZ 
      IZF = 2*IZC - 1
      ICF(1) = CELL_INDEX(LF%NX, IYF  , IZF  )
      ICF(2) = CELL_INDEX(LF%NX, IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
   ENDDO

   !> IOR = 2
   IOR0 = 2
   DO IZC = 1, LC%NZ 
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF  , IYF, IZF  )
         ICF(2) = CELL_INDEX(IXF+1, IYF, IZF  )
         ICF(3) = CELL_INDEX(IXF  , IYF, IZF+1)
         ICF(4) = CELL_INDEX(IXF+1, IYF, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = -2
   IOR0 = -2
   IXF  = 1
   DO IZC = 1, LC%NZ 
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF    , IYF, IZF  )
         ICF(2) = CELL_INDEX(IXF+1  , IYF, IZF  )
         ICF(3) = CELL_INDEX(IXF    , IYF, IZF+1)
         ICF(4) = CELL_INDEX(IXF+1  , IYF, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = 3
   IOR0 = 3
   DO IXC = 1, LC%NX 
      IXF = 2*IXC - 1
      ICF(1) = CELL_INDEX(IXF    , IYF  , 1)
      ICF(2) = CELL_INDEX(IXF+1  , IYF  , 1)
      IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
   ENDDO

   !> IOR = -3
   IOR0 = -3
   DO IXC = 1, LC%NX 
      IXF = 2*IXC - 1
      ICF(1) = CELL_INDEX(IXF  , IYF  , LF%NZ)
      ICF(2) = CELL_INDEX(IXF+1, IYF  , LF%NZ)
      IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
   ENDDO

ELSE

   !> IOR = 1
   IOR0 = 1
   DO IZC = 1, LC%NZ 
      IZF = 2*IZC - 1
      DO IYC = 1, LC%NY 
         IYF = 2*IYC - 1
         ICF(1) = CELL_INDEX(1  , IYF  , IZF  )
         ICF(2) = CELL_INDEX(1  , IYF+1, IZF  )
         ICF(3) = CELL_INDEX(1  , IYF  , IZF+1)
         ICF(4) = CELL_INDEX(1  , IYF+1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = -1
   IOR0 = -1
   DO IZC = 1, LC%NZ 
      IZF = 2*IZC - 1
      DO IYC = 1, LC%NY 
         IYF = 2*IYC - 1
         ICF(1) = CELL_INDEX(LF%NX, IYF  , IZF  )
         ICF(2) = CELL_INDEX(LF%NX, IYF+1, IZF  )
         ICF(3) = CELL_INDEX(LF%NX, IYF  , IZF+1)
         ICF(4) = CELL_INDEX(LF%NX, IYF+1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = 2
   IOR0 = 2
   DO IZC = 1, LC%NZ 
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF  , 1, IZF  )
         ICF(2) = CELL_INDEX(IXF+1, 1, IZF  )
         ICF(3) = CELL_INDEX(IXF  , 1, IZF+1)
         ICF(4) = CELL_INDEX(IXF+1, 1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = -2
   IOR0 = -2
   IXF  = 1
   DO IZC = 1, LC%NZ 
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF    , LF%NY, IZF  )
         ICF(2) = CELL_INDEX(IXF+1  , LF%NY, IZF  )
         ICF(3) = CELL_INDEX(IXF    , LF%NY, IZF+1)
         ICF(4) = CELL_INDEX(IXF+1  , LF%NY, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = 3
   IOR0 = 3
   DO IYC = 1, LC%NY 
      IYF = 2*IYC - 1
      DO IXC = 1, LC%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF    , IYF  , 1)
         ICF(2) = CELL_INDEX(IXF+1  , IYF  , 1)
         ICF(3) = CELL_INDEX(IXF    , IYF+1, 1)
         ICF(4) = CELL_INDEX(IXF+1  , IYF+1, 1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = -3
   IOR0 = -3
   DO IYC = 1, LC%NY 
      IYF = 2*IYC - 1
      DO IXC = 1, LC%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF  , IYF  , LF%NZ)
         ICF(2) = CELL_INDEX(IXF+1, IYF  , LF%NZ)
         ICF(3) = CELL_INDEX(IXF  , IYF+1, LF%NZ)
         ICF(4) = CELL_INDEX(IXF+1, IYF+1, LF%NZ)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

ENDIF

SCARC_COUNT_EXTERNAL_WALL_CELLS = IWC
CALL SCARC_LEAVE_ROUTINE()
RETURN
END FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS

!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWC, NW_INT
INTEGER :: IC, IO
INTEGER :: I, J, K
TYPE (SCARC_CELL_TYPE) , POINTER :: CF, CC
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF, LC
TYPE (SCARC_OBST_TYPE) , POINTER :: OB

CALL SCARC_ENTER_ROUTINE('SCARC_COUNT_INTERNAL_WALL_CELLS')

LF  => SCARC(NM)%LEVEL(NL-1)
LC  => SCARC(NM)%LEVEL(NL)

CF  => SCARC(NM)%LEVEL(NL-1)%CELL
CC  => SCARC(NM)%LEVEL(NL)%CELL

LC%N_OBST = LF%N_OBST                   !> Number of obstructions is the same on all levels

ALLOCATE(LC%OBST(LC%N_OBST), STAT=IERROR)
CALL ChkMemErr('SCARC_COUNT_INTERNAL_WALL_CELLS','OBST',IERROR)

NW_INT = 0
IWC = LC%N_WALL_CELLS_EXT + 1

DO IO = 1, LF%N_OBST

   OB => LC%OBST(IO)                    !> obstruction pointer on coarse level

   OB%I1 = (LF%OBST(IO)%I1+1)/2
   OB%I2 =  LF%OBST(IO)%I2/2

   IF (TWO_D) THEN
      OB%J1 = 0
      OB%J2 = 1
   ELSE
      OB%J1 = (LF%OBST(IO)%J1+1)/2
      OB%J2 =  LF%OBST(IO)%J2/2
   ENDIF

   OB%K1 = (LF%OBST(IO)%K1+1)/2
   OB%K2 =  LF%OBST(IO)%K2/2
   
   !> Analyze IOR = 1
   I = OB%I1
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = CC%DOF(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = CC%DOF(I+1, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 2
   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = CC%DOF(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -2
   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = CC%DOF(I, J+1, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 3
   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = CC%DOF(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -3
   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = CC%DOF(I, J, K+1)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO
         
ENDDO

SCARC_COUNT_INTERNAL_WALL_CELLS = NW_INT
CALL SCARC_LEAVE_ROUTINE()

END FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS



!> -------------------------------------------------------------------------------------------------
!> Count external wall cells on face IOR
!> -------------------------------------------------------------------------------------------------
LOGICAL FUNCTION IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, NCNT, NM, NL)
INTEGER, INTENT(IN) :: IOR0, NM, NL, NCNT
INTEGER, DIMENSION(:), INTENT(IN) :: ICF
INTEGER, POINTER, DIMENSION(:,:), INTENT(IN) :: WALL_INDEX
INTEGER:: I, IWF_LAST, IWF(4)=0
REAL(EB) :: BSUM
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

CALL SCARC_ENTER_ROUTINE('IS_EXTERNAL_WALLCELL')

IS_EXTERNAL_WALLCELL = .FALSE.

L => SCARC(NM)%LEVEL(NL)

DO I = 1, NCNT
   IWF(I) = WALL_INDEX(ICF(I), -IOR0)
ENDDO

BSUM = 0.0_EB
IWF_LAST = 0

DO I = 1, NCNT
   IF (IWF(I)>0) THEN
      BSUM = BSUM + REAL(L%WALL(IWF(I))%BTYPE,EB)
      IWF_LAST = IWF(I)
   ENDIF
ENDDO

IF (IWF_LAST == 0) THEN
   CALL SCARC_LEAVE_ROUTINE()
   RETURN
ENDIF

IF (ABS(BSUM/REAL(NCNT,EB) - REAL(L%WALL(IWF_LAST)%BTYPE,EB)) < 1E-12) THEN
   IS_EXTERNAL_WALLCELL = .TRUE.
   CALL SCARC_LEAVE_ROUTINE()
   RETURN
ELSE
   CALL SCARC_SHUTDOWN('Wrong boundary-type sum for IOR0 = ', 'NONE', IOR0)
ENDIF

END FUNCTION IS_EXTERNAL_WALLCELL


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLCELL_NEIGHBOR(NX1, NX2, NY1, NY2, NZ1, NZ2, IWG, IOR0, NM, NOM, NL)
INTEGER, INTENT(IN) :: NX1, NX2, NY1, NY2, NZ1, NZ2
INTEGER, INTENT(IN) :: IWG, IOR0, NM, NOM, NL
INTEGER :: NOMX, NOMY, NOMZ
INTEGER :: ICG, ICO, ICE, IWL, ICPL, IX, IY, IZ, JL
TYPE (SCARC_LEVEL_TYPE),   POINTER :: L, OL
TYPE (SCARC_MAPPING_TYPE), POINTER :: M, OM

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_WALLCELL_NEIGHBOR')

L  => SCARC(NM)%LEVEL(NL)
M  => L%MAP

OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
OM => OL%MAP

!> store basic information about neighbor and orientation
OL%IOR  = IOR0
OL%NOM  = NOM

ICO  = L%MAP%ICO_PTR
ICE  = L%MAP%ICE_PTR

ICG  = OL%MAP%ICG_PTR
IWL  = OL%MAP%IWL_PTR

OL%NCPLS = OL%NCPL
IF (OL%NCG == OL%NWL) THEN
   OL%NCPLR = OL%NCPL
ELSE IF (OL%NCG == 2*OL%NWL) THEN
   OL%NCPLR =  1
ELSE IF (OL%NWL == 2*OL%NCG) THEN
   OL%NCPLR =  2
ENDIF

!> set neighboring coordinates
L%WALL(IWG)%IXN(1) = NX1
L%WALL(IWG)%IXN(2) = NX2
L%WALL(IWG)%IYN(1) = NY1
L%WALL(IWG)%IYN(2) = NY2
L%WALL(IWG)%IZN(1) = NZ1
L%WALL(IWG)%IZN(2) = NZ2

!> allocate pointer arrays for extended, ghost and neighboring cells
CALL SCARC_ALLOCATE_INT1(L%WALL(IWG)%ICE, 1, OL%NCPL, NSCARC_INIT_UNDEFINED, 'ICE')
CALL SCARC_ALLOCATE_INT1(L%WALL(IWG)%ICG, 1, OL%NCPL, NSCARC_INIT_UNDEFINED, 'ICG')

IWL = IWL + 1
ICO = ICO + 1

IF (OL%SUBDIVISION(1, IOR0) == NSCARC_ZERO_INT) OL%SUBDIVISION(1, IOR0) =  IWL
OL%SUBDIVISION(2, IOR0) = OL%SUBDIVISION(2, IOR0) + 1

NOMX = MESHES(NOM)%IBAR
NOMY = MESHES(NOM)%JBAR
NOMZ = MESHES(NOM)%KBAR

IF (NL > 1) THEN
   DO JL = 2, NL
      NOMX = MESHES(NOM)%IBAR/NL
      IF (.NOT.TWO_D) NOMY = MESHES(NOM)%JBAR/NL
      NOMZ = MESHES(NOM)%KBAR/NL
   ENDDO
ENDIF

!> store information about overlapped cells and set mapping arrays
ICPL = 0
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ICPL = ICPL + 1
         ICG  = ICG  + 1
         ICE  = ICE  + 1

         L%WALL(IWG)%ICE(ICPL) = ICE                    !> number of extended grid cell
         L%WALL(IWG)%ICG(ICPL) = ICG                    !> number of ghost grid cell
   
         M%ICE_TO_IWG(ICE)  = IWG                       !> map extended cell to global wall cell
         M%ICE_TO_IWL(ICE)  = IWL                       !> map extended cell to local wall cell
         M%ICE_TO_ICG(ICE)  = ICG                       !> map extended cell to ghost cell
         M%ICE_TO_ICG(ICE)  = ICG                       !> map extended cell to ghost cell

         OM%ICG_TO_IWG(ICG) = IWG                       !> map ghost cell to global wall cell
         OM%ICG_TO_ICO(ICG) = ICO                       !> map ghost cell to extended grid cell
         OM%ICG_TO_ICE(ICG) = L%WALL(IWG)%ICE(ICPL)     !> map ghost cell to extended grid cell

      ENDDO
   ENDDO
ENDDO

IF (BAMG) THEN
   DO IZ = NZ1, NZ2
      DO IY = NY1, NY2
         DO IX = NX1, NX2
            OM%IWL_TO_ICG(IWL, ICPL) = ICG              !> map extended cell to ghost cell
         ENDDO
      ENDDO
   ENDDO
ENDIF

L%WALL(IWG)%ICO = ICO                                   !> number of overlapping cell
L%WALL(IWG)%IWL = IWL                                   !> number of local wall cell

M%ICO_PTR  = ICO                                        !> store overlapping cell counter
M%ICE_PTR  = ICE                                        !> store extended cell counter

OM%IWL_PTR = IWL                                        !> store local wall cell pointer
OM%ICG_PTR = ICG                                        !> store ghost cell counter

OM%IWL_TO_IWG(IWL) = IWG                                !> map local wall cell to global wall cell
OM%IWL_TO_ICO(IWL) = ICO                                !> map local wall cell to internal grid cell
!OM%IWL_TO_ICW(IWL) = L%WALL(IWG)%ICW                   !> map local wall cell to internal grid cell (AMG only)

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_WALLCELL_NEIGHBOR


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MAPPINGS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_MAPPING_TYPE), POINTER :: M

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_MAPPINGS')

L => SCARC(NM)%LEVEL(NL)
M => L%MAP

!IF (NMESHES > 1) THEN
   CALL SCARC_ALLOCATE_INT1(L%MAP%ICE_TO_IWG, L%NCS+1, L%NCE, NSCARC_INIT_ZERO, 'ICE_TO_IWG')
   CALL SCARC_ALLOCATE_INT1(L%MAP%ICE_TO_IWL, L%NCS+1, L%NCE, NSCARC_INIT_ZERO, 'ICL_TO_IWG')
   CALL SCARC_ALLOCATE_INT1(L%MAP%ICE_TO_ICG, L%NCS+1, L%NCE, NSCARC_INIT_ZERO, 'ICE_TO_IICG')
!ENDIF

IF (BAMG) CALL SCARC_ALLOCATE_INT1(L%AMG%INTERN, 1, L%NCS, NSCARC_INIT_ZERO, 'BDRY')

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_MAPPINGS


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OMAPPINGS(NM, NOM, NL)
INTEGER, INTENT(IN) :: NM, NOM, NL
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: OL
TYPE (SCARC_MAPPING_TYPE), POINTER :: OM

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_OMAPPINGS')

OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
OM => OL%MAP

CALL SCARC_ALLOCATE_INT1(OM%IWL_TO_IWG, 1, OL%NWL, NSCARC_INIT_ZERO, 'IWL_TO_IWG')
!CALL SCARC_ALLOCATE_INT1(OM%IWL_TO_ICW, 1, OL%NWL, NSCARC_INIT_ZERO, 'IWL_TO_ICW')
CALL SCARC_ALLOCATE_INT1(OM%IWL_TO_ICO, 1, OL%NWL, NSCARC_INIT_ZERO, 'IWL_TO_ICO')

CALL SCARC_ALLOCATE_INT1(OM%ICG_TO_IWG, 1, OL%NCG, NSCARC_INIT_ZERO, 'ICG_TO_IWG')
CALL SCARC_ALLOCATE_INT1(OM%ICG_TO_ICO, 1, OL%NCG, NSCARC_INIT_ZERO, 'ICG_TO_ICO')
CALL SCARC_ALLOCATE_INT1(OM%ICG_TO_ICE, 1, OL%NCG, NSCARC_INIT_ZERO, 'ICG_TO_ICE')

IF (BAMG) CALL SCARC_ALLOCATE_INT2(OM%IWL_TO_ICG, 1, OL%NWL, 1, OL%NCPL, NSCARC_INIT_ZERO,'IWL_TO_ICG')

!IF (PRES_ON_WHOLE_DOMAIN.AND.OL%NCS>0) &         !> ACHTUNG: kann verbessert werden??? riesig!
!   CALL SCARC_ALLOCATE_INT1(OM%ICN_TO_ICE, 1, OL%NCS, NSCARC_INIT_ZERO,'ICN_TO_ICE')

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_OMAPPINGS


!> ----------------------------------------------------------------------------------------------------
!> Check divisibility by 2 of a given number of elements (in one grid direction)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_DIVISIBILITY(NN, CDIR)
INTEGER, INTENT(IN) :: NN
CHARACTER(*) , INTENT(IN) :: CDIR
CALL SCARC_ENTER_ROUTINE('SCARC_CHECK_DIVISBILITY')
IF (MOD(NN,2) /= 0) CALL SCARC_SHUTDOWN('Parameter not divisable by 2 for ', CDIR, -999)
CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_CHECK_DIVISIBILITY


!> ----------------------------------------------------------------------------------------------------
!> Set wall cell information on coarse level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)
INTEGER, INTENT(IN)    :: IOR0, NM, NL
INTEGER, INTENT(INOUT) :: IWG
INTEGER:: INBR
TYPE (SCARC_LEVEL_TYPE), POINTER :: LC, LF            !> LEVEL types for coarse and fine 
TYPE (SCARC_FACE_TYPE) , POINTER :: FC, FF            !> FACE types on coarse and fine level

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_FACE_DIMENSIONS')

!> reference coarse and fine LEVEL type
LC => SCARC(NM)%LEVEL(NL)
LF => SCARC(NM)%LEVEL(NL-1)

!> reference coarse and fine FACE type
FC => SCARC(NM)%LEVEL(NL)%FACE(IOR0)
FF => SCARC(NM)%LEVEL(NL-1)%FACE(IOR0)

!> initialize FACE type for coarser mesh
FC%IWG_PTR = IWG

FC%N_NEIGHBORS = FF%N_NEIGHBORS
IF (FC%N_NEIGHBORS /= 0) CALL SCARC_ALLOCATE_INT1(FC%NEIGHBORS, 1, FC%N_NEIGHBORS, NSCARC_INIT_NONE, 'FACE_NEIGHBORS')
DO INBR= 1, FC%N_NEIGHBORS
   FC%NEIGHBORS(INBR) = FF%NEIGHBORS(INBR)
ENDDO

SELECT CASE (ABS(IOR0))

   CASE (1)
      FC%DH => LC%COORD%DXL
      FC%NFX = 1                                              !> no extension in x-direction
      IF (.NOT.TWO_D) THEN                                    !> only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z')              !> number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2
      FC%NFC = LC%NX                                           !> number of cells between opposite faces
   CASE (2)
      FC%DH => LC%COORD%DYL
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X')              !> number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2
      FC%NFY = 1                                              !> no extension in y-direction
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z')              !> number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2
      FC%NFC = LC%NY                                           !> number of cells between opposite faces
   CASE (3)
      FC%DH => LC%COORD%DZL
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X')              !> number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2
      IF (.NOT.TWO_D) THEN                                    !> only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF
      FC%NFZ = 1                                              !> no extension in y-direction
      FC%NFC = LC%NZ                                           !> number of cells between opposite faces
END SELECT

FC%NFW = FC%NFX * FC%NFY * FC%NFZ                             !> get number of wall cells for that face
IWG = IWG + FC%NFW                                            !> increase global wall cell counter

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_FACE_DIMENSIONS


!> ----------------------------------------------------------------------------------------------------
!> Set wall cell information on coarse level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE(IOR0, IWC, IREFINE, NM, NL)
INTEGER, INTENT(INOUT) :: IWC
INTEGER, INTENT(IN) :: NM, NL
INTEGER, INTENT(IN) :: IOR0, IREFINE
INTEGER :: IWF(4) , IBCF(4), NOMF(4)
INTEGER :: I, NCPL
INTEGER :: IX,  IY,  IZ
INTEGER :: NX1, NY1, NZ1
INTEGER :: NX2, NY2, NZ2
INTEGER :: IX1, IY1, IZ1
INTEGER :: IX2, IY2, IZ2
INTEGER :: IDIFF, JDIFF, KDIFF
TYPE (SCARC_WALL_TYPE) , POINTER, DIMENSION(:) :: WC, WF         !> coarse and fine WALL types
TYPE (SCARC_FACE_TYPE) , POINTER               :: FF, FC         !> coarse and fine FACE types
TYPE (SCARC_LEVEL_TYPE), POINTER               :: LC, LF, OLC    !> coarse and fine LEVEL types

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_FACE')

!> reference coarse and fine LEVEL type
LC => SCARC(NM)%LEVEL(NL)
LF => SCARC(NM)%LEVEL(NL-1)

!> reference coarse and fine WALL type
WC => SCARC(NM)%LEVEL(NL)%WALL
WF => SCARC(NM)%LEVEL(NL-1)%WALL

!> reference coarse and fine FACE type
FC => SCARC(NM)%LEVEL(NL)%FACE(IOR0)
FF => SCARC(NM)%LEVEL(NL-1)%FACE(IOR0)

!> set coordinate dimensions for correspoding face
SELECT CASE (ABS(IOR0))
   CASE (1)
      IF (IOR0 > 0) THEN                                           !> set dimensions for wall cell counting
         NX1 = 0
         NX2 = 0
      ELSE
         NX1 = LC%NX+1
         NX2 = LC%NX+1
      ENDIF
      NY1 = 1
      NY2 = LC%NY
      NZ1 = 1
      NZ2 = LC%NZ
   CASE (2)
      NX1 = 1                                                      !> set dimensions for wall cell counting
      NX2 = LC%NX
      IF (IOR0 > 0) THEN
         NY1 = 0
         NY2 = 0
      ELSE
         NY1 = LC%NY+1
         NY2 = LC%NY+1
      ENDIF
      NZ1 = 1
      NZ2 = LC%NZ
   CASE (3)
      NX1 = 1                                                      !> set dimensions for wall cell counting
      NX2 = LC%NX
      NY1 = 1
      NY2 = LC%NY
      IF (IOR0 > 0) THEN
         NZ1 = 0
         NZ2 = 0
      ELSE
         NZ1 =LC%NZ+1
         NZ2 =LC%NZ+1
      ENDIF
END SELECT

!>
!> Loop over all wall cells of face IOR0
!>
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         !> Set orientation of neiboring face, indices of ghost and adjacent cell for coarse IW
         WC(IWC)%IOR = IOR0

         SELECT CASE (IOR0)
            CASE (1)
               WC(IWC)%ICW = (IZ-1)*LC%NX*LC%NY + (IY-1)*LC%NX + IX + 1
            CASE (-1)
               WC(IWC)%ICW = (IZ-1)*LC%NX*LC%NY + (IY-1)*LC%NX + IX - 1
            CASE (2)
               WC(IWC)%ICW = (IZ-1)*LC%NX*LC%NY +  IY   *LC%NX + IX
            CASE (-2)
               WC(IWC)%ICW = (IZ-1)*LC%NX*LC%NY + (IY-2)*LC%NX + IX
            CASE (3)
               WC(IWC)%ICW =  IZ   *LC%NX*LC%NY + (IY-1)*LC%NX + IX
            CASE (-3)
               WC(IWC)%ICW = (IZ-2)*LC%NX*LC%NY + (IY-1)*LC%NX + IX
         END SELECT

         WC(IWC)%IOR = IOR0

         WC(IWC)%IXG = IX
         WC(IWC)%IYG = IY
         WC(IWC)%IZG = IZ

         SELECT CASE (IOR0)
            CASE (1)
               WC(IWC)%IXW = IX+1
               WC(IWC)%IYW = IY
               WC(IWC)%IZW = IZ
            CASE (-1)
               WC(IWC)%IXW = IX-1
               WC(IWC)%IYW = IY
               WC(IWC)%IZW = IZ
            CASE (2)
               WC(IWC)%IXW = IX
               WC(IWC)%IYW = IY+1
               WC(IWC)%IZW = IZ
            CASE (-2)
               WC(IWC)%IXW = IX
               WC(IWC)%IYW = IY-1
               WC(IWC)%IZW = IZ
            CASE (3)
               WC(IWC)%IXW = IX
               WC(IWC)%IYW = IY
               WC(IWC)%IZW = IZ+1
            CASE (-3)
               WC(IWC)%IXW = IX
               WC(IWC)%IYW = IY
               WC(IWC)%IZW = IZ-1
         END SELECT

         !> ------------------------------------------------------------
         !> 2D-version
         !> ------------------------------------------------------------
         IF (TWO_D) THEN

            !> determine fine IW's, which must be merged to one coarse IW
            SELECT CASE (ABS(IOR0))
               CASE ( 1)
                  IWF(1) = FF%IWG_PTR + 2*(IZ-1)
               CASE ( 2)
                  IWF(1) = FF%IWG_PTR + 2*(IZ-1)*LF%NX + 2*(IX - 1)
               CASE ( 3)
                  IWF(1) = FF%IWG_PTR + 2*(IX-1)
            END SELECT
            IWF(2) = IWF(1)+1

            !> set fine cell neighbors (they must be the same for all fine IW's)
            NOMF(1) = WF(IWF(1))%NOM
            NOMF(2) = WF(IWF(2))%NOM
            IF (NOMF(1) /= NOMF(2)) CALL SCARC_SHUTDOWN('Inconsistent neighbors ', 'NONE', NOMF(1))

            WC(IWC)%NOM = NOMF(1)

            !> set corresponding pressure_bc_index on coarser level
            IBCF(1) = WF(IWF(1))%BTYPE
            IBCF(2) = WF(IWF(2))%BTYPE
            IF (IBCF(1) == INTERNAL .OR. IBCF(2) == INTERNAL) THEN
               WC(IWC)%BTYPE = INTERNAL
            ELSE IF (IBCF(1) == DIRICHLET .OR. IBCF(2) == DIRICHLET) THEN
               WC(IWC)%BTYPE = DIRICHLET
            ELSE
               WC(IWC)%BTYPE = NEUMANN
            ENDIF

            !> set corresponding pressure_bc_index on coarser level
            IBCF(1) = WF(IWF(1))%BOUNDARY_TYPE
            IBCF(2) = WF(IWF(2))%BOUNDARY_TYPE
            IF (IBCF(1)==NULL_BOUNDARY .OR. IBCF(2)==NULL_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = NULL_BOUNDARY
            ELSE IF (IBCF(1)==SOLID_BOUNDARY .OR. IBCF(2)==SOLID_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
            ELSE IF (IBCF(1)==OPEN_BOUNDARY .OR. IBCF(2)==OPEN_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = OPEN_BOUNDARY
            ELSE IF (IBCF(1)==INTERPOLATED_BOUNDARY .OR. IBCF(2)==INTERPOLATED_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = INTERPOLATED_BOUNDARY
            ELSE IF (IBCF(1)==MIRROR_BOUNDARY .OR. IBCF(2)==MIRROR_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = MIRROR_BOUNDARY
            ELSE
               CALL SCARC_SHUTDOWN('Wrong boundary type ', 'NONE', IBCF(1))
            ENDIF

            !> in case of an internal boundary set neighboring WALL cells
            IF (NOMF(1) > 0) THEN

               OLC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)

               IY1 = 1
               IY2 = 1
               SELECT CASE (ABS(IOR0))
                  CASE (1)
                     KDIFF = WF(IWF(2))%IZN(1) - WF(IWF(1))%IZN(1)
                     IF (KDIFF == 1) THEN
                        IZ1 = WF(IWF(2))%IZN(2)/2
                        IZ2 = IZ1
                     ELSE IF (KDIFF == 2) THEN
                        IZ1 = WF(IWF(1))%IZN(2)/2
                        IZ2 = WF(IWF(2))%IZN(2)/2
                     ELSE IF (KDIFF == 0) THEN
                        IZ1 = (WF(IWF(1))%IZN(2)+1)/2
                        IZ2 = IZ1
                     ELSE
                        CALL SCARC_SHUTDOWN('Wrong resolution for IOR0 = ', 'NONE', IOR0)
                     ENDIF
                  CASE (3)
                     IDIFF = WF(IWF(2))%IXN(1) - WF(IWF(1))%IXN(1)
                     IF (IDIFF == 1) THEN
                        IX1 = WF(IWF(2))%IXN(2)/2
                        IX2 = IX1
                     ELSE IF (IDIFF == 2) THEN
                        IX1 = WF(IWF(1))%IXN(2)/2
                        IX2 = WF(IWF(2))%IXN(2)/2
                     ELSE IF (IDIFF == 0) THEN
                        IX1 = (WF(IWF(1))%IXN(2)+1)/2
                        IX2 = IX1
                     ELSE
                        CALL SCARC_SHUTDOWN('Wrong resolution for IOR0 = ', 'NONE', IOR0)
                     ENDIF
               END SELECT

               SELECT CASE (IOR0)
                  CASE (1)
                     IX1 = MESHES(NOMF(1))%IBAR/IREFINE
                     IX2 = IX1
                  CASE (-1)
                     IX1 = 1
                     IX2 = 1
                  CASE (3)
                     IZ1 = MESHES(NOMF(1))%KBAR/IREFINE
                     IZ2 = IZ1
                  CASE (-3)
                     IZ1 = 1
                     IZ2 = 1
               END SELECT

               WC(IWC)%IXN(1) = IX1
               WC(IWC)%IYN(1) = 1
               WC(IWC)%IZN(1) = IZ1
               WC(IWC)%IXN(2) = IX2
               WC(IWC)%IYN(2) = 1
               WC(IWC)%IZN(2) = IZ2

               !>!
               !> Allocate and specify ICN and ICE arrays for OC
               !>!
               NCPL = (IZ2-IZ1+1)*(IX2-IX1+1)
               OLC%NCPL = NCPL

               CALL SCARC_SETUP_WALLCELL_NEIGHBOR(IX1, IX2, 1, 1, IZ1, IZ2, IWC, IOR0, NM, NOMF(1), NL)

            ENDIF


         !> ------------------------------------------------------------
         !> 3D-version
         !> ------------------------------------------------------------
         ELSE

            !> determine fine IW's, which must be merged to one coarse IW
            SELECT CASE (ABS(IOR0))
               CASE (1)
                  IWF(1) = FF%IWG_PTR + (2*IZ-2)*LF%NY + 2*IY - 2
                  IWF(3) = FF%IWG_PTR + (2*IZ-1)*LF%NY + 2*IY - 2
               CASE (2)
                  IWF(1) = FF%IWG_PTR + (2*IZ-2)*LF%NX + 2*IX - 2
                  IWF(3) = FF%IWG_PTR + (2*IZ-1)*LF%NX + 2*IX - 2
               CASE (3)
                  IWF(1) = FF%IWG_PTR + (2*IY-2)*LF%NX + 2*IX - 2
                  IWF(3) = FF%IWG_PTR + (2*IY-1)*LF%NX + 2*IX - 2
            END SELECT
            IWF(2) = IWF(1)+1
            IWF(4) = IWF(3)+1

            !> set fine cell neighbors (they must be the same for all fine IW's)
            DO I=1,4
               NOMF(I) = WF(IWF(I))%NOM
            ENDDO

            IF (NOMF(1)/=NOMF(2) .OR. NOMF(1)/=NOMF(3) .OR. NOMF(1)/=NOMF(4)) &
               CALL SCARC_SHUTDOWN('Inconsistent neighbors for IOR0 = ', 'NONE', IOR0)
            WC(IWC)%NOM = NOMF(1)

            !> set corresponding pressure_bc_index on coarser level
            DO I=1,4
               IBCF(I) = WF(IWF(I))%BTYPE
            ENDDO
            IF (IBCF(1)==INTERNAL.OR.IBCF(2)==INTERNAL.OR.&
                IBCF(3)==INTERNAL.OR.IBCF(4)==INTERNAL) THEN
               WC(IWC)%BTYPE =INTERNAL
            ELSE IF (IBCF(1)==DIRICHLET.OR.IBCF(2)==DIRICHLET.OR.&
                     IBCF(3)==DIRICHLET.OR.IBCF(4)==DIRICHLET) THEN
               WC(IWC)%BTYPE =DIRICHLET
            ELSE
               WC(IWC)%BTYPE =NEUMANN
            ENDIF

            !> set corresponding pressure_bc_index on coarser level
            DO I=1,4
               IBCF(I) = WF(IWF(I))%BOUNDARY_TYPE
            ENDDO
            IF (IBCF(1)==NULL_BOUNDARY.OR.IBCF(2)==NULL_BOUNDARY.OR.&
                IBCF(3)==NULL_BOUNDARY.OR.IBCF(4)==NULL_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = NULL_BOUNDARY
            ELSE IF (IBCF(1)==SOLID_BOUNDARY.OR.IBCF(2)==SOLID_BOUNDARY.OR.&
                     IBCF(3)==SOLID_BOUNDARY.OR.IBCF(4)==SOLID_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
            ELSE IF (IBCF(1)==OPEN_BOUNDARY.OR.IBCF(2)==OPEN_BOUNDARY.OR.&
                     IBCF(3)==OPEN_BOUNDARY.OR.IBCF(4)==OPEN_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = OPEN_BOUNDARY
            ELSE IF (IBCF(1)==INTERPOLATED_BOUNDARY.OR.IBCF(2)==INTERPOLATED_BOUNDARY.OR.&
                     IBCF(3)==INTERPOLATED_BOUNDARY.OR.IBCF(4)==INTERPOLATED_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = INTERPOLATED_BOUNDARY
            ELSE IF (IBCF(1)==MIRROR_BOUNDARY.OR.IBCF(2)==MIRROR_BOUNDARY.OR.&
                     IBCF(3)==MIRROR_BOUNDARY.OR.IBCF(4)==MIRROR_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = MIRROR_BOUNDARY
            ELSE
               CALL SCARC_SHUTDOWN('Wrong boundary type ', 'NONE', -999)
            ENDIF

            !> in case of an internal boundary set WALL(10:15,IWC)
            IF (NOMF(1) > 0) THEN

               OLC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)

               SELECT CASE (ABS(IOR0))
                  CASE (1)
                     JDIFF = WF(IWF(2))%IYN(1) - WF(IWF(1))%IYN(1)
                     KDIFF = WF(IWF(3))%IZN(1) - WF(IWF(1))%IZN(1)
                     IF (JDIFF==1 .AND. KDIFF==1) THEN
                        IY1 = WF(IWF(2))%IYN(2)/2
                        IY2 = IY1
                        IZ1 = WF(IWF(3))%IZN(2)/2
                        IZ2 = IZ1
                     ELSE IF (JDIFF==2 .AND. KDIFF==2) THEN
                        IY1 = WF(IWF(1))%IYN(2)/2
                        IY2 = WF(IWF(2))%IYN(2)/2
                        IZ1 = WF(IWF(1))%IZN(2)/2
                        IZ2 = WF(IWF(3))%IZN(2)/2
                     ELSE IF (JDIFF==0 .AND. KDIFF==0) THEN
                        IY1 = WF(IWF(1))%IYN(1)/2
                        IY2 = IY1
                        IZ1 = WF(IWF(1))%IZN(1)/2
                        IZ2 = IZ1
                     ELSE
                        CALL SCARC_SHUTDOWN('Wrong resolutions for IOR0 = ', 'NONE', IOR0)
                     ENDIF
                  CASE (2)
                     IDIFF = WF(IWF(2))%IXN(1) - WF(IWF(1))%IXN(1)
                     KDIFF = WF(IWF(3))%IZN(1) - WF(IWF(1))%IZN(1)
                     IF (IDIFF==1 .AND. KDIFF==1) THEN
                        IX1 = WF(IWF(2))%IXN(2)/2
                        IX2 = IX1
                        IZ1 = WF(IWF(3))%IZN(2)/2
                        IZ2 = IZ1
                     ELSE IF (IDIFF==2 .AND. KDIFF==2) THEN
                        IX1 = WF(IWF(1))%IXN(2)/2
                        IX2 = WF(IWF(2))%IXN(2)/2
                        IZ1 = WF(IWF(1))%IZN(2)/2
                        IZ2 = WF(IWF(3))%IZN(2)/2
                     ELSE IF (IDIFF==0 .AND. KDIFF==0) THEN
                        IX1 = WF(IWF(1))%IXN(2)/2
                        IX2 = IX1
                        IZ1 = WF(IWF(1))%IZN(2)/2
                        IZ2 = IZ1
                     ELSE
                        CALL SCARC_SHUTDOWN('Wrong resolutions for IOR0 = ', 'NONE', IOR0)
                     ENDIF
                  CASE (3)
                     IDIFF = WF(IWF(2))%IXN(1) - WF(IWF(1))%IXN(1)
                     JDIFF = WF(IWF(3))%IYN(1) - WF(IWF(1))%IYN(1)
                     IF (IDIFF==1 .AND. JDIFF==1) THEN
                        IX1 = WF(IWF(2))%IXN(2)/2
                        IX2 = IX1
                        IY1 = WF(IWF(3))%IYN(2)/2
                        IY2 = IY1
                     ELSE IF (IDIFF==2 .AND. JDIFF==2) THEN
                        IX1 = WF(IWF(1))%IXN(2)/2
                        IX2 = WF(IWF(2))%IXN(2)/2
                        IY1 = WF(IWF(1))%IYN(2)/2
                        IY2 = WF(IWF(3))%IYN(2)/2
                     ELSE IF (IDIFF==0 .AND. JDIFF==0) THEN
                        IX1 = WF(IWF(2))%IXN(2)/2
                        IX2 = IX1
                        IY1 = WF(IWF(3))%IYN(2)/2
                        IY2 = IY1
                     ELSE
                        CALL SCARC_SHUTDOWN('Wrong resolutions for IOR0 = ', 'NONE', IOR0)
                     ENDIF
               END SELECT

               SELECT CASE (IOR0)
                  CASE (1)
                     IX1 = MESHES(NOMF(1))%IBAR/IREFINE
                     IX2 = IX1
                  CASE (-1)
                     IX1 = 1
                     IX2 = IX1
                  CASE (2)
                     IY1 = MESHES(NOMF(1))%JBAR/IREFINE
                     IY2 = IY1
                  CASE (-2)
                     IY1 = 1
                     IY2 = IY1
                  CASE (3)
                     IZ1 = MESHES(NOMF(1))%KBAR/IREFINE
                     IZ2 = IZ1
                  CASE (-3)
                     IZ1 = 1
                     IZ2 = IZ1
               END SELECT

               WC(IWC)%IXN(1) = IX1
               WC(IWC)%IYN(1) = IY1
               WC(IWC)%IZN(1) = IZ1
               WC(IWC)%IXN(2) = IX2
               WC(IWC)%IYN(2) = IY2
               WC(IWC)%IZN(2) = IZ2

               !> Allocate and specify ICN and ICE arrays for OLC
               NCPL = (IZ2-IZ1+1)*(IY2-IY1+1)*(IX2-IX1+1)
               OLC%NCPL = NCPL

               CALL SCARC_SETUP_WALLCELL_NEIGHBOR(IX1, IX2, IY1, IY2, IZ1, IZ2, IWC, IOR0, NM, NOMF(1), NL)

            ENDIF
         ENDIF
         IWC = IWC + 1
      ENDDO
   ENDDO
ENDDO

FF%IOFFSET_WALL = FF%IOFFSET_WALL + FF%NFW
CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_FACE


!> -------------------------------------------------------------------------------------------------
!> Set analogues to I_MIN, IMAX, J_MIN, J_MAX, K_MIN, K_MAX on coarse level (only GMG!)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)
INTEGER, INTENT(IN) :: IREFINE, NM, NOM, NL
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IW
LOGICAL :: FOUND
TYPE (OSCARC_TYPE)        , POINTER :: OS
TYPE (SCARC_EXCHANGE_TYPE), POINTER :: OX
TYPE (SCARC_LEVEL_TYPE)   , POINTER :: L, OL

L  => SCARC(NM)%LEVEL(NL)
OS => SCARC(NM)%OSCARC(NOM)
OX => OS%EXCHANGE
OL => OS%LEVEL(NL)

IMIN=0
IMAX=MESHES(NOM)%IBAR/IREFINE+1

IF (TWO_D) THEN
   JMIN=0
   JMAX=2
ELSE
   JMIN=0
   JMAX=MESHES(NOM)%JBAR/IREFINE+1
ENDIF

KMIN=0
KMAX=MESHES(NOM)%KBAR/IREFINE+1

OX%NIC_R = 0
FOUND = .FALSE.

SEARCH_LOOP: DO IW=1, OL%NCG

   !> neighborship structure already known from finest level
   IF (L%WALL(IW)%NOM/=NOM) CYCLE SEARCH_LOOP
   OX%NIC_R = OX%NIC_R + 1
   FOUND = .TRUE.

   SELECT CASE (L%WALL(IW)%IOR)
      CASE ( 1)
         IMIN=MAX(IMIN,L%WALL(NM)%IXN(1)-1)
      CASE (-1)
         IMAX=MIN(IMAX,L%WALL(NM)%IXN(2)+1)
      CASE ( 2)
         JMIN=MAX(JMIN,L%WALL(NM)%IYN(1)-1)
      CASE (-2)
         JMAX=MIN(JMAX,L%WALL(NM)%IYN(2)+1)
      CASE ( 3)
         KMIN=MAX(KMIN,L%WALL(NM)%IZN(1)-1)
      CASE (-3)
         KMAX=MIN(KMAX,L%WALL(NM)%IZN(2)+1)
   END SELECT
ENDDO SEARCH_LOOP

N_EXCHANGE = N_EXCHANGE+1

END SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS

!> ----------------------------------------------------------------------------------------------------
!> Allocate several global structures for data exchange
!> --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NP

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_GLOBALS')

!> Allocate and preset counter and displacement vector for global data exchanges
CALL SCARC_ALLOCATE_INT1 (COUNTS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'COUNTS')
CALL SCARC_ALLOCATE_INT1 (DISPLS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'DISPLS')

CALL SCARC_ALLOCATE_INT1 (LOCAL_INT , 1, NMESHES, NSCARC_INIT_ZERO, 'LOCAL_INT')
CALL SCARC_ALLOCATE_REAL1(LOCAL_REAL, 1, NMESHES, NSCARC_INIT_ZERO, 'LOCAL_REAL')

!> Get number of data to send per process
DO NP = 0, N_MPI_PROCESSES-1
   DO NM = 1, NMESHES
      IF (PROCESS(NM)==NP) COUNTS(NP) = COUNTS(NP) + 1
   ENDDO
ENDDO

!> Get displacements on communication vector for all meshes
DO NP = 1, N_MPI_PROCESSES-1
   DISPLS(NP) = COUNTS(NP-1) + DISPLS(NP-1)
ENDDO

CALL SCARC_SETUP_GLOBALS_STRUCTURED(NLEVEL_MIN)

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_GLOBALS


!> ----------------------------------------------------------------------------------------------------
!> Broadcast a local integer value to all and process exchanged data due to NTYPE
!> corresponding local data is passed in LOCAL_INT
!> --------------------------------------------------------------------------------------------------ss
INTEGER FUNCTION SCARC_BROADCAST_INT(NTYPE)
INTEGER, INTENT(IN) :: NTYPE

IF (N_MPI_PROCESSES > 1) THEN
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,LOCAL_INT,COUNTS,DISPLS,&
                       MPI_INTEGER,MPI_COMM_WORLD,IERROR)
ENDIF

GLOBAL_INT = 0
SELECT CASE (NTYPE)
   CASE (NSCARC_BROADCAST_SUM)
      GLOBAL_INT = SUM(LOCAL_INT(1:NMESHES))
   CASE (NSCARC_BROADCAST_PRODUCT)
      GLOBAL_INT = PRODUCT(LOCAL_INT(1:NMESHES))
   CASE (NSCARC_BROADCAST_FIRST)
      GLOBAL_INT = LOCAL_INT(1)
   CASE (NSCARC_BROADCAST_LAST)
      GLOBAL_INT = LOCAL_INT(NMESHES)
END SELECT

SCARC_BROADCAST_INT = GLOBAL_INT
RETURN
END FUNCTION SCARC_BROADCAST_INT


!> ----------------------------------------------------------------------------------------------------
!> Broadcast a local integer value to all and process exchanged data due to NTYPE
!> corresponding local data is passed in LOCAL_INT
!> --------------------------------------------------------------------------------------------------ss
DOUBLE PRECISION FUNCTION SCARC_BROADCAST_REAL(NTYPE)
INTEGER, INTENT(IN) :: NTYPE

IF (N_MPI_PROCESSES > 1) THEN
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,LOCAL_REAL,COUNTS,DISPLS,&
                       MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
ENDIF

GLOBAL_REAL = 0.0_EB
SELECT CASE (NTYPE)
   CASE (NSCARC_BROADCAST_SUM)
      GLOBAL_REAL = SUM(LOCAL_REAL(1:NMESHES))
   CASE (NSCARC_BROADCAST_MEAN)
      GLOBAL_REAL = SUM(LOCAL_REAL(1:NMESHES))/REAL(NMESHES, EB)
   CASE (NSCARC_BROADCAST_PRODUCT)
      GLOBAL_REAL = PRODUCT(LOCAL_REAL(1:NMESHES))
   CASE (NSCARC_BROADCAST_FIRST)
      GLOBAL_REAL = LOCAL_REAL(1)
   CASE (NSCARC_BROADCAST_LAST)
      GLOBAL_REAL = LOCAL_REAL(NMESHES)
END SELECT

SCARC_BROADCAST_REAL = GLOBAL_REAL
RETURN
END FUNCTION SCARC_BROADCAST_REAL


!> ----------------------------------------------------------------------------------------------------
!> Get number of local and global cells for structured case
!> local values of all meshes are passed by broadcast routine in LOCAL_INT
!> --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_SETUP_GLOBALS_STRUCTURED(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L => SCARC(NM)%LEVEL(NL)
   DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      LOCAL_INT(NM2) = SCARC(NM2)%LEVEL(NL)%NCS
   ENDDO
ENDDO

N_CELLS_GLOBAL(NL) = SCARC_BROADCAST_INT(NSCARC_BROADCAST_SUM)

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L => SCARC(NM)%LEVEL(NL)
   L%CELL%NC_LOCAL(1:NMESHES) = LOCAL_INT(1:NMESHES)
   L%CELL%NC_GLOBAL = N_CELLS_GLOBAL(NL)
   IF (NMESHES > 1) THEN
      DO NM2 = 2, NMESHES
         L%CELL%NC_OFFSET(NM2) = L%CELL%NC_OFFSET(NM2-1) + L%CELL%NC_LOCAL(NM2-1)
      ENDDO
   ENDIF
ENDDO

IF (NL == NLEVEL_MIN) THEN
   DO NM = 1, NMESHES
      SCARC(NM)%N_CELLS = LOCAL_INT(NM)
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_GLOBALS_STRUCTURED

!> -----------------------------------------------------------------------------
!> Get information about global numbers of unknowns for unstructured case
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBALS_UNSTRUCTURED(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

!> For structured discretizations, numbers are the same 
IF (BSTRUCTURED) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM)%LEVEL(NL)
      L%CELL%NC_LOCAL_DOF  = L%CELL%NC_LOCAL
      L%CELL%NC_GLOBAL_DOF = L%CELL%NC_GLOBAL
      L%CELL%NC_OFFSET_DOF = L%CELL%NC_OFFSET
   ENDDO

   N_CELLS_GLOBAL_DOF = N_CELLS_GLOBAL

!> For unstructured discretizations, compute correct numbers base on DOF
ELSE
   MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM)%LEVEL(NL)
      DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         LOCAL_INT(NM2) = SCARC(NM2)%LEVEL(NL)%CELL%NC_LOCAL_DOF(NM2)
      ENDDO
   ENDDO MESHES_LOOP1
   
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,LOCAL_INT,COUNTS,DISPLS,&
                       MPI_INTEGER,MPI_COMM_WORLD,IERROR)
   
   MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM)%LEVEL(NL)
      L%CELL%NC_LOCAL_DOF(1:NMESHES) = LOCAL_INT(1:NMESHES)
      L%CELL%NC_GLOBAL_DOF = SUM(LOCAL_INT(1:NMESHES))
      N_CELLS_GLOBAL_DOF(NL) = L%CELL%NC_GLOBAL_DOF
      IF (NMESHES > 1) THEN
         DO NM2=2,NMESHES
            L%CELL%NC_OFFSET_DOF(NM2) = L%CELL%NC_OFFSET_DOF(NM2-1) + L%CELL%NC_LOCAL_DOF(NM2-1)
         ENDDO
      ENDIF
   ENDDO MESHES_LOOP2

ENDIF

END SUBROUTINE SCARC_SETUP_GLOBALS_UNSTRUCTURED

!> ----------------------------------------------------------------------------------------------------
!> Setup system of equation:
!> Define matrix stencils and initialize matrices and boundary conditions on all needed levels
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM
INTEGER :: NM, NL

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_SYSTEM')

!> ---------------------------------------------------------------------------------------------
!> Setup sizes for system matrices 
!> ---------------------------------------------------------------------------------------------
SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   CASE (NSCARC_METHOD_KRYLOV)
      CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)
      IF (BTWOLEVEL) CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MAX)    !> setup size for coarse grid
      IF (BCGGMG) THEN
         DO NL=NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NL)                    !> setup size for all levels
         ENDDO
      ENDIF

   CASE (NSCARC_METHOD_MULTIGRID)
      SELECT CASE (TYPE_MULTIGRID)
         CASE (NSCARC_MULTIGRID_GEOMETRIC)
            DO NL=NLEVEL_MIN, NLEVEL_MAX
               CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX, NL)                  !> setup size for all levels
            ENDDO
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)
            CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)             !> setup sizes for AMG
      END SELECT

END SELECT SELECT_METHOD


MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   SELECT_SOLVER: SELECT CASE (TYPE_METHOD)

      !> ---------------------------------------------------------------------------------------------
      !> Krylov method (CG/BICG) as main solver, different preconditioners possible
      !> ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_KRYLOV)

         SELECT_PRECON: SELECT CASE (TYPE_PRECON)

            !> ---------------------------------------------------------------------------------------
            !> in case of multigrid as preconditioner:
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_RELAX_MULTIGRID)

               SELECT_PRECON_MG: SELECT CASE (TYPE_MULTIGRID)

                  !> Geometric multigrid:
                  !>    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
                  CASE (NSCARC_MULTIGRID_GEOMETRIC)

                     DO NL = NLEVEL_MIN, NLEVEL_MAX
                        CALL SCARC_SETUP_MATRIX  (NM, NL)
                        CALL SCARC_SETUP_BOUNDARY(NM, NL)
#if defined(WITH_MKL)
                        IF (TYPE_LU_LEVEL(NL) == NSCARC_MKL_LOCAL) CALL SCARC_SETUP_MATRIX_MKL(NM, NL)
                        IF (TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) CALL SCARC_SETUP_MATRIX_MKL(NM, NL)
#endif
                     ENDDO


                  !> Algebraic multigrid:
                  !>    -  use compact storage technique on all levels (no other choise possible!)
                  !>    -  assemble standard n-point-matrix only on finest level
                  !>    -  construct all coarser levels by requested coarsening strategy
                  CASE (NSCARC_MULTIGRID_ALGEBRAIC)

                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               END SELECT SELECT_PRECON_MG

#ifdef WITH_MKL
            !> ---------------------------------------------------------------------------------------
            !> in case of Pardiso as preconditioners and possible coarse grid correction 
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_RELAX_PARDISO)

               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MIN)

               IF (BTWOLEVEL) THEN
                  CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
                  IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MAX)
               ENDIF

            !> ---------------------------------------------------------------------------------------
            !> in case of Cluster_Sparse_Solver as preconditioner 
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_RELAX_CLUSTER)

               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
               CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MIN)
#endif

            !> ---------------------------------------------------------------------------------------
            !> in case of one-level preconditioners (JACOBI/SSOR/FFT):
            !> assemble standard n-point-matrix on finest level with possible coarse grid correction
            !> ---------------------------------------------------------------------------------------
            CASE DEFAULT

               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               IF (BTWOLEVEL) THEN
                  CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
#ifdef WITH_MKL
                  IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MAX)
#endif
               ENDIF


         END SELECT SELECT_PRECON

      !> ---------------------------------------------------------------------------------------------
      !> Multigrid as main solver
      !> ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_MULTIGRID)

         SELECT_MG: SELECT CASE (TYPE_MULTIGRID)

            !> ---------------------------------------------------------------------------------------
            !> Geometric multigrid:
            !>    -  assemble standard n-point-matrix hierarchy on all levels
            !>    -  use MKL coarse grid solver if requested
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_GEOMETRIC)

               DO NL = NLEVEL_MIN, NLEVEL_MAX
                  CALL SCARC_SETUP_MATRIX  (NM, NL)
                  CALL SCARC_SETUP_BOUNDARY(NM, NL)
               ENDDO

#ifdef WITH_MKL
               DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                  IF (TYPE_LU_LEVEL(NL) == NSCARC_MKL_LOCAL   .OR. &
                      TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
                     CALL SCARC_SETUP_MATRIX_MKL(NM, NL)
                  ENDIF
               ENDDO

               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MAX)
#endif

            !> ---------------------------------------------------------------------------------------
            !> Algebraic multigrid:
            !>    -  use compact storage technique (no other choice possible!)
            !>    -  assemble standard n-point-matrix only on finest level
            !>    -  construct all coarser levels later by requested coarsening strategy
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_ALGEBRAIC)

               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         END SELECT SELECT_MG



#ifdef WITH_MKL
      !> ---------------------------------------------------------------------------------------------
      !> MKL-Pardiso method 
      !> ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_LUDECOMP)

         CALL SCARC_SETUP_MATRIX(NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
         CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MIN)
#endif

   END SELECT SELECT_SOLVER

ENDDO MESHES_LOOP

!> ------------------------------------------------------------------------------------------------
!> Exchange matrix stencil information in case of AMG
!> ------------------------------------------------------------------------------------------------
IF (NMESHES>1 .AND. BAMG) THEN

   !> set sizes for matrix stencils on overlapping parts
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_STENCIL, NLEVEL_MIN)
   CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_STENCIL, NLEVEL_MIN)

   !> exchange matrix entries on overlapping parts
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SYSTEM, NLEVEL_MIN)
   CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_SYSTEM, NLEVEL_MIN)

ENDIF

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_SYSTEM


!> ------------------------------------------------------------------------------------------------
!> Allocate matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
!> Use compact storage technique:
!>
!> Compression technique to store sparse matrices, non-zero entries are stored
!> in a 1D-vector B(.), row after row,
!> Each row starts with its diagonal entry followed by the other non-zero entries
!> In order to identify each element, pointer arrays ROW and COL are needed,
!> ROW points to the several diagonal entries in vector B(.),
!> COL points to the columns which non-zero entries in the matrix stencil
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IX, IY, IZ !, IW
INTEGER :: ICNUM, ICTYPE, IWNUM, IWTYPE
INTEGER :: IC, ICI, IP, IOR0
INTEGER :: NAV, NAC, NAR, NSTENCIL
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
INTEGER, POINTER, DIMENSION(:,:)   :: WALL_INDEX
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_MATRIX_TYPE) , POINTER :: A
!#ifdef WITH_MKL
!TYPE (SCARC_MKL_TYPE), POINTER :: MKL
!#endif

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_MATRIX')

L => SCARC(NM)%LEVEL(NL)
A => SCARC(NM)%SYSTEM(NL)%A

IF (TWO_D) THEN
   NSTENCIL = 5
ELSE
   NSTENCIL = 7
ENDIF
NAV = L%NCS*NSTENCIL 
NAC = L%NCS*NSTENCIL 
NAR = L%NCS + 1
CALL SCARC_ALLOCATE_MATRIX(A, NAV, NAC, NAR, NSTENCIL, NSCARC_INIT_ZERO, 'A')

#ifdef WITH_MKL
IF ((TYPE_LUDECOMP == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_LUDECOMP == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
     TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
   CALL SCARC_ALLOCATE_INT1 (A%COLG, 1, NAC, NSCARC_INIT_ZERO, 'A%COLG')
ENDIF
#endif


IF (NL == NLEVEL_MIN) THEN
   CELL_INDEX => MESHES(NM)%CELL_INDEX
   WALL_INDEX => MESHES(NM)%WALL_INDEX
   !NCELLS     =  N_CELLS(NM)             ! ACHTUNG: Kann weg ?
ELSE
   CELL_INDEX => L%CELL%CINDEX
   WALL_INDEX => L%CELL%WINDEX
   !NCELLS     =  L%N_CELLS               ! ACHTUNG: Kann weg ?
ENDIF

IWTYPE=-9999

!> Compute single matrix entries and corresponding row and column pointers
!> Along internal boundaries use placeholders for the neighboring matrix entries
!> which will be communicated in a following step
IF (TWO_D) THEN

   IP  = 1
   IY  = 1
   DO IZ = 1, L%NZ
      DO IX = 1, L%NX

         IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. L%CELL%STATE(IX, IY, IZ)/=NSCARC_DISCRET_GASPHASE) CYCLE

         IWTYPE = 0
         IWNUM  = 0

         IC  = L%CELL%DOF(IX, IY, IZ)
         ICI = CELL_INDEX(IX, IY, IZ)

         CALL SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)

         !> IOR = 1
         IOR0 = 1
         IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
         IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
         IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
         !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            ICNUM  = L%CELL%DOF(IX-1, IY, IZ)
            ICTYPE = L%CELL%STATE(IX-1, IY, IZ)
            CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
         ENDIF

         !> IOR = -1
         IOR0 = -1
         IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
         IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
         IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
         !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            ICNUM  = L%CELL%DOF(IX+1, IY, IZ)
            ICTYPE = L%CELL%STATE(IX+1, IY, IZ)
            CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
         ENDIF

         !> IOR = 3
         IOR0 = 3
         IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
         IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
         IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
         !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            ICNUM  = L%CELL%DOF(IX, IY, IZ-1)
            ICTYPE = L%CELL%STATE(IX, IY, IZ-1)
            CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
         ENDIF

         !> IOR = -3
         IOR0 = -3
         IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
         IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
         IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
         !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            ICNUM  = L%CELL%DOF(IX, IY, IZ+1)
            ICTYPE = L%CELL%STATE(IX, IY, IZ+1)
            CALL SCARC_SETUP_MATRIX_SUBDIAG  (IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
         ENDIF

      ENDDO
   ENDDO

ELSE   ! 3D

   IP  = 1
   DO IZ = 1, L%NZ
      DO IY = 1, L%NY
         DO IX = 1, L%NX
 
            IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. L%CELL%STATE(IX, IY, IZ)/=NSCARC_DISCRET_GASPHASE) CYCLE
 
            IWTYPE = 0
            IWNUM  = 0
 
            IC  = L%CELL%DOF(IX, IY, IZ)
            ICI = CELL_INDEX(IX, IY, IZ)
 
            CALL SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)
 
            !> set subdiagonal entries depending on type of neighboring cell
            !> IOR = 1
            IOR0 = 1
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)          ! get wall index of related cell
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE    ! if there is a boundary, get corresponding type
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN    ! do, if neighbor is internal or on adjacent mesh
               ICNUM  = L%CELL%DOF(IX-1, IY, IZ)
               ICTYPE = L%CELL%STATE(IX-1, IY, IZ)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = -1
            IOR0 = -1
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%CELL%DOF(IX+1, IY, IZ)
               ICTYPE = L%CELL%STATE(IX+1, IY, IZ)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = 2
            IOR0 = 2
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%CELL%DOF(IX, IY-1, IZ)
               ICTYPE = L%CELL%STATE(IX, IY-1, IZ)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IY, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = -2
            IOR0 = -2
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%CELL%DOF(IX, IY+1, IZ)
               ICTYPE = L%CELL%STATE(IX, IY+1, IZ)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IY, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = 3
            IOR0 = 3
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%CELL%DOF(IX, IY, IZ-1)
               ICTYPE = L%CELL%STATE(IX, IY, IZ-1)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = -3
            IOR0 = -3
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%CELL%DOF(IX, IY, IZ+1)
               ICTYPE = L%CELL%STATE(IX, IY, IZ+1)
               CALL SCARC_SETUP_MATRIX_SUBDIAG  (IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
         ENDDO
      ENDDO
   ENDDO

ENDIF
A%ROW(A%NAR) = IP
A%NAV        = IP-1                         !> set correct number of matrix entries
A%NAS        = IP-1                                   
A%NAC        = IP-1

CALL SCARC_REDUCE_MATRIX(A, A%NAV, 'Reduced System-Matrix')

!WRITE(*,*) 'ACHTUNG: COLG muss auch noch reduziert werden!'
CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_MATRIX


!> ------------------------------------------------------------------------------------------------
!> Set main diagonal entry for matrix - full matrix of the global problem
!> In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
!> If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)
INTEGER, INTENT(IN)  :: IC, IX, IY, IZ
INTEGER, INTENT(IN)  :: NM, NL
INTEGER, INTENT(INOUT) :: IP
TYPE (SCARC_LEVEL_TYPE) ,  POINTER :: L
TYPE (SCARC_MATRIX_TYPE) , POINTER :: A
TYPE (SCARC_COORD_TYPE)  , POINTER :: C

A => SCARC(NM)%SYSTEM(NL)%A
L => SCARC(NM)%LEVEL(NL)
C => SCARC(NM)%LEVEL(NL)%COORD

A%VAL(IP) = A%VAL(IP) - 2.0_EB/(C%DXL(IX-1)*C%DXL(IX))
IF (.NOT.TWO_D) A%VAL(IP) = A%VAL(IP) - 2.0_EB/(C%DYL(IY-1)*C%DYL(IY))
A%VAL(IP) = A%VAL(IP) - 2.0_EB/(C%DZL(IZ-1)*C%DZL(IZ))

A%ROW(IC) = IP
A%COL(IP) = IC
      
#ifdef WITH_MKL
IF ((TYPE_LUDECOMP == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_LUDECOMP == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
     TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
   A%COLG(IP) = A%COL(IP) + L%CELL%NC_OFFSET(NM)
ENDIF
#endif
      
IP = IP + 1
END SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG


SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG (IOR0, IC, IP, IPTR, NBR_TYPE, NBR_NUM, NM, NL)
INTEGER, INTENT(IN)  :: IC, IOR0, IPTR, NBR_TYPE, NBR_NUM
INTEGER, INTENT(IN)  :: NM, NL
INTEGER, INTENT(INOUT) :: IP
INTEGER  :: IW
#ifdef WITH_MKL
INTEGER  :: IX, IY, IZ
#endif
REAL(EB) :: DSCAL, DH1, DH2
LOGICAL  :: BINTERNAL
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_MATRIX_TYPE) , POINTER :: A

L => SCARC(NM)%LEVEL(NL)
A => SCARC(NM)%SYSTEM(NL)%A

DH1 = L%FACE(IOR0)%DH(IPTR-1)
DH2 = L%FACE(IOR0)%DH(IPTR)
   
!> set bounds and step sizes depending on orientation of face (lower or upper subdiagonal)
IF (IOR0 > 0) THEN
   BINTERNAL = IPTR > 1
   DSCAL = 2.0_EB/(DH1*(DH1+DH2))
ELSE
   BINTERNAL = IPTR < L%FACE(IOR0)%NFC
   DSCAL = 2.0_EB/(DH2*(DH1+DH2))
ENDIF

!> if IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
IF (BINTERNAL) THEN

   IF (PRES_ON_WHOLE_DOMAIN .OR. NBR_TYPE == NSCARC_DISCRET_GASPHASE) THEN
      A%VAL(IP) = A%VAL(IP) + DSCAL
      A%COL(IP) = NBR_NUM
#ifdef WITH_MKL
      IF ((TYPE_LUDECOMP == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
          (TYPE_LUDECOMP == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
           TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
         A%COLG(IP)= A%COL(IP) + L%CELL%NC_OFFSET(NM)
      ENDIF
#endif
      IP = IP + 1
   ELSE
      WRITE(*,*) 'ACHTUNG: Fehlt! Nichts zu tun?'
   ENDIF
         
!> if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
ELSE IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
   
   IF (HAS_NEIGHBOR(L%CELL%DOF, L%CELL%STATE, L%WALL, IC, IW, L%FACE(IOR0)%IWG_PTR, L%FACE(IOR0)%NFW)) THEN
      A%VAL(IP) = A%VAL(IP) + DSCAL
      A%COL(IP) = L%WALL(IW)%ICE(1)
      !A%COL(IP) = L%WALL(IW)%ICO
         
#ifdef WITH_MKL
      IF ((TYPE_LUDECOMP == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
          (TYPE_LUDECOMP == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
           TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
         IX = L%WALL(IW)%IXG
         IY = L%WALL(IW)%IYG
         IZ = L%WALL(IW)%IZG
         A%COLG(IP) = L%CELL%DOF(IX, IY, IZ) + L%CELL%NC_OFFSET(L%WALL(IW)%NOM)
      ENDIF
#endif
         
      IP = IP + 1
   ENDIF
   
ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG


!> ------------------------------------------------------------------------------------------------
!> Determine if cell IC has a neighbor and, if yes, return corresponding IW-value
!> ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION HAS_NEIGHBOR(DOF0, STATE0, WALL0, IC, IW, IW0, IL0)
TYPE (SCARC_WALL_TYPE), DIMENSION(:), INTENT(IN):: WALL0
INTEGER, DIMENSION(0:,0:,0:), INTENT(IN):: DOF0, STATE0
INTEGER, INTENT(IN)    :: IC, IW0, IL0
INTEGER :: IXW0, IYW0, IZW0
INTEGER :: IXG0, IYG0, IZG0
INTEGER :: IC0
INTEGER, INTENT(INOUT) :: IW

HAS_NEIGHBOR = .FALSE.
SEARCH_WALL_CELLS_LOOP: DO IW = IW0, IW0+IL0-1

  IXW0 = WALL0(IW)%IXW
  IYW0 = WALL0(IW)%IYW
  IZW0 = WALL0(IW)%IZW

  IXG0 = WALL0(IW)%IXG
  IYG0 = WALL0(IW)%IYG
  IZG0 = WALL0(IW)%IZG

  IC0 = DOF0(IXW0, IYW0, IZW0) 
  IF (IC == IC0 .AND. (WALL0(IW)%NOM /= 0.AND.STATE0(IXG0, IYG0, IZG0)/=NSCARC_DISCRET_SOLID)) THEN
     HAS_NEIGHBOR = .TRUE.
     EXIT SEARCH_WALL_CELLS_LOOP
  ENDIF
ENDDO SEARCH_WALL_CELLS_LOOP

RETURN
END FUNCTION HAS_NEIGHBOR


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Build system matrix for MKL solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MKL (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, JC, JC0
INTEGER :: ICS, JCS
INTEGER :: ICOL, JCOL, IAS
INTEGER :: ISYM, JSYM, NSYM
REAL(EB) :: VAL = 0.0_EB, VALS = 0.0_EB, DIFF
LOGICAL  :: BSYM
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A
#ifdef WITH_MKL_FB
TYPE (SCARC_MATRIX_FB_TYPE), POINTER :: ASYM
#else
TYPE (SCARC_MATRIX_TYPE), POINTER :: ASYM
#endif
INTEGER, DIMENSION(:), ALLOCATABLE :: ICOL_AUX, JC_AUX

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_MATRIX_MKL')

L  => SCARC(NM)%LEVEL(NL)
A  => SCARC(NM)%SYSTEM(NL)%A

!> ------------------------------------------------------------------------------------------------
!> Store only symmetric parts of matrix (diagonal and upper part)
!> ------------------------------------------------------------------------------------------------
IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN

   !> First check whether symmetry of system matrix is guaranteed
   DO IC = 1, L%NCS
   
      COLUMN_LOOP: DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
         ICS = A%COL(ICOL)
         VAL = A%VAL(ICOL)
         IF (ICS > IC .AND. ICS <= L%NCS) THEN
            BSYM = .FALSE.
            DO JCOL = A%ROW(ICS)+1, A%ROW(ICS+1)-1
               JCS = A%COL(JCOL)
               IF (JCS == IC) THEN
                  VALS = A%VAL(JCOL)
                  DIFF = ABS(VAL-VALS)
                  IF (ABS(VAL - VALS) < 1E-6) THEN
                     BSYM=.TRUE.
                     CYCLE COLUMN_LOOP
                  ENDIF
               ENDIF
            ENDDO
             IF (.NOT.BSYM) CALL SCARC_SHUTDOWN('System matrix not symmetric ', 'NONE', NM)
         ENDIF
      ENDDO COLUMN_LOOP
   
   ENDDO

   !> Compute number of entries in symmetric matrix
   A%NAS = 0
   DO IC = 1, L%NCS
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         IF (TYPE_LU_LEVEL(NL) == NSCARC_MKL_LOCAL) THEN
            JC = A%COL(ICOL)
            IF (JC >= IC .AND. JC <= L%NCS) A%NAS = A%NAS+1   !really neccesary ???
         ELSE IF (TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
            JC = A%COLG(ICOL)
            IF (JC >= IC + L%CELL%NC_OFFSET(NM)) A%NAS = A%NAS+1
         ELSE
            CALL SCARC_SHUTDOWN('Error in setup matrix, TYPE_LU_LEVEL = ', 'NONE', TYPE_LU_LEVEL(NL))
         ENDIF
      ENDDO
   ENDDO

!> Non-symmetric case: use complete matrix (not only symmetric part) of matrix
ELSE
   A%NAS = A%NAV
ENDIF

!> ------------------------------------------------------------------------------------------------
!> allocate storage for symmetric matrix and its column and row pointers
!> ------------------------------------------------------------------------------------------------
#ifdef WITH_MKL_FB
ASYM => SCARC(NM)%SYSTEM(NL)%ASYM_FB
CALL SCARC_ALLOCATE_MATRIX_FB(ASYM, A%NAS, A%NAS, A%NAR, A%NSTENCIL, NSCARC_INIT_ZERO, 'ASYM')
#else
ASYM => SCARC(NM)%SYSTEM(NL)%ASYM
CALL SCARC_ALLOCATE_MATRIX(ASYM, A%NAS, A%NAS, A%NAR, A%NSTENCIL, NSCARC_INIT_ZERO, 'ASYM')
#endif

IF ((TYPE_LUDECOMP == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_LUDECOMP == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
     TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
   CALL SCARC_ALLOCATE_INT1(ICOL_AUX, 1, A%NSTENCIL, NSCARC_INIT_NONE, 'ICOL_AUX')
   CALL SCARC_ALLOCATE_INT1(JC_AUX  , 1, A%NSTENCIL, NSCARC_INIT_NONE, 'JC_AUX')
ENDIF

!> ------------------------------------------------------------------------------------------------
!> extract symmetric matrix part from usual system matrix
!> ------------------------------------------------------------------------------------------------
IAS = 1
DO IC = 1, L%NCS
   ASYM%ROW(IC) = IAS

   !> blockwise use of MKL solver
   IF (TYPE_LU_LEVEL(NL) == NSCARC_MKL_LOCAL) THEN    !really neccesary ??

      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         JC = A%COL(ICOL)

            IF (JC >= IC .AND. JC <= L%NCS) THEN
               ASYM%COL(IAS) = A%COL(ICOL) 
#ifdef WITH_MKL_FB
               ASYM%VAL(IAS) = REAL(A%VAL(ICOL),FB)
#else
               ASYM%VAL(IAS) = A%VAL(ICOL)
#endif
               IAS = IAS + 1
            ENDIF
      ENDDO

   !> global use of MKL solver
   ELSE IF (TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN   

      !> store indices of all diagonal and upper-diagonal entries
      ICOL_AUX = 0
      JC_AUX   = 99999999
      ISYM = 1
      JC0 = A%COLG(A%ROW(IC))
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         !JC = A%COL(ICOL)
         JC = A%COLG(ICOL)
         IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
            IF (JC >= JC0) THEN
               ICOL_AUX(ISYM) = ICOL
               JC_AUX(ISYM) = JC
               ISYM  = ISYM  + 1
            ENDIF
         ELSE
            ICOL_AUX(ISYM) = ICOL
            JC_AUX(ISYM) = JC
            ISYM  = ISYM  + 1
         ENDIF
      ENDDO
      NSYM = ISYM - 1
   
      !> sort them in increasing order (for the use of DSS-PARDISO functionality)
      JSYM = 1
      SORT_LOOP: DO WHILE (JSYM <= NSYM)
         DO ISYM = 1, NSYM
            JC = JC_AUX(ISYM)
            IF (JC == 99999999) CYCLE
            IF (JC <= MINVAL(JC_AUX)) THEN
               ICOL = ICOL_AUX(ISYM)
#ifdef WITH_MKL_FB
               ASYM%VAL(IAS) = REAL(A%VAL(ICOL),FB)
#else
               ASYM%VAL(IAS) = A%VAL(ICOL)
#endif
               !ASYM%COL(IAS) = A%COL(ICOL)
               ASYM%COL(IAS) = A%COLG(ICOL)
               JC_AUX(ISYM) = 99999999            ! mark entry as already used
               IAS  = IAS  + 1
            ENDIF
         ENDDO
         JSYM = JSYM + 1
      ENDDO SORT_LOOP
   ENDIF
ENDDO

ASYM%ROW(L%NCS+1) = IAS
ASYM%COL = ASYM%COL 

IF ((TYPE_LUDECOMP == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_LUDECOMP == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
     TYPE_LU_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
   DEALLOCATE (ICOL_AUX, STAT=IERROR)
   DEALLOCATE (JC_AUX, STAT=IERROR)
   DEALLOCATE (A%COLG, STAT=IERROR)
ENDIF

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_MATRIX_MKL
#endif

!> ----------------------------------------------------------------------------------------------------
!> Extract overlapping matrix data
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OVERLAPS (NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM
INTEGER :: IW, IWG, IWC, IG, IC, JC, JCA, ICW, ICN, ICE, ICE0, ICG
INTEGER :: IX, IY, IZ
INTEGER :: ICOL, ICCE, ICPL, NCOL
INTEGER :: ICC, ICOLG, ICOLE, IOR0, NCE0
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L, OL, LF, LC, OLF, OLC
TYPE (SCARC_MATRIX_TYPE) , POINTER :: A, OA, PF, OPF


SELECT_TYPE: SELECT CASE (NTYPE)

   !> --------------- set stencil sizes on overlapping parts  -----------------------------
   CASE (NSCARC_MATRIX_STENCIL)

      STENCIL_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%SYSTEM(NL)%A

         STENCIL_EXTERNAL_CELL_LOOP: DO ICE = L%NCS+1, L%NCE

           IW = ABS(L%MAP%ICE_TO_IWG(ICE))
           ICG = L%WALL(IW)%ICG(1)

           NOM = L%WALL(IW)%NOM
           OA => SCARC(NM)%OSCARC(NOM)%SYSTEM(NL)%A

           A%ROW(ICE+1) = A%ROW(ICE) + OA%STENCIL(ICG)

         ENDDO STENCIL_EXTERNAL_CELL_LOOP

      ENDDO STENCIL_MESHES_LOOP

   !> -------------------- after exchange of system matrix ------------------------
   CASE (NSCARC_MATRIX_SYSTEM)

      SYSTEM_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%SYSTEM(NL)%A

         NCE0 = L%NCE+1
         ICE = 1
         DO ICE = L%NCS+1, L%NCE

            IC  = ICE
            IWG = L%MAP%ICE_TO_IWG(ICE)
            IX  = L%WALL(IWG)%IXG
            IY  = L%WALL(IWG)%IYG
            IZ  = L%WALL(IWG)%IZG
            ICN = L%CELL%DOF(IX, IY, IZ)

            NOM = L%WALL(IWG)%NOM
            OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

            ICPL = 0
            DO ICOL = A%ROW(ICE),A%ROW(ICE+1)-1
               JC  = A%COL(ICOL)
               JCA = ABS(A%COL(ICOL))
WRITE(*,*) 'ACHTUNG HIER WEGEN ICN SCHAUEN!'
               IF (JC < 0 .AND. JCA <= OL%NC) THEN
                  !IF (OL%MAP%ICN_TO_ICE(JCA)==0) THEN
                  IF (JC==-1) THEN
!                     OL%MAP%ICN_TO_ICE(JCA) = NCE0
                     L%AMG%MAPPING(2, ICE) = NCE0
                     NCE0 = NCE0 + 1
                  ELSE
!                     L%AMG%MAPPING(2, ICE) = OL%MAP%ICN_TO_ICE(JCA)
                  ENDIF
                  L%AMG%MAPPING(3, ICE) = JC
                  A%COL(ICOL) = L%AMG%MAPPING(2,ICE)
               ELSE
                  A%COL(ICOL) = JC
               ENDIF
            ENDDO
         ENDDO
         !L%NCE2 = NCE0-1              !> needed where ????? OUTDATED ???
      ENDDO SYSTEM_MESHES_LOOP

   !> -------------- after exchange of prolongation matrix ----------------------------
   CASE (NSCARC_MATRIX_PROLONGATION)

      PROLONGATION_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         LF => SCARC(NM)%LEVEL(NL)
         LC => SCARC(NM)%LEVEL(NL+1)
         PF => SCARC(NM)%SYSTEM(NL)%P

         IG   = 0
         IWC  = 0
         ICCE = LF%AMG%NCC
         PROLONGATION_WALL_LOOP1: DO IW = 1, LF%NW

            NOM  = LF%WALL(IW)%NOM
            IF (NOM == 0) CYCLE PROLONGATION_WALL_LOOP1

            OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

            OPF => SCARC(NM)%OSCARC(NOM)%SYSTEM(NL)%P

            IOR0 = LF%WALL(IW)%IOR
            ICW  = LF%WALL(IW)%ICW
            ICC  = LF%AMG%CTYPE(ICW)

            ICE   = LF%WALL(IW)%ICE(1)
            ICG   = LF%WALL(IW)%ICG(1)
            NCOL  = OPF%ROW(ICG +1)-OPF%ROW(ICG)         !> wozu wird das nochmal gebraucht ??

            IF (ICC >= NSCARC_CELL_TYPE_COARSE) THEN
               IWC = IWC + 1
               LC%CELL%WINDEX(ICC, IOR0) = IWC
               LC%WALL(IWC)%NOM = NOM
               LC%WALL(IWC)%ICW = ICC
               LC%WALL(IWC)%IOR = IOR0

            ENDIF

            IF (LF%AMG%CTYPE(ICE) >= NSCARC_CELL_TYPE_COARSE) THEN
               ICCE = ICCE + 1
               LF%AMG%CTYPE(ICE) = ICCE
               !LC%AMG%EXTERN(ICCE) = OCF%STYLE(ICG)
               ICOL = OPF%ROW(ICG)
               IG  = IG + 1
               !OLF%AMG%ICG0 = OLF%AMG%ICG0 + 1
               OLF%AMG%ICG0 = IG
!               LC%MAP%ICN_TO_ICE(OPF%COL(ICOL))   = ICCE
               LC%MAP%ICG_TO_ICE(OLF%AMG%ICG0) = ICCE
            ENDIF

            IF (NL/=NLEVEL_MIN) CYCLE PROLONGATION_WALL_LOOP1

         ENDDO PROLONGATION_WALL_LOOP1

         !> Replace negative cell numbers of neighboring cells for internal wall cells
         !> Note that this only holds true for fine grid cells
         !> (weights of coarse cells don't reach in neighboring mesh)
         PROLONGATION_WALL_LOOP2: DO IW = 1, LF%NW

            NOM  = LF%WALL(IW)%NOM
            IF (NOM == 0) CYCLE PROLONGATION_WALL_LOOP2

            OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

            ICW = LF%WALL(IW)%ICW
            ICC = LF%AMG%CTYPE(ICW)

            IF (ICC < NSCARC_CELL_TYPE_COARSE) THEN

               DO ICOL = PF%ROW(ICW), PF%ROW(ICW+1)-1
                  JC   = PF%COL(ICOL)

                  !> Additionally identify coarse cells from second layer
                  !> adjacent to internal boundary
WRITE(*,*) 'ACHTUNG HIER WEGEN ICN SCHAUEN'
!                  IF (JC < 0) THEN
!                     IC = LF%MAP%ICN_TO_ICE(ABS(JC))
!                     PF%COL(ICOL) = LF%AMG%CTYPE(IC)
!                  ENDIF
               ENDDO
            ENDIF
         ENDDO PROLONGATION_WALL_LOOP2

         PROLONGATION_ECELL_LOOP: DO ICE0 = LF%NCS+1, LF%NCE

            IW  = LF%MAP%ICE_TO_IWG(ICE0)
            NOM = LF%WALL(ABS(IW))%NOM

            OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

            OPF => SCARC(NM)%OSCARC(NOM)%SYSTEM(NL)%P

            ICG = LF%WALL(IW)%ICG(1)
            ICE = LF%WALL(IW)%ICE(1)
            ICOLE = PF%ROW(ICE)

            DO ICOLG = OPF%ROW(ICG), OPF%ROW(ICG+1)-1
               PF%COL(ICOLE) = OPF%COL(ICOLG)
!               PF%COL(ICOLE) = LF%MAP%ICN_TO_ICE(OPF%COL(ICOLG))
               PF%VAL(ICOLE) = OPF%VAL(ICOLG)
               ICOLE = ICOLE + 1

!>              IF (JC > 0 .AND. JC <= OLF%AMG%NCC) THEN
!>                 IC = LC%MAP%ICN_TO_ICE(JC)
!>                 !OPF%COL(ICOLG) = IC
!>                 PF%COL(ICOLE)  = IC
!>              ELSE
!>                 !OPF%COL(ICOLG) = F%CELL%STYLE(ABS(JC))
!>                 PF%COL(ICOLE) = F%CELL%STYLE(ABS(JC))
!>              ENDIF
!>              PF%VAL(ICOLE) = OPF%VAL(ICOLG)
!>              ICOLE = ICOLE + 1
            ENDDO
            PF%ROW(ICE+1) = ICOLE

         ENDDO PROLONGATION_ECELL_LOOP

         LC%NW = IWC

      ENDDO PROLONGATION_MESHES_LOOP

END SELECT SELECT_TYPE

END SUBROUTINE SCARC_SETUP_OVERLAPS


!> ------------------------------------------------------------------------------------------------
!> Set pointer for different structures on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP, NW, NC, ITYPE, ICN, ICE, JC, ICOL, IS=0
REAL(EB) :: DBC
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A
TYPE (SCARC_COORD_TYPE) , POINTER :: C

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_BOUNDARY')

A => SCARC(NM)%SYSTEM(NL)%A
L => SCARC(NM)%LEVEL(NL)
C => SCARC(NM)%LEVEL(NL)%COORD

IF (PRES_ON_WHOLE_DOMAIN) THEN
   NW = L%N_WALL_CELLS_EXT
ELSE
   NW = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
ENDIF

!>
!> If A is a pure Neumann matrix, get neighboring cell indices of communicated stencil legs for condensed system
!> also save values and column indices of last matrix row of last mesh
!>
NO_DIRIC_IF: IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN

   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_INDEX  , NLEVEL_MIN)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALUE, NLEVEL_MIN)

   LAST_CELL_IN_LAST_MESH_IF: IF (NM == NMESHES) THEN

      NC = L%CELL%NC_LOCAL(NMESHES)
      IP = A%ROW(NC)

      !> store column indices and values of diagonal and all off-diagonal entries in last row
      !> index '1' correspondings to main diagonal entry
      IS = IS + 1                                                  
      ICOL = 1
      A%STORE(1)%PTR(ICOL)      = IP
      A%STORE(1)%COL(ICOL)      = A%COL(IP)
      A%STORE(1)%VAL_ORIG(ICOL) = A%VAL(IP)
      A%STORE(1)%VAL_COND(ICOL) = 1.0_EB

      DO IP = A%ROW(NC)+1, A%ROW(NC+1)-1
         ICOL = ICOL + 1
         A%STORE(1)%PTR(ICOL)      = IP
         A%STORE(1)%COL(ICOL)      = A%COL(IP)
         A%STORE(1)%VAL_ORIG(ICOL) = A%VAL(IP)
         A%STORE(1)%VAL_COND(ICOL) = 0.0_EB
      ENDDO
  
      A%STORE(1)%NROW = NC              !> row index of last row
      A%STORE(1)%NCOL = ICOL            !> number of stored columns 

      !> within last mesh: check which other cells have a connection to the last cell;
      !> in each corresponding matrix row store the column index and value of just that matrix entry
      !> for each direction only one value has to be stored
      JC = NC - 1                                             
      DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
         IF (A%COL(IP) == NC) THEN
            IS = IS + 1
            A%STORE(IS)%PTR(1)      = IP
            A%STORE(IS)%COL(1)      = JC
            A%STORE(IS)%VAL_ORIG(1) = A%VAL(IP)
            A%STORE(IS)%VAL_COND(1) = 0.0_EB
            A%STORE(IS)%NROW        = JC          
            A%STORE(IS)%NCOL        = 1          
            EXIT  
         ENDIF 
      ENDDO 

      JC = NC - L%NX
      DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
         IF (A%COL(IP) == NC) THEN
            IS = IS + 1
            A%STORE(IS)%PTR(1)      = IP
            A%STORE(IS)%COL(1)      = JC
            A%STORE(IS)%VAL_ORIG(1) = A%VAL(IP)
            A%STORE(IS)%VAL_COND(1) = 0.0_EB
            A%STORE(IS)%NROW        = JC          
            A%STORE(IS)%NCOL        = 1          
            EXIT  
         ENDIF 
      ENDDO 

      IF (.NOT.TWO_D) THEN
         JC = NC - L%NX * L%NY
         DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
            IF (A%COL(IP) == NC) THEN
               IS = IS + 1
               A%STORE(IS)%PTR(1)      = IP
               A%STORE(IS)%COL(1)      = JC
               A%STORE(IS)%VAL_ORIG(1) = A%VAL(IP)
               A%STORE(IS)%VAL_COND(1) = 0.0_EB
               A%STORE(IS)%NROW        = JC          
               A%STORE(IS)%NCOL        = 1          
               EXIT  
            ENDIF 
         ENDDO 
      ENDIF

   ENDIF LAST_CELL_IN_LAST_MESH_IF

   !> cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled 
   !> with the last cell of last mesh; 
   !> this can be a cell on the opposite side of the own mesh or on a different mesh
   !> if such a cell exists, store corresponding matrix entry
   WALL_CELLS_LOOP1: DO IW = 1, NW

      IOR0 = L%WALL(IW)%IOR
      IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE              

      I    = L%WALL(IW)%IXW
      J    = L%WALL(IW)%IYW
      K    = L%WALL(IW)%IZW

      NOM  = L%WALL(IW)%NOM
      ITYPE = L%CELL%STATE(I,J,K)

      IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.ITYPE/=NSCARC_DISCRET_GASPHASE) CYCLE
      L%WALL(IW)%ICW = L%CELL%DOF(I, J, K)
      IC = L%CELL%DOF(I, J, K)

      PERIODIC_NBR_IF1: IF (NOM == NMESHES) THEN
         ICE = L%WALL(IW)%ICE(1)                            !> adjacent ghost cell number
         ICN = L%MAP%ICE_TO_ICN(ICE)                        !> get column index of neighboring offdiagonal matrix entry
         IF (ICN /= SCARC(NMESHES)%N_CELLS) CYCLE           !> if no relation to last cell in last mesh, cycle
         DO IP = A%ROW(ICN)+1, A%ROW(ICN+1)-1
            IF (A%COL(IP) == ICE) THEN
               IS = IS + 1
               A%STORE(IS)%PTR(1)      = IP
               A%STORE(IS)%COL(1)      = ICN
               A%STORE(IS)%VAL_ORIG(1) = A%VAL(IP)
               A%STORE(IS)%VAL_COND(1) = 0.0_EB
               A%STORE(IS)%NROW        = JC          
               A%STORE(IS)%NCOL        = 1          
               EXIT  
            ENDIF 
         ENDDO 
      ENDIF PERIODIC_NBR_IF1
   ENDDO WALL_CELLS_LOOP1
   A%NSTORE = IS

ENDIF NO_DIRIC_IF

!>
!> Set correct boundary conditions for system matrix
!> Take care of whether the structured or unstructured discretization is used
!>
WALL_CELLS_LOOP2: DO IW = 1, NW

   IOR0 = L%WALL(IW)%IOR
   IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE              !> cycle boundaries in y-direction for 2D-cases

   I    = L%WALL(IW)%IXW
   J    = L%WALL(IW)%IYW
   K    = L%WALL(IW)%IZW

   NOM  = L%WALL(IW)%NOM
   ITYPE = L%CELL%STATE(I,J,K)

   IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.ITYPE/=NSCARC_DISCRET_GASPHASE) CYCLE

   L%WALL(IW)%ICW = L%CELL%DOF(I, J, K)
   IC = L%CELL%DOF(I, J, K)

   SELECT CASE (ABS(IOR0))
      CASE (1)
         DBC= C%DXI2           ! Achtung: Wirklich richtig oder Mittelwert?
      CASE (2)
         DBC= C%DYI2
      CASE (3)
         DBC= C%DZI2
   END SELECT

   !> SPD-matrix with mixture of Dirichlet and Neumann BC's according to the settings of BTYPE 
   IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN

      IP = A%ROW(IC)
      SELECT CASE (L%WALL(IW)%BTYPE)
         CASE (DIRICHLET)
            A%VAL(IP) = A%VAL(IP) - DBC
         CASE (NEUMANN)
            A%VAL(IP) = A%VAL(IP) + DBC
      END SELECT

   !> purely Neumann matrix 
   ELSE
      IP = A%ROW(IC)
      IF (L%WALL(IW)%BTYPE == NEUMANN) A%VAL(IP) = A%VAL(IP) + DBC
   ENDIF

ENDDO WALL_CELLS_LOOP2

!> If there are no Dirichlet BC's transform sytem into condensed one by replacing the
!> matrix entries in last column and row by the stored ones (zeros and one at diaonal position)
SETUP_CONDENSING_IF: IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0 .AND. TYPE_PRECON /= NSCARC_RELAX_FFT) THEN
   DO IS = 1, A%NSTORE
      DO ICOL = 1, A%STORE(IS)%NCOL
         IP = A%STORE(IS)%PTR(ICOL)
         A%VAL(IP) = A%STORE(IS)%VAL_COND(ICOL)
      ENDDO
   ENDDO 
ENDIF SETUP_CONDENSING_IF

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_BOUNDARY


!> ------------------------------------------------------------------------------------------------
!> Setup condensed system in case of periodic or pure Neumann boundary conditions
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CONDENSING (NVECTOR, NL, ITYPE)
INTEGER, INTENT(IN) :: NVECTOR, NL, ITYPE
INTEGER  :: NM, NOM, IW, IWG, IWL, NWL, IFACE, ICN, ICE, ICW, JC, NC, IS
REAL(EB) :: RHS_END
REAL(EB), POINTER, DIMENSION(:) :: VC
TYPE (SCARC_LEVEL_TYPE)   , POINTER :: L, OL
TYPE (SCARC_EXCHANGE_TYPE), POINTER :: OX
TYPE (SCARC_MATRIX_TYPE)  , POINTER :: A

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0 .OR. TYPE_PRECON == NSCARC_RELAX_FFT) RETURN

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_CONDENSING')

!> In last mesh:
!> Subtract B*RHS(end) for internal legs of stencil 
LOCAL_REAL = 0.0_EB
IF (UPPER_MESH_INDEX == NMESHES) THEN

   L  => SCARC(NMESHES)%LEVEL(NL)
   A  => SCARC(NMESHES)%SYSTEM(NL)%A
   VC => POINT_TO_VECTOR(NVECTOR, NMESHES, NL)
   NC =  L%CELL%NC_LOCAL(NMESHES)

   !> process last column entries of all rows except of last one
   !> for those rows only one matrix entry was stored, namely that one which connects to the last cell
   DO IS = 2, A%NSTORE
      JC = A%STORE(IS)%COL(1)
      IF (JC < NC) THEN
         VC(JC) = VC(JC) - A%STORE(IS)%VAL_ORIG(1)*VC(NC)
      ENDIF
   ENDDO

   LOCAL_REAL(NMESHES) = VC(NC)    !> store last entry of RHS
   VC(NC) = 0.0_EB                 !> set last entry of last mesh to zero

ENDIF

IF (ITYPE == 0) RETURN

!> Broadcast last RHS-value of last cell in last mesh to all meshes
RHS_END = SCARC_BROADCAST_REAL(NSCARC_BROADCAST_LAST)
DO NM = 1, NMESHES
   SCARC(NM)%RHS_END = RHS_END
ENDDO


!>
!> Only in case of periodic BC's: 
!> Subtract B*RHS(end) for corresponding entries of all periodic communication partners
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   SNODE = PROCESS(NM)
   RNODE = PROCESS(NMESHES)                         

   OX => SCARC(NM)%OSCARC(NMESHES)%EXCHANGE
   IF (OX%NICMAX_S==0 .OR. OX%NICMAX_R==0) CYCLE             !> if there is no communication with last mesh, cycle

   L  => SCARC(NM)%LEVEL(NL)
   OL => SCARC(NM)%OSCARC(NMESHES)%LEVEL(NL)
   VC => POINT_TO_VECTOR (NVECTOR, NM, NL)

   !> subtract B*RHS(end) at corresponding positions
   DO IFACE = 1, 6                                          !> check if this face has connection to last cell

      IWL = OL%SUBDIVISION(1, FACE_ORIENTATION(IFACE))      !> first wall cell on this face
      NWL = OL%SUBDIVISION(2, FACE_ORIENTATION(IFACE))      !> number of wall cells on this face
   
      DO IW = IWL, IWL + NWL - 1

         IWG = OL%MAP%IWL_TO_IWG(IW)                         !> corresponding global wall cell number
         ICE = L%WALL(IWG)%ICE(1)                            !> adjacent ghost cell number
         ICW = L%WALL(IWG)%ICW                               !> adjacent internal cell number
         NOM = L%WALL(IWG)%NOM                               !> neighbor for that wall cell

         IF (NOM > 0) THEN
            IF (NOM /= NMESHES) CYCLE                           !> only check for common matrix entries with last mesh
            ICN = L%MAP%ICE_TO_ICN(ICE)                         !> get column index of neighboring offdiagonal matrix entry
            IF (ICN /= SCARC(NMESHES)%N_CELLS) CYCLE            !> if no relation to last cell in last mesh, cycle
            VC(ICW) = VC(ICW) - L%MAP%ICE_TO_VAL(ICE)*SCARC(NM)%RHS_END
         ENDIF

      ENDDO
   ENDDO
ENDDO

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_CONDENSING


!> ----------------------------------------------------------------------------------------------------
!> Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_METHODS
INTEGER :: NSTACK

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_METHODS')
SELECT_METHOD: SELECT CASE (TRIM(SCARC_METHOD))

   !> ------------------ Krylov method -------------------------------------
   CASE ('KRYLOV')

      !> Setup basic Krylov solver
      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_KRYLOV
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !> Setup preconditioner for Krylov solver
      NSTACK = NSTACK + 1
      SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)                      

         CASE (NSCARC_RELAX_JACOBI)                                    !> Jacobi preconditioner
            STACK(NSTACK)%SOLVER => PRECON_JACOBI
            CALL SCARC_SETUP_PRECON(NSTACK)

         CASE (NSCARC_RELAX_SSOR)                                      !> SSOR preconditioner
            STACK(NSTACK)%SOLVER => PRECON_SSOR
            CALL SCARC_SETUP_PRECON(NSTACK)

         CASE (NSCARC_RELAX_FFT)                                       !> FFT preconditioner
            STACK(NSTACK)%SOLVER => PRECON_FFT
            CALL SCARC_SETUP_PRECON(NSTACK)
            CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MIN)

         CASE (NSCARC_RELAX_FFT_OVERLAP)                               !> FFT_OVERLAP preconditioner
            STACK(NSTACK)%SOLVER => PRECON_FFT
            CALL SCARC_SETUP_PRECON(NSTACK)
            CALL SCARC_SETUP_FFT_OVERLAP(NLEVEL_MIN, NLEVEL_MIN)

#ifdef WITH_MKL
         CASE (NSCARC_RELAX_PARDISO)                                   !> PARDISO preconditioner
            STACK(NSTACK)%SOLVER => PRECON_PARDISO
            CALL SCARC_SETUP_PRECON(NSTACK)
            CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)

         CASE (NSCARC_RELAX_CLUSTER)                                   !> CLUSTER preconditioner
            STACK(NSTACK)%SOLVER => PRECON_CLUSTER
            CALL SCARC_SETUP_PRECON(NSTACK)
            CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)
#endif

         CASE (NSCARC_RELAX_MULTIGRID)                                 !> Multigrid preconditioner

            STACK(NSTACK)%SOLVER => PRECON_MULTIGRID
            CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_PRECON, NSCARC_SCOPE_TWO, NSTACK, NLEVEL_MIN, NLEVEL_MAX)

            NSTACK = NSTACK + 1
            SELECT CASE (TYPE_SMOOTH)
               CASE (NSCARC_RELAX_JACOBI)                              !> Jacobi smoother
                  STACK(NSTACK)%SOLVER => SMOOTH_JACOBI
                  CALL SCARC_SETUP_SMOOTH(NSTACK)
               CASE (NSCARC_RELAX_SSOR)                                !> SSOR smoother
                  STACK(NSTACK)%SOLVER => SMOOTH_SSOR
                  CALL SCARC_SETUP_SMOOTH(NSTACK)
               CASE (NSCARC_RELAX_FFT)                                 !> FFT smoother
                  STACK(NSTACK)%SOLVER => SMOOTH_FFT
                  CALL SCARC_SETUP_SMOOTH(NSTACK)
                  CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MAX)
               CASE (NSCARC_RELAX_FFT_OVERLAP)                         !> FFT_OVERLAP smoother
                  STACK(NSTACK)%SOLVER => SMOOTH_FFT
                  CALL SCARC_SETUP_SMOOTH(NSTACK)
                  CALL SCARC_SETUP_FFT_OVERLAP(NLEVEL_MIN, NLEVEL_MAX)
#ifdef WITH_MKL
               CASE (NSCARC_RELAX_PARDISO)                             !> PARDISO smoother
                  STACK(NSTACK)%SOLVER => SMOOTH_PARDISO
                  CALL SCARC_SETUP_SMOOTH(NSTACK)
                  CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
#endif
            END SELECT

            NSTACK = NSTACK + 1                                         !> coarse grid solver of MG-preconditioner
            CALL SCARC_SETUP_COARSESOLVER(NSCARC_SCOPE_TWO, NSTACK, NLEVEL_MAX, NLEVEL_MAX)           

      END SELECT SELECT_KRYLOV_PRECON

      !> Twolevel-Krylov: allocate intermediate structures for interpolation and workspace for coarse solver
      IF (BTWOLEVEL) THEN

         CALL SCARC_SETUP_INTERPOLATION(NSCARC_SCOPE_ONE, NLEVEL_MIN+1, NLEVEL_MAX)

         NSTACK = NSTACK + 1
         CALL SCARC_SETUP_COARSESOLVER(NSCARC_SCOPE_ONE, NSTACK, NLEVEL_MAX, NLEVEL_MAX)                     

      ENDIF

   !> ------------------ Multigrid method -------------------------------------
   CASE ('MULTIGRID')

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_MULTIGRID
      CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MAX)

      NSTACK = NSTACK + 1
      SELECT CASE(TYPE_SMOOTH)
         CASE (NSCARC_RELAX_JACOBI)                                   !> Jacobi smoother
            STACK(NSTACK)%SOLVER => SMOOTH_JACOBI
            CALL SCARC_SETUP_SMOOTH(NSTACK)
         CASE (NSCARC_RELAX_SSOR)                                     !> SSOR smoother
            STACK(NSTACK)%SOLVER => SMOOTH_SSOR
            CALL SCARC_SETUP_SMOOTH(NSTACK)
         CASE (NSCARC_RELAX_FFT)                                      !> FFT smoother
            STACK(NSTACK)%SOLVER => SMOOTH_FFT
            CALL SCARC_SETUP_SMOOTH(NSTACK)
            CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MAX)
         CASE (NSCARC_RELAX_FFT_OVERLAP)                              !> FFT_OVERLAP smoother
            STACK(NSTACK)%SOLVER => SMOOTH_FFT
            CALL SCARC_SETUP_SMOOTH(NSTACK)
            CALL SCARC_SETUP_FFT_OVERLAP(NLEVEL_MIN, NLEVEL_MAX)
#ifdef WITH_MKL
         CASE (NSCARC_RELAX_PARDISO)                                  !> PARDISO smoother
            STACK(NSTACK)%SOLVER => SMOOTH_PARDISO
            CALL SCARC_SETUP_SMOOTH(NSTACK)
            CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
         CASE (NSCARC_RELAX_CLUSTER)                                  !> CLUSTER smoother
            STACK(NSTACK)%SOLVER => SMOOTH_CLUSTER
            CALL SCARC_SETUP_SMOOTH(NSTACK)
            CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)
#endif
      END SELECT

      NSTACK = NSTACK + 1                                              !> coarse grid solver
      CALL SCARC_SETUP_COARSESOLVER(NSCARC_SCOPE_ONE, NSTACK, NLEVEL_MAX, NLEVEL_MAX)   


#ifdef WITH_MKL
   !> ------------------ MKL method -------------------------------------
   CASE ('MKL')

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_LUDECOMP
      CALL SCARC_SETUP_LUDECOMP(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      SELECT_MKL: SELECT CASE (TYPE_LUDECOMP)
         CASE (NSCARC_MKL_GLOBAL)
            CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)
         CASE (NSCARC_MKL_LOCAL)
            CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
      END SELECT SELECT_MKL

#endif

END SELECT SELECT_METHOD

!> Store total number of stack entries (used solvers)
N_STACK_TOTAL = NSTACK

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_METHODS

!> ----------------------------------------------------------------------------------------------------
!> Set description pointers to solution vectors related to used scope (main/relax)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POINTERS(BX, BF, BD, BG, BW, BY, BZ, BE, NSTACK)
LOGICAL, INTENT(IN) :: BX, BF, BD, BG, BW, BY, BZ, BE
INTEGER, INTENT(IN) :: NSTACK
TYPE (SCARC_POINTERS_TYPE), POINTER :: PTR
TYPE (SCARC_TYPES_TYPE)   , POINTER :: TYP

PTR => STACK(NSTACK)%SOLVER%POINTERS
TYP => STACK(NSTACK)%SOLVER%TYPES

SELECT CASE (TYP%TYPE_SCOPE)
   CASE (NSCARC_SCOPE_ONE)
      IF (BX) PTR%X = NSCARC_VECTOR_ONE_X
      IF (BF) PTR%F = NSCARC_VECTOR_ONE_F
      IF (BD) PTR%D = NSCARC_VECTOR_ONE_D
      IF (BG) PTR%G = NSCARC_VECTOR_ONE_G
      IF (BW) PTR%W = NSCARC_VECTOR_ONE_W
      IF (BY) PTR%Y = NSCARC_VECTOR_ONE_Y
      IF (BZ) PTR%Z = NSCARC_VECTOR_ONE_Z
      IF (BE) PTR%E = NSCARC_VECTOR_ONE_E
   CASE (NSCARC_SCOPE_TWO)
      IF (BX) PTR%X = NSCARC_VECTOR_TWO_X
      IF (BF) PTR%F = NSCARC_VECTOR_TWO_F
      IF (BD) PTR%D = NSCARC_VECTOR_TWO_D
      IF (BG) PTR%G = NSCARC_VECTOR_TWO_G
      IF (BW) PTR%W = NSCARC_VECTOR_TWO_W
      IF (BY) PTR%Y = NSCARC_VECTOR_TWO_Y
      IF (BZ) PTR%Z = NSCARC_VECTOR_TWO_Z
      IF (BE) PTR%E = NSCARC_VECTOR_TWO_E
END SELECT

END SUBROUTINE SCARC_SETUP_POINTERS

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_VECTORS()
INTEGER :: NM, ISTACK, NL
TYPE (SCARC_LEVEL_TYPE)   , POINTER :: L
TYPE (SCARC_POINTERS_TYPE), POINTER :: PTR
TYPE (SCARC_TYPES_TYPE)   , POINTER :: TYP
TYPE (SCARC_SCOPE_TYPE)   , POINTER :: SCO
TYPE (SCARC_SOLVER_TYPE)  , POINTER :: SOL

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_VECTORS')

DO ISTACK = 1, N_STACK_TOTAL

   SOL => STACK(ISTACK)%SOLVER
   TYP => SOL%TYPES
   PTR => SOL%POINTERS

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      DO NL = TYP%TYPE_NLMIN, TYP%TYPE_NLMAX

         L   => SCARC(NM)%LEVEL(NL)
         SCO => SCARC(NM)%SCOPE(TYP%TYPE_SCOPE, NL)

         CALL SCARC_ALLOCATE_REAL1(SCO%X, 1, L%NCE, NSCARC_INIT_ZERO, 'X')
         CALL SCARC_ALLOCATE_REAL1(SCO%F, 1, L%NCE, NSCARC_INIT_ZERO, 'F')
         CALL SCARC_ALLOCATE_REAL1(SCO%D, 1, L%NCE, NSCARC_INIT_ZERO, 'D')
         CALL SCARC_ALLOCATE_REAL1(SCO%G, 1, L%NCE, NSCARC_INIT_ZERO, 'G')
         CALL SCARC_ALLOCATE_REAL1(SCO%W, 1, L%NCE, NSCARC_INIT_ZERO, 'W')
         CALL SCARC_ALLOCATE_REAL1(SCO%Y, 1, L%NCE, NSCARC_INIT_ZERO, 'Y')
         CALL SCARC_ALLOCATE_REAL1(SCO%Z, 1, L%NCE, NSCARC_INIT_ZERO, 'Z')
         CALL SCARC_ALLOCATE_REAL1(SCO%E, 1, L%NCE, NSCARC_INIT_ZERO, 'Z')

#ifdef WITH_MKL_FB
         CALL SCARC_ALLOCATE_REAL1_FB(SCO%G_FB, 1, L%NCE, NSCARC_INIT_ZERO, 'G_FB')
         CALL SCARC_ALLOCATE_REAL1_FB(SCO%W_FB, 1, L%NCE, NSCARC_INIT_ZERO, 'W_FB')
#endif

      ENDDO
   ENDDO
ENDDO

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_VECTORS

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_KRYLOV(NSOLVER, NSCOPE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTACK, NLMIN, NLMAX
TYPE (SCARC_SOLVER_TYPE)     , POINTER :: SOL
TYPE (SCARC_TYPES_TYPE)      , POINTER :: TYP
TYPE (SCARC_CONVERGENCE_TYPE), POINTER :: CON

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_KRYLOV')

SOL => STACK(NSTACK)%SOLVER

CON => SOL%CONVREQS
TYP => SOL%TYPES

!> Preset types for Krylov method
TYP%TYPE_METHOD    = NSCARC_METHOD_KRYLOV
TYP%TYPE_SOLVER    = NSOLVER
TYP%TYPE_SCOPE     = NSCOPE
TYP%TYPE_NLMIN     = NLMIN
TYP%TYPE_NLMAX     = NLMAX
TYP%TYPE_INTERPOL  = TYPE_INTERPOL
TYP%TYPE_ACCURACY  = TYPE_ACCURACY
TYP%TYPE_PRECISION = TYPE_PRECISION

!> Preset iteration parameters for Krylov method
SELECT CASE(NSOLVER)

   CASE (NSCARC_SOLVER_MAIN)                               !> Used as main solver 

      SOL%CNAME = 'SCARC_MAIN_KRYLOV'

      CON%EPS = SCARC_KRYLOV_ACCURACY
      CON%NIT = SCARC_KRYLOV_ITERATIONS

      TYP%TYPE_RELAX    = TYPE_PRECON                      !> use specified preconditioner
      TYP%TYPE_TWOLEVEL = TYPE_TWOLEVEL                    !> use specified number of levels

   CASE (NSCARC_SOLVER_COARSE)                             !> Used as coarse grid solver 

      SOL%CNAME = 'SCARC_COARSE_KRYLOV'

      CON%EPS = SCARC_COARSE_ACCURACY
      CON%NIT = SCARC_COARSE_ITERATIONS

      TYP%TYPE_RELAX    = NSCARC_RELAX_SSOR                !> only use SSOR-preconditioning for coarse solver
      TYP%TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE             !> only use one level for coarse solver
 
   CASE DEFAULT                                            !> No other choices possible

      CALL SCARC_SHUTDOWN('Error with input parameter ', 'NONE', NSOLVER)

END SELECT

!> Point to solution vectors (in corresponding scope)
CALL SCARC_SETUP_POINTERS(.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE., NSTACK)
CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_KRYLOV

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Geometric Multigrid method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MULTIGRID(NSOLVER, NSCOPE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTACK, NLMIN, NLMAX
TYPE (SCARC_SOLVER_TYPE)     , POINTER :: SOL
TYPE (SCARC_TYPES_TYPE)      , POINTER :: TYP
TYPE (SCARC_CONVERGENCE_TYPE), POINTER :: CON

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_KRYLOV')

SOL => STACK(NSTACK)%SOLVER
CON => SOL%CONVREQS
TYP => SOL%TYPES

!> Preset types for Multigrid method
TYP%TYPE_METHOD    = NSCARC_METHOD_MULTIGRID
TYP%TYPE_SOLVER    = NSOLVER
TYP%TYPE_SCOPE     = NSCOPE 
TYP%TYPE_NLMIN     = NLMIN
TYP%TYPE_NLMAX     = NLMAX
TYP%TYPE_RELAX     = TYPE_SMOOTH
TYP%TYPE_INTERPOL  = TYPE_INTERPOL
TYP%TYPE_ACCURACY  = TYPE_ACCURACY
TYP%TYPE_CYCLING   = TYPE_CYCLING
TYP%TYPE_PRECISION = TYPE_PRECISION

SELECT CASE(NSOLVER)
   CASE (NSCARC_SOLVER_MAIN)                                  !> Used as main solver 
      SOL%CNAME = 'SCARC_MAIN_MULTIGRID'
   CASE (NSCARC_SOLVER_PRECON)                                !> Used as preconditioner in Krylov-method
      SOL%CNAME = 'SCARC_PRECON_MULTIGRID'
   CASE DEFAULT                                               !> No other choices possible
      CALL SCARC_SHUTDOWN('Error with input parameter ', 'NONE', NSOLVER)
END SELECT

!> Preset iteration parameters for Multigrid method
CON%EPS = SCARC_MULTIGRID_ACCURACY
CON%NIT = SCARC_MULTIGRID_ITERATIONS

!> Point to solution vectors (in corresponding scope)
CALL SCARC_SETUP_POINTERS(.TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE., NSTACK)
CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_MULTIGRID

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for MKL-methods
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSESOLVER(NSCOPE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN)    :: NSCOPE, NLMIN, NLMAX
INTEGER, INTENT(INOUT) :: NSTACK

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_COARSESOLVER')

SELECT_COARSE: SELECT CASE (TYPE_COARSE)
   CASE (NSCARC_COARSE_ITERATIVE)
      STACK(NSTACK)%SOLVER => COARSE_KRYLOV
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_COARSE, NSCOPE, NSTACK, NLMIN, NLMAX)
      NSTACK = NSTACK + 1
      TYPE_PRECON = NSCARC_RELAX_SSOR
      STACK(NSTACK)%SOLVER => PRECON_SSOR
      CALL SCARC_SETUP_PRECON(NSTACK)
#ifdef WITH_MKL
   CASE (NSCARC_COARSE_DIRECT)
      IF (N_MPI_PROCESSES > 1) THEN
         STACK(NSTACK)%SOLVER => COARSE_CLUSTER
         CALL SCARC_SETUP_LUDECOMP(NSCARC_SOLVER_COARSE, NSCOPE, NSTACK, NLMIN, NLMAX)
         CALL SCARC_SETUP_CLUSTER(NLMIN, NLMAX)
      ELSE
         STACK(NSTACK)%SOLVER => COARSE_PARDISO
         CALL SCARC_SETUP_LUDECOMP(NSCARC_SOLVER_COARSE, NSCOPE, NSTACK, NLMIN, NLMAX)
         CALL SCARC_SETUP_PARDISO(NLMIN, NLMAX)
      ENDIF
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN('Error with input parameter ', 'NONE', TYPE_COARSE)
END SELECT SELECT_COARSE
CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_COARSESOLVER


#ifdef WITH_MKL
!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for MKL-methods
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LUDECOMP(NSOLVER, NSCOPE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTACK, NLMIN, NLMAX
TYPE (SCARC_SOLVER_TYPE)  , POINTER :: SOL
TYPE (SCARC_TYPES_TYPE)   , POINTER :: TYP

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_LUDECOMP')

SOL => STACK(NSTACK)%SOLVER
SELECT CASE (NSOLVER)
   CASE (NSCARC_SOLVER_MAIN) 
      SOL%CNAME = 'SCARC_MAIN_MKL'
   CASE (NSCARC_SOLVER_COARSE) 
      SOL%CNAME = 'SCARC_COARSE_MKL'
   CASE DEFAULT
      CALL SCARC_SHUTDOWN('Error with input parameter ', 'NONE', NSOLVER)
END SELECT

TYP => SOL%TYPES

!> Preset types for LU-decomposition method
TYP%TYPE_METHOD    = NSCARC_METHOD_LUDECOMP
TYP%TYPE_SOLVER    = NSOLVER
TYP%TYPE_SCOPE     = NSCOPE
TYP%TYPE_NLMIN     = NLMIN
TYP%TYPE_NLMAX     = NLMAX
TYP%TYPE_PRECISION = TYPE_PRECISION

!> Point to solution vectors (in corresponding scope)
CALL SCARC_SETUP_POINTERS(.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_LUDECOMP
#endif

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PRECON(NSTACK)
INTEGER, INTENT(IN) :: NSTACK
TYPE (SCARC_TYPES_TYPE)      , POINTER :: TYP, TYPC
TYPE (SCARC_SOLVER_TYPE)     , POINTER :: SOL, SOLC
TYPE (SCARC_POINTERS_TYPE)   , POINTER :: PTR, PTRC
TYPE (SCARC_CONVERGENCE_TYPE), POINTER :: CON

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_PRECON')

SOL => STACK(NSTACK)%SOLVER
CON => SOL%CONVREQS

SELECT CASE(TYPE_PRECON)
   CASE (NSCARC_RELAX_JACOBI)
      SOL%CNAME = 'SCARC_PRECON_JACOBI'
      CON%EPS   =  SCARC_PRECON_ACCURACY                
      CON%NIT   =  SCARC_PRECON_ITERATIONS
      CON%OMEGA =  SCARC_PRECON_OMEGA
   CASE (NSCARC_RELAX_SSOR)
      SOL%CNAME = 'SCARC_PRECON_SSOR'
      CON%EPS   =  SCARC_PRECON_ACCURACY               
      CON%NIT   =  SCARC_PRECON_ITERATIONS
      CON%OMEGA =  SCARC_PRECON_OMEGA
   CASE (NSCARC_RELAX_FFT)
      SOL%CNAME = 'SCARC_PRECON_FFT'
      CON%EPS   =  SCARC_PRECON_ACCURACY              
      CON%NIT   =  SCARC_PRECON_ITERATIONS
      !CON%OMEGA =  SCARC_PRECON_OMEGA
      CON%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_FFT_OVERLAP)
      SOL%CNAME = 'SCARC_PRECON_FFT_OVERLAP'
      CON%EPS   =  SCARC_PRECON_ACCURACY              
      CON%NIT   =  SCARC_PRECON_ITERATIONS
      !CON%OMEGA =  SCARC_PRECON_OMEGA
      CON%OMEGA = 1.0_EB
#ifdef WITH_MKL
   CASE (NSCARC_RELAX_PARDISO)
      SOL%CNAME = 'SCARC_PRECON_PARDISO'
      CON%EPS   =  SCARC_PRECON_ACCURACY             
      CON%NIT   =  SCARC_PRECON_ITERATIONS
      !CON%OMEGA =  SCARC_PRECON_OMEGA
      CON%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_CLUSTER)
      SOL%CNAME = 'SCARC_PRECON_CLUSTER'
      CON%EPS   =  SCARC_PRECON_ACCURACY            
      CON%NIT   = 1
      CON%OMEGA = 1.0_EB
#endif
   CASE DEFAULT                                 
      CALL SCARC_SHUTDOWN('Error with input parameter ', 'NONE', TYPE_PRECON)
END SELECT

!> Preset iteration parameters for Multigrid method
PTR => SOL%POINTERS
TYP => SOL%TYPES

IF (NSTACK > 1) THEN
   SOLC => STACK(NSTACK-1)%SOLVER                !> point to calling solver
   TYPC => SOLC%TYPES
   PTRC => SOLC%POINTERS
ELSE
   CALL SCARC_SHUTDOWN('Wrong number of solvers in stack ', 'NONE', NSTACK)
ENDIF

!> Preset types for preconditioner
TYP%TYPE_SOLVER    = TYPC%TYPE_SOLVER 
TYP%TYPE_SCOPE     = TYPC%TYPE_SCOPE 
TYP%TYPE_NLMIN     = TYPC%TYPE_NLMIN 
TYP%TYPE_NLMAX     = TYPC%TYPE_NLMAX 
TYP%TYPE_RELAX     = TYPC%TYPE_RELAX 
TYP%TYPE_INTERPOL  = TYPC%TYPE_INTERPOL
TYP%TYPE_ACCURACY  = TYPC%TYPE_ACCURACY
TYP%TYPE_PRECISION = TYPC%TYPE_PRECISION

!> Preset pointers for preconditioner
PTR%X = PTRC%X                                    !> use same pointers as calling Krylov-solver
PTR%F = PTRC%F
PTR%D = PTRC%D
PTR%G = PTRC%G
PTR%W = PTRC%W
PTR%Y = PTRC%Y
PTR%Z = PTRC%Z
PTR%E = PTRC%E

#ifdef WITH_MKL_FB
PTR%X_FB = PTRC%X_FB                              !> use same pointers as calling Krylov-solver
PTR%F_FB = PTRC%F_FB
PTR%D_FB = PTRC%D_FB
PTR%G_FB = PTRC%G_FB
PTR%W_FB = PTRC%W_FB
#endif

CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_PRECON

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SMOOTH(NSTACK)
INTEGER, INTENT(IN) :: NSTACK
TYPE (SCARC_TYPES_TYPE)      , POINTER :: TYP, TYPC
TYPE (SCARC_SOLVER_TYPE)     , POINTER :: SOL, SOLC
TYPE (SCARC_POINTERS_TYPE)   , POINTER :: PTR, PTRC
TYPE (SCARC_CONVERGENCE_TYPE), POINTER :: CON

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_SMOOTH')

SOL => STACK(NSTACK)%SOLVER

SELECT CASE(TYPE_SMOOTH)
   CASE (NSCARC_RELAX_JACOBI)
      SOL%CNAME = 'SCARC_SMOOTH_JACOBI'
   CASE (NSCARC_RELAX_SSOR)
      SOL%CNAME = 'SCARC_SMOOTH_SSOR'
   CASE (NSCARC_RELAX_FFT)
      SOL%CNAME = 'SCARC_SMOOTH_FFT'
   CASE (NSCARC_RELAX_FFT_OVERLAP)
      SOL%CNAME = 'SCARC_SMOOTH_FFT_OVERLAP'
#ifdef WITH_MKL
   CASE (NSCARC_RELAX_PARDISO)
      SOL%CNAME = 'SCARC_SMOOTH_PARDISO'
   CASE (NSCARC_RELAX_CLUSTER)
      SOL%CNAME = 'SCARC_SMOOTH_CLUSTER'
#endif
   CASE DEFAULT                                               !> No other choices possible
      CALL SCARC_SHUTDOWN('Error with input parameter ', 'NONE', TYPE_SMOOTH)
END SELECT

!> Preset iteration parameters for Multigrid method
CON => SOL%CONVREQS
CON%EPS   = SCARC_SMOOTH_ACCURACY                 !> set iteration parameters for preconditioner
CON%NIT   = SCARC_SMOOTH_ITERATIONS
CON%OMEGA = SCARC_SMOOTH_OMEGA

PTR => SOL%POINTERS
TYP => SOL%TYPES

IF (NSTACK > 1) THEN
   SOLC => STACK(NSTACK-1)%SOLVER                !> point to calling solver
   TYPC => SOLC%TYPES
   PTRC => SOLC%POINTERS
ELSE
   CALL SCARC_SHUTDOWN('Wrong number of solvers in stack ' , 'NONE', NSTACK)
ENDIF

!> Preset types for preconditioner
TYP%TYPE_SOLVER    = NSCARC_SOLVER_SMOOTH
TYP%TYPE_SCOPE     = TYPC%TYPE_SCOPE 
TYP%TYPE_NLMIN     = TYPC%TYPE_NLMIN  
TYP%TYPE_NLMAX     = TYPC%TYPE_NLMAX  
TYP%TYPE_RELAX     = TYPC%TYPE_RELAX  
TYP%TYPE_INTERPOL  = TYPC%TYPE_INTERPOL
TYP%TYPE_ACCURACY  = TYPC%TYPE_ACCURACY
TYP%TYPE_PRECISION = TYPC%TYPE_PRECISION

PTR%X = PTRC%X                                    !> use same pointers as calling Krylov-solver
PTR%F = PTRC%F
PTR%D = PTRC%D
PTR%G = PTRC%G
PTR%W = PTRC%W
PTR%Y = PTRC%Y
PTR%Z = PTRC%Z
PTR%E = PTRC%E

#ifdef WITH_MKL_FB
PTR%X_FB = PTRC%X_FB                              !> use same pointers as calling Krylov-solver, single precision
PTR%F_FB = PTRC%F_FB
PTR%D_FB = PTRC%D_FB
PTR%G_FB = PTRC%G_FB
PTR%W_FB = PTRC%W_FB
#endif

CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_SMOOTH

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for additive or multiplicative coarse grid 
!> (corresponding to Schwarz domain decomposition method)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERPOLATION(NSCOPE, NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NSCOPE, NLMIN, NLMAX
INTEGER :: NM, NL
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_SCOPE_TYPE), POINTER :: SCO

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_INTERPOLATION')

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L   => SCARC(NM)%LEVEL(NL)
      SCO => SCARC(NM)%SCOPE(NSCOPE, NL)
      
      CALL SCARC_ALLOCATE_REAL1(SCO%X, 1, L%NCE, NSCARC_INIT_ZERO, 'X')
      CALL SCARC_ALLOCATE_REAL1(SCO%F, 1, L%NCE, NSCARC_INIT_ZERO, 'F')
      CALL SCARC_ALLOCATE_REAL1(SCO%W, 1, L%NCE, NSCARC_INIT_ZERO, 'W')
      CALL SCARC_ALLOCATE_REAL1(SCO%G, 1, L%NCE, NSCARC_INIT_ZERO, 'G')
      CALL SCARC_ALLOCATE_REAL1(SCO%Y, 1, L%NCE, NSCARC_INIT_ZERO, 'Y')
      CALL SCARC_ALLOCATE_REAL1(SCO%Z, 1, L%NCE, NSCARC_INIT_ZERO, 'Z')

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_INTERPOLATION

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for blockwise FFT methods
!> New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFT(NLMIN, NLMAX)
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IERR = 0
TYPE (SCARC_TYPE), POINTER :: S
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_FFT_TYPE)  , POINTER :: F
TYPE (MESH_TYPE)  , POINTER :: M

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_FFT')

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      S => SCARC(NM)
      L => SCARC(NM)%LEVEL(NL)
      F => SCARC(NM)%FFT(NL)
      M => MESHES(NM)
      
      !> Allocate working space for FFT routine
      F%LBC = M%LBC
      F%MBC = M%MBC
      F%NBC = M%NBC

      F%ITRN = L%NX+1
      IF (TWO_D) THEN
         F%JTRN = 1
      ELSE
         F%JTRN = L%NY+1
      ENDIF
      F%KTRN = L%NZ+1

      F%LSAVE = (F%ITRN+1)*F%JTRN*F%KTRN+7*F%ITRN+5*F%JTRN+6*F%KTRN+56
      F%LWORK = (F%ITRN+1)*F%JTRN*F%KTRN

      CALL SCARC_ALLOCATE_REAL1(F%SAVE1, -3, F%LSAVE, NSCARC_INIT_ZERO, 'FFT')
      CALL SCARC_ALLOCATE_REAL1(F%WORK ,  1, F%LWORK, NSCARC_INIT_ZERO, 'FFT')

      !> Allocate stretching vector (set to 1)
      CALL SCARC_ALLOCATE_REAL1(F%HX, 1, L%NX+1, NSCARC_INIT_ONE, 'FFT')

      !> Allocate RHS vector for FFT routine
      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3(F%PRHS, 1, L%NX+1, 1, 1,      1, L%NZ+1, NSCARC_INIT_ZERO, 'FFT')
      ELSE
         CALL SCARC_ALLOCATE_REAL3(F%PRHS, 1, L%NX+1, 1, L%NY+1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'FFT')
      ENDIF
      
      !> Allocate boundary data vector for XS
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(F%BXS, 1, L%NY+1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'BXS')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BXS, 1,      1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'BXS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(F%BXS, 1, L%NY+1, 1, 1, NSCARC_INIT_ZERO, 'BXS')
      ENDIF
      
      !> Allocate boundary data vector for XF
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(F%BXF, 1, L%NY+1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'BXF')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BXF, 1,      1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'BXF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(F%BXF, 1, L%NY+1, 1, 1, NSCARC_INIT_ZERO, 'BXF')
      ENDIF
      
      !> Allocate boundary data vector for YS
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(F%BYS, 1, L%NX+1,1, L%NZ+1, NSCARC_INIT_ZERO, 'BYS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(F%BYS, 1, L%NX+1,1,      1, NSCARC_INIT_ZERO, 'BYS')
      ENDIF
      
      !> Allocate boundary data vector for YF
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(F%BYF, 1, L%NX+1,1, L%NZ+1, NSCARC_INIT_ZERO, 'BYF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(F%BYF, 1, L%NX+1,1,      1, NSCARC_INIT_ZERO, 'BYF')
      ENDIF
      
      !> Allocate boundary data vector for ZS
      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(F%BZS, 1, L%NX+1, 1, L%NY+1, NSCARC_INIT_ZERO, 'BZS')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BZS, 1, L%NX+1, 1,      1, NSCARC_INIT_ZERO, 'BZS')
      
      !> Allocate boundary data vector for ZF
      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(F%BZF, 1, L%NX+1, 1, L%NY+1, NSCARC_INIT_ZERO, 'BZF')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BZF, 1, L%NX+1, 1,      1, NSCARC_INIT_ZERO, 'BZF')
      
      IF (TWO_D) THEN
         CALL H2CZIS(S%XS,S%XF,S%IBAR,F%LBC,S%ZS,S%ZF,S%KBAR,F%NBC,F%HX,F%XLM,F%ITRN,IERR,F%SAVE1)
      ELSE
         CALL H3CZIS(S%XS,S%XF,S%IBAR,F%LBC,S%YS,S%YF,S%JBAR,F%MBC,S%ZS,S%ZF,S%KBAR,F%NBC,&
                     F%HX,F%XLM,F%ITRN,F%JTRN,IERR,F%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_FFT

!> ----------------------------------------------------------------------------------------------------
!> Still experimental code: 
!> Perform FFT-preconditioning with one layer overlap 
!> and use Dirichlet BC at internal (overlapping) boundaries
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFT_OVERLAP(NLMIN, NLMAX)
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IERR = 0
TYPE (SCARC_TYPE), POINTER :: S
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_FFT_TYPE)  , POINTER :: F
TYPE (MESH_TYPE)  , POINTER :: M

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_FFT_OVERLAP')

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      S => SCARC(NM)
      L => SCARC(NM)%LEVEL(NL)
      F => SCARC(NM)%FFT(NL)
      M => MESHES(NM)
      
      !> Allocate working space for FFT routine
      F%LBC = M%LBC
      F%MBC = M%MBC
      F%NBC = M%NBC

      F%IBAR0 = S%IBAR
      F%JBAR0 = S%JBAR
      F%KBAR0 = S%KBAR

      F%XS0 = S%XS
      F%XF0 = S%XF
      F%YS0 = S%YS
      F%YF0 = S%YF
      F%ZS0 = S%ZS
      F%ZF0 = S%ZF

      IF (L%FACE( 1)%N_NEIGHBORS /= 0) THEN
         F%IBAR0 = F%IBAR0 + 1
         F%IS0 = 1
         F%XS0 = F%XS0 - L%COORD%DX
      ENDIF
      IF (L%FACE(-1)%N_NEIGHBORS /= 0) THEN
         F%IBAR0 = F%IBAR0 + 1
         F%IF0 = 1
         F%XF0 = F%XF0 + L%COORD%DX
      ENDIF
      IF (L%FACE( 2)%N_NEIGHBORS /= 0) THEN
         F%JBAR0 = F%JBAR0 + 1
         F%JS0 = 1
         F%YS0 = F%YS0 - L%COORD%DY
      ENDIF
      IF (L%FACE(-2)%N_NEIGHBORS /= 0) THEN
         F%JBAR0 = F%JBAR0 + 1
         F%JF0 = 1
         F%YF0 = F%YF0 + L%COORD%DY
      ENDIF
      IF (L%FACE( 3)%N_NEIGHBORS /= 0) THEN
         F%KBAR0 = F%KBAR0 + 1
         F%KS0 = 1
         F%ZS0 = F%ZS0 - L%COORD%DZ
      ENDIF
      IF (L%FACE(-3)%N_NEIGHBORS /= 0) THEN
         F%KBAR0 = F%KBAR0 + 1
         F%KF0 = 1
         F%ZF0 = F%ZF0 + L%COORD%DZ
      ENDIF
      F%ITRN0 = F%IBAR0 + 1
      IF (TWO_D) THEN
         F%JTRN0 = 1
      ELSE
         F%JTRN0 = F%JBAR0 + 1
      ENDIF
      F%KTRN0 = F%KBAR0 + 1

      F%LSAVE = (F%ITRN0+1)*F%JTRN0*F%KTRN0+7*F%ITRN0+5*F%JTRN0+6*F%KTRN0+56
      F%LWORK = (F%ITRN0+1)*F%JTRN0*F%KTRN0

      CALL SCARC_ALLOCATE_REAL1(F%SAVE1, -3, F%LSAVE, NSCARC_INIT_ZERO, 'FFT')
      CALL SCARC_ALLOCATE_REAL1(F%WORK ,  1, F%LWORK, NSCARC_INIT_ZERO, 'FFT')

      !> Allocate stretching vector (set to 1)
      CALL SCARC_ALLOCATE_REAL1(F%HX, 1, F%ITRN0, NSCARC_INIT_ONE, 'FFT')

      !> Allocate RHS vector for FFT routine
      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3(F%PRHS, 1, F%ITRN0, 1,       1, 1, F%KTRN0, NSCARC_INIT_ZERO, 'FFT')
      ELSE
         CALL SCARC_ALLOCATE_REAL3(F%PRHS, 1, F%ITRN0, 1, F%JTRN0, 1, F%KTRN0, NSCARC_INIT_ZERO, 'FFT')
      ENDIF
      
      !> Allocate boundary data vector for XS
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(F%BXS, 1, F%JTRN0, 1, F%KTRN0, NSCARC_INIT_ZERO, 'BXS')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BXS, 1,       1, 1, F%KTRN0, NSCARC_INIT_ZERO, 'BXS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(F%BXS, 1, F%JTRN0, 1, 1, NSCARC_INIT_ZERO, 'BXS')
      ENDIF
      
      !> Allocate boundary data vector for XF
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(F%BXF, 1, F%JTRN0, 1, F%KTRN0, NSCARC_INIT_ZERO, 'BXF')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BXF, 1,       1, 1, F%KTRN0, NSCARC_INIT_ZERO, 'BXF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(F%BXF, 1, F%JTRN0, 1, 1, NSCARC_INIT_ZERO, 'BXF')
      ENDIF
      
      !> Allocate boundary data vector for YS
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(F%BYS, 1, F%ITRN0,1, F%KTRN0, NSCARC_INIT_ZERO, 'BYS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(F%BYS, 1, F%ITRN0,1,       1, NSCARC_INIT_ZERO, 'BYS')
      ENDIF
      
      !> Allocate boundary data vector for YF
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(F%BYF, 1, F%ITRN0,1, F%KTRN0, NSCARC_INIT_ZERO, 'BYF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(F%BYF, 1, F%ITRN0,1,       1, NSCARC_INIT_ZERO, 'BYF')
      ENDIF
      
      !> Allocate boundary data vector for ZS
      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(F%BZS, 1, F%ITRN0, 1, F%JTRN0, NSCARC_INIT_ZERO, 'BZS')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BZS, 1, F%ITRN0, 1,       1, NSCARC_INIT_ZERO, 'BZS')
      
      !> Allocate boundary data vector for ZF
      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(F%BZF, 1, F%ITRN0, 1, F%JTRN0, NSCARC_INIT_ZERO, 'BZF')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BZF, 1, F%ITRN0, 1,       1, NSCARC_INIT_ZERO, 'BZF')
      
      IF (TWO_D) THEN
         CALL H2CZIS(F%XS0,F%XF0,F%IBAR0,F%LBC,F%ZS0,F%ZF0,F%KBAR0,F%NBC,F%HX,F%XLM,F%ITRN0,IERR,F%SAVE1)
      ELSE
         CALL H3CZIS(F%XS0,F%XF0,F%IBAR0,F%LBC,F%YS0,F%YF0,F%JBAR0,F%MBC,F%ZS0,F%ZF0,F%KBAR0,F%NBC,&
                     F%HX,F%XLM,F%ITRN0,F%JTRN0,IERR,F%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_FFT_OVERLAP



#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CLUSTER(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I !, IC, IP
REAL (EB) :: TNOW
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_MKL_TYPE)  , POINTER :: MKL
#ifdef WITH_MKL_FB
REAL (FB) :: DUMMY(1)=0.0_FB
TYPE (SCARC_MATRIX_FB_TYPE), POINTER :: ASYM
#else
REAL (EB) :: DUMMY(1)=0.0_EB
TYPE (SCARC_MATRIX_TYPE), POINTER :: ASYM
#endif

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_CLUSTER')
TNOW = CURRENT_TIME()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L    => SCARC(NM)%LEVEL(NL)
      MKL  => SCARC(NM)%MKL(NL)

      !> Allocate workspace for parameters needed in MKL-routine
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL_IPARM')
      
      !> Allocate workspace for pointers needed in MKL-routine
      IF (.NOT.ALLOCATED(MKL%CT)) THEN
         ALLOCATE(MKL%CT(64), STAT=IERROR) 
         CALL CHKMEMERR ('SCARC', 'CT', IERROR)
         DO I=1,64
            MKL%CT(I)%DUMMY = 0
         ENDDO
      ENDIF
      
      !> Define corresponding parameters
      !> Note: IPARM-vectory is allocate from 1:64, not from 0:63
      MKL%NRHS   =  1         ! one right hand side
      MKL%MAXFCT =  1         ! one matrix
      MKL%MNUM   =  1         ! number of matrix to be factorized
      MKL%ERROR  =  0         ! initialize error flag
      MKL%MSGLVL =  0         ! do not print statistical information
      IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
         MKL%MTYPE  = -2         ! Matrix type real and symmetric indefinite
      ELSE
         MKL%MTYPE  = 11         ! Matrix type real and non-symmetric
      ENDIF
      MKL%IPARM(1)  =  1      ! supply own parameters
      MKL%IPARM(2)  =  3      ! supply own parameters
      MKL%IPARM(4)  =  0      ! supply own parameters
      MKL%IPARM(5)  =  0      ! supply own parameters
      MKL%IPARM(6)  =  0      ! write solution to x
      MKL%IPARM(8)  =  2      ! automatic setting of iterative refinement steps
      MKL%IPARM(10) = 13      ! pivoting perturbation
      MKL%IPARM(11) =  1      ! pivoting perturbation
      MKL%IPARM(13) =  1      ! pivoting perturbation
      MKL%IPARM(14) =  0      ! pivoting perturbation
      MKL%IPARM(18) =  0      ! pivoting perturbation
      MKL%IPARM(19) =  0      ! pivoting perturbation
      MKL%IPARM(20) =  0      ! pivoting perturbation
      MKL%IPARM(21) =  1      ! Bunch-Kaufman pivoting which is default in case of IPARM(0)=0
      MKL%IPARM(24) =  0      ! Bunch-Kaufman pivoting which is default in case of IPARM(0)=0
      MKL%IPARM(27) =  1      ! use matrix checker
      MKL%IPARM(40) = 2                                             ! provide matrix in distributed format
      MKL%IPARM(41) = L%CELL%NC_OFFSET(NM) + 1                      ! first global cell number for mesh NM
      MKL%IPARM(42) = L%CELL%NC_OFFSET(NM) + L%CELL%NC_LOCAL(NM)    ! last global cell number for mesh NM
      !MKL%IPARM(39) = 2                                            ! provide matrix in distributed format
      !MKL%IPARM(40) = L%CELL%NC_OFFSET(NM)+1                       ! first global cell number for mesh NM
      !MKL%IPARM(41) = L%CELL%NC_OFFSET(NM)+L%CELL%NC_LOCAL(NM)     ! last global cell number for mesh NM
    
#ifdef WITH_MKL_FB
      ASYM => SCARC(NM)%SYSTEM(NL)%ASYM_FB
      MKL%IPARM(28)=1         ! single precision

      ! perform only reordering and symbolic factorization
      MKL%PHASE = 11
      CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%CELL%NC_GLOBAL, &
                                   ASYM%VAL, ASYM%ROW, ASYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
                     
      ! perform only factorization
      MKL%PHASE = 22 
      CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%CELL%NC_GLOBAL, &
                                   ASYM%VAL, ASYM%ROW, ASYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
         
#else
      ASYM => SCARC(NM)%SYSTEM(NL)%ASYM
      MKL%IPARM(28)=0         ! double precision

      ! perform only reordering and symbolic factorization
      MKL%PHASE = 11
      CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%CELL%NC_GLOBAL, &
                                   ASYM%VAL, ASYM%ROW, ASYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
                     
      ! perform only factorization
      MKL%PHASE = 22 
      CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%CELL%NC_GLOBAL, &
                                   ASYM%VAL, ASYM%ROW, ASYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
               
#endif

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

CALL SCARC_LEAVE_ROUTINE()
TSETUP(MYID+1)%CLUSTER=TSETUP(MYID+1)%CLUSTER+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SETUP_CLUSTER

!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PARDISO(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I, IDUMMY(1)=0 
REAL (EB) :: TNOW
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_MKL_TYPE)  , POINTER :: MKL
#ifdef WITH_MKL_FB
REAL (FB) :: DUMMY(1)=0.0_FB
TYPE (SCARC_MATRIX_FB_TYPE), POINTER :: ASYM
#else
REAL (EB) :: DUMMY(1)=0.0_EB
TYPE (SCARC_MATRIX_TYPE), POINTER :: ASYM
#endif

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_PARDISO')
TNOW = CURRENT_TIME()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L    => SCARC(NM)%LEVEL(NL)
      MKL  => SCARC(NM)%MKL(NL)
      
      !> Allocate workspace for parameters needed in MKL-routine
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL_IPARM')
      
      !> Allocate workspace for pointers needed in MKL-routine
      IF (.NOT.ALLOCATED(MKL%PT)) THEN
         ALLOCATE(MKL%PT(64), STAT=IERROR) 
         CALL CHKMEMERR ('SCARC', 'PT', IERROR)
         DO I=1,64
            MKL%PT(I)%DUMMY = 0
         ENDDO
      ENDIF
      
      !> Define corresponding parameters 
      !> Note: IPARM-vectory is allocate from 1:64, not from 0:63
      MKL%NRHS   = 1
      MKL%MAXFCT = 1
      MKL%MNUM   = 1
      
      MKL%IPARM(1)  =  1      ! no solver default
      MKL%IPARM(2)  =  3      ! parallel (OpenMP) version of the nested dissection algorithm
      MKL%IPARM(4)  =  0      ! factorization computed as required by phase
      MKL%IPARM(5)  =  0      ! user permutation ignored
      MKL%IPARM(6)  =  0      ! write solution on x
      MKL%IPARM(8)  =  2      ! numbers of iterative refinement steps
      MKL%IPARM(10) = 13      ! perturb the pivot elements with 1E-13
      MKL%IPARM(11) =  0      ! disable scaling (default for SPD)
      MKL%IPARM(13) =  0      ! disable matching
      MKL%IPARM(18) = -1      ! Output: number of nonzeros in the factor LU
      MKL%IPARM(19) = -1      ! Output: number of floating points operations
      MKL%IPARM(20) =  1      ! Output: Numbers of CG Iterations
      MKL%IPARM(27) =  1      ! use matrix checker
      MKL%IPARM(37) =  0      ! matrix storage in CSR-format
      
      MKL%ERROR  =  0         ! initialize error flag
      MKL%MSGLVL =  0         ! do not print statistical information
      MKL%MTYPE  = -2         ! Matrix type real non-symmetric
      
#ifdef WITH_MKL_FB
      ASYM => SCARC(NM)%SYSTEM(NL)%ASYM_FB
      MKL%IPARM(28)=1         ! single precision

      ! perform only reordering and symbolic factorization
      MKL%PHASE = 11
      CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NCS, ASYM%VAL, ASYM%ROW, ASYM%COL, &
                     IDUMMY, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
                     
      ! perform only Factorization
      MKL%PHASE = 22 
      CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NCS, ASYM%VAL, ASYM%ROW, ASYM%COL, &
                     IDUMMY, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)

#else
      ASYM => SCARC(NM)%SYSTEM(NL)%ASYM
      MKL%IPARM(28)=0         ! double precision

      ! perform only reordering and symbolic factorization
      MKL%PHASE = 11
      CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NCS, ASYM%VAL, ASYM%ROW, ASYM%COL, &
                     IDUMMY, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
                     
      ! perform only Factorization
      MKL%PHASE = 22 
      CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NCS, ASYM%VAL, ASYM%ROW, ASYM%COL, &
                     IDUMMY, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)

#endif

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP
      
CALL SCARC_LEAVE_ROUTINE()
TSETUP(MYID+1)%PARDISO=TSETUP(MYID+1)%PARDISO+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SETUP_PARDISO
#endif


!> ------------------------------------------------------------------------------------------------
!> Set sizes for transfer matrices
!> ------------------------------------------------------------------------------------------------
SUBROUTINE  SCARC_SETUP_SIZES(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: IW, IC, IG, ICP, ICCE
INTEGER :: NM, NOM
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE)   , POINTER :: L, OL, LF, LC, OLF, OLC
TYPE (SCARC_CELL_TYPE)    , POINTER :: CC, CF
TYPE (SCARC_MATRIX_TYPE)  , POINTER :: A, AF, OA
TYPE (SCARC_MAPPING_TYPE) , POINTER :: MC
TYPE (SCARC_EXCHANGE_TYPE), POINTER :: OX

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_SIZES')

SELECT CASE (NTYPE)

   !> --------------------------------------------------------------------------------------------------
   !> Define sizes for system matrix A (including extended regions related to overlaps)
   !> --------------------------------------------------------------------------------------------------
   CASE (NSCARC_SIZE_MATRIX)

      LEVEL_SYSTEM_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%SYSTEM(NL)%A

         IF (TWO_D) THEN
            A%NSTENCIL = 5
         ELSE
            A%NSTENCIL = 7
         ENDIF

         A%NAV = L%NCS * A%NSTENCIL
         A%NAC = A%NAV 
         A%NAR = L%NCS + 1
         A%NAS = L%NCS * A%NSTENCIL

         !> Determine sizes of overlapped parts for later communication with corresponding neighbors
         DO IW = 1, L%NW
            NOM = L%WALL(IW)%NOM
            IF (NOM /= 0) THEN
               OA => SCARC(NM)%OSCARC(NOM)%SYSTEM(NL)%A
               OA%NAV = OA%NAV + A%NSTENCIL
            ENDIF
         ENDDO

      ENDDO LEVEL_SYSTEM_MESHES_LOOP

      !> Exchange sizes for AMG matrices and vectors
      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SIZE, NL)

      !> Determine extended sizes for extended prolongation and restriction matrix (AMG only)
      IF (BAMG) THEN
         LEVEL_SYSTEM_MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
            L => SCARC(NM)%LEVEL(NL)
            A => SCARC(NM)%SYSTEM(NL)%A
            A%NAE = A%NAV
   
            CALL SCARC_ALLOCATE_INT2(L%AMG%MAPPING, 1, 3, L%NCS+1, 2*L%NCE, NSCARC_INIT_ZERO, 'AMG%MAP')
   
            LEVEL_OTHER_MESHES_LOOP2: DO NOM = 1, NMESHES
   
               !IF (NOM == NM) CYCLE LEVEL_OTHER_MESHES_LOOP2
               OS => SCARC(NM)%OSCARC(NOM)
               OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
               OA => SCARC(NM)%OSCARC(NOM)%SYSTEM(NL)%A
               OX => OS%EXCHANGE
               IF (OX%NICMAX_S==0 .AND. OX%NICMAX_R==0) CYCLE LEVEL_OTHER_MESHES_LOOP2
   
               OA => OS%SYSTEM(NL)%A
               A%NAE  = A%NAE  + OA%NAV
   
               CALL SCARC_ALLOCATE_INT1(OA%STENCIL, 1, OL%NCG, NSCARC_INIT_ZERO, 'OA_STENCIL')
   
            ENDDO LEVEL_OTHER_MESHES_LOOP2
   
         ENDDO LEVEL_SYSTEM_MESHES_LOOP2
      ENDIF

   !> -------------------------------------------------------------------------------------------
   !> Define sizes for transfer matrices P and R (including extended regions related to overlaps)
   !> -------------------------------------------------------------------------------------------
   CASE (NSCARC_SIZE_TRANSFER)

      TRANSFER_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         LF => SCARC(NM)%LEVEL(NL)
         LC => SCARC(NM)%LEVEL(NL+1)

         AF => SCARC(NM)%SYSTEM(NL)%A

         CF => SCARC(NM)%LEVEL(NL)%CELL
         CC => SCARC(NM)%LEVEL(NL+1)%CELL

         !> Determine dimensions of restriction and prolongation matrices in own mesh
         LF%AMG%NCF  = 0
         LF%AMG%NCC  = 0
         LF%AMG%NP   = 0
         LF%AMG%NR   = 0
         LF%AMG%NCCI = 0

         DO IC = 1, LF%NC
            IF (LF%AMG%CTYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
               LF%AMG%NCC = LF%AMG%NCC + 1
               LF%AMG%NP  = LF%AMG%NP  + 1
               LF%AMG%NR  = LF%AMG%NP
    !>         F%CELL%STYLE(IC)  = F%AMG%NCC
            ELSE
               LF%AMG%NCF = LF%AMG%NCF + 1
               LF%AMG%NP  = LF%AMG%NP  + AF%ROW(IC+1)-AF%ROW(IC) - 1
               LF%AMG%NR  = LF%AMG%NP
            ENDIF
         ENDDO

    !>   DF%AMG%NCCE = F%AMG%NCC
    !>   F%AMG%NCFE = F%AMG%NCF
    !>   DO IC = F%NC+1, F%NCE
    !>      IF (F%CELL%STYLE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
    !>         F%AMG%NCCE = F%AMG%NCCE + 1
    !>         F%CELL%STYLE(IC)  = F%AMG%NCCE
    !>      ELSE
    !>         F%AMG%NCFE = F%AMG%NCFE + 1
    !>      ENDIF
    !>   ENDDO

         LF%AMG%NCCI = LF%AMG%NCC
         LF%AMG%NCCE = LF%AMG%NCC

         !> Determine number of coarse and fine cells and check correctness of computation
         IF (LF%AMG%NCC + LF%AMG%NCF /= LF%NC) &
            CALL SCARC_SHUTDOWN('Inconsistent number of coarse and fine cells ', 'NONE', LF%AMG%NCC)

         !> define variables for overlapping parts
         TRANSFER_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

            !IF (NOM == NM) CYCLE TRANSFER_OTHER_MESHES_LOOP
            OS  => SCARC(NM)%OSCARC(NOM)
            OLF => OS%LEVEL(NL)
            OX  => OS%EXCHANGE
            IF (OX%NICMAX_S==0 .AND. OX%NICMAX_R==0) CYCLE TRANSFER_OTHER_MESHES_LOOP

            OLF%AMG%NP   = 0
            OLF%AMG%NR   = 0
            OLF%AMG%NPS  = 0
            OLF%AMG%NRS  = 0
            OLF%AMG%NCC  = 0
            OLF%AMG%NCF  = 0
            OLF%AMG%NCCS = 0
            OLF%AMG%NCFS = 0
            OLF%AMG%ICG0 = 0

         ENDDO TRANSFER_OTHER_MESHES_LOOP

         !> Determine sizes of overlapped parts for later communication with corresponding neighbors
         ICCE = LF%AMG%NCC
         LF%NCW = 0

         DO IW = 1, LF%NW
            NOM = LF%WALL(IW)%NOM
            IF (NOM /= 0) THEN

               OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

               !WRITE(*,*) 'ACHTUNG: HIER WALL%ICN aendern auf OFI!'

               IC = LF%WALL(IW)%ICW
               !ICE = LF%WALL(IW)%ICE(1)

               IF (LF%AMG%CTYPE(IC)>0) LF%NCW = LF%NCW + 1

               !IF (CF%STYLE(ICE) >= NSCARC_CELL_TYPE_COARSE) THEN
               IF (LF%AMG%CTYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
                  OLF%AMG%NPS  = OLF%AMG%NPS  + 1
                  OLF%AMG%ICG0 = OLF%AMG%ICG0  + 1
                  OLF%AMG%NRS  = OLF%AMG%NPS
                  OLF%AMG%NCCS = OLF%AMG%NCCS + 1
               ELSE
                  OLF%AMG%NPS  = OLF%AMG%NPS  + AF%ROW(IC+1)-AF%ROW(IC) - 1
                  OLF%AMG%NRS  = OLF%AMG%NPS
                  OLF%AMG%NCFS = OLF%AMG%NCFS + 1
               ENDIF

               OLF%AMG%NCC=OLF%AMG%NCCS
               OLF%AMG%NCF=OLF%AMG%NCFS

            ENDIF
         ENDDO

         !>!
         !> Determine new numbering for coarse cells in interior of mesh
         !>!
         ICP   = 0
         DO IC = 1, LF%NCE
            IF (LF%AMG%CTYPE(IC) == NSCARC_CELL_TYPE_COARSE) THEN
               ICP = ICP + 1
               LF%AMG%CTYPE(IC) = ICP
            ENDIF
         ENDDO

      ENDDO TRANSFER_MESHES_LOOP

      !> Exchange sizes for AMG matrices and vectors
      IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_TRANSFER_SIZE, NL)

      !> Determine extended sizes for extended prolongation and restriction matrix
      TRANSFER_MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         LF => SCARC(NM)%LEVEL(NL)
         LC => SCARC(NM)%LEVEL(NL+1)
         MC => LC%MAP

         LF%AMG%NPE  = LF%AMG%NP
         LF%AMG%NRE  = LF%AMG%NR
         LF%AMG%NCCE = LF%AMG%NCC
         LF%AMG%NCFE = LF%AMG%NCF

         TRANSFER_OTHER_MESHES_LOOP2: DO NOM = 1, NMESHES

            !IF (NOM == NM) CYCLE TRANSFER_OTHER_MESHES_LOOP2

            OS => SCARC(NM)%OSCARC(NOM)
            OX => OS%EXCHANGE
            IF (OX%NICMAX_S==0 .AND. OX%NICMAX_R==0) CYCLE TRANSFER_OTHER_MESHES_LOOP2

            OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

            OLC%NC  = OLF%AMG%NCCS
            OLC%NCG = OLF%AMG%ICG0

            OLF%AMG%NCC = OLF%AMG%NCCS

            WRITE(*,*) 'ACHTUNG: STIMMT DAS MIT OLF%AMG%NCC?'

            CALL SCARC_ALLOCATE_INT1(LC%MAP%ICG_TO_ICE, 1, OLC%NC   , NSCARC_INIT_ZERO, 'ICG_TO_ICG')
!            CALL SCARC_ALLOCATE_INT1(LC%MAP%ICN_TO_ICE, 1, OLF%AMG%NCCI, NSCARC_INIT_ZERO, 'ICN_TO_ICE')

            DO IG = 1, OLC%NC
               LC%MAP%ICG_TO_ICE(IG) = LF%AMG%NCCE + IG
            ENDDO

            LF%AMG%NPE  = LF%AMG%NPE  + OLF%AMG%NP
            LF%AMG%NRE  = LF%AMG%NRE  + OLF%AMG%NR
            LF%AMG%NCCE = LF%AMG%NCCE + OLF%AMG%NCC
            LF%AMG%NCFE = LF%AMG%NCFE + OLF%AMG%NCF

         ENDDO TRANSFER_OTHER_MESHES_LOOP2
      ENDDO TRANSFER_MESHES_LOOP2

      !> Determine extended sizes for extended prolongation and restriction matrix
      TRANSFER_MESHES_LOOP3: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         LF => SCARC(NM)%LEVEL(NL)
         LC => SCARC(NM)%LEVEL(NL+1)

         CF => SCARC(NM)%LEVEL(NL)%CELL
         CC => SCARC(NM)%LEVEL(NL+1)%CELL

         LC%NC  = LF%AMG%NCC
         LC%NCE = LF%AMG%NCCE
         LC%NW  = LF%NCW

         ALLOCATE (LC%WALL(1:LC%NW), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'LC%WALL', IERROR)

         CALL SCARC_ALLOCATE_INT1 (LC%AMG%CTYPE    , 1, LC%NCE, NSCARC_INIT_UNDEFINED, 'CELL_TYPE')
         CALL SCARC_ALLOCATE_REAL1(LC%AMG%MEASURE  , 1, LC%NCE, NSCARC_INIT_UNDEFINED, 'MEASURE')
         CALL SCARC_ALLOCATE_INT1 (LC%AMG%INTERN, 1, LC%NC , NSCARC_INIT_UNDEFINED, 'BDRY')

         CALL SCARC_ALLOCATE_INT2 (CC%WINDEX    , 1, LC%NC, -3, 3, NSCARC_INIT_ZERO, 'WALL')

         IF (LC%NCE > LC%NC) THEN
            CALL SCARC_ALLOCATE_INT1(LC%MAP%ICE_TO_IWG, LC%NC+1, LC%NCE, NSCARC_INIT_NONE, 'ICE_TO_IWG')
            CALL SCARC_ALLOCATE_INT1(LC%AMG%EXTERN    , LC%NC+1, LC%NCE, NSCARC_INIT_NONE, 'EXTERNAL')
         ENDIF

      ENDDO TRANSFER_MESHES_LOOP3

END SELECT

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_SIZES



!> ------------------------------------------------------------------------------------------------
!> Initialize global 3D-solver methods (cg/mg)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AMG
INTEGER :: NL, NM, NOM
TYPE (SCARC_LEVEL_TYPE)   , POINTER :: L, LF, LC, OL
TYPE (SCARC_MATRIX_TYPE)  , POINTER :: A, S, ST, P, R, OP, OR
TYPE (SCARC_EXCHANGE_TYPE), POINTER :: OX
TYPE (SCARC_COORD_TYPE)   , POINTER :: CC
TYPE (SCARC_AMG_TYPE)     , POINTER :: AMG
TYPE (OSCARC_TYPE)        , POINTER :: OS

IF (.NOT.BAMG) RETURN
CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_AMG')

!> Determine number of multigrid levels
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX-1
!LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MIN

   !> ---------------------------------------------------------------------------------------------
   !> Determine coarser meshes corresponding to requested coarsening strategy
   !>  --- allocate necessary arrays
   !>  --- setup measures of single cells
   !>  --- setup CELL_TYPEs of single cells
   !>  --- setup sizes of transformation matrices (prolongation/restriction)
   !> ---------------------------------------------------------------------------------------------
   !IF (NL == NLEVEL_MIN) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      L   => SCARC(NM)%LEVEL(NL)
      AMG => SCARC(NM)%LEVEL(NL)%AMG

      A   => SCARC(NM)%SYSTEM(NL)%A
      S   => SCARC(NM)%SYSTEM(NL)%S
      ST  => SCARC(NM)%SYSTEM(NL)%ST

      CALL SCARC_ALLOCATE_REAL1(AMG%MEASURE, 1, L%NCE, NSCARC_INIT_UNDEFINED, 'MEASURE')
      CALL SCARC_ALLOCATE_INT1 (AMG%CTYPE  , 1, L%NCE, NSCARC_INIT_UNDEFINED, 'CTYPE')
      CALL SCARC_ALLOCATE_INT1 (AMG%GRAPH  , 1, L%NCE, NSCARC_INIT_UNDEFINED, 'GRAPH')
      CALL SCARC_ALLOCATE_INT1 (AMG%MARKER , 1, L%NCE, NSCARC_INIT_UNDEFINED, 'MARKER')

      CALL SCARC_COPY_MATRIX(A, S, 'S')
      CALL SCARC_ALLOCATE_MATRIX(ST, A%NAC+50, A%NAR+50, A%NAC+50, A%NSTENCIL, NSCARC_INIT_ZERO, 'ST')

      OTHER_MESH_INDEX1: DO NOM = 1, NMESHES

         !IF (NOM == NM) CYCLE OTHER_MESH_INDEX1

         OS  => SCARC(NM)%OSCARC(NOM)
         OX  => OS%EXCHANGE
         IF (OX%NICMAX_S==0 .AND. OX%NICMAX_R==0) CYCLE OTHER_MESH_INDEX1

         OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

         CALL SCARC_ALLOCATE_REAL1(OL%AMG%MEASURE, 1, OL%NCG, NSCARC_INIT_UNDEFINED, 'OMEASURE')
         CALL SCARC_ALLOCATE_INT1 (OL%AMG%CTYPE  , 1, OL%NCG, NSCARC_INIT_UNDEFINED, 'OCELL_TYPE')
         CALL SCARC_ALLOCATE_INT1 (OL%AMG%GRAPH  , 1, OL%NCG, NSCARC_INIT_UNDEFINED, 'OCELL_TYPE')

      ENDDO OTHER_MESH_INDEX1

   ENDDO

   !ENDIF

   CALL SCARC_SETUP_STRENGTH_MATRIX (NSCARC_COARSENING_RS3, NL)

   !> Then set measures and CELL_TYPEs on internal cells due to chosen coarsening strategy
   SELECT CASE (TYPE_COARSENING)
      CASE (NSCARC_COARSENING_RS3)
         CALL SCARC_SETUP_COLORING (NSCARC_COARSENING_RS3, NL)
      CASE (NSCARC_COARSENING_FALGOUT)
         CALL SCARC_SETUP_COLORING (NSCARC_COARSENING_RS3, NL)
         CALL SCARC_SETUP_COLORING (NSCARC_COARSENING_FALGOUT, NL)
   END SELECT

   !> Set dimensions for coarser mesh and define sizes of prolongation and restriction matrices
   CALL SCARC_SETUP_SIZES (NSCARC_SIZE_TRANSFER, NL)

   !CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE,  NL)

   !> --------------------------------------------------------------------------------------------
   !> Allocate and define grid transfer matrices
   !>  NCE  : number of extended cells in fine grid
   !>  NCCE : number of extended cells in coarse grid
   !>  NCW  : number of wall cells in coarse grid
   !> --------------------------------------------------------------------------------------------
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      L => SCARC(NM)%LEVEL(NL)
      A => SCARC(NM)%SYSTEM(NL)%A
      P => SCARC(NM)%SYSTEM(NL)%P
      R => SCARC(NM)%SYSTEM(NL)%R

      !> allocate prolongation and restriction matrix 
      CALL SCARC_ALLOCATE_MATRIX(P, L%AMG%NPE, L%AMG%NPE, L%NCE     , -1, NSCARC_INIT_ZERO, 'P')
      CALL SCARC_ALLOCATE_MATRIX(R, L%AMG%NRE, L%AMG%NRE, L%AMG%NCCE, -1, NSCARC_INIT_ZERO, 'P')

      !> allocate auxiliary tag arrays to mark relevant positions in A and P
      CALL SCARC_ALLOCATE_INT1(A%TAG, 1, L%NCE , NSCARC_INIT_ZERO, 'ATAG')
      CALL SCARC_ALLOCATE_INT1(P%TAG, 1, L%AMG%NCCE, NSCARC_INIT_ZERO, 'PTAG')

      OTHER_MESH_INDEX2: DO NOM = 1, NMESHES

         !IF (NOM == NM) CYCLE OTHER_MESH_INDEX2
         OS => SCARC(NM)%OSCARC(NOM)
         OX => OS%EXCHANGE

         IF (OX%NICMAX_S==0 .AND. OX%NICMAX_R==0) CYCLE OTHER_MESH_INDEX2
         OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
         OP => SCARC(NM)%OSCARC(NOM)%SYSTEM(NL)%P
         OR => SCARC(NM)%OSCARC(NOM)%SYSTEM(NL)%R

         CALL SCARC_ALLOCATE_MATRIX(OP, OL%AMG%NP+10, OL%AMG%NP+10, OL%NC+10, -1, NSCARC_INIT_ZERO, 'OP')
         CALL SCARC_ALLOCATE_MATRIX(OR, OL%AMG%NP+10, OL%AMG%NP+10, OL%NC+10, -1, NSCARC_INIT_ZERO, 'OR')   !

      ENDDO OTHER_MESH_INDEX2
   ENDDO

   !>
   !> determine prolongation and restriction matrix
   !> set corresponding overlap information between neighboring meshes for coarser grid level
   !>
   CALL SCARC_SETUP_PROLONGATION(NL)
   CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_MATRIX_TRANSPOSE(P%VAL, P%ROW, P%COL, &
                               R%VAL, R%ROW, R%COL, L%NCE, L%AMG%NCCE )

   !> -----------------------------------------------------------------------------------------
   !> Allocate coarse grid matrix including pointer arrays
   !> Note: number of cells on coarse level corresponds to number of c-points on fine level
   !> Compute coarse grid matrix by multiplication with restriction and prolongation matrix:
   !>  A_coarse := R * A_fine * P
   !> -----------------------------------------------------------------------------------------
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LF  => SCARC(NM)%LEVEL(NL)               !> Pointer to fine level
      LC  => SCARC(NM)%LEVEL(NL+1)             !> Pointer to coarse level

      CC => SCARC(NM)%LEVEL(NL+1)%COORD

      LC%NC  = LF%AMG%NCC
      LC%NCE = LF%AMG%NCCE
      LF%NW  = LF%NCW               ! ACHTUNG: HIER C%NW ???

      !> Allocate WALL-information for coarser grids (own and neighboring grids)
      CALL SCARC_ALLOCATE_INT1(LC%MAP%ICG_TO_ICE, LC%NC+1, LC%NCE, NSCARC_INIT_ZERO, 'ICG_TO_ICE')

   ENDDO
   IF (NMESHES>1) CALL SCARC_SETUP_WALLS_AMG (NL)          !ACHTUNG: HIER NOCHMAL CHECKEN ----

   CALL SCARC_SETUP_SUBDIVISION_AMG(NL)                       !HIER NOCHMAL CHECKEN ----

   CALL SCARC_TRANSFER_MATRIX (NL)
   CALL SCARC_SETUP_SYSTEM_AMG (NL)

   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SUBDIAG, NL+1)

!>! ACHTUNG : Deallocate auxiliary arrays ----!>!
   IF (LC%NC <= 4) THEN
      NLEVEL_MAX = NL + 1
      EXIT LEVEL_LOOP
   ENDIF

ENDDO LEVEL_LOOP

!DO NL = NLEVEL_MIN, NLEVEL_MAX
!>  DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!>     LF => SCARC(NM)%LEVEL(NL)
!>     DEALLOCATE(LF%AMG%MEASURE)
!>     DEALLOCATE(CF%STYLE)
!>  ENDDO
!ENDDO
CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_AMG



!> ----------------------------------------------------------------------------------------------------
!> Store subdivision information
!> temporarily not used
!> ----------------------------------------------------------------------------------------------------
!SUBROUTINE SCARC_SETUP_COORDINATES_AMG(NM, NL)
!INTEGER, INTENT(IN) :: NM, NL
!INTEGER IXC, IYC, IZC
!TYPE (SCARC_LEVEL_TYPE), POINTER :: LC , LF
!REAL(EB), DIMENSION(:) , POINTER :: CXCOR, CYCOR, CZCOR
!REAL(EB), DIMENSION(:) , POINTER :: FXCOR, FYCOR, FZCOR
!
!LF => SCARC(NM)%LEVEL(NL)
!LC => SCARC(NM)%LEVEL(NL+1)
!
!IF (NL == NLEVEL_MIN) THEN
!   FXCOR => MESHES(NM)%X
!   FYCOR => MESHES(NM)%Y
!   FZCOR => MESHES(NM)%Z
!ELSE
!   FXCOR => LF%COORD%XCOR
!   FYCOR => LF%COORD%YCOR
!   FZCOR => LF%COORD%ZCOR
!ENDIF
!
!CXCOR => LC%COORD%XCOR
!CYCOR => LC%COORD%YCOR
!CZCOR => LC%COORD%ZCOR
!
!IF (TWO_D) THEN
!   DO IXC = 0, LC%NX-1
!      CXCOR(IXC) = FXCOR(2*IXC)
!   ENDDO
!   CXCOR(IXC)  = FXCOR(LF%NX)
!   DO IZC = 0, LC%NZ-1
!      CZCOR(IZC) = FZCOR(2*IZC)
!   ENDDO
!   CZCOR(IZC) = FZCOR(LF%NZ)
!ELSE
!   DO IZC = 1, LC%NZ
!      DO IYC = 1, LC%NY
!         DO IXC = 1, LC%NX
!            WRITE(*,*) 'ACHTUNG: SETUP_COORDINATES_AMG 3D: STILL MISSING'
!         ENDDO
!      ENDDO
!   ENDDO
!ENDIF
!
!END SUBROUTINE SCARC_SETUP_COORDINATES_AMG

!> ----------------------------------------------------------------------------------------------------
!> Store subdivision information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIVISION(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IOR_LAST, NOM, INBR
INTEGER :: NEIGHBORS(20,-3:3)
TYPE (SCARC_LEVEL_TYPE), POINTER :: LC

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_SUBDIVISION')

IOR_LAST    = 0
NEIGHBORS   = 0

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LC => SCARC(NM)%LEVEL(NL)
   LC%SUBDIVISION = 0

   WALL_CELLS_LOOP: DO IW = 1, LC%NW

      IOR0 = LC%WALL(IW)%IOR

      IF (IOR_LAST /= IOR0) LC%SUBDIVISION(1,IOR0) = IW
      LC%SUBDIVISION(2,IOR0) = LC%SUBDIVISION(2,IOR0) + 1

      NOM= LC%WALL(IW)%NOM

      IF (NOM /= 0) THEN
         NEIGHBOR_LOOP: DO INBR = 1, 20
            IF (NOM == NEIGHBORS(INBR, IOR0)) THEN
               EXIT NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) /= 0) THEN
               CYCLE NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) == 0) THEN
               NEIGHBORS(INBR, IOR0) = NOM
               LC%SUBDIVISION(3,IOR0) = LC%SUBDIVISION(3,IOR0) + 1
               EXIT NEIGHBOR_LOOP
            ELSE
               CALL SCARC_SHUTDOWN('More than 20 neighbors at one face not allowed ', 'NONE', -999)
            ENDIF
         ENDDO NEIGHBOR_LOOP
      ENDIF

      IOR_LAST = IOR0

   ENDDO WALL_CELLS_LOOP
ENDDO MESHES_LOOP1

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_SUBDIVISION


!> ----------------------------------------------------------------------------------------------------
!> Store subdivision information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIVISION_AMG(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IOR_LAST, NOM, INBR
INTEGER :: NEIGHBORS(20,-3:3)
TYPE (SCARC_LEVEL_TYPE), POINTER :: LC

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_SUBDIVISION_AMG')

IOR_LAST    = 0
NEIGHBORS   = 0

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LC => SCARC(NM)%LEVEL(NL+1)
   LC%SUBDIVISION = 0

   WALL_CELLS_LOOP: DO IW = 1, LC%NW

      IOR0 = LC%WALL(IW)%IOR

      IF (IOR_LAST /= IOR0) LC%SUBDIVISION(1,IOR0) = IW
      LC%SUBDIVISION(2,IOR0) = LC%SUBDIVISION(2,IOR0) + 1

      NOM= LC%WALL(IW)%NOM

      IF (NOM /= 0) THEN
         NEIGHBOR_LOOP: DO INBR = 1, 20
            IF (NOM == NEIGHBORS(INBR, IOR0)) THEN
               EXIT NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) /= 0) THEN
               CYCLE NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) == 0) THEN
               NEIGHBORS(INBR, IOR0) = NOM
               LC%SUBDIVISION(3,IOR0) = LC%SUBDIVISION(3,IOR0) + 1
               EXIT NEIGHBOR_LOOP
            ELSE
               CALL SCARC_SHUTDOWN('More than 20 neighbors at one face not allowed ', 'NONE', -999)
            ENDIF
         ENDDO NEIGHBOR_LOOP
      ENDIF

      IOR_LAST = IOR0

   ENDDO WALL_CELLS_LOOP
ENDDO MESHES_LOOP1

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_SUBDIVISION_AMG



!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLS_AMG(NL)
INTEGER , INTENT(IN) :: NL
INTEGER :: NM, NOM, ICOL, IW, ICW, NCPL, JC
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE)   , POINTER :: LC, OLC
TYPE (SCARC_MATRIX_TYPE)  , POINTER :: AC
TYPE (SCARC_EXCHANGE_TYPE), POINTER :: OX

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_WALLS_AMG')

!> -------------------------------------------------------------------------
!> Loop over all boundary cells IW of fine grid
!> Get corresponding adjacent and ghost cell
!> -------------------------------------------------------------------------
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LC => SCARC(NM)%LEVEL(NL+1)
   AC => SCARC(NM)%SYSTEM(NL+1)%A

   DO IW = 1, LC%NW

      NOM = LC%WALL(IW)%NOM
      OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

      ICW = LC%WALL(IW)%ICW

      WRITE(*,*) 'L: HIER WALL%ICN aendernauf OFI!'

!>    LC%WALL(IW)%NCPL = NCPL
!>    CALL SCARC_ALLOCATE_INT1(LC%WALL(IW)%ICN, 1, NCPL, NSCARC_INIT_NONE, 'ICN')
!>    CALL SCARC_ALLOCATE_INT1(LC%WALL(IW)%ICE, 1, NCPL, NSCARC_INIT_NONE, 'ICE')
!>    CALL SCARC_ALLOCATE_INT1(LC%WALL(IW)%ICG, 1, NCPL, NSCARC_INIT_NONE, 'ICG')

      NCPL = 0
      DO ICOL = AC%ROW(ICW)+1, AC%ROW(ICW+1)-1
         JC = AC%COL(ICOL)
         IF (JC > LC%NC) THEN
            NCPL = NCPL + 1
!>           C%WALL(IW)%ICE(NCPL) = JC
!>           C%WALL(IW)%ICN(NCPL) = C%AMG%EXTERN(JC)
         ENDIF
      ENDDO

    ENDDO

ENDDO MESHES_LOOP

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MESHINFO, NL+1)

MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      !IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
      OS => SCARC(NM)%OSCARC(NOM)
      OX => OS%EXCHANGE

      IF (OX%NICMAX_S==0 .AND. OX%NICMAX_R==0) CYCLE OTHER_MESHES_LOOP
      OLC =>SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

      ALLOCATE(OLC%WALL(1:OLC%NW), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_WALLS_AMG','OLC%WALL',IERROR)

   ENDDO OTHER_MESHES_LOOP
ENDDO MESHES_LOOP2

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WALLINFO, NL+1)

CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_WALLS_AMG


!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION STRONGLY_COUPLED(A, IC, ICOL)
TYPE (SCARC_MATRIX_TYPE), POINTER, INTENT(INOUT) :: A
INTEGER , INTENT(IN) :: IC, ICOL
INTEGER :: JCOL
REAL(EB) :: AMG_TOL, VAL_MAX

AMG_TOL = 0.25_EB
VAL_MAX = 0.00_EB

DO  JCOL= A%ROW(IC)+1, A%ROW(IC+1)-1
   IF (JCOL /= ICOL) VAL_MAX = MAX(VAL_MAX, A%VAL(JCOL))
ENDDO

IF (A%VAL(ICOL) >= AMG_TOL * VAL_MAX) THEN
   STRONGLY_COUPLED=.TRUE.
ELSE
   STRONGLY_COUPLED=.FALSE.
ENDIF
RETURN

END FUNCTION STRONGLY_COUPLED

!> ------------------------------------------------------------------------------------------------
!> Determine measure of cells corresponding to requested coarsening type
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STRENGTH_MATRIX(NTYPE, NL)
INTEGER, PARAMETER  :: NCOUPLINGS = 20
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER  :: NM, IC, JC, ICOL, IDIAG, IS , IPOS
REAL(EB) :: ADIAG, ACOL, ROW_SCALE, ROW_SUM, MAX_ROW_SUM = 0.9_EB, THRESHOLD = 0.25_EB
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, S, ST

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_STRENGTH_MATRIX')

!> Only dummy (NTYPE really used ?)
IC = NTYPE

STRENGTH_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L  => SCARC(NM)%LEVEL(NL)
   A  => SCARC(NM)%SYSTEM(NL)%A
   S  => SCARC(NM)%SYSTEM(NL)%S
   ST => SCARC(NM)%SYSTEM(NL)%ST

   STRENGTH_CELL_LOOP: DO IC = 1, L%NC

      IDIAG = A%ROW(IC)
      ADIAG = A%VAL(IDIAG)

      !> get row sum and scaling factor
      ROW_SCALE = 0.0_EB
      ROW_SUM   = ADIAG

      DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
         JC = A%COL(ICOL)
         IF (JC <= L%NC) THEN
            ACOL = A%VAL(ICOL)
            ROW_SCALE = MAX(ROW_SCALE, ACOL)
            ROW_SUM   = ROW_SUM + ACOL
         ENDIF
      ENDDO

      !> get row entries of strength matrix S
      ROW_SUM = ABS(ROW_SUM/ADIAG)
      S%COL(IDIAG) = -1

      IF ((ROW_SUM > MAX_ROW_SUM) .AND. (MAX_ROW_SUM < 1.0_EB)) THEN
         !> set all dependencies to be weak
         DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
            JC = A%COL(ICOL)
            IF (JC <= L%NC) S%COL(ICOL) = -1
         ENDDO
      ELSE
         !> set dependencies to be weak related to threshold
         DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
            JC = A%COL(ICOL)
            IF (JC <= L%NC) THEN
               IF (A%VAL(ICOL) <= THRESHOLD * ROW_SCALE) S%COL(ICOL) = -1
            !ELSE
            !>   S%COL(ICOL) = -1
            ENDIF
          ENDDO
      ENDIF

   ENDDO STRENGTH_CELL_LOOP

   !> Compress strength matrix
   IS = 1
   STRENGTH_CELL_LOOP2: DO IC = 1, L%NC
      S%ROW(IC) = IS
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         !IF (S%COL(ICOL) > -1.AND.S%COL(ICOL)<=L%NC) THEN
         IF (S%COL(ICOL) > -1) THEN
            S%COL(IS) = S%COL(ICOL)
            IS = IS + 1
         ENDIF
      ENDDO
   ENDDO STRENGTH_CELL_LOOP2
   S%ROW(L%NC+1) = IS

   DO IC = 1, L%NCE+1
      ST%ROW(IC) = 0
   ENDDO

   IS = S%ROW(L%NC+1)-1
   DO ICOL = 1, IS
      ST%ROW(S%COL(ICOL)+1) = ST%ROW(S%COL(ICOL)+1) + 1
   ENDDO
   ST%ROW(1) = 1

   DO IC = 1, L%NCE
      ST%ROW(IC+1)= ST%ROW(IC+1) + ST%ROW(IC)
   ENDDO
   DO IC = 1, L%NC
      DO ICOL = S%ROW(IC), S%ROW(IC+1)-1
         IPOS = S%COL(ICOL)
         ST%COL(ST%ROW(IPOS)) = IC
         ST%ROW(IPOS) = ST%ROW(IPOS) + 1
      ENDDO
   ENDDO
   DO IC = L%NCE+1, 2, -1
      ST%ROW(IC) = ST%ROW(IC-1)
   ENDDO
   ST%ROW(1) = 1

ENDDO STRENGTH_MESHES_LOOP

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_STRENGTH_MATRIX


!> ------------------------------------------------------------------------------------------------
!> Determine measure of cells corresponding to requested coarsening type
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COLORING(NTYPE, NL)
INTEGER, PARAMETER  :: NCOUPLINGS = 20
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, IC, IC0, JC, KC, ICOL, JCOL, IRAND, REMAINING_CELLS, ICG, IW
INTEGER :: IGRAPH, IGRAPHE, IGRAPH_GLOBAL, ICYCLE
INTEGER :: FCELL = -1, ZCELL =-2, ICT2 = -1, ICT
LOGICAL :: BEMPTY=.FALSE., BDEBUG_NONEMPTY=.FALSE., BNEIGHBOR, BREAK = .TRUE., BCYCLE = .TRUE.
REAL(EB) :: RAND_NUM, MEASURE, NEW_MEASURE
REAL(EB) :: MEASURE_MAX, EPS
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_MATRIX_TYPE) , POINTER :: S, ST

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_COLORING')

EPS = 1.0E-12
MEASURE_MAX = 0.0_EB

!> Select coarsening strategy
SELECT CASE (NTYPE)

   !> ---------------------------------------------------------------------------------------------
   !> RS3-coarsening:
   !>      - Original Ruge-Stuben method with parallel postprocessing
   !>      - Produces good C/F splittings but is inherently serial.
   !>      - May produce AMG hierarchies with relatively high operator complexities.
   !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_RS3)

      RS3_MESH_INDEX: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L  => SCARC(NM)%LEVEL(NL)
         S  => SCARC(NM)%SYSTEM(NL)%S
         ST => SCARC(NM)%SYSTEM(NL)%ST

         REMAINING_CELLS = 0

         !> Currently the measures are computed as row sums of ST (number of influences for IC)
         RS3_MEASURE_LOOP0: DO IC = 1, L%NC
            L%AMG%MEASURE(IC) = ST%ROW(IC+1)-ST%ROW(IC)
         ENDDO RS3_MEASURE_LOOP0

         !> Subdivide in fine and coarse grid cells
         RS3_MEASURE_LOOP1: DO IC = 1, L%NC

            IF (S%ROW(IC+1)-S%ROW(IC) == 0) THEN
               L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_FINE0
               L%AMG%MEASURE(IC) = 0.0_EB
             !> IF (AGGRESSIVE2) L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE0
            ELSE
               L%AMG%CTYPE(IC) = NSCARC_UNDEFINED_INT
               REMAINING_CELLS   = REMAINING_CELLS + 1
            ENDIF
         ENDDO RS3_MEASURE_LOOP1

         RS3_MEASURE_LOOP2: DO IC = 1, L%NC
            MEASURE = L%AMG%MEASURE(IC)
            IF (L%AMG%CTYPE(IC) /= NSCARC_CELL_TYPE_FINE0 .AND. L%AMG%CTYPE(IC) /= NSCARC_CELL_TYPE_COARSE0) THEN
               IF (L%AMG%MEASURE(IC) > 0.0_EB) THEN
                  WRITE(*,*) 'ACHTUNG: KORREKT?'
                  !L%AMG%MEASURE(IC) = L%AMG%MEASURE(IC)
               ELSE
                  IF (L%AMG%MEASURE(IC) < 0.0_EB) WRITE(*,*) 'SCARC_SETUP_MEASURE: Negative measure !>'
                  L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_FINE
                  DO ICOL = S%ROW(IC), S%ROW(IC+1)-1
                     JC = S%COL(ICOL)
                     IF (L%AMG%CTYPE(JC) /= NSCARC_CELL_TYPE_COARSE0 .AND. L%AMG%CTYPE(JC) /= NSCARC_CELL_TYPE_FINE0) THEN
                        IF (JC < IC) THEN
                           NEW_MEASURE = L%AMG%MEASURE(JC)
                           IF (NEW_MEASURE > 0.0_EB) L%AMG%MEASURE(JC) = 0.0_EB
                           NEW_MEASURE = L%AMG%MEASURE(JC)+1
                           L%AMG%MEASURE(JC) = NEW_MEASURE
                        ELSE
                           NEW_MEASURE = L%AMG%MEASURE(JC)+1
                        ENDIF
                     ENDIF
                  ENDDO
                  REMAINING_CELLS = REMAINING_CELLS - 1
               ENDIF
            ENDIF
         ENDDO RS3_MEASURE_LOOP2

         RS3_CYCLE_LOOP: DO WHILE (REMAINING_CELLS > 0)

            !> get maximum (remaining) measure for all cells
            MEASURE_MAX = MAXVAL(L%AMG%MEASURE(1:L%NC))
            IF (MEASURE_MAX <= EPS) EXIT RS3_CYCLE_LOOP

            RS3_CELL_LOOP: DO IC = 1, L%NC

             !> Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, L%AMG%MEASURE(IC))) THEN

                  L%AMG%MEASURE(IC)  = NSCARC_MEASURE_NONE
                  L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE
                  REMAINING_CELLS = REMAINING_CELLS - 1

                  !> Determine set of fine cells
                  DO ICOL = ST%ROW(IC), ST%ROW(IC+1)-1

                   !> IF JC hasn't been marked yet, set it to be a fine cell which is no longer measured
                     JC = ST%COL(ICOL)

                     IF (JC > L%NC) CYCLE
                     IF (L%AMG%CTYPE(JC) == NSCARC_UNDEFINED_INT) THEN

                        L%AMG%MEASURE(JC)  = NSCARC_MEASURE_NONE
                        L%AMG%CTYPE(JC) = NSCARC_CELL_TYPE_FINE
                        REMAINING_CELLS = REMAINING_CELLS - 1

                        !>  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                        DO JCOL = S%ROW(JC), S%ROW(JC+1)-1
                           KC = S%COL(JCOL)
                           IF (L%AMG%CTYPE(KC)==NSCARC_UNDEFINED_INT) THEN
                              L%AMG%MEASURE(KC) = L%AMG%MEASURE(KC) + 1.0_EB
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
                  EXIT RS3_CELL_LOOP
               ENDIF
            ENDDO RS3_CELL_LOOP

            DO ICOL = S%ROW(IC), S%ROW(IC+1)-1
               JC = S%COL(ICOL)
               IF (JC > L%NC) CYCLE
               IF (L%AMG%CTYPE(JC) == NSCARC_UNDEFINED_INT) THEN
                  MEASURE = L%AMG%MEASURE(JC) - 1
                  L%AMG%MEASURE(JC) = MEASURE
                  IF (MEASURE > 0.0_EB) THEN
                     L%AMG%MEASURE(JC) = MEASURE
                  ELSE
                     L%AMG%CTYPE(JC) = NSCARC_CELL_TYPE_FINE
                     REMAINING_CELLS = REMAINING_CELLS - 1
                     DO JCOL = S%ROW(JC), S%ROW(JC+1)-1
                        KC = S%COL(JCOL)
                        IF (L%AMG%CTYPE(KC)==NSCARC_UNDEFINED_INT) THEN
                           L%AMG%MEASURE(KC) = L%AMG%MEASURE(KC) + 1.0_EB
                           MEASURE_MAX = MAX(MEASURE_MAX, L%AMG%MEASURE(KC))
                        ENDIF
                     ENDDO
                  ENDIF

               ENDIF
            ENDDO

         ENDDO RS3_CYCLE_LOOP
         L%NCW = 0

         DO IC = 1, L%NC
            IF (L%AMG%CTYPE(IC) == NSCARC_CELL_TYPE_COARSE0) L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE
         ENDDO

         !> exchange information with neighboring meshes
         IF (NMESHES > 1) THEN

            DO IC = 1, L%NCE
               L%AMG%GRAPH(IC) = -1
            ENDDO

            IC0 = 1
            DO IC = 1, L%NC
               IF (ICT2 /= IC) ICT = -1
               IF (L%AMG%CTYPE(IC) == NSCARC_CELL_TYPE_FINE) THEN

                  DO ICOL = S%ROW(IC), S%ROW(IC+1)-1
                     JC = S%COL(ICOL)
                     IF (JC <= L%NC) THEN
                        IF (L%AMG%CTYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                           L%AMG%GRAPH(JC) = IC
                        ENDIF
                     ENDIF
                  ENDDO

                  !> Hier fehlt noch die Abfrage nach Nachbarmesh !>!

                  DO ICOL = S%ROW(IC), S%ROW(IC+1)-1
                     JC = S%COL(ICOL)
                     IF (JC <= L%NC) THEN
                     IF (L%AMG%CTYPE(JC) == NSCARC_CELL_TYPE_FINE) THEN

                        !> DIESER PART WIRD ERST FUER ANDERE VERFEINERUNGEN AKTIV
                        !> ACHTUNG: DANN NOCHMAL UEBERPRUEFEN!
                        BEMPTY = .TRUE.
                        DO JCOL = S%ROW(JC), S%ROW(JC+1)-1
                           KC = S%COL(JCOL)
                           IF (L%AMG%GRAPH(KC) == IC) THEN
                              BEMPTY = .FALSE.
                              EXIT
                           ENDIF
                           IF (BEMPTY) THEN
                              IF (BDEBUG_NONEMPTY) THEN
                                 L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE
                                 IF (ICT > -1) THEN
                                    L%AMG%CTYPE(ICT) = NSCARC_CELL_TYPE_FINE
                                    ICT = -1
                                 ENDIF
                                 !> Hier fehlt noch Nachbaranteil
                                 BDEBUG_NONEMPTY = .FALSE.
                                 EXIT
                              ELSE
                                 ICT  = JC
                                 ICT2 = IC
                                 L%AMG%CTYPE(JC) = NSCARC_CELL_TYPE_COARSE
                                 BDEBUG_NONEMPTY = .FALSE.
                                 EXIT
                                 IC0 = IC0 - 1
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDIF
                     ENDIF
                  ENDDO
               ENDIF
               IC0 = IC0 + 1 !> Achtung, nochmal pruefen, eventuell IC0 als Index verwenden?
            ENDDO
         ENDIF
      ENDDO RS3_MESH_INDEX

      !> Exchange CELL_TYPE-data along internal boundaries
      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)

      !> Third pass: identify fine cells along boundary and get their coarse neighbors
      RS3_MESH_INDEX2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         DO IC = 1, L%NC
            L%AMG%GRAPH(IC) = -1
         ENDDO
      ENDDO RS3_MESH_INDEX2


!> ----------------------------------------------------------------------------------------
!> Falgout coarsening
!> ----------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_FALGOUT)

      FCELL = ZCELL
      FALGOUT_MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         S => SCARC(NM)%SYSTEM(NL)%S

         !> Reset the measures as row sums of ST (number of influences for IC)
         !> plus random number to make them unique
         FALGOUT_CELL_LOOP1: DO IC = 1, L%NC
            RAND_NUM = 0.01_EB
            DO IRAND = 1, 5
               CALL RANDOM_NUMBER(RAND_NUM)
               RAND_NUM = RAND_NUM + RAND_NUM/10**(IRAND-1)
               L%AMG%MEASURE(IC) = L%AMG%MEASURE(IC) + RAND_NUM
            ENDDO
            L%AMG%MEASURE(IC) = S%ROW(IC+1)-S%ROW(IC) + RAND_NUM
         ENDDO FALGOUT_CELL_LOOP1

      ENDDO FALGOUT_MESHES_LOOP1

      !> Initial exchange of measure array
      IF (NMESHES > 1) THEN
         TYPE_VECTOR = NSCARC_VECTOR_MEASURE
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR , NL)
      ENDIF

      FALGOUT_INIT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)

         !> reset CELL_TYPE for cells with neighbors in other meshes
         FALGOUT_INTERNAL_CELL_LOOP1: DO IC = 1, L%NC
            BNEIGHBOR = .FALSE.
            DO ICOL = S%ROW(IC), S%ROW(IC+1)-1
               JC = S%COL(ICOL)
               IF (JC > L%NC) THEN
                  BNEIGHBOR = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (BNEIGHBOR .OR. L%AMG%CTYPE(IC) == NSCARC_CELL_TYPE_FINE) L%AMG%CTYPE(IC)=NSCARC_UNDEFINED_INT
         ENDDO FALGOUT_INTERNAL_CELL_LOOP1

         !> initialize GRAPH and reset CELL_TYPE on ghost cells
         FALGOUT_EXTENDED_CELL_LOOP: DO IC = L%NC+1, L%NCE
            L%AMG%GRAPH(IC) = IC
            L%AMG%CTYPE(IC) = NSCARC_UNDEFINED_INT
         ENDDO FALGOUT_EXTENDED_CELL_LOOP

         !> reset CELL_TYPE on internal wall cells
         DO IW = 1, L%NW
            IF (L%WALL(IW)%NOM /= 0) THEN
               IC = L%WALL(IW)%ICW
               L%AMG%CTYPE(IC) = NSCARC_UNDEFINED_INT
            ENDIF
         ENDDO

         !> reset CELL_TYPE on internal fine cells
         ICG = 1
         FALGOUT_INTERNAL_CELL_LOOP2: DO IC = 1, L%NC
            IF (L%AMG%CTYPE(IC)<NSCARC_UNDEFINED_INT) THEN
               L%AMG%CTYPE(IC)=NSCARC_UNDEFINED_INT
            ENDIF
            IF (L%AMG%CTYPE(IC) == NSCARC_CELL_TYPE_ZPNT) THEN
               IF (L%AMG%MEASURE(IC) >= 1.0_EB .OR. (S%ROW(IC+1)-S%ROW(IC)) > 0) THEN
                  L%AMG%CTYPE(IC) = NSCARC_UNDEFINED_INT
                  L%AMG%GRAPH(ICG) = IC
                  ICG = ICG + 1
               ELSE
                  L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_FINE
               ENDIF
            ELSE IF (L%AMG%CTYPE(IC) == NSCARC_CELL_TYPE_SFPNT) THEN
               L%AMG%MEASURE(IC) = NSCARC_MEASURE_NONE
            ELSE
               L%AMG%GRAPH(ICG) = IC
               ICG = ICG + 1
            ENDIF
         ENDDO FALGOUT_INTERNAL_CELL_LOOP2

         IGRAPH  = ICG-1
         IGRAPHE = L%NCE-L%NC

      ENDDO FALGOUT_INIT_LOOP

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         FALGOUT_EXTERNAL_LOOP: DO IC = L%NC+1, L%NCE
            L%AMG%MEASURE(IC) = NSCARC_MEASURE_NONE
         ENDDO FALGOUT_EXTERNAL_LOOP
      ENDDO

      !> Coloring loop until all cells are coarse or fine
      FALGOUT_COLORING_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         ICYCLE = 0
         FALGOUT_GRAPH_LOOP: DO WHILE (BCYCLE)

            IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MEASURE_ADD, NL)
            !IF (NMESHES > 1) THEN
            !   TYPE_VECTOR = NSCARC_VECTOR_MEASURE
            !   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)
            !ENDIF

            IF (ICYCLE > 0) THEN

               ICG = 1
               FALGOUT_GRAPH_LOOP1: DO WHILE (ICG <= IGRAPH)

                  IC = L%AMG%GRAPH(ICG)

                  !> if cell isn't marked as coarse cell yet and has measure
                  !> less than 1, mark it as fine cell
                  !> take care that all dependencies have been taken into account
                  IF (L%AMG%CTYPE(IC) /= NSCARC_CELL_TYPE_COARSE .AND. &
                      L%AMG%MEASURE(IC)  <  NSCARC_MEASURE_ONE) THEN
                     L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_FINE
                     DO ICOL = S%ROW(IC), S%ROW(IC+1)-1
                        JC = S%COL(ICOL)
                        IF (JC < 0) CYCLE
                        L%AMG%CTYPE(IC) = NSCARC_UNDEFINED_INT
                     ENDDO
                  ENDIF

                  !> if cell is already marked as fine or coarse, set its measure to zero
                  !> and extract it from the graph (put it at the end of the graph
                  !> array and decrease number of relevant graph entries)
                  IF (L%AMG%CTYPE(IC) /= NSCARC_UNDEFINED_INT) THEN
                     L%AMG%MEASURE(IC) = NSCARC_MEASURE_NONE
                     L%AMG%GRAPH(ICG) = L%AMG%GRAPH(IGRAPH)
                     L%AMG%GRAPH(IGRAPH) = IC
                     IGRAPH = IGRAPH - 1
                     ICG = ICG - 1
                  ENDIF

                  ICG = ICG + 1
               ENDDO FALGOUT_GRAPH_LOOP1

            ENDIF

            !> exchanges measure and CELL_TYPEs of neighboring cells
            IF (NMESHES > 1) THEN
               TYPE_VECTOR = NSCARC_VECTOR_MEASURE
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)
            ENDIF

            !> Get global number of relevant graph entries (was ist im seriellen Fall??)
            IF (NMESHES>1 .AND. N_MPI_PROCESSES>1) THEN
               WRITE(*,*) 'ACHTUNG, HIER CHECKEN WEGEN N_PROCESSES'
               CALL MPI_ALLREDUCE (IGRAPH, IGRAPH_GLOBAL, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERROR)
            ENDIF

            !> If all cells have been encountered for, leave loop
            IF (IGRAPH_GLOBAL == 0) EXIT FALGOUT_GRAPH_LOOP

            !> search for an independent set of points with maximal measure
            IF (ICYCLE > 0) CALL SCARC_SETUP_GRAPHSET(IGRAPH, NM, NL)
            ICYCLE = ICYCLE + 1

            !> exchanges CELL_TYPEs of neighboring cells (with new information from graphset)
            IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)

            ICG = L%NC+1

            FALGOUT_GRAPH_LOOP2: DO WHILE (ICG <= L%NC+IGRAPHE)
               IC = L%AMG%GRAPH(ICG)

               IF (IC < 0) CYCLE
               IF (L%AMG%CTYPE(IC) < NSCARC_UNDEFINED_INT) THEN
                  L%AMG%GRAPH(ICG) = L%AMG%GRAPH(L%NC+IGRAPHE)
                  L%AMG%GRAPH(L%NC+IGRAPHE) = IC
                  IGRAPHE = IGRAPHE - 1
               ENDIF
               ICG = ICG + 1
            ENDDO FALGOUT_GRAPH_LOOP2

            L%AMG%MEASURE(L%NC+1:L%NCE) = NSCARC_MEASURE_NONE

            FALGOUT_GRAPH_LOOP3: DO ICG = 1, IGRAPH

               IC = L%AMG%GRAPH(ICG)

               !> Coarse cells don't interpolate from influencing neighbors
               IF (L%AMG%CTYPE(IC) > NSCARC_UNDEFINED_INT) THEN

                  L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE

                  DO ICOL = S%ROW(IC), S%ROW(IC+1)-1
                     JC = S%COL(ICOL)

                     !> remove edge from S and decrement measures of unmarked neighbors
                     IF (JC > -1) THEN
                        S%COL(ICOL) = - S%COL(ICOL) - 1
                        IF (L%AMG%CTYPE(JC) == NSCARC_UNDEFINED_INT) L%AMG%MEASURE(JC) = L%AMG%MEASURE(JC) - 1
                     ENDIF
                  ENDDO

               ELSE

                  !> dependencies which are already marked
                  DO ICOL = S%ROW(IC), S%ROW(IC+1)-1

                     JC = S%COL(ICOL)
                     IF (JC < 0) JC = -JC - 1

                     !> remove edge from S and temporarily reset CELL_TYPE
                     IF (L%AMG%CTYPE(JC) > NSCARC_UNDEFINED_INT) THEN
                        IF (S%COL(ICOL) > -1) S%COL(ICOL) = -S%COL(ICOL) - 1
                        L%AMG%CTYPE(JC) = NSCARC_CELL_TYPE_CPNT
                     ELSE IF (L%AMG%CTYPE(JC) == NSCARC_CELL_TYPE_SFPNT) THEN    !> necessary ??
                        IF (S%COL(ICOL) > -1) S%COL(ICOL) = -S%COL(ICOL) - 1
                     ENDIF
                  ENDDO

                  !> dependencies which aren't marked yet
                  DO ICOL = S%ROW(IC), S%ROW(IC+1)-1

                     JC = S%COL(ICOL)
                     IF (JC > -1 .AND. JC<=L%NC) THEN

                        BREAK = .TRUE.

                        !> check if there are common C-points
                        DO JCOL = S%ROW(JC), S%ROW(JC+1)-1
                           KC = S%COL(JCOL)
                           IF (KC < 0) KC = -KC - 1
                           IF (KC <= L%NC) THEN
                              !> remove edge from S and update measure
                              IF (L%AMG%CTYPE(KC) == NSCARC_CELL_TYPE_CPNT) THEN
                                 S%COL(ICOL) = - S%COL(ICOL) - 1
                                 L%AMG%MEASURE(JC) = L%AMG%MEASURE(JC) - 1
                                 BREAK = .FALSE.
                                 EXIT
                              ENDIF
                           ENDIF
                        ENDDO

                        IF (BREAK) THEN
                           DO JCOL = S%ROW(JC), S%ROW(JC+1)-1
                              KC = S%COL(JCOL)
                              IF (KC < 0) KC = -KC - 1
                              IF (KC > L%NC) THEN
                                 !> remove edge from S and update measure
                                 IF (L%AMG%CTYPE(KC) == NSCARC_CELL_TYPE_CPNT) THEN
                                    S%COL(ICOL) = - S%COL(ICOL) - 1
                                    L%AMG%MEASURE(JC) = L%AMG%MEASURE(JC) - 1
                                    BREAK = .FALSE.
                                    EXIT
                                 ENDIF
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF

             !> reset CELL_TYPEs
               DO ICOL = S%ROW(IC), S%ROW(IC+1)-1
                  JC = S%COL(ICOL)
                  IF (JC < 1) JC = -JC - 1
                  IF (L%AMG%CTYPE(JC) == NSCARC_CELL_TYPE_CPNT) L%AMG%CTYPE(JC) = NSCARC_CELL_TYPE_COARSE
               ENDDO

            ENDDO FALGOUT_GRAPH_LOOP3

         ENDDO FALGOUT_GRAPH_LOOP

       !> reset S-matrix
         DO ICOL = 1, S%ROW(L%NC+1)
            IF (S%COL(ICOL) < 0) S%COL(ICOL) = -S%COL(ICOL)-1
         ENDDO

      ENDDO FALGOUT_COLORING_LOOP

END SELECT

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_COLORING


!> ----------------------------------------------------------------------------------------
!> Select an independent set from the graph
!> ----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRAPHSET(IGRAPH, NM, NL)
INTEGER, INTENT(IN):: IGRAPH, NM, NL
INTEGER :: IG, IC, ICOL, JC
TYPE (SCARC_LEVEL_TYPE),  POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_GRAPHSET')

L => SCARC(NM)%LEVEL(NL)
A => SCARC(NM)%SYSTEM(NL)%A

!> First mark every cell from the 'internal' graphset with measure bigger
!> than one as a coarse cell
DO IG = 1, IGRAPH
   IC = L%AMG%GRAPH(IG)
   IF (L%AMG%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE
   ENDIF
ENDDO

!> Do the same with cells on the overlapping areas
DO IG = L%NC+1, L%NCE
   IC = L%AMG%GRAPH(IG)
   IF (IC < 0) CYCLE
   IF (L%AMG%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE
   ENDIF
ENDDO

!> remove nodes from the initial independent set depending on their measure
!> For each cell consider every connection in stencil and set cell with
!> 1G/smaller measure to Zero-Type
DO IG = 1, IGRAPH
   IC = L%AMG%GRAPH(IG)

   IF (IC < 0) CYCLE

   IF (L%AMG%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
         JC = A%COL(ICOL)
         IF (JC < 0) CYCLE
         IF (L%AMG%MEASURE(JC) > NSCARC_MEASURE_ONE) THEN
            IF (L%AMG%MEASURE(IC) > L%AMG%MEASURE(JC)) THEN
               L%AMG%CTYPE(JC) = NSCARC_UNDEFINED_INT
            ELSE IF (L%AMG%MEASURE(JC) > L%AMG%MEASURE(IC)) THEN
               L%AMG%CTYPE(IC) = NSCARC_UNDEFINED_INT
            ENDIF
         ENDIF
      ENDDO
   ENDIF
ENDDO

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_GRAPHSET

!> ------------------------------------------------------------------------------------------------
!> Numbering of single patches
!> temporarily not used
!> ------------------------------------------------------------------------------------------------
!SUBROUTINE SCARC_PATCH_CELL_TYPES(CELL_TYPE, NX, NZ, NX1, NX2, NY1, NY2, IZ, INCRX, INCRY, INCRZ)
!INTEGER, INTENT(IN):: NX, NZ, NX1, NX2, NY1, NY2, INCRX, INCRY, INCRZ, IZ
!INTEGER, DIMENSION(:), INTENT(OUT) :: CELL_TYPE
!INTEGER :: IX, INX, INY, INZ, PX(3), PY(3), PZ(3), PC(3,3), INX0, INZ0
!
!!> Only dummy (check, if variables are still needed!)
!IX = NY1
!IX = NY2
!IX = INCRY
!
!INX0 = 1
!INZ0 = 1
!INY  = 1
!PY   = 1
!DO IX = NX1, NX2, INCRX
!
!   DO INZ = 1, INCRZ
!
!      PZ(INZ) = IZ + INZ - 1
!      IF (NZ>3.AND.(PZ(INZ)==1 .OR. PZ(INZ)== NZ)) INZ0 = INZ
!
!      DO INX = 1, INCRX
!
!         PX(INX) = IX + INX - 1
!         IF (NX>3.AND.(PX(INX)==1 .OR. PX(INX)== NX)) INX0 = INX
!
!         PC(INX,INZ) = (PZ(INZ)-1)*NX + PX(INX)
!         CELL_TYPE(PC(INX,INZ)) = NSCARC_CELL_TYPE_SFINE
!      ENDDO
!   ENDDO
!   CELL_TYPE(PC(INX0,INZ0)) = NSCARC_CELL_TYPE_COARSE
!
!ENDDO
!
!END SUBROUTINE SCARC_PATCH_CELL_TYPES
!
!
!!> ------------------------------------------------------------------------------------------------
!!> Setup CELL_TYPEs of mesh corresponding to requested coarsening strategy
!!> temporarily not used
!!> ------------------------------------------------------------------------------------------------
!SUBROUTINE SCARC_SETUP_CELL_TYPES(NTYPE, NL)
!INTEGER, PARAMETER  :: NCOUPLINGS = 20, NCYC_MAX=1000
!INTEGER, INTENT(IN) :: NTYPE, NL
!INTEGER  :: NM, REMAINING_CELLS, MEASURE_TYPE
!INTEGER  :: IA, IC, JC, IW, KC, IG, ICOL, JCOL, ILOOP, IGRAPH, ICOUNT, IROW, JROW
!REAL(EB) :: MEASURE_MAX, EPS
!LOGICAL  :: BREMOVE=.FALSE.
!TYPE (SCARC_LEVEL_TYPE),  POINTER :: L
!TYPE (SCARC_MATRIX_TYPE), POINTER :: A, S
!
!CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_CELL_TYPES')
!
!EPS = 1.0E-12
!MEASURE_MAX = 0.0_EB
!
!
!!> Define CELL_TYPEs for corresponding coarsening strategy
!SELECT CASE (NTYPE)
!
!   !> ---------------------------------------------------------------------------------------------
!   !> standard Rue-Stueben coarsening
!   !>  - identify first point IC with maximum measure and set it to be a C-point
!   !>  - identify surrounding strongly connected points JC and set them to be F-points
!   !>  - increase measures of strongly connected neighbours KC of JC
!   !>  - decrease measures of strongly connected neighbours of IC
!   !> ---------------------------------------------------------------------------------------------
!
!   CASE (NSCARC_COARSENING_RS3)
!
!      MEASURE_TYPE = 1
!      RS3_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!
!         L => SCARC(NM)%LEVEL(NL)
!         A => SCARC(NM)%SYSTEM(NL)%A
!
!         REMAINING_CELLS = L%NC
!
!         RS3_CYCLE_LOOP: DO WHILE (REMAINING_CELLS > 0)
!
!            !> get maximum (remaining) measure for all cells
!            MEASURE_MAX = MAXVAL(L%AMG%MEASURE(1:L%NC))
!            IF (MEASURE_MAX <= EPS) EXIT RS3_CYCLE_LOOP
!
!            RS3_CELL_LOOP: DO IC = 1, L%NC
!
!               !> Take first cell with maximum measure as next coarse cell
!               IF (MATCH(MEASURE_MAX, L%AMG%MEASURE(IC))) THEN
!
!                  L%AMG%MEASURE(IC)  = NSCARC_MEASURE_NONE
!                  L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE
!                  REMAINING_CELLS = REMAINING_CELLS - 1
!
!                  !> Determine set of fine cells
!                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
!
!                     IF (STRONGLY_COUPLED(A, IC, ICOL)) THEN
!
!                        !> IF JC hasn't been marked yet, set it to be a fine cell which is no longer measured
!                        JC = A%COL(ICOL)
!                        IF (L%AMG%CTYPE(JC) == NSCARC_UNDEFINED_INT) THEN
!
!                           L%AMG%MEASURE(JC)  = NSCARC_MEASURE_NONE
!                           L%AMG%CTYPE(JC) = NSCARC_CELL_TYPE_FINE
!                           REMAINING_CELLS = REMAINING_CELLS - 1
!
!                           !>  increase measures of cells KC adjacent to fine cells JC based on strong couplings
!                           DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
!                              IF (STRONGLY_COUPLED(A, JC, JCOL)) THEN
!                                 KC = A%COL(JCOL)
!                                 IF (KC /= 0) THEN
!                                    IF (KC /= IC .AND. (L%AMG%CTYPE(KC)==NSCARC_UNDEFINED_INT)) THEN
!                                       L%AMG%MEASURE(KC) = L%AMG%MEASURE(KC) + 1.0_EB
!                                       MEASURE_MAX = MAX(MEASURE_MAX, L%AMG%MEASURE(KC))
!                                    ENDIF
!                                 ENDIF
!                              ENDIF
!                           ENDDO
!
!                        ENDIF
!                     ENDIF
!                  ENDDO
!                  EXIT RS3_CELL_LOOP
!               ENDIF
!            ENDDO RS3_CELL_LOOP
!         ENDDO RS3_CYCLE_LOOP
!         L%NCW = 0
!      ENDDO RS3_LOOP
!
!   !> ---------------------------------------------------------------------------------------------
!   !> set CELL_TYPEs for Falgout coarsening
!   !> Care: ZPOINT is equal to FINE in this case !>!
!   !> ---------------------------------------------------------------------------------------------
!   CASE (NSCARC_COARSENING_FALGOUT)
!
!      FALGOUT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!
!         L => SCARC(NM)%LEVEL(NL)
!         A => SCARC(NM)%SYSTEM(NL)%A
!         S => SCARC(NM)%SYSTEM(NL)%S
!
!         !> copy column pointers for matrix
!         DO IA = 1, A%NAC
!            S%COL(IA) = A%COL(IA)
!         ENDDO
!
!         !> this part may be changed for other coarsening strategies 
!         ICOUNT = 1
!         FALGOUT_CELL_LOOP: DO IC = 1, L%NC
!            IF (L%AMG%CTYPE(IC) == NSCARC_CELL_TYPE_FINE) THEN
!               IF (L%AMG%MEASURE(IC) >= 1.0_EB .AND. A%ROW(IC+1)-A%ROW(IC)+1 > 0) THEN
!                  L%AMG%CTYPE(IC)  = NSCARC_UNDEFINED_INT
!                  L%AMG%GRAPH(ICOUNT) = IC
!                  ICOUNT = ICOUNT + 1
!                  IF (ICOUNT > L%NC) &
!                     CALL SCARC_SHUTDOWN('Graph size must be increased ', 'NONE', ICOUNT)
!               ELSE
!                  L%AMG%CTYPE(IC)  = NSCARC_CELL_TYPE_FINE
!               ENDIF
!            ELSEIF (L%AMG%CTYPE(IC) == NSCARC_CELL_TYPE_FINE0) THEN
!               L%AMG%MEASURE(IC) = NSCARC_MEASURE_NONE
!            ELSE
!               L%AMG%GRAPH(ICOUNT) = IC
!               ICOUNT = ICOUNT + 1
!            ENDIF
!         ENDDO FALGOUT_CELL_LOOP
!
!         IGRAPH = ICOUNT
!         ILOOP  = 0
!         FALGOUT_COARSENING_LOOP: DO
!
!            !> Care: Is graph_size always L%NC ?
!            IF (ILOOP > 0) THEN
!
!               DO IG = 1, IGRAPH
!
!                  IC = L%AMG%GRAPH(IG)
!                  JC = IG
!
!                  !> make IC a fine cell and look for all dependencies
!                  IF ((L%AMG%CTYPE(IC) /= NSCARC_CELL_TYPE_COARSE) .AND. (L%AMG%MEASURE(IC) < 1)) THEN
!                     L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_FINE
!                     DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
!                        IF (S%COL(IROW) > -1) L%AMG%CTYPE(IC) = NSCARC_UNDEFINED_INT
!                     ENDDO
!                  ENDIF
!
!                  !> remove cell from graph
!                  IF (L%AMG%CTYPE(IC) /= NSCARC_UNDEFINED_INT) THEN
!                     L%AMG%MEASURE(JC) = NSCARC_MEASURE_NONE
!                     IGRAPH = IGRAPH - 1
!                     L%AMG%GRAPH(JC) = L%AMG%GRAPH(IGRAPH)
!                     L%AMG%GRAPH(IGRAPH) = IC
!                     JC = JC - 1               !> correct ???
!                  ENDIF
!
!               ENDDO
!
!            ENDIF
!            ILOOP = ILOOP + 1
!
!            FALGOUT_GRAPH_LOOP: DO IG = 1, IGRAPH
!
!               IC = L%AMG%GRAPH(IG)
!
!               !> C-points are not influenced by neighbors (no interpolation)
!               IF (L%AMG%CTYPE(IC) > 0) THEN
!
!                  L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE            !> define it as C-point
!
!                  DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
!                     JC = S%COL(IROW)
!                     IF (JC > -1) THEN
!                        S%COL(IROW) = -S%COL(IROW) - 1  !> remove edge
!                        IF (L%AMG%CTYPE(JC) == NSCARC_UNDEFINED_INT) THEN
!                           L%AMG%MEASURE(JC) = L%AMG%MEASURE(JC) - 1.0_EB   !> decrement measure of unmarked neigh.
!                        ENDIF
!                     ENDIF
!                  ENDDO
!
!               ELSE
!
!                  !> dependencies which have already been marked
!                  DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
!                     JC = S%COL(IROW)
!                     IF (JC < 0) JC = -JC-1
!                     IF (L%AMG%CTYPE(JC) > 0) THEN
!                        IF (S%COL(IROW) > -1) THEN
!                           S%COL(IROW) = -S%COL(IROW) -1   !> remove edge
!                           L%AMG%CTYPE(JC)    = NSCARC_CELL_TYPE_COMMON   !> temporarily set CELL_TYPE to common
!                        ENDIF
!                     ELSEIF (L%AMG%CTYPE(JC) == NSCARC_CELL_TYPE_FINE0) THEN
!                        IF (S%COL(IROW) > -1) THEN
!                           S%COL(IROW) = -S%COL(IROW) -1   !> remove edge
!                        ENDIF
!                     ENDIF
!                  ENDDO
!
!                  !> dependencies which haven't been marked yet
!                  DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
!                     IF (S%COL(IROW) > -1) THEN
!                        BREMOVE = .TRUE.
!                        JC = S%COL(IROW)
!                        DO JROW = A%ROW(JC)+1, A%ROW(JC+1)-1
!                           KC = S%COL(JROW)
!                           IF (KC < 0) KC = -KC-1                        !> check for all dependencies !>
!                           IF (L%AMG%CTYPE(KC) == NSCARC_CELL_TYPE_COMMON) THEN
!                              S%COL(IROW) = -S%COL(IROW)-1
!                              L%AMG%MEASURE(JC) = L%AMG%MEASURE(JC) - 1.0_EB
!                              BREMOVE = .FALSE.
!                              EXIT
!                           ENDIF
!                        ENDDO
!                     ENDIF
!                  ENDDO
!
!               ENDIF
!
!               !> reset CELL_TYPES
!               DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
!                  JC = S%COL(IROW)
!                  IF (JC < 0) JC = -JC-1
!                  IF (L%AMG%CTYPE(JC) == NSCARC_CELL_TYPE_COMMON) L%AMG%CTYPE(JC)=NSCARC_CELL_TYPE_COARSE
!               ENDDO
!
!            ENDDO FALGOUT_GRAPH_LOOP
!         ENDDO FALGOUT_COARSENING_LOOP
!
!      ENDDO FALGOUT_LOOP
!
!   !> ---------------------------------------------------------------------------------------------
!   !> first set CELL_TYPEs for cells in internal boundary layers (adjacent to neighboring meshes)
!   !> ---------------------------------------------------------------------------------------------
!   CASE (NSCARC_COARSENING_BDRY)
!
!      BDRY_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!
!         L => SCARC(NM)%LEVEL(NL)
!
!         !> First define coarse cells on boundary layer(adjacent to neighboring meshes)
!         BDRY_CYCLE_LOOP: DO
!
!            !>!
!            !> get maximum (remaining) measure for all cells
!            !>!
!            MEASURE_MAX = 0.0_EB
!            BDRY_LOOP1: DO IW = 1, L%NW
!               IF (L%WALL(IW)%NOM == 0) CYCLE BDRY_LOOP1
!               IC = L%WALL(IW)%ICW
!               MEASURE_MAX = MAX(MEASURE_MAX, L%AMG%MEASURE(IC))
!            ENDDO BDRY_LOOP1
!            IF (MEASURE_MAX <= EPS) EXIT BDRY_CYCLE_LOOP
!
!            BDRY_LOOP2: DO IW = 1, L%NW
!
!               IF (L%WALL(IW)%NOM == 0) CYCLE BDRY_LOOP2
!               IC = L%WALL(IW)%ICW
!
!               !> Take first cell with maximum measure as next coarse cell
!               IF (MATCH(MEASURE_MAX, L%AMG%MEASURE(IC))) THEN
!
!                  L%AMG%MEASURE(IC)  = 0.0_EB
!                  L%AMG%CTYPE(IC) = NSCARC_CELL_TYPE_COARSE
!
!                  !> Determine set of fine cells
!                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
!
!                   !> JC is set to be a fine cell which is no longer measured
!                     JC = A%COL(ICOL)
!
!                     L%AMG%MEASURE(JC)  = 0.0_EB
!                     L%AMG%CTYPE(JC) = NSCARC_CELL_TYPE_SFINE
!
!                     !>  increase measures of cells KC adjacent to fine cells JC based on strong couplings
!                     DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
!                        KC = A%COL(JCOL)
!                        IF (KC > 0 .AND. KC /= IC .AND. (L%AMG%CTYPE(KC)==NSCARC_UNDEFINED_INT)) THEN
!                           L%AMG%MEASURE(KC) = L%AMG%MEASURE(KC) + 1.0_EB
!                           MEASURE_MAX = MAX(MEASURE_MAX, L%AMG%MEASURE(KC))
!                        ENDIF
!                     ENDDO
!
!                  ENDDO
!
!                  EXIT BDRY_LOOP2
!               ENDIF
!
!            ENDDO BDRY_LOOP2
!         ENDDO BDRY_CYCLE_LOOP
!      ENDDO BDRY_LOOP
!
!END SELECT
!
!!> -----------------------------------------------------------------------------------------------------
!!> Exchange CELL_TYPEs along internal boundaries
!!> -----------------------------------------------------------------------------------------------------
!IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)
!
!CALL SCARC_LEAVE_ROUTINE()
!END SUBROUTINE SCARC_SETUP_CELL_TYPES


!>> ------------------------------------------------------------------------------------------------
!>> Set cell types for A1 coarsening
!>> ------------------------------------------------------------------------------------------------
!SUBROUTINE SET_TYPES_A1(MEASURE, CELL_TYPE, JC, OFFSET)
!REAL(EB), DIMENSION(:), INTENT(OUT) :: MEASURE
!INTEGER , DIMENSION(:), INTENT(OUT) :: CELL_TYPE
!INTEGER , INTENT(IN) :: JC, OFFSET
!
!MEASURE(JC)         = 0.0_EB
!MEASURE(JC-OFFSET)  = 0.0_EB
!MEASURE(JC+OFFSET)  = 0.0_EB
!
!CELL_TYPE(JC)        = NSCARC_CELL_TYPE_SFINE
!CELL_TYPE(JC-OFFSET) = NSCARC_CELL_TYPE_WFINE
!CELL_TYPE(JC+OFFSET) = NSCARC_CELL_TYPE_WFINE
!
!END SUBROUTINE SET_TYPES_A1


!> ------------------------------------------------------------------------------------------------
!> Get maximum measure in internal of mesh (cells along boundaries not regarded)
!> ------------------------------------------------------------------------------------------------
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


!> ------------------------------------------------------------------------------------------------
!> Get cell number of internal neighbor for ghost cell IW
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION GET_WALL_CELL(WALL, IW, NX, NY)
TYPE (SCARC_WALL_TYPE), DIMENSION(:), INTENT(IN) :: WALL
INTEGER, INTENT(IN) :: IW, NX, NY
INTEGER :: IC

IF (TWO_D) THEN
   IF (ABS(WALL(IW)%NOM) == 2) THEN
      IC = -1
   ELSE
      IC = (WALL(IW)%IZW-1)*NX + WALL(IW)%IXW
   ENDIF
ELSE
   IC = (WALL(IW)%IZW-1)*NX*NY + (WALL(IW)%IYW-1)*NX + WALL(IW)%IXW
ENDIF

GET_WALL_CELL = IC
RETURN

END FUNCTION GET_WALL_CELL


!> ------------------------------------------------------------------------------------------------
!> Perform interpolation to coarse level
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PROLONGATION(NL)
INTEGER , INTENT(IN) :: NL
INTEGER  :: NM, IP, IC, JC, KC, ICG, ICOL, JCOL, KCOL, ICA, SIGN0
INTEGER  :: IW, IW0, IW2, IDIAG, JDIAG, JCO, ICOL0, ICOL_FIRST, ICOL_LAST
REAL(EB) :: SUM_COUPLED, SUM_CPOINTS, SCAL, SUM_COARSE, SUM_DIAG
REAL(EB) :: DATA_SUM, DATA_DIAG, DATA_INTERPOL
REAL(EB) :: VALUES(20), WEIGHTS(20)
INTEGER  :: NEIGHBOR(20), NWEIGHTS, NWEIGHTS2
INTEGER  :: COARSE_CELL(20), COARSE_INDEX(20)
TYPE (SCARC_LEVEL_TYPE),  POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, P, R

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_PROLONGATION')

SELECT_INTERPOLATION: SELECT CASE (TYPE_INTERPOL)

   !> ------------------------------------------------------------------------------------------
   !> Standard interpolation
   !> ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_STANDARD)

      MESHES_LOOP_STANDARD: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_SHUTDOWN('Standard interpolation not yet implemented', 'NONE', -999)
      ENDDO MESHES_LOOP_STANDARD

   !> ------------------------------------------------------------------------------------------
   !> Classical interpolation
   !> ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_CLASSICAL)

      MESHES_LOOP_CLASSICAL: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%SYSTEM(NL)%A
         P => SCARC(NM)%SYSTEM(NL)%P
         R => SCARC(NM)%SYSTEM(NL)%R

         IP = 1
         INTERNAL_CELL_LOOP_CLASSICAL: DO IC = 1, L%NC

            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0

            !>!
            !> If IC is a coarse cell, its value is taken
            !>!
            IF (L%AMG%CTYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               NEIGHBOR(1)= L%AMG%CTYPE(IC)
               WEIGHTS(1) = 1.0_EB
               NWEIGHTS   = 1

            !>!
            !> If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
            !>!
            ELSE

             !> Get main diagonal entry a_ii for that fine cell
               IDIAG = A%ROW(IC)
               SUM_DIAG = A%VAL(IDIAG)

               !> First search for all neighboring coarse grid cells (store them in NEIGHBOR)
               IW = 1
               DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                  JC = A%COL(ICOL)
                  IF (L%AMG%CTYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     NEIGHBOR(IW) = JC
                     WEIGHTS(IW)  = -A%VAL(ICOL)
                     IW = IW + 1
                  ENDIF
               ENDDO
               NWEIGHTS = IW - 1

               !> Then search for the strongly and weakly coupled fine grid cells
               DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1

                  IW2 = 1
                  SUM_COARSE = 0.0_EB

                  JC = A%COL(ICOL)

                  SELECT CASE (L%AMG%CTYPE(JC))

                     CASE (NSCARC_CELL_TYPE_SFINE)

                        !> search for couplings KC of the strongly coupled JC which belong to the
                        !> coarse interpolatory set of IC
                        DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
                           KC = A%COL(JCOL)

                           NWEIGHTS_SFINE_LOOP: DO IW = 1, NWEIGHTS
                              IF (KC == NEIGHBOR(IW) ) THEN
                                 COARSE_CELL (IW2) = JCOL
                                 COARSE_INDEX(IW2) = IW
                                 SUM_COARSE = SUM_COARSE + A%VAL(JCOL)
                                 IW2 = IW2 + 1
                                 EXIT NWEIGHTS_SFINE_LOOP
                              ENDIF
                           ENDDO NWEIGHTS_SFINE_LOOP
                        ENDDO
                        NWEIGHTS2 = IW2 - 1

                        DO IW2 = 1, NWEIGHTS2
                           JCOL = COARSE_CELL(IW2)
                           IW   = COARSE_INDEX(IW2)
                           WEIGHTS(IW) = WEIGHTS(IW) - A%VAL(ICOL)*A%VAL(JCOL)/REAL(SUM_COARSE,EB)
                        ENDDO

                     CASE (NSCARC_CELL_TYPE_WFINE)

                        SUM_DIAG = SUM_DIAG + A%VAL(ICOL)

                  END SELECT
               ENDDO

               DO IW = 1, NWEIGHTS
                  NEIGHBOR(IW) = L%AMG%CTYPE(NEIGHBOR(IW))
                  WEIGHTS(IW)  = WEIGHTS(IW)/SUM_DIAG
               ENDDO

            ENDIF

            !>
            !> Define corresponding entry in interpolation matrix P by means of the upper weights
            !>!
            P%ROW(IC) = IP
            DO IW = 1, NWEIGHTS
               IF  (NEIGHBOR(IW) /= -1) THEN
                  P%COL(IP) = NEIGHBOR(IW)
                  P%VAL(IP)     = WEIGHTS(IW)
                  IP = IP +1
               ENDIF
            ENDDO
            P%ROW(L%NCE+1) = IP

         ENDDO INTERNAL_CELL_LOOP_CLASSICAL

      ENDDO MESHES_LOOP_CLASSICAL

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)

   !> ------------------------------------------------------------------------------------------
   !> Classical interpolation2
   !> ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_CLASSICAL2)

      MESHES_LOOP_CLASSICAL2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%SYSTEM(NL)%A
         P => SCARC(NM)%SYSTEM(NL)%P

         !> Build Interpolation
         ICOL0 = 1
         INTERNAL_CELL_LOOP2_CLASSICAL2: DO IC = 1, L%NC

            !> For a C-point IC, use identity and set mapping
            IF (L%AMG%CTYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               P%ROW(IC) = ICOL0
               P%COL(ICOL0) = L%AMG%CTYPE(IC)
               P%VAL(ICOL0) = 1.0_EB

               ICOL0 = ICOL0 + 1

            !> For a F-point IC, build interpolation
            ELSE

               P%ROW(IC)= ICOL0               !> diagonal part of P
               ICOL_FIRST = ICOL0

               DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                  JC = A%COL(ICOL)

                  !> If JC is a C-point, initialize interpolation to zero and set column number in P%COL
                  IF (L%AMG%CTYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     L%AMG%MARKER(JC) = ICOL0
                     P%COL(ICOL0) = L%AMG%CTYPE(JC)
                     P%VAL(ICOL0) = 0.0_EB
                     ICOL0 = ICOL0 + 1
                  !> If JC is a F-point, set it to be a strong F-point with relevant connections
                  ELSE IF (L%AMG%CTYPE(JC) /= NSCARC_CELL_TYPE_WFINE) THEN
                     L%AMG%MARKER(JC) = NSCARC_CELL_TYPE_SFINE
                  ENDIF

               ENDDO
               ICOL_LAST = ICOL0 - 1

               !> consider the i-th row of A, start with diagonal part
               IDIAG = A%ROW(IC)
               DATA_DIAG = A%VAL(IDIAG)

               DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                  JC = A%COL(ICOL)

                  !> If JC is a strongly coupled C-point to IC, the value a_(IC, JC) must be
                  !> incorporated to the interpolation weight
                  IF (L%AMG%MARKER(JC) >= ICOL_FIRST) THEN
                     JCOL = L%AMG%MARKER(JC)
                     P%VAL(JCOL) = P%VAL(JCOL) + A%VAL(ICOL)

                  !> If JC is a strongly coupled F-point to IC, the value a_(IC, JC) must be distributed
                  !> to C-points which are strongly coupled to IC (no distribution to the diagonal part)
                  ELSE IF (L%AMG%MARKER(JC) >= NSCARC_CELL_TYPE_SFINE) THEN
                     DATA_SUM  = 0.0_EB
                     SIGN0 = 1
                     IF (A%VAL(IDIAG) > 0) SIGN0=-1

                     !> search JC in the row of A and get sum of the couplings to C-points strongly coupled to IC
                     DO JCOL = A%ROW(JC), A%ROW(JC+1)-1
                        KC = A%COL(JCOL)
                        IF (L%AMG%MARKER(KC) > ICOL_FIRST .AND. SIGN0*A%VAL(JCOL) >0) THEN
                           DATA_SUM = DATA_SUM + A%VAL(JCOL)
                        ENDIF
                     ENDDO

                     IF (DATA_SUM .NE. 0) THEN

                        DATA_INTERPOL = A%VAL(ICOL)/DATA_SUM

                        !> loop over row of A for JC and spread data
                        DO JCOL = A%ROW(JC), A%ROW(JC+1)-1
                           KC = A%COL(JCOL)
                           IF (L%AMG%MARKER(KC) >= ICOL_FIRST .AND. SIGN0*A%VAL(JCOL) > 0) THEN
                              KCOL = L%AMG%MARKER(KC)
                              P%VAL(KCOL) = P%VAL(KCOL) + DATA_INTERPOL * A%VAL(JCOL)
                           ENDIF
                        ENDDO


                     ENDIF

                  !> If JC is a weakly coupled F-point to IC, the value a_(IC, JC) must be
                  !> distributed to the diagonal
                  ELSE IF (L%AMG%CTYPE(JC) .NE. NSCARC_CELL_TYPE_WFINE) THEN
                     DATA_DIAG = DATA_DIAG + A%VAL(ICOL)
                  ENDIF

               ENDDO

               IF (DATA_DIAG == 0.0_EB) THEN
                  WRITE(*,*) "SCARC_SETUP_PROLONGATION: WARNING - zero diagonal data !"
                  DO ICOL = ICOL_FIRST, ICOL_LAST
                     P%VAL(ICOL) = 0.0_EB
                  ENDDO
               ELSE
                  DO ICOL = ICOL_FIRST, ICOL_LAST
                     P%VAL(ICOL) = - P%VAL(ICOL)/DATA_DIAG
                  ENDDO
               ENDIF
            ENDIF

         ENDDO INTERNAL_CELL_LOOP2_CLASSICAL2
         P%ROW(IC) = ICOL0

      ENDDO MESHES_LOOP_CLASSICAL2

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_PROL, NL)

   !> ---------------------------------------------------------------------------------------------
   !> Direct interpolation
   !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_DIRECT)

      MESHES_LOOP_DIRECT: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%SYSTEM(NL)%A
         P => SCARC(NM)%SYSTEM(NL)%P

         IP = 1
         INTERNAL_CELL_LOOP: DO IC = 1, L%NC

            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0

            !>!
            !> If IC is a coarse cell, its value is taken
            !>!
            IF (L%AMG%CTYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               NEIGHBOR(1)= L%AMG%CTYPE(IC)
               NEIGHBOR(1)= IC
               WEIGHTS(1) = 1.0_EB
               NWEIGHTS   = 1

            !>!
            !> If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
            !>!
            ELSE

               !> Get main diagonal entry a_ii for that fine cell
               IDIAG = A%ROW(IC)

               !> Select type of fine cell (weakly/strongly coupled)
               SELECT_FPOINT_TYPE: SELECT CASE(L%AMG%CTYPE(IC))

                  !>!
                  !> Strongly coupled fine cell IC
                  !> approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                  !> Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                  !>!
                  CASE (NSCARC_CELL_TYPE_SFINE)

                     !> Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                     SUM_COUPLED = 0.0_EB
                     SUM_CPOINTS = 0.0_EB

                     DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                        JC = A%COL(ICOL)
                        SUM_COUPLED = SUM_COUPLED + A%VAL(ICOL)
                        IF (L%AMG%CTYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + A%VAL(ICOL)
                     ENDDO

                     SCAL = - SUM_COUPLED/SUM_CPOINTS

                     !> for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                     IW = 1
                     DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                        JC = A%COL(ICOL)
                        IF (L%AMG%CTYPE(JC) > 0 ) THEN
                           NEIGHBOR(IW) = L%AMG%CTYPE(JC)
                           NEIGHBOR(IW) = JC
                           WEIGHTS(IW)  = SCAL * A%VAL(ICOL)/A%VAL(IDIAG)
                           IW = IW +1
                        ENDIF
                     ENDDO
                     NWEIGHTS = IW - 1

                  !>!
                  !> Weakly coupled fine cell IC:
                  !> Determine strongly coupled fine cells JC surrounding IC and, in turn, replace
                  !> each of them by a mean value of their surrounding strongly coupled coarse cells
                  !>!
                  CASE (NSCARC_CELL_TYPE_WFINE)

                     IW = 1                                             !> weights counter
                     DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1         !> loop over couplings of weakly coupled fine cell

                        !> regard subdiagonal matrix entries a_ij of weakly coupled fine cell JC to IC
                        JC = A%COL(ICOL)

                        !> Find all surrounding (coupled) points of JC and compute scaling factor
                        !> compute scaling factor SCAL = a_ij/a_ii * [sum_(k in N_j) a_jk]/[sum_(l in C_j) ajl]
                        SUM_COUPLED = 0.0_EB
                        SUM_CPOINTS = 0.0_EB

                        DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
                           KC = A%COL(JCOL)
                           SUM_COUPLED = SUM_COUPLED + A%VAL(JCOL)
                           IF (L%AMG%CTYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + A%VAL(JCOL)
                        ENDDO

                        IF (SUM_CPOINTS == 0.0_EB) &
                           CALL SCARC_SHUTDOWN('Error in interpolation routine for cell ', 'NONE', IC)
                        SCAL =  A%VAL(ICOL)/A%VAL(IDIAG) * SUM_COUPLED/SUM_CPOINTS

                        !> Get diagonal matrix a_jj for point JC
                        JDIAG = A%ROW(JC)

                        !> Compute interpolation weights for all strong coarse cell couplings KC of JC
                        !> note that a coarse cell KC may be considered several times for different JC's
                        DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
                           KC = A%COL(JCOL)
                           IF (L%AMG%CTYPE(KC) > 0) THEN
                             NEIGHBOR(IW) = L%AMG%CTYPE(KC)
                             NEIGHBOR(IW) = KC
                             WEIGHTS(IW)  = SCAL * A%VAL(JCOL)/A%VAL(JDIAG)
                             IW = IW +1
                           ENDIF
                        ENDDO

                     ENDDO
                     NWEIGHTS = IW - 1

                     !> make weights unique (add weights for multiple coarse cells)
                     DO IW0 = 1, NWEIGHTS
                        DO IW = IW0+1, NWEIGHTS
                           IF  (NEIGHBOR(IW0) /= -1 .AND. NEIGHBOR(IW) == NEIGHBOR(IW0)) THEN
                              WEIGHTS(IW0) = WEIGHTS(IW0) + WEIGHTS(IW)
                              WEIGHTS(IW)  = 0.0_EB
                              NEIGHBOR(IW) = -1
                           ENDIF
                        ENDDO
                     ENDDO

               END SELECT SELECT_FPOINT_TYPE

            ENDIF


            !>
            !> Define corresponding entry in interpolation matrix P by means of the upper weights
            !>!
            P%ROW(IC) = IP
            DO IW = 1, NWEIGHTS
               IF  (NEIGHBOR(IW) /= -1) THEN
                  P%COL(IP) = NEIGHBOR(IW)
                  P%VAL(IP)     = WEIGHTS(IW)
                  IP = IP +1
               ENDIF
            ENDDO

         ENDDO INTERNAL_CELL_LOOP

         P%ROW(L%NC+1) = IP
         !L%NROW0 = IP                !> really neede - OUTDATED ?

         !> Determine new coarse number and store it in CELL_TYPE
         INTERNAL_CELL_LOOP2: DO IC = 1, L%NC
            DO ICOL = P%ROW(IC), P%ROW(IC+1)-1
              JC = P%COL(ICOL)
              IF (JC <= L%NC) THEN
                 P%COL(ICOL) = L%AMG%CTYPE(JC)
              ENDIF
            ENDDO
         ENDDO INTERNAL_CELL_LOOP2


         !> In the multi-mesh case determine the new coarse cell numbers for neighboring cells
         IF (NMESHES > 1) THEN
            GHOST_CELL_LOOP: DO IW = 1, L%NW
               IF (L%WALL(IW)%NOM==0) CYCLE GHOST_CELL_LOOP
               ICG = L%WALL(IW)%ICG(1)
               ICA = L%WALL(IW)%ICW
               DO ICOL = P%ROW(ICA), P%ROW(ICA+1)-1
                 JC = P%COL(ICOL)
                 IF (JC > L%NC) THEN
                    JCO = L%CELL%DOF(L%WALL(IW)%IXG,L%WALL(IW)%IYG,L%WALL(IW)%IZG)  ! ACHTUNG: SCHAUEN!
                    P%COL(ICOL) = - JCO
                 ENDIF
               ENDDO
            ENDDO GHOST_CELL_LOOP
         ENDIF

      ENDDO MESHES_LOOP_DIRECT

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)

 !> ---------------------------------------------------------------------------------------------
 !> Direct interpolation with special treatment of boundary layers
 !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_DIRECT_BDRY)

      MESHES_LOOP_DIRECT_BDRY: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%SYSTEM(NL)%A
         P => SCARC(NM)%SYSTEM(NL)%P

         IP = 1
         INTERNAL_CELL_BDRY_LOOP: DO IC = 1, L%NC

            DIRECT_BDRY_IF: IF (L%AMG%INTERN(IC)/=0) THEN

               VALUES   = 0.0_EB
               WEIGHTS  = 0.0_EB
               NEIGHBOR = 0

               !>!
               !> If IC is a coarse cell, its value is taken
               !>!
               IF (L%AMG%CTYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

                  NEIGHBOR(1)= L%AMG%CTYPE(IC)
                  NEIGHBOR(1)= IC
                  WEIGHTS(1) = 1.0_EB
                  NWEIGHTS   = 1


               !>!
               !> If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
               !>!
               ELSE


                  !> Get main diagonal entry a_ii for that fine cell
                  IDIAG = A%ROW(IC)

                  !> Strongly coupled fine cell IC
                  !> approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                  !> Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                  !> Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                  SUM_COUPLED = 0.0_EB
                  SUM_CPOINTS = 0.0_EB

                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     JC = A%COL(ICOL)
                     SUM_COUPLED = SUM_COUPLED + A%VAL(ICOL)
                     IF (JC > L%NC) CYCLE
                     IF (L%AMG%INTERN(JC)/=0.AND.L%AMG%CTYPE(JC)>0) SUM_CPOINTS = SUM_CPOINTS + A%VAL(ICOL)
                  ENDDO

                  SCAL = - SUM_COUPLED/SUM_CPOINTS

                  !> for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                  IW = 1
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     JC = A%COL(ICOL)
                     IF (JC > L%NC) CYCLE
                     IF (L%AMG%INTERN(JC)/=0.AND.L%AMG%CTYPE(JC) > 0 ) THEN
                        NEIGHBOR(IW) = L%AMG%CTYPE(JC)
                        NEIGHBOR(IW) = JC
                        WEIGHTS(IW)  = SCAL * A%VAL(ICOL)/A%VAL(IDIAG)
                        IW = IW +1
                     ENDIF
                  ENDDO
                  NWEIGHTS = IW - 1

               ENDIF


             !> ----------------------------!>> real internal cell ------------------------
            ELSE

               VALUES   = 0.0_EB
               WEIGHTS  = 0.0_EB
               NEIGHBOR = 0

               !>!
               !> If IC is a coarse cell, its value is taken
               !>!
               IF (L%AMG%CTYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

                  NEIGHBOR(1)= L%AMG%CTYPE(IC)
                  NEIGHBOR(1)= IC
                  WEIGHTS(1) = 1.0_EB
                  NWEIGHTS   = 1

               !>!
               !> If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
               !>!
               ELSE

                  !> Get main diagonal entry a_ii for that fine cell
                  IDIAG = A%ROW(IC)

                  !> Select type of fine cell (weakly/strongly coupled)
                  SELECT_FPOINT_BDRY_TYPE: SELECT CASE(L%AMG%CTYPE(IC))

                     !>!
                     !> Strongly coupled fine cell IC
                     !> approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                     !> Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                     !>!
                     CASE (NSCARC_CELL_TYPE_SFINE)

                        !> Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                        SUM_COUPLED = 0.0_EB
                        SUM_CPOINTS = 0.0_EB

                        DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                           JC = A%COL(ICOL)
                           SUM_COUPLED = SUM_COUPLED + A%VAL(ICOL)
                           IF (L%AMG%CTYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + A%VAL(ICOL)
                        ENDDO

                        SCAL = - SUM_COUPLED/SUM_CPOINTS

                        !> for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                        IW = 1
                        DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                           JC = A%COL(ICOL)
                           IF (L%AMG%CTYPE(JC) > 0 ) THEN
                              NEIGHBOR(IW) = L%AMG%CTYPE(JC)
                              NEIGHBOR(IW) = JC
                              WEIGHTS(IW)  = SCAL * A%VAL(ICOL)/A%VAL(IDIAG)
                              IW = IW +1
                           ENDIF
                        ENDDO
                        NWEIGHTS = IW - 1

                     !>!
                     !> Weakly coupled fine cell IC:
                     !> Determine strongly coupled fine cells JC surrounding IC and, in turn, replace
                     !> each of them by a mean value of their surrounding strongly coupled coarse cells
                     !>!
                     CASE (NSCARC_CELL_TYPE_WFINE)

                        IW = 1                                             !> weights counter
                        DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1         !> loop over couplings of weakly coupled fine cell

                           !> regard subdiagonal matrix entries a_ij of weakly coupled fine cell JC to IC
                           JC = A%COL(ICOL)

                           !> Find all surrounding (coupled) points of JC and compute scaling factor
                           !> compute scaling factor SCAL = a_ij/a_ii * [sum_(k in N_j) a_jk]/[sum_(l in C_j) ajl]
                           SUM_COUPLED = 0.0_EB
                           SUM_CPOINTS = 0.0_EB

                           DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
                              KC = A%COL(JCOL)
                              SUM_COUPLED = SUM_COUPLED + A%VAL(JCOL)
                              IF (L%AMG%CTYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + A%VAL(JCOL)
                           ENDDO

                           IF (SUM_CPOINTS == 0.0_EB) &
                              CALL SCARC_SHUTDOWN('Error in interpolation routine for cell ', 'NONE', IC)
                           SCAL =  A%VAL(ICOL)/A%VAL(IDIAG) * SUM_COUPLED/SUM_CPOINTS

                           !> Get diagonal matrix a_jj for point JC
                           JDIAG = A%ROW(JC)

                           !> Compute interpolation weights for all strong coarse cell couplings KC of JC
                           !> note that a coarse cell KC may be considered several times for different JC's
                           DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
                              KC = A%COL(JCOL)
                              IF (L%AMG%CTYPE(KC) > 0) THEN
                                NEIGHBOR(IW) = L%AMG%CTYPE(KC)
                                NEIGHBOR(IW) = KC
                                WEIGHTS(IW)  = SCAL * A%VAL(JCOL)/A%VAL(JDIAG)
                                IW = IW +1
                              ENDIF
                           ENDDO

                        ENDDO
                        NWEIGHTS = IW - 1

                        !> make weights unique (add weights for multiple coarse cells)
                        DO IW0 = 1, NWEIGHTS
                           DO IW = IW0+1, NWEIGHTS
                              IF  (NEIGHBOR(IW0) /= -1 .AND. NEIGHBOR(IW) == NEIGHBOR(IW0)) THEN
                                 WEIGHTS(IW0) = WEIGHTS(IW0) + WEIGHTS(IW)
                                 WEIGHTS(IW)  = 0.0_EB
                                 NEIGHBOR(IW) = -1
                              ENDIF
                           ENDDO
                        ENDDO

                  END SELECT SELECT_FPOINT_BDRY_TYPE

               ENDIF


            ENDIF DIRECT_BDRY_IF
            !>
            !> Define corresponding entry in interpolation matrix P by means of the upper weights
            !>!
            P%ROW(IC) = IP
            DO IW = 1, NWEIGHTS
               IF  (NEIGHBOR(IW) /= -1) THEN
                  P%COL(IP) = NEIGHBOR(IW)
                  P%VAL(IP)     = WEIGHTS(IW)
                  IP = IP +1
               ENDIF
            ENDDO

         ENDDO INTERNAL_CELL_BDRY_LOOP

         P%ROW(L%NC+1) = IP
         !L%NROW0 = IP                !> really neede - OUTDATED ?

         !> Determine new coarse number and store it in CELL_TYPE
         INTERNAL_CELL_BDRY_LOOP2: DO IC = 1, L%NC
            DO ICOL = P%ROW(IC), P%ROW(IC+1)-1
              JC = P%COL(ICOL)
              IF (JC <= L%NC) THEN
                 P%COL(ICOL) = L%AMG%CTYPE(JC)
              ENDIF
            ENDDO
         ENDDO INTERNAL_CELL_BDRY_LOOP2


      ENDDO MESHES_LOOP_DIRECT_BDRY

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)

END SELECT SELECT_INTERPOLATION

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_PROLONGATION

!> ------------------------------------------------------------------------------------------------
!> Build transpose R of prolongation matrix P
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATRIX_TRANSPOSE(M, M_ROW, M_COL, MT, MT_ROW, MT_COL, NCE, NCCE)
INTEGER,  DIMENSION(:)  , INTENT(IN)  :: M_ROW , M_COL
INTEGER,  DIMENSION(:)  , INTENT(OUT) :: MT_ROW, MT_COL
REAL(EB), DIMENSION(:)  , INTENT(IN)  :: M
REAL(EB), DIMENSION(:)  , INTENT(OUT) :: MT
INTEGER, INTENT(IN)  :: NCE, NCCE
INTEGER :: IC, JC, ICOL, JCOL

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_MATRIX_TRANSPOSE')

!> identify the number of non-zero entries in every column of M (corresponds to a row in MT)
!> And store it in the MT_ROW-array (caution: starts from the second position)
MT_ROW(1) = 1
DO ICOL = 1, M_ROW(NCE+1)-1
   IC = M_COL(ICOL)
   MT_ROW(IC+1) = MT_ROW(IC+1) + 1
ENDDO

!> shift it to the first position while summing up the entries for M_ROW
DO IC = 2, NCCE+1
   MT_ROW(IC) = MT_ROW(IC) + MT_ROW(IC-1)
ENDDO

!> set correct entries for MT
DO IC = 1, NCE
   DO ICOL = M_ROW(IC), M_ROW(IC+1)-1
      JC   = M_COL(ICOL)
      JCOL = MT_ROW(JC)
      MT_COL(JCOL) = IC
      MT(JCOL)     = M(ICOL)
      MT_ROW(JC)   = MT_ROW(JC) + 1
   ENDDO
ENDDO

DO IC = NCCE, 2, -1
   MT_ROW(IC) = MT_ROW(IC-1)
ENDDO
MT_ROW(1)=1

END SUBROUTINE SCARC_MATRIX_TRANSPOSE


!> ------------------------------------------------------------------------------------------------
!> Build transpose R of prolongation matrix P
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM_AMG(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM
INTEGER :: IROW, IROW_INIT=0, IROW_SAVE
INTEGER :: ICOL1, ICOL2, ICOL3, IC0, IC1, IC2, IC3
LOGICAL :: BSQUARE = .TRUE.
REAL(EB) :: R_VALUE, RA_VALUE, RAP_VALUE
TYPE (SCARC_LEVEL_TYPE) , POINTER :: LF, LC
TYPE (SCARC_MATRIX_TYPE), POINTER :: AF, AC, PF, PC, RF, RC

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_SYSTEM_AMG')

MESHES_LOOP_SYSTEM_AMG: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LF => SCARC(NM)%LEVEL(NL)                     !> fine level
   LC => SCARC(NM)%LEVEL(NL+1)                   !> coarse level

   AF => SCARC(NM)%SYSTEM(NL)%A
   AC => SCARC(NM)%SYSTEM(NL+1)%A

   PF => SCARC(NM)%SYSTEM(NL)%P
   PC => SCARC(NM)%SYSTEM(NL+1)%P

   RF => SCARC(NM)%SYSTEM(NL)%R
   RC => SCARC(NM)%SYSTEM(NL+1)%R

   IROW = 1

!> 
!> Still experimental code:
!> First look at all exterior entries
!> 
!>  IF (LC%NW > 0) THEN
!>
!>     LOOP_NBR_CELLS: DO IC0 = 1, L%AMG%NCC
!>
!>        IF (L%AMG%INTERN(IC) < 0) CYCLE LOOP_NBR_CELLS
!>
!>        IROW_INIT = IROW
!>
!>      !> loop over all entries in row IC0 of R
!>        LOOP_NBR_R_ENTRIES: DO ICOL1 = RF%ROW(IC0), RF%ROW(IC0+1)-1
!>
!>           IC1 = RF%COL(ICOL1)
!>           IF (BDEBUG_LESS) WRITE(MSG%LU_DEBUG,*) '--------- ICOL1 =', ICOL1,'   : IC1=',IC1
!>
!>         !> loop over all entries in row IC1 of A
!>           LOOP_NBR_A_ENTRIES: DO ICOL2 = A%ROW(IC1), A%ROW(IC1+1)-1
!>
!>              IC2 = AF%COL(ICOL2)
!>              IF (BDEBUG_LESS) WRITE(MSG%LU_DEBUG,*) '   ------ ICOL2 =', ICOL2,'   : IC2=',IC2
!>
!>              IF (AF%TAG(IC2) /= IC0) THEN
!>
!>               !> mark IC2 as already considered
!>                 AF%TAG(IC2) = IC0
!>                 IF (BDEBUG_LESS) WRITE(MSG%LU_DEBUG,*) '          A%TAG(',IC2,') =', AF%TAG(IC2)
!>
!>               !> loop over all entries in row IC2 of P
!>                 LOOP_NBR_P_ENTRIES: DO ICOL3 = PF%ROW(IC2), PF%ROW(IC2+1)-1
!>
!>                    IC3 = PF%COL(ICOL3)
!>                    IF (BDEBUG_LESS) WRITE(MSG%LU_DEBUG,*) '    !> ---- ICOL3 =', ICOL3,'   : IC3=',IC3
!>                    !IF (BDEBUG_LESS) WRITE(MSG%LU_DEBUG,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
!>                  !>   '          P_TAG(',IC3,') =', PF%TAG(IC3),': IROW_INIT=',IROW_INIT,' : IROW=',IROW
!>
!>                  !> verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
!>                    IF (PF%TAG(IC3) < IROW_INIT) THEN
!>                       PF%TAG(IC3) = IROW
!>                       IF (BDEBUG_LESS) WRITE(MSG%LU_DEBUG,*) '          P_TAG(',IC3,') =', PF%TAG(IC3)
!>                       IROW = IROW + 1
!>                       IF (BDEBUG_LESS) WRITE(MSG%LU_DEBUG,*) '          IROW =', IROW
!>                    ENDIF
!>                 ENDDO LOOP_NBR_P_ENTRIES
!>
!>              ENDIF
!>
!>           ENDDO LOOP_NBR_A_ENTRIES
!>        ENDDO LOOP_NBR_R_ENTRIES
!>
!>      !> Store counters
!>        IROW_SAVE = IROW
!>
!>
!>     ENDDO LOOP_NBR_CELLS
!>
!>  ENDIF


   !> Determine size of matrix RAP
   !> loop over interior c-cells
   LOOP1_OWN_CELLS: DO IC0 = 1, LF%AMG%NCCE

      IROW_INIT = IROW

      IF (BSQUARE) THEN
         PF%TAG(IC0) = IROW
         IROW = IROW + 1
      ENDIF

      !> loop over all entries in row IC0 of R
      LOOP1_R_ENTRIES: DO ICOL1 = RF%ROW(IC0), RF%ROW(IC0+1)-1
         IC1 = RF%COL(ICOL1)

         !> loop over all entries in row IC1 of A
         LOOP1_A_ENTRIES: DO ICOL2 = AF%ROW(IC1), AF%ROW(IC1+1)-1
            IC2 = AF%COL(ICOL2)

            IF (AF%TAG(IC2) /= IC0) THEN
               AF%TAG(IC2) = IC0                           !> mark IC2 as already considered

               !> loop over all entries in row IC2 of P
               LOOP1_P_ENTRIES: DO ICOL3 = PF%ROW(IC2), PF%ROW(IC2+1)-1

                  IC3 = PF%COL(ICOL3)

                  !> verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
                  IF (IC3 >0) THEN   !> ONLY TEMPORARILY
                  IF (PF%TAG(IC3) < IROW_INIT) THEN
                     PF%TAG(IC3) = IROW
                     IROW = IROW + 1
                  ENDIF
                  ENDIF   !> ONLY TEMPORARILY
               ENDDO LOOP1_P_ENTRIES

            ENDIF

         ENDDO LOOP1_A_ENTRIES
      ENDDO LOOP1_R_ENTRIES

      !> Store counters
      IROW_SAVE = IROW

   ENDDO LOOP1_OWN_CELLS

   !> Determine size of matrix RAP  == LC%A
   AC%NAV = IROW
   AC%NAE = AC%NAV + 100  !> ONLY TEMPORARILY

   CALL SCARC_ALLOCATE_MATRIX(AC, AC%NAE+1, AC%NAE+1, LC%NCE+1, -1, NSCARC_INIT_ZERO, 'AC')
   AC%ROW(LF%AMG%NCC+1) = IROW

   AF%TAG = 0
   PF%TAG = 0

   IROW = 1

   !> loop over interior c-cells
   LOOP2_C_CELLS: DO IC0 = 1, LF%AMG%NCC
      IROW_INIT = IROW
      AC%ROW(IC0) = IROW_INIT

      IF (BSQUARE) THEN
         PF%TAG(IC0) = IROW
         AC%COL(IROW) = IC0
         AC%VAL(IROW) = 0.0_EB
         IROW = IROW + 1
      ENDIF

      !> loop over all entries in row IC0 of R
      LOOP2_R_ENTRIES: DO ICOL1 = RF%ROW(IC0), RF%ROW(IC0+1)-1

         IC1 = RF%COL(ICOL1)
         R_VALUE = RF%VAL(ICOL1)

         !> loop over all entries in row IC1 of A
         LOOP2_A_ENTRIES: DO ICOL2 = AF%ROW(IC1), AF%ROW(IC1+1)-1

            IC2 = AF%COL(ICOL2)
            RA_VALUE = R_VALUE * AF%VAL(ICOL2)

            !> Hasn't cell IC2 been considered before? (new values for RAP can only be done for unmarked cells)
            IF (AF%TAG(IC2) /= IC0) THEN

               !> mark IC2 as already considered
               AF%TAG(IC2) = IC0

               !> loop over all entries in row IC2 of P
               LOOP2_P_ENTRIES1: DO ICOL3 = PF%ROW(IC2), PF%ROW(IC2+1)-1

                  IC3 = PF%COL(ICOL3)
                  RAP_VALUE = RA_VALUE * PF%VAL(ICOL3)

                  !> verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
                  IF (IC3 >0) THEN     !> ONLY TEMPORARILY

                  IF (PF%TAG(IC3) < IROW_INIT) THEN
                     PF%TAG(IC3)  = IROW
                     AC%COL(IROW) = PF%COL(ICOL3)
                     AC%VAL(IROW) = RAP_VALUE
                     IROW = IROW + 1
                  ELSE
                     AC%VAL(PF%TAG(IC3)) = AC%VAL(PF%TAG(IC3)) + RAP_VALUE
                  ENDIF

                  ENDIF  !> ONLY TEMPORARILY
               ENDDO LOOP2_P_ENTRIES1

            !> or has it been already considered
            ELSE

               !> loop over all entries in row IC2 of P
               LOOP2_P_ENTRIES2: DO ICOL3 = PF%ROW(IC2), PF%ROW(IC2+1)-1
                  IC3 = PF%COL(ICOL3)
                  RAP_VALUE = RA_VALUE * PF%VAL(ICOL3)
                  AC%VAL(PF%TAG(IC3)) = AC%VAL(PF%TAG(IC3)) + RAP_VALUE
               ENDDO LOOP2_P_ENTRIES2
            ENDIF

         ENDDO LOOP2_A_ENTRIES
      ENDDO LOOP2_R_ENTRIES

      !> Store counters
      IROW_SAVE = IROW

   ENDDO LOOP2_C_CELLS
ENDDO MESHES_LOOP_SYSTEM_AMG

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_SYSTEM_AMG

!> ------------------------------------------------------------------------------------------------
!> Define restriction matrix R (currently transpose of prolongation matrix P)
!>  - In spite of the additinal need for the storing of R, this is done to save computational time
!>  - during the later matrix transfer operations
!> temporarily not used
!> ------------------------------------------------------------------------------------------------
!SUBROUTINE SCARC_SETUP_RESTRICTION(NL)
!INTEGER , INTENT(IN) :: NL
!INTEGER  :: NM, NOM, ICP, ICP2, IC, IC2, IG, IG0, ICC
!LOGICAL  :: BFIRST
!TYPE (OSCARC_TYPE), POINTER :: OS
!TYPE (SCARC_LEVEL_TYPE)   , POINTER :: L, OL
!TYPE (SCARC_MATRIX_TYPE)  , POINTER :: P, R, OP, OR
!TYPE (SCARC_EXCHANGE_TYPE), POINTER :: OX
!
!CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_RESTRICTION')
!
!SELECT CASE (TYPE_INTERPOL)
!
! !> ------------------------------------------------------------------------------------------
! !> Classical
! !> ------------------------------------------------------------------------------------------
!   CASE (NSCARC_INTERPOL_CLASSICAL2)
!
!      MESHES_LOOP_CLASSICAL2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!         L => SCARC(NM)%LEVEL(NL)
!         P => SCARC(NM)%SYSTEM(NL)%P
!         R => SCARC(NM)%SYSTEM(NL)%R
!         CALL SCARC_MATRIX_TRANSPOSE(P%VAL, P%ROW, P%COL, R%VAL, R%ROW, R%COL, L%NCE, L%AMG%NCCE)
!      ENDDO MESHES_LOOP_CLASSICAL2
!
!
! !> ------------------------------------------------------------------------------------------
! !> DEFAULT
! !> ------------------------------------------------------------------------------------------
!   CASE DEFAULT
!
!      MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!
!         L => SCARC(NM)%LEVEL(NL)
!         R => SCARC(NM)%SYSTEM(NL)%R
!         P => SCARC(NM)%SYSTEM(NL)%P
!
!         IC2 = 1
!         DO ICP = 1, L%AMG%NCCE
!            BFIRST = .TRUE.
!            DO IC = 1, L%NCE
!
!               ROW_LOOP: DO ICP2 = P%ROW(IC),P%ROW(IC+1)-1
!                  IF (P%COL(ICP2) == ICP) THEN
!                     R%VAL(IC2) = P%VAL(ICP2)
!                     IF (BFIRST) THEN
!                        R%ROW(ICP) = IC2
!                     ENDIF
!                     R%COL(IC2) = IC
!                     IC2 = IC2 + 1
!                     BFIRST = .FALSE.
!                     EXIT ROW_LOOP
!                  ENDIF
!               ENDDO ROW_LOOP
!
!            ENDDO
!
!         ENDDO
!         R%ROW(L%AMG%NCCE+1)=IC2
!
!         OTHER_MESHES_LOOP: DO NOM = 1, NMESHES
!
!            !IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
!            OS => SCARC(NM)%OSCARC(NOM)
!            OX => OS%EXCHANGE
!            IF (OX%NICMAX_S==0 .AND. OX%NICMAX_R==0) CYCLE OTHER_MESHES_LOOP
!
!            OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
!            OP => SCARC(NM)%OSCARC(NOM)%SYSTEM(NL)%P
!            OR => SCARC(NM)%OSCARC(NOM)%SYSTEM(NL)%R
!
!            IC2 = 1
!            ICC = 0
!            OTHER_GHOSTCELL_LOOP: DO IG0 = 1, OL%NCG
!
!               ICP = OL%AMG%CTYPE(IG0)
!               IF (ICP < NSCARC_CELL_TYPE_COARSE) CYCLE OTHER_GHOSTCELL_LOOP
!               ICC = ICC + 1
!
!               BFIRST = .TRUE.
!               DO IG = 1, OL%NCG
!
!                  OTHER_ROW_LOOP: DO ICP2 = OP%ROW(IG),OP%ROW(IG+1)-1
!                     IF (OP%COL(ICP2) == ICC) THEN
!                        OR%VAL(IC2) = OP%VAL(ICP2)
!                        IF (BFIRST) THEN
!                           OR%ROW(ICC) = IC2
!                        ENDIF
!                        OR%COL(IC2) = IG
!                        IC2 = IC2 + 1
!                        BFIRST = .FALSE.
!                        EXIT OTHER_ROW_LOOP
!                     ENDIF
!                  ENDDO OTHER_ROW_LOOP
!
!               ENDDO
!
!            ENDDO OTHER_GHOSTCELL_LOOP
!            OR%ROW(ICC+1)=IC2
!
!         ENDDO OTHER_MESHES_LOOP
!
!      ENDDO MESHES_LOOP
!
!END SELECT
!CALL SCARC_LEAVE_ROUTINE()
!
!END SUBROUTINE SCARC_SETUP_RESTRICTION


!> ------------------------------------------------------------------------------------------------
!> Compute coarser matrix A_coarse by use of interpolation matrices:
!>
!>      A_coarse = I_fine^coarse * A_fine * I_coarse_fine
!>
!> Note the different storage techniques (compact storage technique for
!> transfer matrices and coarse matrix,  for finest matrix)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_TRANSFER_MATRIX(NL)
INTEGER , INTENT(IN) :: NL
INTEGER :: NM, IC, ICO2, ICP1, ICP2, IP, ICOL, IDIAG
REAL (EB), ALLOCATABLE, DIMENSION(:) :: VAL
!REAL(EB):: MAUX1(50,50)=0.0_EB!, MAUX2(30,30)=0.0_EB
REAL(EB) :: AUX1, AUX2, PW !, MATRIX(16,16)
TYPE (SCARC_LEVEL_TYPE) , POINTER :: LF, LC
TYPE (SCARC_MATRIX_TYPE), POINTER :: AF, AC

WRITE(*,*) 'ACHTUNG, MUSS ALTERNATIV ZUR FOLGEROUTINE GENUTZT WERDEN'
CALL SCARC_ENTER_ROUTINE('SCARC_TRANSFER_MATRIX')

SELECT CASE (TYPE_COARSENING)

 !> ---------------------------------------------------------------------------------
 !> Default case
 !> ---------------------------------------------------------------------------------
   CASE DEFAULT

      MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         LF => SCARC(NM)%LEVEL(NL)             !> pointer to fine level
         LC => SCARC(NM)%LEVEL(NL+1)           !> pointer to coarse level

         AF => SCARC(NM)%SYSTEM(NL)%A
         AC => SCARC(NM)%SYSTEM(NL+1)%A

         IF (NMESHES==1) THEN
            CALL SCARC_ALLOCATE_REAL1(VAL, 1, LF%NC , NSCARC_INIT_ZERO, 'VAL')
         ELSE
            CALL SCARC_ALLOCATE_REAL1(VAL, 1, LF%NCE, NSCARC_INIT_ZERO, 'VAL')
         ENDIF

         IP = 1
         ICP1_LOOP: DO ICP1 = 1, LF%AMG%NCC


            !> ---------------------------------------------------------------------------------
            !> First: loop over the internal coarse cells
            !> ---------------------------------------------------------------------------------
            ICP2_LOOP1: DO ICP2 = 1, LF%AMG%NCCE

               AUX2 = 0.0_EB

               DO IC = 1, LF%NC

                  IDIAG = AF%ROW(IC)
                  PW    = P_WEIGHT(IC, ICP2, NM, NL)
                  AUX1  = AF%VAL(IDIAG) * PW

                  COUPLINGS_LOOP: DO ICOL = AF%ROW(IC)+1, AF%ROW(IC+1)-1
                     ICO2 = AF%COL(ICOL)
                     PW = P_WEIGHT(ICO2, ICP2, NM, NL)
                     AUX1 = AUX1 + AF%VAL(ICOL) * PW
                  ENDDO COUPLINGS_LOOP

                  AUX2 = AUX2 + P_WEIGHT(IC, ICP1, NM, NL) * AUX1

               ENDDO
               VAL(ICP2) = AUX2

            ENDDO ICP2_LOOP1

            !> analyze new matrix line and store it corresponding to compact storage technique:
            !> (diagonal entry first)
            AC%VAL(IP)       = VAL(ICP1)
            AC%ROW(ICP1) = IP
            AC%COL(IP)   = ICP1

            IP  = IP + 1
            ICP2_LOOP2: DO ICP2 = 1, LF%AMG%NCCE
            !DO ICP2 = 1, LF%AMG%NCC
               IF (ICP2 /= ICP1 .AND. ABS(VAL(ICP2)) >= 1.0E-12_EB) THEN
                  AC%VAL(IP) = VAL(ICP2)
                  AC%COL(IP) = ICP2
                  IP  = IP + 1
               ENDIF
            ENDDO ICP2_LOOP2

            VAL = 0.0_EB
         ENDDO ICP1_LOOP

         AC%ROW(LC%NC+1) = IP
         AC%NAV = IP - 1

         DEALLOCATE(VAL)

      ENDDO MESHES_LOOP

END SELECT

IF (NMESHES > 1) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_MATRIX_SUBDIAG, NL)

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_TRANSFER_MATRIX


!> ------------------------------------------------------------------------------------------------
!> Determine if corresponding coarse point contributes a non-zero interpolation weight or not
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION P_WEIGHT(IC, ICP, NM, NL)
INTEGER, INTENT(IN):: IC, ICP, NM, NL
INTEGER :: ICOL
REAL(EB) :: VAL
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: P

L => SCARC(NM)%LEVEL(NL)
P => SCARC(NM)%SYSTEM(NL)%P

VAL = 0.0_EB
!IF (IC <= L%NC) THEN
   P_WEIGHT_LOOP: DO ICOL = P%ROW(IC), P%ROW(IC+1)-1
      IF (ICOL > 0) THEN
        IF (P%COL(ICOL) /= ICP) CYCLE P_WEIGHT_LOOP
      ENDIF
      VAL = P%VAL(ICOL)
   ENDDO P_WEIGHT_LOOP
!ENDIF

P_WEIGHT = VAL
RETURN
END FUNCTION P_WEIGHT

!> ------------------------------------------------------------------------------------------------
!> Determine if corresponding coarse point contributes a non-zero interpolation weight or not
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION R_WEIGHT(IC, ICP, NM, NL)
INTEGER, INTENT(IN):: IC, ICP, NM, NL
INTEGER :: ICOL
REAL(EB) :: VAL
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: R

L => SCARC(NM)%LEVEL(NL)
R => SCARC(NM)%SYSTEM(NL)%R

VAL = 0.0_EB
IF (ICP <= L%AMG%NCC) THEN
   R_WEIGHT_LOOP: DO ICOL = R%ROW(ICP), R%ROW(ICP+1)-1
      IF (ICOL > 0) THEN
        IF (R%COL(ICOL) /= IC) CYCLE R_WEIGHT_LOOP
      ENDIF
      VAL = R%VAL(ICOL)
   ENDDO R_WEIGHT_LOOP
ELSE
  WRITE(*,*)
ENDIF

R_WEIGHT = VAL
RETURN
END FUNCTION R_WEIGHT


!> -----------------------------------------------------------------------------------------------
!> Interface for the call of different ScaRC-solvers
!> -----------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SOLVER
REAL (EB) :: TNOW

CALL SCARC_ENTER_ROUTINE ('SCARC_SOLVER')
TNOW = CURRENT_TIME()

ITE_PRES = ITE_PRES + 1

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !> ---------------- Krylov method (CG/BICG) --------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV: SELECT CASE (TYPE_KRYLOV)
         CASE (NSCARC_KRYLOV_CG)
            CALL SCARC_METHOD_CG  (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
         CASE (NSCARC_KRYLOV_BICG)
            CALL SCARC_METHOD_BICG(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
      END SELECT SELECT_KRYLOV

   !> ---------------- Multigrid method ---------------------------------------
   CASE (NSCARC_METHOD_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)

   !> ---------------- MKL method ---------------------------------------------
#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LUDECOMP)

      SELECT_MKL: SELECT CASE (TYPE_LUDECOMP)
         CASE (NSCARC_MKL_GLOBAL) 
            CALL SCARC_METHOD_CLUSTER(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
         CASE (NSCARC_MKL_LOCAL) 
            CALL SCARC_METHOD_PARDISO(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
      END SELECT SELECT_MKL
#endif

END SELECT SELECT_METHOD

IF (STOP_STATUS==SETUP_STOP) RETURN    

CALL SCARC_LEAVE_ROUTINE()
T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
TSTEP(MYID+1)%SOLVER=MAX(TSTEP(MYID+1)%SOLVER,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%SOLVER =TSUM(MYID+1)%SOLVER+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SOLVER

!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for compact storage technique
!> ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_SCOPE(NTYPE, NM, NL)
TYPE (SCARC_SCOPE_TYPE), POINTER :: POINT_TO_SCOPE
INTEGER, INTENT(IN):: NTYPE, NM, NL

SELECT CASE (NTYPE)
   CASE (NSCARC_SOLVER_MAIN)
      POINT_TO_SCOPE => SCARC(NM)%SCOPE(NSCARC_SCOPE_ONE, NL)
   CASE (NSCARC_SOLVER_PRECON)
      POINT_TO_SCOPE => SCARC(NM)%SCOPE(NSCARC_SCOPE_TWO, NL)
   CASE (NSCARC_SOLVER_SMOOTH)
      IF (BCGGMG) THEN
         POINT_TO_SCOPE => SCARC(NM)%SCOPE(NSCARC_SCOPE_TWO, NL)
      ELSE
         POINT_TO_SCOPE => SCARC(NM)%SCOPE(NSCARC_SCOPE_ONE, NL)
      ENDIF
   CASE (NSCARC_SOLVER_COARSE)
      IF (BCGGMG) THEN
         POINT_TO_SCOPE => SCARC(NM)%SCOPE(NSCARC_SCOPE_TWO, NL)
      ELSE
         POINT_TO_SCOPE => SCARC(NM)%SCOPE(NSCARC_SCOPE_ONE, NL)
      ENDIF
END SELECT
RETURN

END FUNCTION POINT_TO_SCOPE

!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for compact storage technique
!> ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_VECTOR(NVECTOR, NM, NL)
REAL(EB), POINTER, DIMENSION(:) :: POINT_TO_VECTOR
INTEGER, INTENT(IN):: NVECTOR, NM, NL
TYPE (SCARC_TYPE), POINTER :: S

S => SCARC(NM)
SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_ONE_X)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_ONE, NL)%X
   CASE (NSCARC_VECTOR_ONE_F)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_ONE, NL)%F
   CASE (NSCARC_VECTOR_ONE_D)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_ONE, NL)%D
   CASE (NSCARC_VECTOR_ONE_G)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_ONE, NL)%G
   CASE (NSCARC_VECTOR_ONE_W)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_ONE, NL)%W
   CASE (NSCARC_VECTOR_ONE_Y)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_ONE, NL)%Y
   CASE (NSCARC_VECTOR_ONE_Z)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_ONE, NL)%Z
   CASE (NSCARC_VECTOR_ONE_E)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_ONE, NL)%E
   CASE (NSCARC_VECTOR_TWO_X)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_TWO, NL)%X
   CASE (NSCARC_VECTOR_TWO_F)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_TWO, NL)%F
   CASE (NSCARC_VECTOR_TWO_D)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_TWO, NL)%D
   CASE (NSCARC_VECTOR_TWO_G)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_TWO, NL)%G
   CASE (NSCARC_VECTOR_TWO_W)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_TWO, NL)%W
   CASE (NSCARC_VECTOR_TWO_Y)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_TWO, NL)%Y
   CASE (NSCARC_VECTOR_TWO_Z)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_TWO, NL)%Z
   CASE (NSCARC_VECTOR_TWO_E)
      POINT_TO_VECTOR => S%SCOPE(NSCARC_SCOPE_TWO, NL)%E
END SELECT

RETURN
END FUNCTION POINT_TO_VECTOR

#ifdef WITH_MKL_FB
!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for compact storage technique
!> ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_VECTOR_FB(NVECTOR, NM, NL)
REAL(FB), POINTER, DIMENSION(:) :: POINT_TO_VECTOR_FB
INTEGER, INTENT(IN):: NVECTOR, NM, NL
TYPE (SCARC_TYPE), POINTER :: S

S => SCARC(NM)
SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_ONE_X)
      POINT_TO_VECTOR_FB => S%SCOPE(NSCARC_SCOPE_ONE, NL)%X_FB
   CASE (NSCARC_VECTOR_ONE_F)
      POINT_TO_VECTOR_FB => S%SCOPE(NSCARC_SCOPE_ONE, NL)%F_FB
   CASE (NSCARC_VECTOR_ONE_D)
      POINT_TO_VECTOR_FB => S%SCOPE(NSCARC_SCOPE_ONE, NL)%D_FB
   CASE (NSCARC_VECTOR_ONE_G)
      POINT_TO_VECTOR_FB => S%SCOPE(NSCARC_SCOPE_ONE, NL)%G_FB
   CASE (NSCARC_VECTOR_ONE_W)
      POINT_TO_VECTOR_FB => S%SCOPE(NSCARC_SCOPE_ONE, NL)%W_FB
END SELECT

RETURN
END FUNCTION POINT_TO_VECTOR_FB
#endif

!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for banded storage technique
!> ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_HVECTOR(NVECTOR, NM)
REAL(EB), POINTER, DIMENSION(:,:,:) :: POINT_TO_HVECTOR
INTEGER, INTENT(IN):: NVECTOR, NM

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_H)
      POINT_TO_HVECTOR => MESHES(NM)%H
   CASE (NSCARC_VECTOR_HS)
      POINT_TO_HVECTOR => MESHES(NM)%HS
END SELECT

RETURN
END FUNCTION POINT_TO_HVECTOR

!> ------------------------------------------------------------------------------------------------
!> Compute global matrix-vector product (including data exchange along internal boundaries)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATVEC_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER, INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), POINTER, DIMENSION(:) :: V1, V2
REAL(EB) :: TNOW
INTEGER , POINTER :: NC, NCE
INTEGER :: NM, IC, JC, ICOL
TYPE (SCARC_MATRIX_TYPE), POINTER:: A

TNOW = CURRENT_TIME()
CALL SCARC_ENTER_ROUTINE('SCARC_MATVEC_PRODUCT')

!>    
!> Exchange internal boundary values of vector1 such that the ghost values contain the corresponding
!> overlapped values of adjacent neighbor
!>
TYPE_VECTOR = NVECTOR1
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

!> 
!> Perform global matrix-vector product:
!> Note: - matrix already contains subdiagonal values from neighbor along internal boundaries
!>       - if vector1 contains neighboring values, then correct values of global matvec are achieved
!> 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   V1 => POINT_TO_VECTOR (NVECTOR1, NM, NL)
   V2 => POINT_TO_VECTOR (NVECTOR2, NM, NL)

   NC  => SCARC(NM)%LEVEL(NL)%NCS                                     !> number of cells
   NCE => SCARC(NM)%LEVEL(NL)%NCE                                     !> number of extended cells

   A => SCARC(NM)%SYSTEM(NL)%A                                        !> system matrix

   DO IC = 1, NC

      ICOL = A%ROW(IC)                                                !> diagonal entry
      JC   = A%COL(ICOL)
      V2(IC) = A%VAL(ICOL)* V1(JC)

      DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1                            !> subdiagonal entries
         JC = A%COL(ICOL)
         V2(IC) =  V2(IC) + A%VAL(ICOL)* V1(JC)
!         WRITE(MSG%LU_DEBUG,'(A,3i4,3e14.6)') 'IC, ICOL, JC, V1(IC), V2(IC), A(ICOL):', &
!                          IC, ICOL, JC, V1(IC), V2(IC), A%VAL(ICOL)

      ENDDO
   ENDDO
ENDDO

CALL SCARC_LEAVE_ROUTINE()
TSTEP(MYID+1)%MATVEC=MAX(TSTEP(MYID+1)%MATVEC,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%MATVEC =TSUM(MYID+1)%MATVEC+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_MATVEC_PRODUCT


!> ------------------------------------------------------------------------------------------------
!> Compute global scalarproductt (including global data exchange)
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_SCALAR_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER, INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), DIMENSION(:)    , POINTER ::  V1, V2
REAL(EB) :: TNOW
INTEGER , POINTER :: NC
INTEGER  :: NM, NL0
#if defined(WITH_MKL)
REAL(EB) :: DDOT
EXTERNAL :: DDOT
#else
INTEGER :: IC
#endif

TNOW = CURRENT_TIME()
LOCAL_REAL = 0.0_EB

!> Compute local scalar products on single meshes und process group
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   V1 => POINT_TO_VECTOR (NVECTOR1, NM, NL)
   V2 => POINT_TO_VECTOR (NVECTOR2, NM, NL)

   NC => SCARC(NM)%LEVEL(NL)%NCS

#if defined(WITH_MKL)
   LOCAL_REAL(NM) = DDOT(NC, V1, 1, V2, 1)
#else   
   LOCAL_REAL(NM) = 0.0_EB
   DO IC = 1, NC
      LOCAL_REAL(NM) = LOCAL_REAL(NM) + V1(IC) * V2(IC)
   ENDDO
#endif

ENDDO

!> Compute global scalar product as sum of local scalar products
NL0 = NL
GLOBAL_REAL = 0.0_EB

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,LOCAL_REAL,COUNTS,DISPLS,&
                       MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
ENDIF
GLOBAL_REAL = SUM(LOCAL_REAL(1:NMESHES))

SCARC_SCALAR_PRODUCT = GLOBAL_REAL

TSTEP(MYID+1)%SCALPROD=MAX(TSTEP(MYID+1)%SCALPROD,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%SCALPROD =MAX(TSUM(MYID+1)%SCALPROD,CURRENT_TIME()-TNOW)

RETURN
END FUNCTION SCARC_SCALAR_PRODUCT


!> ------------------------------------------------------------------------------------------------
!> Compute global L2-norm (including global data exchange)
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_L2NORM(NVECTOR1, NL)
INTEGER, INTENT(IN):: NVECTOR1, NL
REAL(EB) :: TNOW
TNOW = CURRENT_TIME()

GLOBAL_REAL = SCARC_SCALAR_PRODUCT(NVECTOR1, NVECTOR1, NL)
GLOBAL_REAL = SQRT (GLOBAL_REAL)

SCARC_L2NORM = GLOBAL_REAL

TSTEP(MYID+1)%L2NORM=MAX(TSTEP(MYID+1)%L2NORM,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%L2NORM =TSUM(MYID+1)%L2NORM+CURRENT_TIME()-TNOW

RETURN
END FUNCTION SCARC_L2NORM


!> ------------------------------------------------------------------------------------------------
!> Compute linear combination of two vectors for banded storage technique
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_SUM(NVECTOR1, NVECTOR2, SCAL1, SCAL2, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), INTENT(IN):: SCAL1, SCAL2
REAL(EB), DIMENSION(:)    , POINTER ::  V1, V2
INTEGER  :: NM
#if defined(WITH_MKL)
INTEGER, POINTER  :: NC
EXTERNAL :: DAXPBY
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   V1 => POINT_TO_VECTOR(NVECTOR1, NM, NL)
   V2 => POINT_TO_VECTOR(NVECTOR2, NM, NL)
#if defined(WITH_MKL)
   NC => SCARC(NM)%LEVEL(NL)%NCS
   CALL DAXPBY(NC, SCAL1, V1, 1, SCAL2, V2, 1)
#else
   V2 = SCAL1 * V1 + SCAL2 * V2
#endif
ENDDO

END SUBROUTINE SCARC_VECTOR_SUM


!> ------------------------------------------------------------------------------------------------
!> Copy one integer array to another (possible scaled)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_COPY_INT(X, Y, SCAL, NLEN)
INTEGER,  INTENT(IN) , DIMENSION(:) :: X
INTEGER,  INTENT(OUT), DIMENSION(:) :: Y
INTEGER,  INTENT(IN) :: SCAL, NLEN
INTEGER :: ILEN

IF (SCAL == 1.0_EB) THEN
   DO ILEN = 1, NLEN
      Y(ILEN) = X(ILEN)
   ENDDO
ELSE IF (SCAL == 0.0_EB) THEN
   DO ILEN = 1, NLEN
      Y(ILEN) = 0
   ENDDO
ELSE
   DO ILEN = 1, NLEN
      Y(ILEN) = SCAL * X(ILEN)
   ENDDO
ENDIF

END SUBROUTINE SCARC_COPY_INT


!> ------------------------------------------------------------------------------------------------
!> Copy one real array to another (possible scaled)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_COPY_REAL(X, Y, SCAL, NLEN)
REAL(EB), INTENT(IN) , DIMENSION(:) :: X
REAL(EB), INTENT(OUT), DIMENSION(:) :: Y
REAL(EB), INTENT(IN) :: SCAL
INTEGER,  INTENT(IN) :: NLEN
INTEGER :: ILEN

IF (SCAL == 1.0_EB) THEN
   DO ILEN = 1, NLEN
      Y(ILEN) = X(ILEN)
   ENDDO
ELSE IF (SCAL == 0.0_EB) THEN
   DO ILEN = 1, NLEN
      Y(ILEN) = 0.0_EB
   ENDDO
ELSE
   DO ILEN = 1, NLEN
      Y(ILEN) = SCAL * X(ILEN)
   ENDDO
ENDIF

END SUBROUTINE SCARC_COPY_REAL


!> ------------------------------------------------------------------------------------------------
!> Define vector2 to be a scaled copy of vector 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_COPY(NVECTOR1, NVECTOR2, SCAL1, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), INTENT(IN):: SCAL1
REAL(EB), DIMENSION(:), POINTER ::  V1, V2
INTEGER  :: NM

#if defined(WITH_MKL)
INTEGER , POINTER :: NC
EXTERNAL :: DCOPY, DSCAL
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   V1 => POINT_TO_VECTOR(NVECTOR1, NM, NL)
   V2 => POINT_TO_VECTOR(NVECTOR2, NM, NL)

#if defined(WITH_MKL)
   NC => SCARC(NM)%LEVEL(NL)%NCS
   CALL DCOPY(NC, V1, 1, V2, 1)
   CALL DSCAL(NC, SCAL1, V2, 1)
#else
   V2 = SCAL1 * V1
#endif

ENDDO

END SUBROUTINE SCARC_VECTOR_COPY


!> ------------------------------------------------------------------------------------------------
!> Clear vector
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_CLEAR(NVECTOR, NL)
INTEGER , INTENT(IN):: NVECTOR, NL
REAL(EB), DIMENSION(:), POINTER ::  VC
INTEGER  :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   VC => POINT_TO_VECTOR(NVECTOR, NM, NL)
   VC =  0.0_EB
ENDDO

END SUBROUTINE SCARC_VECTOR_CLEAR


!> ------------------------------------------------------------------------------------------------
!> Preset vector with specified value
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_INIT (V, VAL, NL)
INTEGER, INTENT(IN):: V, NL
REAL (EB), INTENT(IN) :: VAL
REAL (EB), POINTER, DIMENSION(:) :: VC
INTEGER :: IC, NM, I, J, K
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L => SCARC(NM)%LEVEL(NL)
   VC => POINT_TO_VECTOR (V, NM, NL)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%CELL%STATE(I,J,K) /= NSCARC_DISCRET_GASPHASE) CYCLE
            IC = L%CELL%DOF(I,J,K)
            VC(IC) = VAL
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_VECTOR_INIT

!> ------------------------------------------------------------------------------------------------
!> Perform preconditioning
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_BLOCK_SOLVER (NV1, NV2, NSTACK, NPARENT, NL)
USE POIS, ONLY: H2CZSS, H3CZSS
INTEGER, INTENT(IN):: NV1, NV2, NSTACK, NPARENT, NL
INTEGER , POINTER:: NC => NULL()
INTEGER  :: NM, IC, IW, I, J, K, ICOL !, IOR0, IROW, JCOL
INTEGER  :: IXW, IYW, IZW, ICW, ICE
REAL(EB) :: AUX, OMEGA_SSOR=1.5_EB, VAL
REAL(EB), DIMENSION(:), POINTER ::  V1, V2
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_WALL_TYPE)  , POINTER :: WC
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A
TYPE (SCARC_FFT_TYPE)   , POINTER :: F
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: MKL
#ifdef WITH_MKL_FB
REAL(FB), DIMENSION(:), POINTER ::  V1_FB, V2_FB
TYPE (SCARC_MATRIX_FB_TYPE), POINTER :: ASYM_FB
#else
TYPE (SCARC_MATRIX_TYPE), POINTER :: ASYM
#endif
#endif

REAL (EB) :: TNOW

TNOW = CURRENT_TIME()
CALL SCARC_ENTER_ROUTINE('SCARC_BLOCK_SOLVER')

SELECT CASE (STACK(NSTACK-1)%SOLVER%TYPES%TYPE_RELAX)

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise Jacobi
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_JACOBI)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         NC => SCARC(NM)%LEVEL(NL)%NCS
         A  => SCARC(NM)%SYSTEM(NL)%A
         V2 => POINT_TO_VECTOR(NV2, NM, NL)

         DO IC = 1, NC
            V2(IC) = V2(IC) / A%VAL(A%ROW(IC))
         ENDDO

      ENDDO

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise SSOR
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_SSOR)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         NC => SCARC(NM)%LEVEL(NL)%NCS
         A  => SCARC(NM)%SYSTEM(NL)%A
         V2 => POINT_TO_VECTOR(NV2, NM, NL)

         !> forward SOR step
         FORWARD_CELL_LOOP: DO IC = 1, NC
            AUX = 0.0_EB
            LOWER_DIAG_LOOP: DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
               IF (A%COL(ICOL) >= IC) CYCLE LOWER_DIAG_LOOP
               IF (A%COL(ICOL) <= NC) AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))
            ENDDO LOWER_DIAG_LOOP
            V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / A%VAL(A%ROW(IC))
         ENDDO FORWARD_CELL_LOOP

         !> backward SOR step
         BACKWARD_CELL_LOOP: DO IC = NC-1, 1, -1
            AUX = 0.0_EB
            UPPER_DIAG_LOOP: DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
               IF (A%COL(ICOL) <= IC) CYCLE UPPER_DIAG_LOOP
               IF (A%COL(ICOL) <= NC) AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))
            ENDDO UPPER_DIAG_LOOP
            V2(IC) = V2(IC) - AUX * OMEGA_SSOR / A%VAL(A%ROW(IC))
         ENDDO BACKWARD_CELL_LOOP

      ENDDO

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise Geometric Multigrid 
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID (NSTACK, NPARENT, NLEVEL_MIN)


   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise FFT based on Crayfishpak
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_FFT)

      TYPE_VECTOR = NV1
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
        
         M  => MESHES(NM)

         V1  => POINT_TO_VECTOR(NV1, NM, NL)
         V2  => POINT_TO_VECTOR(NV2, NM, NL)

         L => SCARC(NM)%LEVEL(NL)
         F => SCARC(NM)%FFT(NL)

         DO K = 1, L%NZ
            DO J = 1, L%NY
               DO I = 1, L%NX
                  IC = (K-1) * L%NX * L%NY + (J-1) * L%NX + I
                  F%PRHS(I, J, K) = V1(IC)
               ENDDO
            ENDDO
         ENDDO

         DO IW = 1, L%N_WALL_CELLS_EXT

            WC => L%WALL(IW)

            IXW = L%WALL(IW)%IXW
            IYW = L%WALL(IW)%IYW
            IZW = L%WALL(IW)%IZW

            ICW = L%CELL%DOF(IXW,IYW,IZW)

            !IF (L%WALL(IW)%NOM /= 0) THEN
            !   ICE = L%WALL(IW)%ICE(1)
            !   VAL = 0.5_EB*(V1(ICW) + V1(ICE))
            !ELSE
            !   VAL = 0.0_EB
            !ENDIF

            !> Use zero BC's for the moment
            VAL = 0.0_EB

            SELECT CASE(WC%IOR)
               CASE( 1)
                  F%BXS(IYW,IZW) = VAL
               CASE(-1)
                  F%BXF(IYW,IZW) = VAL
               CASE( 2)
                  F%BYS(IXW,IZW) = VAL
               CASE(-2)
                  F%BYF(IXW,IZW) = VAL
               CASE( 3)
                  F%BZS(IXW,IYW) = VAL
               CASE(-3)
                  F%BZF(IXW,IYW) = VAL
            END SELECT

         ENDDO

         IF (TWO_D) THEN
            CALL H2CZSS (F%BXS, F%BXF, F%BZS, F%BZF, L%NX+1, &
                         F%PRHS, F%POIS_PTB, F%SAVE1, F%WORK, F%HX)
         ELSE
            CALL H3CZSS (F%BXS, F%BXF, F%BYS, F%BYF, F%BZS, F%BZF, L%NX+1, L%NY+1, &
                         F%PRHS, F%POIS_PTB, F%SAVE1, F%WORK, F%HX)
         ENDIF

         DO K = 1, L%NZ
            DO J = 1, L%NY
               DO I = 1, L%NX
                  IC = (K-1) * L%NX * L%NY + (J-1) * L%NX + I
                  V2(IC) = F%PRHS(I, J, K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise FFT based on Crayfishpak
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_FFT_OVERLAP)

      TYPE_VECTOR = NV1
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
        
         M  => MESHES(NM)

         V1  => POINT_TO_VECTOR(NV1, NM, NL)
         V2  => POINT_TO_VECTOR(NV2, NM, NL)

         L => SCARC(NM)%LEVEL(NL)
         F => SCARC(NM)%FFT(NL)

         DO K = 1, L%NZ
            DO J = 1, L%NY
               DO I = 1, L%NX
                  IC = (K-1) * L%NX * L%NY + (J-1) * L%NX + I
                  F%PRHS(I+F%IS0, J+F%JS0, K+F%KS0) = V1(IC)
               ENDDO
            ENDDO
         ENDDO

         DO IW = 1, L%N_WALL_CELLS_EXT

            WC => L%WALL(IW)

            IXW = WC%IXW
            IYW = WC%IYW
            IZW = WC%IZW

            IF (WC%NOM /= 0) THEN
               ICE = L%WALL(IW)%ICE(1)
               SELECT CASE(WC%IOR)
                  CASE ( 1)
                     F%PRHS(      1  , IYW+F%JS0, IZW+F%KS0) = V1(ICE)
                  CASE (-1)
                     F%PRHS(F%IBAR0  , IYW+F%JS0, IZW+F%KS0) = V1(ICE)
                  CASE ( 2)
                     F%PRHS(IXW+F%IS0,         1, IZW+F%KS0) = V1(ICE)
                  CASE (-2)
                     F%PRHS(IXW+F%IS0,   F%JBAR0, IZW+F%KS0) = V1(ICE)
                  CASE ( 3)
                     F%PRHS(IXW+F%IS0, IYW+F%JS0,         1) = V1(ICE)
                  CASE (-3)
                     F%PRHS(IXW+F%IS0, IYW+F%JS0,   F%KBAR0) = V1(ICE)
               END SELECT
            ENDIF
         ENDDO

         !> Proof of concept for simple test case - defining diagonal entries as mean values
         !IF (MYID == 0) THEN
         !   F%PRHS(F%IBAR0,1,F%KBAR0) = 0.5_EB*(F%PRHS(F%IBAR0,1,F%KBAR0-1) + F%PRHS(F%IBAR0-1,1,F%KBAR0))
         !ELSE IF (MYID == 1) THEN
         !   F%PRHS(1,1,F%KBAR0) = 0.5_EB*(F%PRHS(1,1,F%KBAR0-1) + F%PRHS(2,1,F%KBAR0))
         !ELSE IF (MYID == 2) THEN
         !   F%PRHS(F%IBAR0,1,1) = 0.5_EB*(F%PRHS(F%IBAR0-1,1,1) + F%PRHS(F%IBAR0,1,2))
         !ELSE IF (MYID == 3) THEN
         !   F%PRHS(1,1,1) = 0.5_EB*(F%PRHS(1,1,2) + F%PRHS(2,1,1))
         !ENDIF

         F%BXS = 0.0_EB
         F%BXF = 0.0_EB
         F%BYS = 0.0_EB
         F%BYF = 0.0_EB
         F%BZS = 0.0_EB
         F%BZF = 0.0_EB

         IF (TWO_D) THEN
            CALL H2CZSS (F%BXS, F%BXF, F%BZS, F%BZF, F%IBAR0+1, &
                         F%PRHS, F%POIS_PTB, F%SAVE1, F%WORK, F%HX)
         ELSE
            CALL H3CZSS (F%BXS, F%BXF, F%BYS, F%BYF, F%BZS, F%BZF, F%IBAR0+1, F%JBAR0+1, &
                         F%PRHS, F%POIS_PTB, F%SAVE1, F%WORK, F%HX)
         ENDIF

         DO K = 1, L%NZ
            DO J = 1, L%NY
               DO I = 1, L%NX
                  IC = (K-1) * L%NX * L%NY + (J-1) * L%NX + I
                  V2(IC) = F%PRHS(I+F%IS0, J+F%JS0, K+F%KS0)
               ENDDO
            ENDDO
         ENDDO

         DO IW = 1, L%N_WALL_CELLS_EXT

            WC => L%WALL(IW)

            IXW = WC%IXW
            IYW = WC%IYW
            IZW = WC%IZW

            !> don't do it temporarily
            !IF (WC%NOM /= 0) THEN
            IF (WC%NOM == 123456) THEN
               SELECT CASE(WC%IOR)
                  CASE ( 1)
                     V2(ICE) = F%PRHS(    1, IYW, IZW) 
                  CASE (-1)
                     V2(ICE) = F%PRHS(F%IBAR0, IYW, IZW) 
                  CASE ( 2)
                     V2(ICE) = F%PRHS(IXW,     1, IZW) 
                  CASE (-2)
                     V2(ICE) = F%PRHS(IXW, F%JBAR0, IZW) 
                  CASE ( 3)
                     V2(ICE) = F%PRHS(IXW, IYW,     1) 
                  CASE (-3)
                     V2(ICE) = F%PRHS(IXW, IYW, F%KBAR0) 
               END SELECT
            ENDIF
         ENDDO
      ENDDO

#ifdef WITH_MKL
   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by Cluster Sparse Solver from MKL
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_CLUSTER)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      
         L   => SCARC(NM)%LEVEL(NL)
         MKL => SCARC(NM)%MKL(NL)
         NC  => SCARC(NM)%LEVEL(NL)%NCS

         MKL%PHASE  = 33                            ! only solving

         V1 => POINT_TO_VECTOR (NV1, NM, NL)
         V2 => POINT_TO_VECTOR (NV2, NM, NL)

#ifdef WITH_MKL_FB

         ASYM_FB => SCARC(NM)%SYSTEM(NL)%ASYM_FB

         V1_FB => POINT_TO_VECTOR_FB (NV1, NM, NL)
         V2_FB => POINT_TO_VECTOR_FB (NV2, NM, NL)

         V1_FB(1:NC) = REAL(V1(1:NC), FB)
         V2_FB(1:NC) = 0.0_FB

         CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%CELL%NC_GLOBAL, &
                                      ASYM_FB%VAL, ASYM_FB%ROW, ASYM_FB%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)

         V2(1:NC) = REAL(V2_FB(1:NC), EB)

#else

         ASYM => SCARC(NM)%SYSTEM(NL)%ASYM
         CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%CELL%NC_GLOBAL, &
                                      ASYM%VAL, ASYM%ROW, ASYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)
#endif
         IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN('The following MKL-error was detected ', 'NONE', MKL%ERROR)
      
      ENDDO

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by Pardiso Solver from MKL
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_PARDISO)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      
         L   => SCARC(NM)%LEVEL(NL)
         MKL => SCARC(NM)%MKL(NL)
         NC  => SCARC(NM)%LEVEL(NL)%NCS

         MKL%PHASE  = 33                            ! only solving

         V1 => POINT_TO_VECTOR (NV1, NM, NL)
         V2 => POINT_TO_VECTOR (NV2, NM, NL)

#ifdef WITH_MKL_FB

         ASYM_FB => SCARC(NM)%SYSTEM(NL)%ASYM_FB

         V1_FB => POINT_TO_VECTOR_FB (NV1, NM, NL)
         V2_FB => POINT_TO_VECTOR_FB (NV2, NM, NL)

         V1_FB(1:NC) = REAL(V1(1:NC), FB)
         V2_FB(1:NC) = 0.0_FB

         CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, NC, ASYM_FB%VAL, ASYM_FB%ROW, ASYM_FB%COL, &
                        MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

         V2(1:NC) = REAL(V2_FB(1:NC), EB)
#else

         ASYM => SCARC(NM)%SYSTEM(NL)%ASYM

         V1 => POINT_TO_VECTOR (NV1, NM, NL)
         V2 => POINT_TO_VECTOR (NV2, NM, NL)

         CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, NC, ASYM%VAL, ASYM%ROW, ASYM%COL, &
                        MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1, V2, MKL%ERROR)

#endif
         IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN('The following MKL-error was detected ', 'NONE', MKL%ERROR)
      
      ENDDO

#endif

END SELECT

TSTEP(MYID+1)%PRECON=MAX(TSTEP(MYID+1)%PRECON,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%PRECON =TSUM(MYID+1)%PRECON+CURRENT_TIME()-TNOW

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_BLOCK_SOLVER


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CLUSTER(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NL
REAL (EB) :: TNOW
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MKL_TYPE)   , POINTER :: MKL
REAL(EB), POINTER, DIMENSION(:) :: V1, V2
#ifdef WITH_MKL_FB
REAL(FB), POINTER, DIMENSION(:) :: V2_FB, V1_FB
TYPE (SCARC_MATRIX_FB_TYPE), POINTER :: ASYM_FB
#else
TYPE (SCARC_MATRIX_TYPE), POINTER :: ASYM
#endif

CALL SCARC_ENTER_ROUTINE('SCARC_METHOD_CLUSTER')
TNOW = CURRENT_TIME()

NL = NLEVEL

CALL SCARC_SETUP_SOLVER(NSTACK, NPARENT)
CALL SCARC_SETUP_WORKSPACE(NSTACK, NLEVEL)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L    => SCARC(NM)%LEVEL(NL)
   MKL  => SCARC(NM)%MKL(NL)                      

   MKL%PHASE  = 33                                !> only solving

   V1 => POINT_TO_VECTOR (F, NM, NL)
   V2 => POINT_TO_VECTOR (X, NM, NL)

#ifdef WITH_MKL_FB

   ASYM_FB => SCARC(NM)%SYSTEM(NL)%ASYM_FB           !> point to symmetric system matrix - single precision

   V1_FB => POINT_TO_VECTOR_FB (F, NM, NL)
   V2_FB => POINT_TO_VECTOR_FB (X, NM, NL)

   V1_FB = REAL(V1, FB)
   V2_FB = 0.0_FB

   CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%CELL%NC_GLOBAL, &
                                ASYM_FB%VAL, ASYM_FB%ROW, ASYM_FB%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)
   V2 = REAL(V2_FB, EB)

#else

   ASYM => SCARC(NM)%SYSTEM(NL)%ASYM              !> point to symmetric system matrix - double precision

   V1 => POINT_TO_VECTOR (F, NM, NL)

   CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%CELL%NC_GLOBAL, &
                                ASYM%VAL, ASYM%ROW, ASYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)

#endif

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN('The following MKL-error was detected ', 'NONE', MKL%ERROR)

ENDDO MESHES_LOOP

TYPE_VECTOR = X
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NSTACK, NPARENT)
CALL SCARC_LEAVE_ROUTINE()

TSTEP(MYID+1)%CLUSTER=MAX(TSTEP(MYID+1)%CLUSTER,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%CLUSTER =TSUM(MYID+1)%CLUSTER+CURRENT_TIME()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_CLUSTER
#endif


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_PARDISO(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NL !, IP, IC
REAL (EB) :: TNOW
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MKL_TYPE)   , POINTER :: MKL
REAL(EB), POINTER, DIMENSION(:) :: V1, V2
#ifdef WITH_MKL_FB
REAL(FB), POINTER, DIMENSION(:) :: V2_FB, V1_FB
TYPE (SCARC_MATRIX_FB_TYPE), POINTER :: ASYM_FB
#else
TYPE (SCARC_MATRIX_TYPE), POINTER :: ASYM
#endif

CALL SCARC_ENTER_ROUTINE('SCARC_METHOD_PARDISO')
TNOW = CURRENT_TIME()

NL = NLEVEL

CALL SCARC_SETUP_SOLVER(NSTACK, NPARENT)
CALL SCARC_SETUP_WORKSPACE(NSTACK, NLEVEL)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L   => SCARC(NM)%LEVEL(NL)
   MKL => SCARC(NM)%MKL(NL)
   MKL%PHASE  = 33         ! only solving

   V1 => POINT_TO_VECTOR (F, NM, NL)
   V2 => POINT_TO_VECTOR (X, NM, NL)

#ifdef WITH_MKL_FB

   ASYM_FB => SCARC(NM)%SYSTEM(NL)%ASYM_FB

   V1_FB => POINT_TO_VECTOR_FB (F, NM, NL)
   V2_FB => POINT_TO_VECTOR_FB (X, NM, NL)

   V1_FB = REAL(V1, FB)
   V2_FB = 0.0_FB

   CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NCS, &
                  ASYM_FB%VAL, ASYM_FB%ROW, ASYM_FB%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                  MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

   V2 = REAL(V2_FB, EB)

#else

   ASYM => SCARC(NM)%SYSTEM(NL)%ASYM

   V1 => POINT_TO_VECTOR (F, NM, NL)
   V2 => POINT_TO_VECTOR (X, NM, NL)

   V2 = 0.0_EB

   CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NCS, &
                  ASYM%VAL, ASYM%ROW, ASYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                  MKL%MSGLVL, V1, V2, MKL%ERROR)

#endif

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN('The following MKL-error was detected ', 'NONE', MKL%ERROR)

ENDDO MESHES_LOOP

!CALL SCARC_VECTOR_SUM(X, E, 1.0_EB, -1.0_EB, NLEVEL)
!ERR = SCARC_L2NORM (E, NLEVEL)                                        
!CALL SCARC_DUMP_RESIDUAL(0, NLEVEL)

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NSTACK, NPARENT)
CALL SCARC_LEAVE_ROUTINE()

TSTEP(MYID+1)%PARDISO=MAX(TSTEP(MYID+1)%PARDISO,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%PARDISO =TSUM(MYID+1)%PARDISO+CURRENT_TIME()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_PARDISO
#endif

!> ------------------------------------------------------------------------------------------------
!> Increase corresponding iteration count (just for visualization of convergence behavior)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INCREASE_ITERATION_COUNTS(ITE)
INTEGER, INTENT(IN) :: ITE

SELECT CASE (TYPE_SOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      SELECT CASE (TYPE_METHOD)
         CASE (NSCARC_METHOD_KRYLOV)
            ITE_CG = ITE
         CASE (NSCARC_METHOD_MULTIGRID)
            ITE_MG = ITE
         CASE (NSCARC_METHOD_LUDECOMP)
            ITE_LU = ITE
      END SELECT
   CASE (NSCARC_SOLVER_SMOOTH)
      ITE_SMOOTH = ITE
   CASE (NSCARC_SOLVER_COARSE)
      ITE_COARSE = ITE
END SELECT
ITE_TOTAL = ITE_TOTAL + 1

END SUBROUTINE SCARC_INCREASE_ITERATION_COUNTS

!> ------------------------------------------------------------------------------------------------
!> Perform global CG-method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CG(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER   :: NSTATE, NS, NL, NP
REAL (EB) :: SIGMA0, SIGMA1, ALPHA0, GAMMA0
REAL (EB) :: TNOW

CALL SCARC_ENTER_ROUTINE('SCARC_METHOD_CG')
TNOW = CURRENT_TIME()

NS = NSTACK
NL = NLEVEL
NP = NPARENT

!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Get parameters for current scope (note: NL denotes the finest level)
!>   - Get right hand side vector and clear solution vectors
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL)

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_VECTOR_INIT (X, 0.0_EB, NL)                    ! set x to zero
   CALL SCARC_FILTER_MEANVALUE(F, NL)                        ! filter out mean value of F
   CALL SCARC_SETUP_CONDENSING (F, NL, 1)                    ! setup condensed system
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Compute initial residual and perform initial preconditioning
!> ------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (X, W, NL)                         !>  W := A*X
CALL SCARC_VECTOR_SUM     (F, W, -1.0_EB, 1.0_EB, NL)        !>  W := W - F

RES = SCARC_L2NORM (W, NL)                                   !>  RESIN := ||W||
RESIN = RES

ITE = 0
NSTATE = SCARC_CONVERGENCE_STATE (0, NL)                     !>  RES < TOL already ??

IF (NSTATE /= NSCARC_STATE_CONV0) THEN                       !>  if no convergence yet, start precon
   CALL SCARC_PRECONDITIONER(NS, NS, NL)
   SIGMA0 = SCARC_SCALAR_PRODUCT(W, G, NL)                   !>  SIGMA0 := (W,G)
   CALL SCARC_VECTOR_COPY (G, D, -1.0_EB, NL)                !>  D := -G
ELSE
   NIT=0                                                     !>  if already convergence, don't iterate
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Perform conjugate gradient looping
!> ------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   CALL SCARC_MATVEC_PRODUCT (D, Y, NL)                      !>  Y := A*D

   ALPHA0 = SCARC_SCALAR_PRODUCT (D, Y, NL)                  !>  ALPHA0 := (D,Y)
   ALPHA0 = SIGMA0/ALPHA0

   CALL SCARC_VECTOR_SUM (D, X, ALPHA0, 1.0_EB, NL)          !>  X := ALPHA0*D + X
   CALL SCARC_VECTOR_SUM (Y, W, ALPHA0, 1.0_EB, NL)          !>  W := ALPHA0*Y + W

   RES = SCARC_L2NORM (W, NL)                                !>  RES := ||W||
   NSTATE = SCARC_CONVERGENCE_STATE (0, NL)                  !>  RES < TOL ??
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP

   CALL SCARC_PRECONDITIONER(NS, NS, NL)

   SIGMA1 = SCARC_SCALAR_PRODUCT (W, G, NL)                  !>  SIGMA1 := (W,G)
   GAMMA0 = SIGMA1/SIGMA0
   SIGMA0 = SIGMA1

   CALL SCARC_VECTOR_SUM (G, D, -1.0_EB, GAMMA0, NL)         !>  D := -G + GAMMA0*D

ENDDO CG_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of CG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(NSTATE)

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_RESTORE_LAST_CELL(X, NL)
   CALL SCARC_FILTER_MEANVALUE(X, NL)
ENDIF

IF (TYPE_SCOPE == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS  (NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS (NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NS, NP)
CALL SCARC_LEAVE_ROUTINE()

TSTEP(MYID+1)%KRYLOV=MAX(TSTEP(MYID+1)%KRYLOV,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%KRYLOV =TSUM(MYID+1)%KRYLOV+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_METHOD_CG




!> ------------------------------------------------------------------------------------------------
!> Perform global BICGstab-method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_BICG(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER   :: NS, NL, NP
INTEGER   :: NSTATE
REAL (EB) :: ALPHA0, ALPHA1, ALPHA2, RHO0, RHO1, DTHETA, DBETA
REAL (EB) :: TNOW

CALL SCARC_ENTER_ROUTINE('SCARC_METHOD_BICG')
TNOW = CURRENT_TIME()

NS = NSTACK
NL = NLEVEL
NP = NPARENT

!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Save SETTINGS (in case that subsequent solvers with different settings are called)
!>   - Define parameters for current scope (note: NL denotes the finest level)
!>   - Initialize solution, right hand side vector
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL)

!> ------------------------------------------------------------------------------------------------
!> Compute initial defect and perform (double) initial preconditioning
!> ------------------------------------------------------------------------------------------------
ALPHA0 = 1.0_EB
RHO0   = 1.0_EB
RHO1   = 0.0_EB
DBETA  = 0.0_EB
DTHETA = 1.0_EB

!CALL SCARC_VECTOR_COPY    (F, W, 1.0_EB, NL)                         !>  W := F
!CALL SCARC_BLOCK_SOLVER (W, W, NS, NP, NL)                           !>  W := PRECON(W)
CALL SCARC_MATVEC_PRODUCT (X, W, NL)                                  !>  W := A*X
CALL SCARC_VECTOR_SUM     (F, W, 1.0_EB, -1.0_EB, NL)                 !>  W := F - W
WRITE(*,*) 'ACHTUNG: HIER RICHTIGE VEKTOREN CHECKEN! EHEMALS W,W'
CALL SCARC_PRECONDITIONER(NS, NS, NL)

ITE = 0
RESIN = SCARC_L2NORM (W, NL)                                 !>  RESIN := ||W||
NSTATE = SCARC_CONVERGENCE_STATE (0, NL)                     !>  RES < TOL already ??

IF (NSTATE /= NSCARC_STATE_CONV0) THEN                       !>  if no convergence yet, start precon
   CALL SCARC_VECTOR_COPY (W, G, 1.0_EB, NL)                 !>  G := W
ELSE
   NIT = 0
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Perform bi-conjugate gradient looping:
!> ------------------------------------------------------------------------------------------------
BICG_LOOP: DO ITE = 1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   RHO1  = SCARC_SCALAR_PRODUCT (G, W, NL)                            !> RHO1 := (G,W)
   DBETA = (RHO1*DTHETA)/(RHO0*ALPHA0)
   RHO0  = RHO1

   CALL SCARC_VECTOR_SUM (W, Z, 1.0_EB       , DBETA , NL)            !> Z := W + DBETA*Z
   CALL SCARC_VECTOR_SUM (Y, Z, -DBETA*ALPHA0, 1.0_EB, NL)            !> Z := -DBETA*ALPHA0*Y + Z
   CALL SCARC_MATVEC_PRODUCT (Z, Y, NL)                               !> Y := A*Z
   WRITE(*,*) 'ACHTUNG: HIER RICHTIGEN VEKTOREN CHECKEN! EHEMALS Y,Y'
   CALL SCARC_PRECONDITIONER(NS, NS, NL)

   DTHETA = SCARC_SCALAR_PRODUCT (G, Y, NL)                           !> DTHETA := (G,Y)
   DTHETA = RHO1/DTHETA

   CALL SCARC_VECTOR_SUM (Y, W, -DTHETA, 1.0_EB, NL)                  !> W := -DTHETA*Y + W
   CALL SCARC_MATVEC_PRODUCT (W, D, NL)                               !> D := A*W
   WRITE(*,*) 'ACHTUNG: HIER RICHTIGEN VEKTOREN CHECKEN!, EHEMALS D,D'
   CALL SCARC_PRECONDITIONER(NS, NS, NL)

   ALPHA1 = SCARC_SCALAR_PRODUCT (D, W, NL)                           !> ALPHA1 := (D,W)
   ALPHA2 = SCARC_SCALAR_PRODUCT (D, D, NL)                           !> ALPHA2 := (D,D)
   ALPHA0 = ALPHA1/ALPHA2

   CALL SCARC_VECTOR_SUM (Z, X,  DTHETA, 1.0_EB, NL)                  !> X :=  DTHETA*Z + X
   CALL SCARC_VECTOR_SUM (W, X,  ALPHA0, 1.0_EB, NL)                  !> X :=  ALPHA0*W + X
   CALL SCARC_VECTOR_SUM (D, W, -ALPHA0, 1.0_EB, NL)                  !> W := -ALPHA0*D + W

   RES = SCARC_L2NORM (W, NL)                                         !> RES := ||W||

   NSTATE = SCARC_CONVERGENCE_STATE(0, NL)                            !> RES < TOL ???
IF (BVERBOSE_LESS.AND.MYID==0) &
   WRITE(0,'(a,i3,a,e14.5,a,e14.5)') ' BICG-Iteration  =',ITE,': Residuum=',SCARC_RESIDUAL
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT BICG_LOOP

ENDDO BICG_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of BICG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(NSTATE)

IF (TYPE_SCOPE == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS(NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NS, NP)
CALL SCARC_LEAVE_ROUTINE()

TSTEP(MYID+1)%KRYLOV=MAX(TSTEP(MYID+1)%KRYLOV,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%KRYLOV =TSUM(MYID+1)%KRYLOV+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_METHOD_BICG

!> -----------------------------------------------------------------------------------------------
!> Preconditioning method
!> -----------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRECONDITIONER(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER :: NS, NL, NP, NL0

CALL SCARC_ENTER_ROUTINE('SCARC_PRECONDITIONER')

NS = NSTACK
NL = NLEVEL
NP = NPARENT

SELECT CASE (TYPE_TWOLEVEL)

   !> --------------------------------------------------------------------
   !> classical one-level preconditioning
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_NONE)

      CALL SCARC_VECTOR_COPY (W, G, 1.0_EB, NL)                   !>  G := W 
      CALL SCARC_BLOCK_SOLVER (W, G, NS+1, NP, NL)                !>  G := PRECON(W)

   !> --------------------------------------------------------------------
   !> additive two-level preconditioning 
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_ADD)

      CALL SCARC_VECTOR_COPY (W, F, 1.0_EB, NL)                   !>  G := W
      DO NL0 = NL, NLEVEL_MAX-1
         CALL SCARC_RESTRICTION (F, F, NL0, NL0+1)                !>  F_coarse := rest(W_fine)
      ENDDO
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  X_coarse := A_coarse^{-1}(F_coarse)
      CALL SCARC_VECTOR_COPY (X, Z, 1.0_EB, NLEVEL_MAX)           !>  G := W
      DO NL0 = NLEVEL_MAX-1, NL, -1
         CALL SCARC_PROLONGATION(Z, Z, NL0+1, NL0)                !>  Z := prol(X_coarse)
      ENDDO
      CALL SCARC_VECTOR_COPY (W, G, 1.0_EB, NL)                   !>  G := W
      CALL SCARC_BLOCK_SOLVER (W, G, NS+1, NP, NL)                !>  G := PRECON(W)

      CALL SCARC_VECTOR_SUM (Z, G, 1.0_EB, 1.0_EB, NL)            !>  G := Z + G 

   !> --------------------------------------------------------------------
   !> multiplicative two-level preconditioning (coarse first, fine second)
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_MUL)

      CALL SCARC_VECTOR_COPY (W, F, 1.0_EB, NL)                   !>  G := W
      DO NL0 = NL, NLEVEL_MAX-1
         CALL SCARC_RESTRICTION (F, F, NL0, NL0+1)                !>  F_coarse := rest(W_fine)
      ENDDO

      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  X_coarse := A_coarse^{-1}(F_coarse)
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !>  G := W
      DO NL0 = NLEVEL_MAX-1, NL, -1
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !>  G := prol(X_coarse)
      ENDDO
      CALL SCARC_MATVEC_PRODUCT (Y, Z, NL)                        !>  Z := A_fine*G

      CALL SCARC_VECTOR_SUM (W, Z, 1.0_EB, -1.0_EB, NL)           !>  Z := W - Z
      CALL SCARC_VECTOR_COPY (Z, G, 1.0_EB, NL)                   !>  G := Z
      CALL SCARC_BLOCK_SOLVER (Z, G, NS+1, NP, NL)                !>  G := PRECON(Z)
      CALL SCARC_VECTOR_SUM (Y, G, 1.0_EB, 1.0_EB, NL)            !>  Z := W - Z

   !> --------------------------------------------------------------------
   !> multiplicative two-level preconditioning (fine first, coarse second)
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_MUL2)

      CALL SCARC_VECTOR_COPY (W, G, 1.0_EB, NL)                   !>  G := W
      CALL SCARC_BLOCK_SOLVER (W, G, NS+1, NP, NL)                !>  G := PRECON(W)
      CALL SCARC_MATVEC_PRODUCT (G, Z, NL)                        !>  Z := A_fine*G

      CALL SCARC_VECTOR_SUM (W, Z, 1.0_EB, -1.0_EB, NL)           !>  Z := W - Z

      CALL SCARC_RESTRICTION (Z, F, NL, NL+1)                     !>  F_coarse := rest(W_fine)
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  X_coarse := A_coarse^{-1}(F_coarse)
      CALL SCARC_PROLONGATION (X, Z, NL+1, NL)                    !>  G := prol(X_coarse)
      CALL SCARC_VECTOR_SUM (Z, G, 1.0_EB, 1.0_EB, NL)            !>  Z := W - Z


   !> --------------------------------------------------------------------
   !> only coarse grid preconditioner
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_COARSE)
      CALL SCARC_VECTOR_COPY (W, F, 1.0_EB, NL)                   !>  G := W
      DO NL0 = NL, NLEVEL_MAX-1
         CALL SCARC_RESTRICTION (F, F, NL0, NL0+1)                !>  F_coarse := rest(W_fine)
      ENDDO

      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  X_coarse := A_coarse^{-1}(F_coarse)
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !>  G := W
      DO NL0 = NLEVEL_MAX-1, NL, -1
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !>  G := prol(X_coarse)
      ENDDO
      CALL SCARC_VECTOR_COPY (Y, G, 1.0_EB, NL)                   !>  G := W

END SELECT 
CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_PRECONDITIONER


!> ------------------------------------------------------------------------------------------------
!> Call requested coarse grid solver (iterative/direct)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_COARSE(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL

SELECT CASE (TYPE_COARSE)
   CASE (NSCARC_COARSE_ITERATIVE)
      CALL SCARC_METHOD_CG (NSTACK, NPARENT, NLEVEL)
   CASE (NSCARC_COARSE_DIRECT)
#if defined(WITH_MKL)
      IF (N_MPI_PROCESSES > 1) THEN
         CALL SCARC_METHOD_CLUSTER (NSTACK, NPARENT, NLEVEL)
      ELSE
         CALL SCARC_METHOD_PARDISO (NSTACK, NPARENT, NLEVEL)
      ENDIF
#else
      WRITE(*,*) 'SCARC_METHOD_DIRECT not working yet '
#endif
END SELECT

END SUBROUTINE SCARC_METHOD_COARSE


!> ------------------------------------------------------------------------------------------------
!> Perform geometric multigrid method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MULTIGRID(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER   :: NS, NL, NP
INTEGER   :: NSTATE, ICYCLE
REAL (EB) :: TNOW, TNOW_COARSE

CALL SCARC_ENTER_ROUTINE('SCARC_METHOD_MULTIGRID')
TNOW = CURRENT_TIME()

NS = NSTACK
NL = NLEVEL
NP = NPARENT
!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Save SETTINGS (in case that subsequent solvers with different settings are called)
!>   - Define parameters for current scope (note: NL denotes the finest level)
!>   - Initialize solution, right hand side vector
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL)

!> ------------------------------------------------------------------------------------------------
!> Compute initial defect:  RESIN := || F - A*X ||
!>   - Initialize cycle counts for MG-iteration
!>   - Perform initial matrix-vector product on finest level
!>   - calculate norm of initial residual on finest level
!> ------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (X, D, NL)                                  !>  D := A*X
CALL SCARC_VECTOR_SUM (F, D, 1.0_EB, -1.0_EB, NL)                     !>  D := F - D

ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_SETUP, NL)

ITE = 0
RESIN = SCARC_L2NORM (D, NL)                                          !>  RESIN := ||D||
NSTATE = SCARC_CONVERGENCE_STATE (0, NL)                              !>  RES < TOL already ??

!> ------------------------------------------------------------------------------------------------
!> Perform multigrid-looping (start each iteration on finest level)
!> ------------------------------------------------------------------------------------------------
MULTIGRID_LOOP: DO ITE = 1, NIT

   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLING_EXIT)

      !> presmoothing  (smoothing/restriction till coarsest level is reached)
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_PRESMOOTH, NS+1, NS, NL)         !> D_fine   := smooth(defect)
         CALL SCARC_RESTRICTION (D, F, NL, NL+1)                              !> F_coarse := rest(D_fine)
         CALL SCARC_VECTOR_CLEAR (X, NL+1)                                    !> X_coarse := 0.0
         NL = NL + 1                                                          !> set coarser level
      ENDDO PRESMOOTHING_LOOP

      !> coarse grid solver
      TNOW_COARSE = CURRENT_TIME()
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)                          !> X_coarse := exact_sol(.)
      TSTEP(MYID+1)%COARSE=MAX(TSTEP(MYID+1)%COARSE,CURRENT_TIME()-TNOW_COARSE)
      TSUM(MYID+1)%COARSE =TSUM(MYID+1)%COARSE+CURRENT_TIME()-TNOW_COARSE

      !> postsmoothing (smoothing/restriction till finest level is reached again)
      POSTSMOOTHING_LOOP: DO WHILE (NL > NLEVEL_MIN)
         NL=NL-1
         CALL SCARC_PROLONGATION (X, D, NL+1, NL)                             !> D_fine := prol(X_coarse)
         CALL SCARC_VECTOR_SUM (D, X, 1.0_EB, 1.0_EB, NL)                     !> X_fine := D_fine + X_fine
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_POSTSMOOTH, NS+1, NS, NL)        !> D_fine := smooth(defect)
         ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_PROCEED, NL)           !> perform requested cycle
         IF (ICYCLE /= NSCARC_CYCLING_POSTSMOOTH) CYCLE CYCLE_LOOP
      ENDDO POSTSMOOTHING_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) CALL SCARC_SHUTDOWN('Wrong level for multigrid method ', 'NONE', NL)

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

 !> ---------------------------------------------------------------------------------------------
 !> Compute norm of new residual on finest level and  leave loop correspondingly
 !> ---------------------------------------------------------------------------------------------
   CALL SCARC_MATVEC_PRODUCT (X, D, NL)                                       !> D := A*X
   CALL SCARC_VECTOR_SUM (F, D, 1.0_EB, -1.0_EB, NL)                          !> D := F - D

   RES = SCARC_L2NORM (D, NL)                                                 !> RES := ||D||
   NSTATE = SCARC_CONVERGENCE_STATE(0, NL)                                    !> convergence ?

   IF (BVERBOSE_LESS.AND.MYID==0) &
      WRITE(LU_OUTPUT,'(a,i3,a,e14.5,a,e14.5)') '       SCARC_MG-Iteration  =',ITE,': Residuum=',SCARC_RESIDUAL
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT MULTIGRID_LOOP

ENDDO MULTIGRID_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of MG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries (consistency!)
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(NSTATE)

SELECT CASE (TYPE_SCOPE)
   CASE (NSCARC_SOLVER_MAIN)
      CALL SCARC_UPDATE_PRESSURE_MAINCELLS(NLEVEL_MIN)
      CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
   CASE (NSCARC_SOLVER_PRECON)
      CALL SCARC_UPDATE_PRECONDITIONER(NLEVEL_MIN)
END SELECT

CALL SCARC_RELEASE_SOLVER(NS, NP)
CALL SCARC_LEAVE_ROUTINE()

TSTEP(MYID+1)%MULTIGRID=MAX(TSTEP(MYID+1)%MULTIGRID,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%MULTIGRID =TSUM(MYID+1)%MULTIGRID+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_METHOD_MULTIGRID


!> ------------------------------------------------------------------------------------------------
!> Control multigrid cycling (F/V/W)
!> Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CYCLING_CONTROL(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NL0, ICYCLE
TYPE (SCARC_CYCLING_TYPE), POINTER :: CYC

SELECT CASE (NTYPE)

 !> ---------------------------------------------------------------------------------------------
 !> initialize cycle counts at beginning of multigrid method
 !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLING_SETUP)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CYC => SCARC(NM)%CYCLING(NL)
         CYC%COUNTER(2)=1

         DO NL0 = NLEVEL_MIN+1, NLEVEL_MAX - 1
            CYC => SCARC(NM)%CYCLING(NL0)
            IF (TYPE_CYCLING==NSCARC_CYCLING_F) THEN
               CYC%COUNTER(2)=2
            ELSE
               CYC%COUNTER(2)=TYPE_CYCLING
            ENDIF
         ENDDO
      ENDDO

      ICYCLE = NSCARC_CYCLING_NEXT

 !> ---------------------------------------------------------------------------------------------
 !> reset cycle counts at beginning of each new multigrid iteration
 !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLING_RESET)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         DO NL0 = NLEVEL_MIN, NLEVEL_MAX
            CYC => SCARC(NM)%CYCLING(NL0)
            CYC%COUNTER(1)=CYC%COUNTER(2)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLING_NEXT

 !> ---------------------------------------------------------------------------------------------
 !> determine where to proceed with cycling
 !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLING_PROCEED)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CYC => SCARC(NM)%CYCLING(NL)
         CYC%COUNTER(1)=CYC%COUNTER(1)-1

         IF (CYC%COUNTER(1)==0) THEN
            IF (TYPE_CYCLING==NSCARC_CYCLING_F) THEN
               CYC%COUNTER(1)=1
            ELSE
               CYC%COUNTER(1)=CYC%COUNTER(2)
            ENDIF
            IF (NL == NLEVEL_MIN) THEN
               ICYCLE = NSCARC_CYCLING_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLING_POSTSMOOTH
            ENDIF
         ELSE
            IF (NL == NLEVEL_MIN) THEN
               ICYCLE = NSCARC_CYCLING_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLING_NEXT
            ENDIF
         ENDIF

      ENDDO

END SELECT

SCARC_CYCLING_CONTROL = ICYCLE
RETURN

END FUNCTION SCARC_CYCLING_CONTROL


!> ------------------------------------------------------------------------------------------------
!> Perform smoothing
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SMOOTHER(NTYPE, NSTACK, NPARENT, NLEVEL)
INTEGER , INTENT(IN) :: NTYPE, NSTACK, NPARENT, NLEVEL
INTEGER :: NS, NL, NP, NSTATE=0
REAL(EB):: TNOW
LOGICAL :: BMATVEC, BL2NORM

CALL SCARC_ENTER_ROUTINE('SCARC_SMOOTHER')
TNOW = CURRENT_TIME()

NS = NSTACK
NP = NPARENT
NL = NLEVEL

!> ------------------------------------------------------------------------------------------------
!> Initialization
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_SOLVER(NS, NP)
BL2NORM  = .TRUE.
IF (NTYPE == NSCARC_CYCLING_PRESMOOTH.AND.NL==1) THEN
   BMATVEC = .FALSE.
ELSE
   BMATVEC = .TRUE.
ENDIF
!> only temporarily for debugging purposes
BMATVEC = .TRUE.
BL2NORM = .TRUE.


!> ------------------------------------------------------------------------------------------------
!> Calculate initial defect on level NL (only if BMATVEC = .TRUE.)
!> Because initial vector is set to zero, this defect corresponds to F
!> ------------------------------------------------------------------------------------------------
IF (BMATVEC) THEN
   CALL SCARC_MATVEC_PRODUCT (X, D, NL)                                  !>  D := A*X
   CALL SCARC_VECTOR_SUM (F, D, 1.0_EB, -1.0_EB, NL)                     !>  D := F - D
ENDIF

IF (BL2NORM.AND.BMATVEC) THEN
   RESIN = SCARC_L2NORM (D, NL)                                          !>  RESIN := ||D||
ELSE
   RESIN = SCARC_RESIDUAL
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Smoothing loop
!> ------------------------------------------------------------------------------------------------
SMOOTH_LOOP: DO ITE=1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

#ifdef WITH_MKL
   IF (TYPE_SMOOTH == NSCARC_RELAX_PARDISO .OR. TYPE_SMOOTH == NSCARC_RELAX_CLUSTER) THEN
      CALL SCARC_VECTOR_COPY(D, Z, 1.0_EB, NL)
      CALL SCARC_BLOCK_SOLVER (Z, D, NS, NP, NL)                        !>  D := PRECON (D)
   ELSE
      CALL SCARC_BLOCK_SOLVER (D, D, NS, NP, NL)                        !>  D := PRECON (D)
   ENDIF
#else
   CALL SCARC_BLOCK_SOLVER (D, D, NS, NP, NL)                           !>  D := PRECON (D)
#endif

   CALL SCARC_VECTOR_SUM      (D, X, OMEGA, 1.0_EB, NL)                 !>  X := OMEGA*D + X
   CALL SCARC_MATVEC_PRODUCT  (X, D, NL)                                !>  D := A*X

   CALL SCARC_VECTOR_SUM      (F, D, 1.0_EB, -1.0_EB, NL)               !>  D := F - D

   IF (BL2NORM.OR.ITE==NIT) THEN
      RES = SCARC_L2NORM (D, NL)                                        !>  RES := ||D||
      NSTATE = SCARC_CONVERGENCE_STATE(NTYPE, NL)  
      IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP              !>  RES < TOL ?
   ENDIF

ENDDO SMOOTH_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> ------------------------------------------------------------------------------------------------
!ITE = ITE -1
!IF (BL2NORM) CALL SCARC_CONVERGENCE_RATE(NS, NSTATE)

CALL SCARC_RELEASE_SOLVER(NS, NP)
CALL SCARC_LEAVE_ROUTINE()

TSTEP(MYID+1)%SMOOTH=MAX(TSTEP(MYID+1)%SMOOTH,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%SMOOTH =TSUM(MYID+1)%SMOOTH+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SMOOTHER


!> ------------------------------------------------------------------------------------------------
!> Setup environement in every solver CALL (i.e. set pointers to used vectors) related to NSTACK
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SOLVER(NSTACK, NPARENT)
INTEGER, INTENT(IN) :: NSTACK, NPARENT
TYPE (SCARC_SOLVER_TYPE),      POINTER :: SOL
TYPE (SCARC_CONVERGENCE_TYPE), POINTER :: CON, CONP
TYPE (SCARC_TYPES_TYPE)      , POINTER :: TYP
TYPE (SCARC_POINTERS_TYPE)   , POINTER :: PTR

SOL => STACK(NSTACK)%SOLVER

CON => SOL%CONVREQS
TYP => SOL%TYPES
PTR => SOL%POINTERS

CNAME = SOL%CNAME

NIT   = CON%NIT
EPS   = CON%EPS
OMEGA = CON%OMEGA

!> if not first solver in stack, restore last iteration parameters of predecessor
IF (NPARENT > 0) THEN
   CONP => STACK(NPARENT)%SOLVER%CONVREQS
   CONP%ITE   = ITE
   CONP%RES   = RES
   CONP%RESIN = RESIN
   CONP%ERR   = ERR  
   CONP%CAPPA = CAPPA
ENDIF

TYPE_PARENT   = NPARENT
TYPE_METHOD   = TYP%TYPE_METHOD
TYPE_SOLVER   = TYP%TYPE_SOLVER
TYPE_SCOPE    = TYP%TYPE_SCOPE
TYPE_RELAX    = TYP%TYPE_RELAX
TYPE_INTERPOL = TYP%TYPE_INTERPOL
TYPE_TWOLEVEL = TYP%TYPE_TWOLEVEL
TYPE_CYCLING  = TYP%TYPE_CYCLING
TYPE_ACCURACY = TYP%TYPE_ACCURACY

X = PTR%X
F = PTR%F
D = PTR%D
G = PTR%G
W = PTR%W
Y = PTR%Y
Z = PTR%Z
E = PTR%E

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) ITE_TOTAL = 0

END SUBROUTINE SCARC_SETUP_SOLVER


!> ------------------------------------------------------------------------------------------------
!> Reset settings of calling CURRENT-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELEASE_SOLVER(NSTACK, NPARENT)
INTEGER, INTENT(IN)  :: NSTACK, NPARENT
TYPE (SCARC_SOLVER_TYPE),      POINTER :: PAR
TYPE (SCARC_CONVERGENCE_TYPE), POINTER :: CON
TYPE (SCARC_POINTERS_TYPE)   , POINTER :: PTR
TYPE (SCARC_TYPES_TYPE)      , POINTER :: TYP

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   SCARC_CAPPA      = CAPPA
   SCARC_RESIDUAL   = RES
   SCARC_ITERATIONS = ITE
   !CALL SCARC_CLOSE_CSV_FILE()
ENDIF

CON => STACK(NSTACK)%SOLVER%CONVREQS

CON%RESIN = RESIN
CON%RES   = RES
CON%ITE   = ITE
CON%ERR   = ERR  

!> if not first solver in stack, reset environment of parent (calling) routine
IF (NPARENT > 0) THEN

   PAR  => STACK(NPARENT)%SOLVER
   CNAME = PAR%CNAME

   CON  => PAR%CONVREQS
   PTR  => PAR%POINTERS
   TYP  => PAR%TYPES

   ITE   = CON%ITE
   NIT   = CON%NIT
   EPS   = CON%EPS
   RESIN = CON%RESIN
   RES   = CON%RES
   OMEGA = CON%OMEGA
   CAPPA = CON%CAPPA

   TYPE_PARENT   = TYP%TYPE_PARENT
   TYPE_METHOD   = TYP%TYPE_METHOD
   TYPE_SOLVER   = TYP%TYPE_SOLVER
   TYPE_SCOPE    = TYP%TYPE_SCOPE
   TYPE_RELAX    = TYP%TYPE_RELAX
   TYPE_INTERPOL = TYP%TYPE_INTERPOL
   TYPE_TWOLEVEL = TYP%TYPE_TWOLEVEL
   TYPE_CYCLING  = TYP%TYPE_CYCLING
   TYPE_ACCURACY = TYP%TYPE_ACCURACY

   X = PTR%X
   F = PTR%F
   D = PTR%D
   G = PTR%G
   W = PTR%W
   Y = PTR%Y
   Z = PTR%Z
   E = PTR%E

ENDIF

END SUBROUTINE SCARC_RELEASE_SOLVER

!> ----------------------------------------------------------------------------------------------------
!> Dump residual information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_RESIDUAL(ISM, NL)
INTEGER, INTENT(IN) :: ISM, NL
IF (BCSV_NONE) RETURN
IF (MYID == 0) THEN
   IF (ITE_TOTAL == 0) THEN
      WRITE(MSG%LU_CSV,1000) ITE_PRES, ITE_TOTAL, ITE_CG, ITE_MG, ITE_LU, ITE_COARSE, ITE_SMOOTH, &
                             ISM, NL, TYPE_SOLVER, RESIN
   ELSE
      WRITE(MSG%LU_CSV,1000) ITE_PRES, ITE_TOTAL, ITE_CG, ITE_MG, ITE_LU, ITE_COARSE, ITE_SMOOTH, &
                             ISM, NL, TYPE_SOLVER, RES
   ENDIF
ENDIF
!1000 FORMAT(I10,',',i4,',',i7,',',i7,',',i7,',',i4,',',i3,',',i7,',',i7,',',E20.12,',',E20.12)
1000 FORMAT(I8,',',I8,',',I8,',',I8,',',I8,',',I8,',',I4,',',I4,',',I4,',',I4,',',E20.12)
END SUBROUTINE SCARC_DUMP_RESIDUAL

!> ----------------------------------------------------------------------------------------------------
!> Dump residual information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_CAPPA
IF (BCSV_NONE) RETURN
IF (MYID == 0) WRITE(MSG%LU_CSV,1000) -1, ITE_TOTAL, -1, -1, -1, -1, -1, -1, -1, -1, CAPPA
1000 FORMAT(I8,',',I8,',',I8,',',I8,',',I8,',',I8,',',I4,',',I4,',',I4,',',I4,',',E20.12)
END SUBROUTINE SCARC_DUMP_CAPPA

!> ----------------------------------------------------------------------------------------------------
!> Set initial solution corresponding to boundary data in BXS, BXF, ...
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WORKSPACE(ISTACK, NL)
INTEGER, INTENT(IN) :: ISTACK, NL
INTEGER :: NM, IW, IOR0, I, J, K, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: PRHS, HP
REAL(EB) :: VAL
TYPE (MESH_TYPE)        , POINTER :: M
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_WALL_TYPE)  , POINTER :: WC
TYPE (SCARC_SCOPE_TYPE) , POINTER :: SC
TYPE (SCARC_SOLVER_TYPE), POINTER :: SOL
TYPE (SCARC_TYPES_TYPE) , POINTER :: TYP

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_WORKSPACE')

SOL => STACK(ISTACK)%SOLVER
TYP => SOL%TYPES

SELECT CASE (TYP%TYPE_SOLVER)

   !> --------------- IF used as main solver use values from pres-routine as initialization 
   CASE (NSCARC_SOLVER_MAIN) 

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
  
         M  => MESHES(NM)
         L  => SCARC(NM)%LEVEL(NL)
         SC => SCARC(NM)%SCOPE(TYP%TYPE_SCOPE, NL)
      
         PRHS => M%PRHS
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF
 
         !> get right hand side (PRHS from pres.f90) and initial vector (H or HS from last time step)
         DO K = 1, M%KBAR
            DO J = 1, M%JBAR
               DO I = 1, M%IBAR
                  IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%CELL%STATE(I,J,K) /= NSCARC_DISCRET_GASPHASE) CYCLE
                  IC = L%CELL%DOF(I,J,K)
                  SC%F(IC) = PRHS(I, J, K)                 !> use right hand side from pres-routine
                  SC%X(IC) = HP(I, J, K)                   !> use last iterate as initial solution
               ENDDO
            ENDDO
         ENDDO

         LEVEL_WALL_CELLS_LOOP: DO IW = 1, M%N_EXTERNAL_WALL_CELLS

            WC => L%WALL(IW)

            I    = WC%IXW
            J    = WC%IYW
            K    = WC%IZW

            IF (TWO_D .AND. J /= 1) CALL SCARC_SHUTDOWN('Wrong index for J = ', 'NONE', J)
            IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. L%CELL%STATE(I,J,K) /= NSCARC_DISCRET_GASPHASE) CYCLE

            IOR0 = WC%IOR
            IC   = L%CELL%DOF(I,J,K)

            !> Dirichlet BC's:
            !> these are based on the settings in BTYPE
            !> in the structured case this corresponds to the face-wise settings according to the FFT
            !> (this allows to use local FFT's as preconditioners)
            !> in the unstructured case only open boundary cells lead to Dirichlet BC's
            IF_DIRICHLET: IF (WC%BTYPE == DIRICHLET) THEN

               SELECT CASE (IOR0)
                  CASE (1)
                     VAL = - 2.0_EB * L%COORD%DXI2 * M%BXS(J,K) 
                  CASE (-1)
                     VAL = - 2.0_EB * L%COORD%DXI2 * M%BXF(J,K) 
                  CASE (2)
                     VAL = - 2.0_EB * L%COORD%DYI2 * M%BYS(I,K)
                  CASE (-2)
                     VAL = - 2.0_EB * L%COORD%DYI2 * M%BYF(I,K) 
                  CASE (3)
                     VAL = - 2.0_EB * L%COORD%DZI2 * M%BZS(I,J)
                  CASE (-3)
                     VAL = - 2.0_EB * L%COORD%DZI2 * M%BZF(I,J)
               END SELECT

               SC%F(IC) = SC%F(IC) + VAL

            ENDIF IF_DIRICHLET

            !> Neumann BC's:
            !> Note for the unstructured case only:
            !> Here, the matrix also contains Neumann BC's for those cells which have a
            !> PRESSURE_BC_INDEX == DIRICHLET but are NOT open; these cells must be excluded below,
            !> because BXS, BXF, ... contain the Dirichlet information from pres.f90 there;
            !> excluding them corresponds to a homogeneous Neumann condition for these cells
            IF_NEUMANN: IF (WC%BTYPE == NEUMANN) THEN

               IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. M%WALL(IW)%PRESSURE_BC_INDEX /= NEUMANN) CYCLE

               SELECT CASE (IOR0)
                  CASE (1)
                     VAL =   L%COORD%DXI * M%BXS(J,K)        
                  CASE (-1)
                     VAL = - L%COORD%DXI * M%BXF(J,K)       
                  CASE (2)
                     VAL =   L%COORD%DYI * M%BYS(I,K)      
                  CASE (-2)
                     VAL = - L%COORD%DYI * M%BYF(I,K)     
                  CASE (3)
                     VAL =   L%COORD%DZI * M%BZS(I,J)    
                  CASE (-3)
                     VAL = - L%COORD%DZI * M%BZF(I,J)   
               END SELECT

               SC%F(IC) = SC%F(IC) + VAL

            ENDIF IF_NEUMANN

         ENDDO LEVEL_WALL_CELLS_LOOP
      ENDDO
   
      !> In case of a Krylov method clear overlapping parts of auxiliary vectors
      IF (BCG.OR.BTWOLEVEL) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            SC => SCARC(NM)%SCOPE(TYP%TYPE_SCOPE, NL)
            SC%D(L%NCS+1:L%NCE) = 0.0_EB
            SC%G(L%NCS+1:L%NCE) = 0.0_EB
            SC%Y(L%NCS+1:L%NCE) = 0.0_EB
            SC%W(L%NCS+1:L%NCE) = 0.0_EB
            SC%Z(L%NCS+1:L%NCE) = 0.0_EB
         ENDDO
      ENDIF

      !> In case of a multigrid method as main solver clear
      !> overlapping parts of auxiliary vectors and coarse grid solver vectors
      IF (BGMG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            SC => SCARC(NM)%SCOPE(TYP%TYPE_SCOPE, NL)
            SC%D(L%NCS+1:L%NCE) = 0.0_EB
            SC%Z(L%NCS+1:L%NCE) = 0.0_EB
         ENDDO
      ENDIF

      !> In case of pure Neumann or periodic BCs, broadcast RHS(end) from last mesh
      !> to all and store it on all meshes
      IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
         IF (UPPER_MESH_INDEX == NMESHES) THEN
            L  => SCARC(NMESHES)%LEVEL(NL)
            SC => SCARC(NMESHES)%SCOPE(TYP%TYPE_SCOPE, NL)
            LOCAL_REAL = SC%F(L%NC)
         ELSE
            LOCAL_REAL = 0.0_EB
         ENDIF
         GLOBAL_REAL = SCARC_BROADCAST_REAL(NSCARC_BROADCAST_LAST)
         DO NM = 1, NMESHES
            SCARC(NM)%RHS_END = GLOBAL_REAL
         ENDDO
      ENDIF


   !> --------------- If MG is used as Krylov-preconditioner, vector G of main Krylov is the RHS for MG 
   CASE (NSCARC_SOLVER_PRECON)

      IF (BCGGMG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            SC => SCARC(NM)%SCOPE(TYP%TYPE_SCOPE, NL)
            SC%F = SCARC(NM)%SCOPE(NSCARC_SCOPE_ONE, NL)%G
            SC%D(L%NCS+1:L%NCE) = 0.0_EB
            SC%Z(L%NCS+1:L%NCE) = 0.0_EB
         ENDDO
      ENDIF

   !> --------------- If used as coarse grid solver start with zero initialization 
   !> ACHTUNG: Löschen der Vektoren nur fuer iterativn Coarse-Solver noetig !
   CASE (NSCARC_SOLVER_COARSE)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         SC => SCARC(NM)%SCOPE(TYP%TYPE_SCOPE, NL)
         SC%X = 0.0_EB
         SC%D = 0.0_EB
         SC%G = 0.0_EB
         SC%W = 0.0_EB
         SC%Y = 0.0_EB
         SC%Z = 0.0_EB
      ENDDO         

END SELECT
CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_WORKSPACE

!> ------------------------------------------------------------------------------------------------
!> Check if solver converges or diverges and print out residual information
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CONVERGENCE_STATE(ISM, NL)
INTEGER, INTENT(IN) :: NL, ISM
INTEGER :: NSTATE

NSTATE = NSCARC_STATE_PROCEED

IF (BVERBOSE_LESS.AND.MYID==0) WRITE(LU_OUTPUT,    1000) TRIM(CNAME), NL, ITE, RES
IF (BDEBUG_LESS)               WRITE(MSG%LU_DEBUG, 1000) TRIM(CNAME), NL, ITE, RES

CALL SCARC_DUMP_RESIDUAL(ISM, NL)

SELECT CASE (TYPE_ACCURACY)
   CASE (NSCARC_ACCURACY_RELATIVE)
      IF (RES <= RESIN*EPS)                    NSTATE = NSCARC_STATE_CONV
      IF (RES <= NSCARC_THRESHOLD_CONVERGENCE) NSTATE = NSCARC_STATE_CONV
   CASE (NSCARC_ACCURACY_ABSOLUTE)
      !IF (RES <= EPS .AND. RES <= RESIN*1.0E-2) THEN
      IF (RES <= EPS .AND. RES <= RESIN) THEN
         IF (ITE == 0) THEN
            NSTATE = NSCARC_STATE_CONV0
         ELSE
            NSTATE = NSCARC_STATE_CONV
         ENDIF
      ENDIF
END SELECT
IF (RES > NSCARC_THRESHOLD_DIVGERGENCE) NSTATE = NSCARC_STATE_DIVG

SCARC_CONVERGENCE_STATE = NSTATE
RETURN

1000 FORMAT (5X,A30,': level=',i4,': #ite= ',i4,': res =',e25.16)
END FUNCTION SCARC_CONVERGENCE_STATE


!> ------------------------------------------------------------------------------------------------
!> Compute convergence rate and print out residual information for final loop
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_RATE(NSTATE)
INTEGER, INTENT(IN) :: NSTATE

IF (NSTATE == NSCARC_STATE_DIVG) THEN
   ITE    = - 1
   CAPPA  = 1.0_EB
ELSE
   IF (NSTATE == NSCARC_STATE_CONV0) THEN
     ITE= 0
   ELSE IF (NSTATE == NSCARC_STATE_CONV) THEN
     ITE= ITE
   ELSE
     ITE= ITE-1
   ENDIF
   IF (RESIN >= TWO_EPSILON_EB) THEN
      IF (ITE== 0) THEN
         !CAPPA = (RES/RESIN)
         CAPPA = 0.0_EB
      ELSE
         IF (NSTATE == NSCARC_STATE_CONV0) THEN
            CAPPA = 0.0E0
         ELSE
            CAPPA = (RES/RESIN) ** (1.0_EB/ITE)
         ENDIF
      ENDIF
   ELSE
      CAPPA = 0.0_EB
   ENDIF
ENDIF

CALL SCARC_DUMP_CAPPA()

IF (BDEBUG_LESS) WRITE(MSG%LU_DEBUG,2000) CNAME, ITE, CAPPA
IF (MYID==0.AND.BVERBOSE_LESS) THEN
   IF (TRIM(CNAME) /=  'SCARC_COARSE_CG') WRITE(LU_OUTPUT,2000) CNAME, ITE, CAPPA
ENDIF

2000 FORMAT (/,7X,A25,': iterations: ',i6,':  convergence rate =',e14.6,/)
END SUBROUTINE SCARC_CONVERGENCE_RATE


!> ------------------------------------------------------------------------------------------------
!> Perform restriction from finer to coarser grid in multigrid method
!>    - 'VF' corresponds to vector on fine   grid
!>    - 'VC' corresponds to vector on coarse grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTRICTION (NVECTORF, NVECTORC, NLF, NLC)
INTEGER, INTENT(IN) :: NVECTORF, NVECTORC, NLF, NLC
INTEGER :: NM, ICOL, IC
REAL(EB), POINTER, DIMENSION(:) :: VC, VF
INTEGER , POINTER :: NXC, NYC, NZC, NCC
INTEGER  :: NXF, NYF, NZF
INTEGER  :: IXF, IYF, IZF, ICF(8)=0, ICFB(-2:2,-2:2)=0
INTEGER  :: IXC, IYC, IZC, ICC
REAL(EB) :: AUX
TYPE (SCARC_LEVEL_TYPE) , POINTER :: LC, LF
TYPE (SCARC_CELL_TYPE)  , POINTER :: CC, CF
TYPE (SCARC_MATRIX_TYPE), POINTER :: RF


CALL SCARC_ENTER_ROUTINE('SCARC_RESTRICTION')
IF (BAMG) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_VECTOR, NLF)

!>
!> ------------------ Twolevel-CG or Geometric multigrid (as main solver or preconditioner) --------------
!>
IF (BMULTILEVEL) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LC => SCARC(NM)%LEVEL(NLC)                 !> pointer to coarse grid
      LF => SCARC(NM)%LEVEL(NLF)                 !> pointer to fine grid

      CC => SCARC(NM)%LEVEL(NLC)%CELL            !> pointer to coarse grid
      CF => SCARC(NM)%LEVEL(NLF)%CELL            !> pointer to fine grid

      NXC => LC%NX
      NYC => LC%NY
      NZC => LC%NZ

      NXF = 2*NXC
      NYF = 2*NYC
      NZF = 2*NZC

      VF => POINT_TO_VECTOR(NVECTORF, NM, NLF)
      VC => POINT_TO_VECTOR(NVECTORC, NM, NLC)

      IF (TWO_D) THEN

         SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)
  
            CASE (NSCARC_INTERPOL_CONSTANT)

               DO IZC = 1, NZC
                  DO IXC = 1, NXC

                     IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. CC%STATE(IXC,1,IZC)/=NSCARC_DISCRET_GASPHASE) CYCLE

                     IXF = 2*IXC
                     IZF = 2*IZC

                     ICC = CC%DOF(IXC, 1, IZC) 

                     ICF(1) = CF%DOF(IXF-1, 1, IZF-1) 
                     ICF(2) = CF%DOF(IXF-1, 1, IZF  ) 
                     ICF(3) = CF%DOF(IXF  , 1, IZF-1) 
                     ICF(4) = CF%DOF(IXF  , 1, IZF  ) 

                     VC(ICC) = 0.25_EB * (  VF(ICF(1)) &
                                          + VF(ICF(2)) &
                                          + VF(ICF(3)) &
                                          + VF(ICF(4)) )
                  ENDDO
               ENDDO

            CASE (NSCARC_INTERPOL_BILINEAR)
    
               VC=0.0_EB
      
               DO IZC = 1, NZC
                  DO IXC = 1, NXC
                     
                     IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. CC%STATE(IXC,1,IZC)/=NSCARC_DISCRET_GASPHASE) CYCLE
         
                     IXF = 2*IXC
                     IZF = 2*IZC
         
                     ICC = CC%DOF(IXC, 1, IZC) 
         
                     ICFB(-2,-2) = CF%DOF(IXF-2, 1, IZF-2) 
                     ICFB(-1,-2) = CF%DOF(IXF-1, 1, IZF-2) 
                     ICFB( 1,-2) = CF%DOF(IXF  , 1, IZF-2) 
                     ICFB( 2,-2) = CF%DOF(IXF+1, 1, IZF-2) 
      
                     ICFB(-2,-1) = CF%DOF(IXF-2, 1, IZF-1) 
                     ICFB(-1,-1) = CF%DOF(IXF-1, 1, IZF-1) 
                     ICFB( 1,-1) = CF%DOF(IXF  , 1, IZF-1) 
                     ICFB( 2,-1) = CF%DOF(IXF+1, 1, IZF-1) 
      
                     ICFB(-2, 1) = CF%DOF(IXF-2, 1, IZF) 
                     ICFB(-1, 1) = CF%DOF(IXF-1, 1, IZF) 
                     ICFB( 1, 1) = CF%DOF(IXF  , 1, IZF) 
                     ICFB( 2, 1) = CF%DOF(IXF+1, 1, IZF) 
      
                     ICFB(-2, 2) = CF%DOF(IXF-2, 1, IZF+1) 
                     ICFB(-1, 2) = CF%DOF(IXF-1, 1, IZF+1) 
                     ICFB( 1, 2) = CF%DOF(IXF  , 1, IZF+1) 
                     ICFB( 2, 2) = CF%DOF(IXF+1, 1, IZF+1) 
         
                     IF (IXC==1.AND.IZC==1) THEN
                        VC(ICC) = SCALR*( &
                                     W4 *VF(ICFB(-1, 2)) + W3 *VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) + &
                                     W12*VF(ICFB(-1, 1)) + W9 *VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) + &
                                     W16*VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) + W4*VF(ICFB(2,-1)) )
                     ELSE IF (IXC==NXC.AND.IZC==  1) THEN
                        VC(ICC) = SCALR*( &
                                     W1 *VF(ICFB(-2, 2)) + W3 *VF(ICFB(-1, 2)) + W4 *VF(ICFB(1, 2)) + &
                                     W3 *VF(ICFB(-2, 1)) + W9 *VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) + &
                                     W4 *VF(ICFB(-2,-1)) + W12*VF(ICFB(-1,-1)) + W16*VF(ICFB(1,-1)) )
                     ELSE IF (IXC==  1.AND.IZC==NZC) THEN
                        VC(ICC) = SCALR*( &
                                     W16*VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) + W4*VF(ICFB(2, 1)) + &
                                     W12*VF(ICFB(-1,-1)) + W9 *VF(ICFB(1,-1)) + W3*VF(ICFB(2,-1)) + &
                                     W4 *VF(ICFB(-1,-2)) + W3 *VF(ICFB(1,-2)) + W1*VF(ICFB(2,-2)) )
                     ELSE IF (IXC==NXC.AND.IZC==NZC) THEN
                        VC(ICC) = SCALR*( &
                                     W4 *VF(ICFB(-2, 1)) + W12*VF(ICFB(-1, 1)) + W16*VF(ICFB(1, 1)) + &
                                     W3 *VF(ICFB(-2,-1)) + W9 *VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) + &
                                     W1 *VF(ICFB(-2,-2)) + W3 *VF(ICFB(-1,-2)) + W4 *VF(ICFB(1,-2)) )
                     ELSE IF (IZC==  1) THEN
                        VC(ICC) = SCALR*( &
                                     W1*VF(ICFB(-2, 2)) + W3 *VF(ICFB(-1, 2)) + W3 *VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) + &
                                     W3*VF(ICFB(-2, 1)) + W9 *VF(ICFB(-1, 1)) + W9 *VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) + &
                                     W4*VF(ICFB(-2,-1)) + W12*VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) + W4*VF(ICFB(2,-1)) )
                     ELSE IF (IZC==NZC) THEN
                        VC(ICC) = SCALR*( &
                                     W4*VF(ICFB(-2, 1)) + W12*VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) + W4*VF(ICFB(2, 1)) + &
                                     W3*VF(ICFB(-2,-1)) + W9 *VF(ICFB(-1,-1)) + W9 *VF(ICFB(1,-1)) + W3*VF(ICFB(2,-1)) + &
                                     W1*VF(ICFB(-2,-2)) + W3 *VF(ICFB(-1,-2)) + W3 *VF(ICFB(1,-2)) + W1*VF(ICFB(2,-2)) )
                     ELSE IF (IXC==  1) THEN
                        VC(ICC) = SCALR*( &
                                     W4 *VF(ICFB(-1, 2)) + W3*VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) +&
                                     W12*VF(ICFB(-1, 1)) + W9*VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) +&
                                     W12*VF(ICFB(-1,-1)) + W9*VF(ICFB(1,-1)) + W3*VF(ICFB(2,-1)) +&
                                     W4 *VF(ICFB(-1,-2)) + W3*VF(ICFB(1,-2)) + W1*VF(ICFB(2,-2)) )
                     ELSE IF (IXC==NXC) THEN
                        VC(ICC) = SCALR*( &
                                     W1*VF(ICFB(-2, 2)) + W3*VF(ICFB(-1, 2)) + W4 *VF(ICFB(1, 2)) + &
                                     W3*VF(ICFB(-2, 1)) + W9*VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) +&
                                     W3*VF(ICFB(-2,-1)) + W9*VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) +&
                                     W1*VF(ICFB(-2,-2)) + W3*VF(ICFB(-1,-2)) + W4 *VF(ICFB(1,-2)) )
                     ELSE 
                        VC(ICC) = SCALR*( &
                                     W1*VF(ICFB(-2,-2)) + W3*VF(ICFB(-1,-2)) + W3*VF(ICFB(1,-2)) + W1*VF(ICFB(2,-2)) +&
                                     W3*VF(ICFB(-2,-1)) + W9*VF(ICFB(-1,-1)) + W9*VF(ICFB(1,-1)) + W3*VF(ICFB(2,-1)) +&
                                     W3*VF(ICFB(-2, 1)) + W9*VF(ICFB(-1, 1)) + W9*VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) +&
                                     W1*VF(ICFB(-2, 2)) + W3*VF(ICFB(-1, 2)) + W3*VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) )
                     ENDIF
                  ENDDO
               ENDDO

         END SELECT SELECT_INTERPOL
      
      ELSE

         DO IZC = 1, NZC
            DO IYC = 1, NYC
               DO IXC = 1, NXC

                  IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. CC%STATE(IXC,IYC,IZC)/=NSCARC_DISCRET_GASPHASE) CYCLE

                  IXF = 2*IXC
                  IYF = 2*IYC
                  IZF = 2*IZC

                  ICC = CC%DOF(IXC, IYC, IZC) 

                  ICF(1) = CF%DOF(IXF-1, IYF-1, IZF-1) 
                  ICF(2) = CF%DOF(IXF-1, IYF-1, IZF  ) 
                  ICF(3) = CF%DOF(IXF-1, IYF  , IZF-1) 
                  ICF(4) = CF%DOF(IXF-1, IYF  , IZF  ) 
                  ICF(5) = CF%DOF(IXF  , IYF-1, IZF-1) 
                  ICF(6) = CF%DOF(IXF  , IYF-1, IZF  ) 
                  ICF(7) = CF%DOF(IXF  , IYF  , IZF-1) 
                  ICF(8) = CF%DOF(IXF  , IYF  , IZF  ) 

                  VC(ICC) = 0.125_EB * (  VF(ICF(1)) &
                                        + VF(ICF(2)) &
                                        + VF(ICF(3)) &
                                        + VF(ICF(4)) &
                                        + VF(ICF(5)) &
                                        + VF(ICF(6)) &
                                        + VF(ICF(7)) &
                                        + VF(ICF(8)) )

               ENDDO
            ENDDO
         ENDDO

      ENDIF

   ENDDO

!>
!> ------------------ Algebraic multigrid -------------------------------------------
!>
ELSE IF (BAMG) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LC => SCARC(NM)%LEVEL(NLC)                   !> pointer to coarse grid
      LF => SCARC(NM)%LEVEL(NLF)                   !> pointer to fine grid

      RF => SCARC(NM)%SYSTEM(NLF)%R            !> pointer to restriction matrix of fine grid

      VC => POINT_TO_VECTOR(NVECTORC, NM, NLC)    !> pointer to coarse vector
      VF => POINT_TO_VECTOR(NVECTORF, NM, NLF)    !> pointer to fine vector

      NCC => LC%NCS

      DO ICC = 1, NCC
         AUX = 0.0_EB
         DO ICOL = RF%ROW(ICC), RF%ROW(ICC+1)-1
            IC = RF%COL(ICOL)
            AUX = AUX + VF(IC) * RF%VAL(ICOL)
         ENDDO
         VC(ICC) = AUX
      ENDDO
   ENDDO

ENDIF

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_RESTRICTION


!> ------------------------------------------------------------------------------------------------
!> Perform prolongation from coarser to finer grid 
!>    - 'VC' corresponds to coarser grid
!>    - 'VF' corresponds to finer   grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PROLONGATION (NVECTORC, NVECTORF, NLC, NLF)
INTEGER, INTENT(IN) :: NVECTORC, NVECTORF, NLC, NLF
INTEGER :: NM, ICOL, IC, I
REAL(EB), POINTER, DIMENSION(:) :: VC, VF
INTEGER , POINTER :: NXC, NYC, NZC, NCF, NCEF
INTEGER  :: NXF, NYF, NZF
INTEGER  :: IXF, IYF, IZF, ICF(8)=0, ICFB(-1:1,-1:1)=0
INTEGER  :: IXC, IYC, IZC, ICC
REAL(EB) :: AUX, SCAL
TYPE (SCARC_LEVEL_TYPE) , POINTER :: LC, LF
TYPE (SCARC_CELL_TYPE)  , POINTER :: CC, CF
TYPE (SCARC_MATRIX_TYPE), POINTER :: PF

CALL SCARC_ENTER_ROUTINE('SCARC_PROLONGATION')

!>
!> ------------------ Twolevel CG or Geometric Multigrid -------------------------------------------
!>
IF (BMULTILEVEL) THEN

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         LC => SCARC(NM)%LEVEL(NLC)
         LF => SCARC(NM)%LEVEL(NLF)

         CC => SCARC(NM)%LEVEL(NLC)%CELL
         CF => SCARC(NM)%LEVEL(NLF)%CELL

         NXC => LC%NX
         NYC => LC%NY
         NZC => LC%NZ

         NXF = 2*NXC
         NYF = 2*NYC
         NZF = 2*NZC

         VC => POINT_TO_VECTOR(NVECTORC, NM, NLC)
         VF => POINT_TO_VECTOR(NVECTORF, NM, NLF)

         IF (TWO_D) THEN

            SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)

               CASE (NSCARC_INTERPOL_CONSTANT)

                  DO IZC = 1, NZC
                     DO IXC = 1, NXC
      
                        IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.CC%STATE(IXC,1,IZC)/=NSCARC_DISCRET_GASPHASE) CYCLE
      
                        IXF = 2*IXC
                        IYF = 1
                        IZF = 2*IZC
      
                        ICC = CC%DOF(IXC, 1, IZC) 
      
                        ICF(1) = CF%DOF(IXF-1, 1, IZF-1) 
                        ICF(2) = CF%DOF(IXF-1, 1, IZF  ) 
                        ICF(3) = CF%DOF(IXF  , 1, IZF-1) 
                        ICF(4) = CF%DOF(IXF  , 1, IZF  ) 
      
                        DO I = 1, 4
                           VF(ICF(I)) = VC(ICC)
                        ENDDO
                     ENDDO
                  ENDDO

               CASE (NSCARC_INTERPOL_BILINEAR)

                  DO IZC = 1, NZC
                     DO IXC = 1, NXC
         
                        IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.CC%STATE(IXC,1,IZC)/=NSCARC_DISCRET_GASPHASE) CYCLE
         
                        IXF = 2*IXC
                        IZF = 2*IZC
         
                        ICC = CC%DOF(IXC, 1, IZC) 
         
                        ICFB(-1,-1) = CF%DOF(IXF-1, 1, IZF-1) 
                        ICFB(-1, 1) = CF%DOF(IXF-1, 1, IZF  ) 
                        ICFB( 1,-1) = CF%DOF(IXF  , 1, IZF-1) 
                        ICFB( 1, 1) = CF%DOF(IXF  , 1, IZF  ) 
         
                        IF (IXC==1.AND.IZC==1) THEN
                           VF(ICFB(-1,-1)) = VC(ICC)
                           VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+NXC))
                           VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                           VF(ICFB( 1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+NXC)+W1*VC(ICC+NXC+1))
                        ELSE IF (IXC==1 .AND. IZC==NZC) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-NXC))
                           VF(ICFB(-1, 1)) = VC(ICC)
                           VF(ICFB( 1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-NXC)+W1*VC(ICC-NXC+1))
                           VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                        ELSE IF (IXC==NXC .AND. IZC==1) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                           VF(ICFB(-1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+NXC)+W1*VC(ICC+NXC-1))
                           VF(ICFB( 1,-1)) = VC(ICC)
                           VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+NXC))
                        ELSE IF (IXC==NXC .AND. IZC==NZC) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-NXC)+W1*VC(ICC-NXC-1))
                           VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                           VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-NXC))
                           VF(ICFB( 1, 1)) = VC(ICC)
                        ELSE IF (IZC==1) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                           VF(ICFB(-1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+NXC)+W1*VC(ICC+NXC-1))
                           VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                           VF(ICFB( 1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+NXC)+W1*VC(ICC+NXC+1))
                        ELSE IF (IZC==NZC) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-NXC)+W1*VC(ICC-NXC-1))
                           VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                           VF(ICFB( 1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-NXC)+W1*VC(ICC-NXC+1))
                           VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                        ELSE IF (IXC==1) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-NXC))
                           VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+NXC))
                           VF(ICFB( 1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-NXC)+W1*VC(ICC-NXC+1))
                           VF(ICFB( 1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+NXC)+W1*VC(ICC+NXC+1))
                        ELSE IF (IXC==NXC) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-NXC)+W1*VC(ICC-NXC-1))
                           VF(ICFB(-1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+NXC)+W1*VC(ICC+NXC-1))
                           VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-NXC))
                           VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+NXC))
                        ELSE
                           VF(ICFB(-1,-1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-NXC)+W1*VC(ICC-NXC-1))
                           VF(ICFB(-1, 1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+NXC)+W1*VC(ICC+NXC-1))
                           VF(ICFB( 1,-1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-NXC)+W1*VC(ICC-NXC+1))
                           VF(ICFB( 1, 1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+NXC)+W1*VC(ICC+NXC+1))
                        ENDIF
                     ENDDO
                  ENDDO

            END SELECT SELECT_INTERPOL

         ELSE

            ! Note: 3D-bilinear case is still missing
            DO IZC = 1, NZC
               DO IYC = 1, NYC
                  DO IXC = 1, NXC

                     IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.CC%STATE(IXC,IYC,IZC)/=NSCARC_DISCRET_GASPHASE) CYCLE

                     IXF = 2*IXC
                     IYF = 2*IYC
                     IZF = 2*IZC

                     ICC = CC%DOF(IXC, IYC, IZC) 

                     ICF(1) = CF%DOF(IXF-1, IYF-1, IZF-1) 
                     ICF(2) = CF%DOF(IXF-1, IYF-1, IZF  ) 
                     ICF(3) = CF%DOF(IXF-1, IYF  , IZF-1) 
                     ICF(4) = CF%DOF(IXF-1, IYF  , IZF  ) 
                     ICF(5) = CF%DOF(IXF  , IYF-1, IZF-1) 
                     ICF(6) = CF%DOF(IXF  , IYF-1, IZF  ) 
                     ICF(7) = CF%DOF(IXF  , IYF  , IZF-1) 
                     ICF(8) = CF%DOF(IXF  , IYF  , IZF  ) 

                     DO I = 1, 8
                        VF(ICF(I)) = VC(ICC)
                     ENDDO

                  ENDDO
               ENDDO
            ENDDO

         ENDIF

      ENDDO

!>
!> ------------------ Algebraic multigrid -------------------------------------------
!>
ELSE IF (BAMG) THEN

   SCAL = 1.0_EB                                              ! why ?

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LC => SCARC(NM)%LEVEL(NLC)
      LF => SCARC(NM)%LEVEL(NLF)

      PF => SCARC(NM)%SYSTEM(NLF)%P

      VC => POINT_TO_VECTOR(NVECTORC, NM, NLC)
      VF => POINT_TO_VECTOR(NVECTORF, NM, NLF)

      NCF  => LF%NCS
      NCEF => LF%NCE

      DO IC = 1, NCF
         AUX = 0.0_EB
         DO ICOL = PF%ROW(IC), PF%ROW(IC+1)-1
            ICC = PF%COL(ICOL)
            AUX = VC(ICC) * PF%VAL(ICOL)
         ENDDO
         VF(IC) = SCAL*AUX
      ENDDO
      DO IC = NCF+1, NCEF
         VF(IC) = 0.0_EB
      ENDDO
   ENDDO

ENDIF
CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_PROLONGATION

!> ------------------------------------------------------------------------------------------------
!> Copy final solution from GMG (as preconditioner) to corresponding vector of CG (as main solver)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRECONDITIONER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   SCARC(NM)%SCOPE(NSCARC_SCOPE_ONE, NL)%G = SCARC(NM)%SCOPE(NSCARC_SCOPE_TWO, NL)%X
ENDDO
END SUBROUTINE SCARC_UPDATE_PRECONDITIONER

!> ------------------------------------------------------------------------------------------------
!> Finalize data
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRESSURE_MAINCELLS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, I, J, K, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
TYPE (SCARC_TYPE)      , POINTER :: S
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_SCOPE_TYPE), POINTER :: SC

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S  => SCARC(NM)
   L  => S%LEVEL(NL)
   SC => S%SCOPE(NSCARC_SCOPE_ONE, NL)

   IF (PREDICTOR) THEN
      HP => MESHES(NM)%H
   ELSE
      HP => MESHES(NM)%HS
   ENDIF

   DO K = 1, S%KBAR
      DO J = 1, S%JBAR
         DO I = 1, S%IBAR
            IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. L%CELL%STATE(I,J,K) /= NSCARC_DISCRET_GASPHASE) CYCLE
            IC = L%CELL%DOF(I,J,K)
            HP(I, J, K) = SC%X(IC)
         ENDDO
      ENDDO
   ENDDO

ENDDO

END SUBROUTINE SCARC_UPDATE_PRESSURE_MAINCELLS


!> ------------------------------------------------------------------------------------------------
!> Set correct boundary values at external and internal boundaries
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRESSURE_GHOSTCELLS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IXG, IYG, IZG, IXW, IYW, IZW
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_COORD_TYPE), POINTER :: C

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   L => SCARC(NM)%LEVEL(NL)
   C => SCARC(NM)%LEVEL(NL)%COORD

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

   !> compute ghost cell values
   WALL_CELLS_LOOP: DO IW = 1, L%N_WALL_CELLS_EXT

      IXG = L%WALL(IW)%IXG
      IYG = L%WALL(IW)%IYG
      IZG = L%WALL(IW)%IZG

      IXW = L%WALL(IW)%IXW
      IYW = L%WALL(IW)%IYW
      IZW = L%WALL(IW)%IZW

      IOR0 = L%WALL(IW)%IOR

      SELECT CASE (IOR0)
         CASE ( 1)
            IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXG,IYW,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXS(IYW,IZW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) - C%DX *M%BXS(IYW,IZW)
            ENDIF
         CASE (-1)
            IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXG,IYW,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXF(IYW,IZW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) + C%DX *M%BXF(IYW,IZW)
            ENDIF
         CASE ( 2)
            IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYS(IXW,IZW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) - C%DY *M%BYS(IXW,IZW)
            ENDIF
         CASE (-2)
            IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYF(IXW,IZW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) + C%DY *M%BYF(IXW,IZW)
            ENDIF
         CASE ( 3)
            IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZS(IXW,IYW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) - C%DZ *M%BZS(IXW,IYW)
            ENDIF
         CASE (-3)
         IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZF(IXW,IYW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) + C%DZ *M%BZF(IXW,IYW)
            ENDIF
      END SELECT
   ENDDO WALL_CELLS_LOOP

ENDDO

!> -----------------------------------------------------------------------------------------------
!> Perform data exchange to achieve consistency of ghost values along internal boundaries
!> -----------------------------------------------------------------------------------------------
CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_PRESSURE, NL)

END SUBROUTINE SCARC_UPDATE_PRESSURE_GHOSTCELLS


!> ------------------------------------------------------------------------------------------------
!>  Perform data exchange corresponding to requested exchange type (CALL receive and send-routines)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE (NTYPE, NL)
INTEGER, INTENT(IN):: NTYPE, NL

N_REQ = 0
TYPE_EXCHANGE = NTYPE

CALL SCARC_EXCHANGE_RECEIVE(NL)
CALL SCARC_EXCHANGE_SEND(NL)

END SUBROUTINE SCARC_EXCHANGE


!> ------------------------------------------------------------------------------------------------
!>  Receive data from neighbors (corresponds to POST_RECEIVES)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_RECEIVE (NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM
TYPE (SCARC_TYPE)         , POINTER ::  S
TYPE (OSCARC_TYPE)        , POINTER ::  OS
TYPE (SCARC_LEVEL_TYPE)   , POINTER ::  OL
TYPE (SCARC_EXCHANGE_TYPE), POINTER ::  OX

CALL SCARC_ENTER_ROUTINE('SCARC_EXCHANGE_RECEIVE')

RECEIVE_MESH_INDEX: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   RECEIVE_OMESH_INDEX: DO NOM = 1, NMESHES

      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_OMESH_INDEX

      S  => SCARC(NM)
      OS => SCARC(NM)%OSCARC(NOM)

      OX => OS%EXCHANGE
      IF (OX%NICMAX_S==0 .AND. OX%NICMAX_R==0) CYCLE RECEIVE_OMESH_INDEX

      OL => OS%LEVEL(NL)

      SELECT_EXCHANGE_TYPE: SELECT CASE (TYPE_EXCHANGE)


         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring step size along internal boundary
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OX%RECV_INT_BASIC(1),1,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring wall data
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALLINFO)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OX%RECV_INT(1),SIZE(OX%RECV_INT),MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring wall data
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_DISCRETIZATION)

            N_REQ = N_REQ+1
            OX%RECV_INT = 0
            CALL MPI_IRECV(OX%RECV_INT(1),SIZE(OX%RECV_INT),MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring grid dimensions
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MESHINFO)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OX%RECV_INT_BASIC(1),8,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring step size along internal boundary
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WIDTHINFO)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OX%RECV_REAL_BASIC(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange number of neighboring cells for AMG method (compact type only)
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_TRANSFER_SIZE)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OX%RECV_INT_BASIC(1),9,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange number of neighboring cells for AMG method (compact type only)
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OX%RECV_INT_BASIC(1),1,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange neighboring CELL_TYPE or CELL_INDEX
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_TYPE, NSCARC_EXCHANGE_CELL_INDEX)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OX%RECV_INT(1),SIZE(OX%RECV_INT),MPI_INTEGER,&
                           SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Perform exchanges for
         !>    - internal values for matrix-vector multiplication
         !>    - internal boundary values
         !>    - internal subdiagonal matrix values
         !>    - internal subdiagonal or ghost matrix values
         !>    - internal measure/CELL_TYPE values
         !> ---------------------------------------------------------------------------------------
         CASE DEFAULT

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OX%RECV_REAL(1),SIZE(OX%RECV_REAL),MPI_DOUBLE_PRECISION,&
                           SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

      END SELECT SELECT_EXCHANGE_TYPE
   ENDDO RECEIVE_OMESH_INDEX
ENDDO RECEIVE_MESH_INDEX
CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_EXCHANGE_RECEIVE


!> ------------------------------------------------------------------------------------------------
!> Send data to neighbors (corresponds to MESH_EXCHANGE)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_SEND (NL)
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, NOM, NCOL, NROW, NWL
INTEGER  :: IW, IWL, IWG, IWW, IW0, IROW, ICOL, IPTR, ICPL, ICELL_TYPE
INTEGER  :: IOR_NBR, IOR_OWN, IC, JC, ICC, JCC, ICE, ICO, ICN, ICG, ICW, IFACE
INTEGER  :: IX, IY, IZ
INTEGER  :: IXW, IYW, IZW
INTEGER  :: IXG, IYG, IZG
!INTEGER  :: IXN, IYN, IZN
INTEGER  :: I, J, K, LL
REAL(EB) :: ZSUM
INTEGER , POINTER, DIMENSION(:)     ::  RECV_INT, RECV_INT_BASIC
REAL(EB), POINTER, DIMENSION(:)     ::  RECV_REAL, RECV_REAL_BASIC
REAL(EB), POINTER, DIMENSION(:)     ::  VECTOR
REAL(EB), POINTER, DIMENSION(:,:,:) :: HVECTOR
TYPE (SCARC_TYPE)         , POINTER :: S
TYPE (OSCARC_TYPE)        , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE)   , POINTER :: L, OL
TYPE (SCARC_MATRIX_TYPE)  , POINTER :: A, P, R, OA, OP, OR
TYPE (SCARC_MAPPING_TYPE) , POINTER :: M, OM
TYPE (SCARC_EXCHANGE_TYPE), POINTER :: OX
TYPE (SCARC_COORD_TYPE)   , POINTER :: C, OC

CALL SCARC_ENTER_ROUTINE('SCARC_EXCHANGE_SEND')

!> ------------------------------------------------------------------------------------------------
!> Collect data for sending corresponding to requested exchange type
!> ------------------------------------------------------------------------------------------------
MESH_PACK_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S   => SCARC(NM)
   L   => SCARC(NM)%LEVEL(NL)
   M   => L%MAP
   C   => L%COORD

   OMESH_PACK_LOOP: DO NOM = 1, NMESHES

      OS => SCARC(NM)%OSCARC(NOM)

      OX => OS%EXCHANGE
      IF (OX%NICMAX_S == 0 .AND. OX%NICMAX_R == 0) CYCLE OMESH_PACK_LOOP

      OL => OS%LEVEL(NL)
      OC => OL%COORD
      OM => OL%MAP

      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      OMESH_PACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_BASIC: pack neighboring sizes for exchange
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC)

            OX%SEND_INT_BASIC(1)=OL%NWL

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OX%SEND_INT_BASIC(1),1,MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_WIDTHINFO: pack neighboring grid resolution information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WIDTHINFO)

            SELECT CASE(OL%IOR)
               CASE (1)
                  OC%DH = C%DXL(0)
               CASE (-1)
                  OC%DH = C%DXL(L%NX)
                CASE (2)
                  OC%DH = C%DYL(0)
               CASE (-2)
                  OC%DH = C%DYL(L%NY)
               CASE (3)
                  OC%DH = C%DZL(0)
               CASE (-3)
                  OC%DH = C%DZL(L%NZ)
            END SELECT
            OX%SEND_REAL_BASIC(1) = OC%DH

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL_BASIC(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_WALLINFO: pack neighboring wallinfo information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALLINFO)

            IPTR=1
            DO IWL = 1, OL%NWL
               IWG = OM%IWL_TO_IWG(IWL)
               OX%SEND_INT(IPTR   ) = L%WALL(IWG)%IXG
               OX%SEND_INT(IPTR+ 1) = L%WALL(IWG)%IYG
               OX%SEND_INT(IPTR+ 2) = L%WALL(IWG)%IZG
               OX%SEND_INT(IPTR+ 3) = L%WALL(IWG)%IXW
               OX%SEND_INT(IPTR+ 4) = L%WALL(IWG)%IYW
               OX%SEND_INT(IPTR+ 5) = L%WALL(IWG)%IZW
               OX%SEND_INT(IPTR+ 6) = L%WALL(IWG)%IXN(1)
               OX%SEND_INT(IPTR+ 7) = L%WALL(IWG)%IXN(2)
               OX%SEND_INT(IPTR+ 8) = L%WALL(IWG)%IYN(1)
               OX%SEND_INT(IPTR+ 9) = L%WALL(IWG)%IYN(2)
               OX%SEND_INT(IPTR+10) = L%WALL(IWG)%IZN(1)
               OX%SEND_INT(IPTR+11) = L%WALL(IWG)%IZN(2)
               OX%SEND_INT(IPTR+12) = L%WALL(IWG)%NOM
               IPTR = IPTR + 13
               DO ICPL=1,OL%NCPLS
                  OX%SEND_INT(IPTR)=L%WALL(IWG)%ICE(ICPL)
                  IPTR = IPTR + 1
               ENDDO
            ENDDO

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OX%SEND_INT(1),SIZE(OX%SEND_INT),MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_DISCRETIZATION: pack neighboring DISCRET information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_DISCRETIZATION)

            IPTR=1
            OX%SEND_INT=0
            DO IWL = 1, OL%NWL
               IWG = OM%IWL_TO_IWG(IWL)
               OX%SEND_INT(IPTR  ) = L%CELL%STATE(L%WALL(IWG)%IXW,L%WALL(IWG)%IYW,L%WALL(IWG)%IZW)
               OX%SEND_INT(IPTR+1) = L%CELL%DOF(L%WALL(IWG)%IXW,L%WALL(IWG)%IYW,L%WALL(IWG)%IZW)
               IPTR = IPTR + 2
            ENDDO

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OX%SEND_INT(1),SIZE(OX%SEND_INT),MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MESHINFO: pack neighboring mesh information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MESHINFO)

            OX%SEND_INT_BASIC(1)=L%NX
            OX%SEND_INT_BASIC(2)=L%NY
            OX%SEND_INT_BASIC(3)=L%NZ
            OX%SEND_INT_BASIC(4)=L%NC
            OX%SEND_INT_BASIC(5)=L%NCS
            OX%SEND_INT_BASIC(6)=L%NW
            OX%SEND_INT_BASIC(7)=L%N_WALL_CELLS_EXT
            OX%SEND_INT_BASIC(8)=L%N_WALL_CELLS_INT

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OX%SEND_INT_BASIC(1),8,MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_PRESSURE: pack overlapping parts of a H or HS
         ! (pressure vectors for PREDICTOR and CORRECTOR
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_PRESSURE)

            IF (PREDICTOR) THEN
               HVECTOR => POINT_TO_HVECTOR (NSCARC_VECTOR_H , NM)
            ELSE
               HVECTOR => POINT_TO_HVECTOR (NSCARC_VECTOR_HS, NM)
            ENDIF

            LL  = 1
            DO ICG=1, OL%NCG
               IWG = OM%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXW
               IY  = L%WALL(IWG)%IYW
               IZ  = L%WALL(IWG)%IZW
               OX%SEND_REAL(LL) = HVECTOR(IX,IY,IZ)
               LL = LL+1
            ENDDO

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL(1), SIZE(OX%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_VECTOR: pack overlapping parts of a given vector
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_VECTOR)

            VECTOR => POINT_TO_VECTOR(TYPE_VECTOR, NM, NL)

            LL = 1
            DO ICG= 1, OL%NCG
               !ZSUM = 0.0_EB
               IWG = OM%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXW
               IY  = L%WALL(IWG)%IYW
               IZ  = L%WALL(IWG)%IZW

               !ZSUM = VECTOR(L%CELL%DOF(IX, IY, IZ))
               !OX%SEND_REAL(LL) = ZSUM/REAL(OL%NCPLR,EB)

               IF (PRES_ON_WHOLE_DOMAIN.OR.L%CELL%STATE(IX, IY, IZ) == NSCARC_DISCRET_GASPHASE) THEN
                  OX%SEND_REAL(LL) = VECTOR(L%CELL%DOF(IX, IY, IZ))
               ELSE
                  OX%SEND_REAL(LL) = NSCARC_HUGE_REAL
               ENDIF

               LL = LL + 1
            ENDDO 

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL(1), SIZE(OX%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_CELL_INDEX: pack cell indices of communication partners
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_INDEX)

            OX%SEND_INT = 0
            LL = 1
            DO ICG= 1, OL%NCG
               IWG = OM%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXW
               IY  = L%WALL(IWG)%IYW
               IZ  = L%WALL(IWG)%IZW

               IF (PRES_ON_WHOLE_DOMAIN.OR.L%CELL%STATE(IX, IY, IZ) == NSCARC_DISCRET_GASPHASE) THEN
                  OX%SEND_INT(LL) = L%CELL%DOF(IX, IY, IZ)
               ELSE
                  OX%SEND_INT(LL) = NSCARC_HUGE_INT
               ENDIF

               LL = LL + 1
            ENDDO 

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_INT(1), SIZE(OX%SEND_INT), MPI_INTEGER, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_VALUE: pack cell indices of communication partners
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_VALUE)

            A => SCARC(NM)%SYSTEM(NL)%A

            LL = 1
            DO ICG= 1, OL%NCG
               IWG = OM%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXW
               IY  = L%WALL(IWG)%IYW
               IZ  = L%WALL(IWG)%IZW

               IF (PRES_ON_WHOLE_DOMAIN.OR.L%CELL%STATE(IX, IY, IZ) == NSCARC_DISCRET_GASPHASE) THEN
                  IC = L%CELL%DOF(IX, IY, IZ)
                  ICE = OM%ICG_TO_ICE(ICG)
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     IF (A%COL(ICOL) == ICE) THEN
                       OX%SEND_REAL(LL) = A%VAL(ICOL)
                       EXIT
                     ENDIF
                  ENDDO
               ELSE
                  OX%SEND_REAL(LL) = NSCARC_HUGE_REAL
               ENDIF

               LL = LL + 1
            ENDDO 

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL(1), SIZE(OX%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MEASURE_ADD: pack neighboring measure information (AMG only) with adding
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MEASURE_ADD)

            LL  = 1
            PACK_MEASURE_ADD: DO IWL=1, OL%NWL
               ICO = OM%IWL_TO_ICO(IWL)
               OX%SEND_REAL(LL) = L%AMG%MEASURE(ICO)
               LL = LL+1
            ENDDO PACK_MEASURE_ADD

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL(1), SIZE(OX%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> Send data along internal boundaries corresponding to requested exchange type
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_TYPE)

            LL  = 1
            PACK_CELL_TYPE: DO IWL=1, OL%NWL
               ICW = OM%IWL_TO_ICW(IWL)
               OX%SEND_INT(LL) = L%AMG%CTYPE(ICW)
               LL = LL+1
            ENDDO PACK_CELL_TYPE

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_INT(1), SIZE(OX%SEND_INT), MPI_INTEGER, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_SIZE: Send sizes of overlapping matrix areas
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            OA => OS%SYSTEM(NL)%A
            OX%SEND_INT_BASIC(1) = OA%NAV

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OX%SEND_INT_BASIC(1),1,MPI_INTEGER,SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_SUBDIAG: pack subdiagonal entries of system matrix on overlapping parts
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG)

            LL = 1
            IWW = 0

            WRITE(*,*) '!>!>! EXCHANGE_MATRIX_SUBDIAG, PACK: Achtung, hier nochmal checken  AMG !>!>'
            DO ICG=1,OL%NCG
               IOR_OWN = ABS(OL%WALL(ICG)%IOR)
               IWG = OM%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXG
               IY  = L%WALL(IWG)%IYG
               IZ  = L%WALL(IWG)%IZG
            ENDDO 

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL(1), SIZE(OX%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_STENCIL: pack matrix stencil information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_STENCIL)

            A => SCARC(NM)%SYSTEM(NL)%A

            LL = 1
            OX%SEND_REAL = 0.0_EB

            DO ICG=1, OL%NCG
  
               WRITE(*,*) '!>!>! ACHTUNG: EXCHANGE_MATRIX_STENCIL, PACK:  nochmal checken  AMG !>!>'
 
               IWG = OM%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXG
               IY  = L%WALL(IWG)%IYG
               IZ  = L%WALL(IWG)%IZG
               DO ICPL = 1, OL%NCPLR
                  IC = L%CELL%DOF(IX, IY, IZ)
                  NCOL = 0
                  DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
                     JC = A%COL(ICOL)
                     IF (JC > L%NC) THEN
                        IW0 = M%ICE_TO_IWG(JC)
                        IF (L%WALL(IW0)%NOM == NOM) NCOL = NCOL + 1
                     ELSE
                        NCOL = NCOL + 1
                     ENDIF
                  ENDDO
                  OX%SEND_REAL(LL) = REAL(NCOL)
                  LL = LL + 1
               ENDDO
            ENDDO 

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL(1), SIZE(OX%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_SYSTEM: pack overlapping parts of system matrix
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SYSTEM)

            A => SCARC(NM)%SYSTEM(NL)%A

            LL = 1
            OX%SEND_REAL = 0.0_EB

            WRITE(*,*) '!>!>! EXCHANGE_MATRIX_SYSTEM, PACK: Achtung, hier nochmal checken  AMG !>!>'
            !> Pack first cell layer
            PACK_MATRIX_SYSTEM: DO ICG=1, OL%NCG
               IF (OL%WALL(ICG)%NOM/=NM) CYCLE PACK_MATRIX_SYSTEM
               IWG = OM%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXG
               IY  = L%WALL(IWG)%IYG
               IZ  = L%WALL(IWG)%IZG

               DO ICPL = 1, OL%NCPLR
                  IC = L%CELL%DOF(IX, IY, IZ)

                  MATRIX_COLUMN_LOOP: DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
                     JC = A%COL(ICOL)
                     IF (JC > L%NC) THEN
                        IWL = M%ICE_TO_IWL(JC)
                        IF (OL%WALL(IWL)%NOM /= NM) CYCLE MATRIX_COLUMN_LOOP
                        OX%SEND_REAL(LL) = - REAL(M%ICE_TO_IWG(JC),EB)
                        OX%SEND_REAL(LL) = - REAL(JC)
                     ELSE
                        OX%SEND_REAL(LL) =   REAL(JC,EB)
                     ENDIF
                     OX%SEND_REAL(LL+1) = A%VAL(ICOL)
                     LL = LL + 2
                  ENDDO MATRIX_COLUMN_LOOP

               ENDDO

            ENDDO PACK_MATRIX_SYSTEM

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL(1), SIZE(OX%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_PROL: pack overlapping parts of prolongation matrix
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_PROL)

            LL = 1
            P => SCARC(NM)%SYSTEM(NL)%P

            WRITE(*,*) '!>!>! EXCHANGE_MATRIX_PROL, PACK: Achtung, hier nochmal checken  AMG !>!>'

            PACK_MATRIX_PROL: DO ICG=1, OL%NCG
               IF (OL%WALL(ICG)%NOM/=NM) CYCLE PACK_MATRIX_PROL
               IWG = OM%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXG
               IY  = L%WALL(IWG)%IYG
               IZ  = L%WALL(IWG)%IZG

               DO ICPL = 1, OL%NCPLR
                  IC = L%CELL%DOF(IX, IY, IZ)

                  OX%SEND_REAL(LL)   = REAL(L%AMG%CTYPE(IC),EB)
                  OX%SEND_REAL(LL+1) = REAL(P%ROW(IC+1)-P%ROW(IC),EB)

                  LL = LL + 2

                  DO ICOL = P%ROW(IC), P%ROW(IC+1)-1
                     JC = P%COL(ICOL)
                     IF (JC >= L%NC) THEN
                        OX%SEND_REAL(LL) = - REAL(JC,EB)
                     ELSE
                        OX%SEND_REAL(LL) =   REAL(JC,EB)
                     ENDIF
                     OX%SEND_REAL(LL+1) = P%VAL(ICOL)
                     LL = LL + 2
                  ENDDO

               ENDDO

            ENDDO PACK_MATRIX_PROL

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL(1), SIZE(OX%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_REST: pack overlapping parts of restriction matrix
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_REST)

            WRITE(*,*) '!>!>! EXCHANGE_MATRIX_REST, PACK: Achtung, hier nochmal checken  AMG !>!>'
            LL = 1
            R => SCARC(NM)%SYSTEM(NL)%R

            PACK_MATRIX_REST: DO ICG=1, OL%NCG
               IF (OL%WALL(IWG)%NOM/=NM) CYCLE PACK_MATRIX_REST
               
                  WRITE(*,*) 'ACHTUNG: MATRIX_REST: Wegen DISCRET schauen!>'

                  IX = L%WALL(IWG)%IXG
                  IY = L%WALL(IWG)%IYG
                  IZ = L%WALL(IWG)%IZG
                  IC = L%CELL%DOF(IX, IY, IZ)

                  ICC = L%AMG%CTYPE(IC)
                  IF (ICC > 0) THEN
                     OX%SEND_REAL(LL)   = REAL(ICC)
                     OX%SEND_REAL(LL+1) = REAL(R%ROW(ICC+1)-R%ROW(ICC),EB)

                     LL = LL + 2

                     DO ICOL = R%ROW(ICC), R%ROW(ICC+1)-1
                        JCC = R%COL(ICOL)
                        IF (JCC >= L%AMG%NCC) THEN
                           OX%SEND_REAL(LL) = - REAL(JCC,EB)
                        ELSE
                           OX%SEND_REAL(LL) =   REAL(JCC,EB)
                        ENDIF
                        OX%SEND_REAL(LL+1) = R%VAL(ICOL)
                        LL = LL + 2
                     ENDDO

                  ENDIF

               !ENDDO REST_FIRST_LAYER_LOOP

            ENDDO PACK_MATRIX_REST

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OX%SEND_REAL(1), SIZE(OX%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_TRANSER_SIZE: pack sizes of transfer matrices (AMG only)
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_TRANSFER_SIZE)

            OX%SEND_INT_BASIC(1) = L%NC
            OX%SEND_INT_BASIC(2) = L%NW
            OX%SEND_INT_BASIC(3) = L%NCE
            OX%SEND_INT_BASIC(4) = L%NW
            OX%SEND_INT_BASIC(5) = L%AMG%NCCI
            OX%SEND_INT_BASIC(6) = OL%AMG%NPS
            OX%SEND_INT_BASIC(7) = OL%AMG%NRS
            OX%SEND_INT_BASIC(8) = OL%AMG%NCCS
            OX%SEND_INT_BASIC(9) = OL%AMG%NCFS

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OX%SEND_INT_BASIC(1),9,MPI_INTEGER,SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF


      END SELECT OMESH_PACK_SELECT
   ENDDO OMESH_PACK_LOOP
ENDDO MESH_PACK_LOOP


!> ------------------------------------------------------------------------------------------------
!> Information from Mesh NM is received by Mesh NOM  (NOM receiver, NM sender)
!> ------------------------------------------------------------------------------------------------
IF (N_MPI_PROCESSES>1.AND.N_REQ/=0) &
   CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),MPI_STATUSES_IGNORE,IERROR)

!> ------------------------------------------------------------------------------------------------
!> Extract communication data from corresponding RECEIVE-buffers
!> ------------------------------------------------------------------------------------------------
MESH_UNPACK_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S   => SCARC(NM)
   L   => SCARC(NM)%LEVEL(NL)
   M   => L%MAP
   C   => L%COORD

   OMESH_UNPACK_LOOP: DO NOM=1,NMESHES

      SNODE  = PROCESS(NM)
      RNODE  = PROCESS(NOM)

      OS => SCARC(NM)%OSCARC(NOM)
      OX => OS%EXCHANGE

      OMESH_UNPACK_IF: IF (OX%NICMAX_S/=0 .AND. OX%NICMAX_R/=0) THEN

         OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
         OC => OL%COORD
         OM => OL%MAP

         IF (RNODE/=SNODE) THEN
            RECV_INT        => SCARC(NM)%OSCARC(NOM)%EXCHANGE%RECV_INT
            RECV_INT_BASIC  => SCARC(NM)%OSCARC(NOM)%EXCHANGE%RECV_INT_BASIC
            RECV_REAL       => SCARC(NM)%OSCARC(NOM)%EXCHANGE%RECV_REAL
            RECV_REAL_BASIC => SCARC(NM)%OSCARC(NOM)%EXCHANGE%RECV_REAL_BASIC
         ELSE
            RECV_INT        => SCARC(NOM)%OSCARC(NM)%EXCHANGE%SEND_INT
            RECV_INT_BASIC  => SCARC(NOM)%OSCARC(NM)%EXCHANGE%SEND_INT_BASIC
            RECV_REAL       => SCARC(NOM)%OSCARC(NM)%EXCHANGE%SEND_REAL
            RECV_REAL_BASIC => SCARC(NOM)%OSCARC(NM)%EXCHANGE%SEND_REAL_BASIC
         ENDIF

         OMESH_UNPACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_BASIC: unpack basic neighboring sizes for exchange
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_BASIC)

               OL%NCG  = RECV_INT_BASIC(1)

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WIDTHINFO: unpack neighboring step size information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_WIDTHINFO)

               OC%DH = RECV_REAL_BASIC(1)

               SELECT CASE (OL%IOR)
                  CASE ( 1)
                     C%DXL(0)    = 0.5_EB*(OC%DH + C%DXL(0))
                  CASE (-1)
                     C%DXL(L%NX) = 0.5_EB*(OC%DH + C%DXL(L%NX))
                  CASE ( 2)
                     C%DYL(0)    = 0.5_EB*(OC%DH + C%DYL(0))
                  CASE (-2)
                     C%DYL(L%NY) = 0.5_EB*(OC%DH + C%DYL(L%NY))
                  CASE ( 3)
                     C%DZL(0)    = 0.5_EB*(OC%DH + C%DZL(0))
                  CASE (-3)
                     C%DZL(L%NZ) = 0.5_EB*(OC%DH + C%DZL(L%NZ))
               END SELECT

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WALLINFO: unpack neighboring wall information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_WALLINFO)
               IPTR=1
               DO ICG = 1, OL%NCG
                  OL%WALL(ICG)%IXG    = RECV_INT(IPTR    )
                  OL%WALL(ICG)%IYG    = RECV_INT(IPTR + 1)
                  OL%WALL(ICG)%IZG    = RECV_INT(IPTR + 2)
                  OL%WALL(ICG)%IXW    = RECV_INT(IPTR + 3)
                  OL%WALL(ICG)%IYW    = RECV_INT(IPTR + 4)
                  OL%WALL(ICG)%IZW    = RECV_INT(IPTR + 5)
                  OL%WALL(ICG)%IXN(1) = RECV_INT(IPTR + 6)
                  OL%WALL(ICG)%IXN(2) = RECV_INT(IPTR + 7)
                  OL%WALL(ICG)%IYN(1) = RECV_INT(IPTR + 8)
                  OL%WALL(ICG)%IYN(2) = RECV_INT(IPTR + 9)
                  OL%WALL(ICG)%IZN(1) = RECV_INT(IPTR +10)
                  OL%WALL(ICG)%IZN(2) = RECV_INT(IPTR +11)
                  OL%WALL(ICG)%NOM    = RECV_INT(IPTR +12)
                  IPTR = IPTR + 13
                  ALLOCATE (OL%WALL(ICG)%ICE(OL%NCPLR))
                  DO ICPL=1,OL%NCPLR
                     OL%WALL(ICG)%ICE(ICPL) = RECV_INT(IPTR)
                     IPTR = IPTR + 1
                  ENDDO
               ENDDO

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_DISCRETIZATION: unpack neighboring DISCRET information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_DISCRETIZATION)
               IPTR=1
               DO ICG = 1, OL%NCG

                  OL%WALL(ICG)%STATE = RECV_INT(IPTR  )
                  OL%WALL(ICG)%DOF   = RECV_INT(IPTR+1)
                  IPTR = IPTR + 2

                  IWG = OM%ICG_TO_IWG(ICG)

                  IXG = L%WALL(IWG)%IXG
                  IYG = L%WALL(IWG)%IYG
                  IZG = L%WALL(IWG)%IZG

                  IXW= L%WALL(IWG)%IXW
                  IYW= L%WALL(IWG)%IYW
                  IZW= L%WALL(IWG)%IZW

                  !IF (OL%WALL(ICG)%STATE == NSCARC_DISCRET_GASPHASE.AND.L%CELL%STATE(IXW,IYW,IZW)/=NSCARC_DISCRET_SOLID) THEN
                  IF (OL%WALL(ICG)%STATE == NSCARC_DISCRET_GASPHASE) THEN
                     L%CELL%STATE(IXG, IYG, IZG) = OL%WALL(ICG)%STATE
                     L%CELL%DOF(IXG, IYG, IZG)   = OL%WALL(ICG)%DOF
                  ENDIF

               ENDDO


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WALLINFO: unpack neighboring mesh information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MESHINFO)

               OL%NX  = RECV_INT_BASIC(1)
               OL%NY  = RECV_INT_BASIC(2)
               OL%NZ  = RECV_INT_BASIC(3)
               OL%NC  = RECV_INT_BASIC(4)
               OL%NCS = RECV_INT_BASIC(5)
               OL%NW  = RECV_INT_BASIC(6)
               OL%N_WALL_CELLS_EXT = RECV_INT_BASIC(7)
               OL%N_WALL_CELLS_INT = RECV_INT_BASIC(8)

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WALLINFO: unpack neighboring matrix size information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

               OA => OS%SYSTEM(NL)%A
               OA%NAV = RECV_INT_BASIC(1)

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WALLINFO: unpack neighboring transfer size information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_TRANSFER_SIZE)

               OL%NC  =  RECV_INT_BASIC(1)
               OL%NW  =  RECV_INT_BASIC(2)
               OL%NCE =  RECV_INT_BASIC(3)
               OL%NW  =  RECV_INT_BASIC(4)

               OL%AMG%NCCI =  RECV_INT_BASIC(5)
               OL%AMG%NP   =  RECV_INT_BASIC(6)
               OL%AMG%NR   =  RECV_INT_BASIC(7)
               OL%AMG%NCC  =  RECV_INT_BASIC(8)
               OL%AMG%NCF  =  RECV_INT_BASIC(9)

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_PRESSURE: unpack overlapping parts of a H or HS
            ! (pressure vectors for PREDICTOR and CORRECTOR
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_PRESSURE)

               IF (PREDICTOR) THEN
                  HVECTOR => POINT_TO_HVECTOR (NSCARC_VECTOR_H, NM)
               ELSE
                  HVECTOR => POINT_TO_HVECTOR (NSCARC_VECTOR_HS, NM)
               ENDIF

               LL = 1
               DO IFACE = 1, 6
                  IOR_NBR = FACE_ORIENTATION(IFACE)
                  IWL = OL%SUBDIVISION(1, -IOR_NBR)
                  NWL = OL%SUBDIVISION(2, -IOR_NBR)

                  UNPACK_PRESSURE: DO IW = IWL, IWL + NWL - 1
   
                     ZSUM=0.0_EB
                     IWG = OM%IWL_TO_IWG(IW)
                     DO ICPL = 1, OL%NCPL
                        ZSUM=ZSUM+RECV_REAL(LL)
                        LL = LL+1
                     ENDDO

                     I=L%WALL(IWG)%IXG
                     J=L%WALL(IWG)%IYG
                     K=L%WALL(IWG)%IZG
   
                     HVECTOR(I, J, K) = ZSUM/REAL(OL%NCPL,EB)
   
                  ENDDO UNPACK_PRESSURE
               ENDDO

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_VECTOR: unpack overlapping parts of a given vector
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_VECTOR)

               VECTOR => POINT_TO_VECTOR (TYPE_VECTOR, NM, NL)

               LL = 1
               DO IFACE = 1, 6
                  IOR_NBR = FACE_ORIENTATION(IFACE)
                  IWL = OL%SUBDIVISION(1, -IOR_NBR)
                  NWL = OL%SUBDIVISION(2, -IOR_NBR)

                  DO IW = IWL, IWL + NWL - 1
                     IWG = OM%IWL_TO_IWG(IW)
   
                     ICE     = L%WALL(IWG)%ICE(1)
                     IOR_OWN = L%WALL(IWG)%IOR

                     VECTOR(ICE) = RECV_REAL(LL)
   
                     LL = LL + 1
                  ENDDO 
               ENDDO

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_CELL_INDEX: unpack cell indices of communication partners
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_CELL_INDEX)

               LL = 1
               DO IFACE = 1, 6
                  IOR_NBR = FACE_ORIENTATION(IFACE)
                  IWL = OL%SUBDIVISION(1, -IOR_NBR)
                  NWL = OL%SUBDIVISION(2, -IOR_NBR)

                  DO IW = IWL, IWL + NWL - 1
                     IWG = OM%IWL_TO_IWG(IW)
                     ICE     = L%WALL(IWG)%ICE(1)
                     IOR_OWN = L%WALL(IWG)%IOR

                     L%MAP%ICE_TO_ICN(ICE) = RECV_INT(LL)
                     LL = LL + 1
                  ENDDO 
               ENDDO

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_VALUE: unpack offdiagonal entry of matrix for condensed system
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_VALUE)

               LL = 1
               DO IFACE = 1, 6
                  IOR_NBR = FACE_ORIENTATION(IFACE)
                  IWL = OL%SUBDIVISION(1, -IOR_NBR)
                  NWL = OL%SUBDIVISION(2, -IOR_NBR)

                  DO IW = IWL, IWL + NWL - 1
                     IWG = OM%IWL_TO_IWG(IW)
                     ICE     = L%WALL(IWG)%ICE(1)
                     IOR_OWN = L%WALL(IWG)%IOR

                     L%MAP%ICE_TO_VAL(ICE) = RECV_REAL(LL)

                     LL = LL + 1
                  ENDDO 
               ENDDO


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MEASURE_ADD: unpack neighboring measure information (AMG only) with adding
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MEASURE_ADD)

               LL = 1
               DO IWL = 1, OL%NWL
                  ICW = OM%IWL_TO_ICW(IWL)
                  L%AMG%MEASURE(ICW) = L%AMG%MEASURE(ICW)  + RECV_REAL(LL)
                  LL = LL + 1
               ENDDO 

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_CELL_TYPE: unpack neighboring CELL_TYPE information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_CELL_TYPE)

               LL = 1
               DO IWL = 1, OL%NWL
                  DO ICPL = 1, OL%NCPL
                     ICG = OM%IWL_TO_ICG(IWL, ICPL)
                     OL%AMG%CTYPE(ICG) = RECV_INT(LL)
                     LL = LL+1
                  ENDDO
               ENDDO 


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_SUBDIAG: unpack subdiagonal information from !neighboring system matrix
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG)

               A => SCARC(NM)%SYSTEM(NL)%A

               LL = 1
               UNPACK_MATRIX_SUBDIAG: DO IWL = 1, OL%NWL

                  ZSUM=0.0_EB
                  ! ACHTUNG WEGEN NCPL SCHAUEN 
                  IWG = OM%IWL_TO_IWG(IWL)
                  IX = L%WALL(IWG)%IXG
                  IY = L%WALL(IWG)%IYG
                  IZ = L%WALL(IWG)%IZG
                  IC = L%CELL%DOF(IX, IY, IZ)
                  ZSUM=ZSUM+RECV_REAL(LL)
                  LL = LL+1

                  IC=L%WALL(IWG)%ICW

                  IROW = A%ROW(IC)
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     IF (IW == -A%COL(ICOL)) THEN
                        A%VAL(ICOL)     = ZSUM/REAL(L%WALL(IW)%NCPL,EB)
                        A%COL(ICOL) = L%WALL(IW)%ICE(1)
                     ENDIF
                  ENDDO

               ENDDO UNPACK_MATRIX_SUBDIAG

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_STENCIL: unpack stencil information from neighboring system matrix
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_STENCIL)

               OA => OS%SYSTEM(NL)%A

               LL = 1
               UNPACK_MATRIX_STENCIL: DO IWL = 1, OL%NWL
                  !IWG = OM%IWL_TO_IWG(IWL)
                  !ICO = OM%IWL_TO_ICO(IWL)
                  !ICO = L%WALL(IWG)%ICO

                  OA%STENCIL(IWL) = NINT(RECV_REAL(LL))
                  LL = LL + 1

                ENDDO UNPACK_MATRIX_STENCIL

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_SYSTEM: extract neighboring system matrix information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SYSTEM)

               A => SCARC(NM)%SYSTEM(NL)%A

               LL = 1
               UNPACK_MATRIX_SYSTEM: DO IWL = 1, OL%NWL

                  IWG = OM%IWL_TO_IWG(IWL)

                  ICE = L%WALL(IWG)%ICE(1)
                  IX  = L%WALL(IWG)%IXG
                  IY  = L%WALL(IWG)%IYG
                  IZ  = L%WALL(IWG)%IZG
                  ICN = L%CELL%DOF(IX, IY, IZ)
                  ICG = L%WALL(IWG)%ICG(1)
                  ICW = L%WALL(IWG)%ICW

                  DO ICOL = A%ROW(ICE), A%ROW(ICE+1)-1
                     A%COL(ICOL)= NINT(RECV_REAL(LL))
                     A%VAL(ICOL)    = RECV_REAL(LL+1)
                     LL = LL + 2
                  ENDDO

                ENDDO UNPACK_MATRIX_SYSTEM


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_PROL: Extract neighboring prolongation matrix information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_PROL)

               OP => OS%SYSTEM(NL)%P

               LL = 1
               NROW=1
               ICC = 1
               UNPACK_MATRIX_PROL: DO IWL = 1, OL%NWL

                  IWG = OM%IWL_TO_IWG(IWL)

                  ICG = L%WALL(IWG)%ICG(1)
                  ICE = L%WALL(IWG)%ICE(1)
                  IX  = L%WALL(IWG)%IXG
                  IY  = L%WALL(IWG)%IYG
                  IZ  = L%WALL(IWG)%IZG
                  ICN = L%CELL%DOF(IX, IY, IZ)
                  ICELL_TYPE = NINT(RECV_REAL(LL))
                  NCOL       = NINT(RECV_REAL(LL+1))

                  OL%AMG%CTYPE(ICG) = ICELL_TYPE
                  OP%ROW(ICG) = NROW

                  LL = LL + 2

                  DO ICOL = NROW, NROW+NCOL-1
                     OP%COL(ICOL) = NINT(RECV_REAL(LL))
                     OP%VAL(ICOL) = RECV_REAL(LL+1)
                     LL = LL + 2
                  ENDDO
                  NROW = NROW + NCOL
                  OP%ROW(ICG+1) = NROW

               ENDDO UNPACK_MATRIX_PROL


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_REST: exchange neighboring restriction matrix information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_REST)

               OR => OS%SYSTEM(NL)%R

               LL = 1
               NROW=1
               ICC = 1
               UNPACK_MATRIX_REST: DO

                  IW = NINT(RECV_REAL(LL))
                  LL = LL + 1
                  IF (IW==-888.OR.IW==-999) EXIT UNPACK_MATRIX_REST

                  WRITE(*,*) 'ACHTUNG, MUSS WIEDER AUSGEARBEITET WERDEN!'

                  !ICG = L%WALL(IW)%ICG(1)
                  !ICE = L%WALL(IW)%ICE(1)
                  !ICN = L%WALL(IW)%ICN(1)
                  !ICELL_TYPE = NINT(RECV_REAL(LL))
                  !NCOL      = NINT(RECV_REAL(LL+1))
                 !>
                 !> OL%AMG%CTYPE(ICG) = ICELL_TYPE
                 !> OP%ROW(ICG) = NROW
         !>
         !>        LL = LL + 2
         !>
         !>        DO ICOL = NROW, NROW+NCOL-1
         !>           OP%COL(ICOL)= NINT(RECV_REAL(LL))
         !>           OP%VAL(ICOL)    = RECV_REAL(LL+1)
         !>           LL = LL + 2
         !>        ENDDO
         !>        NROW = NROW + NCOL
         !>        OP%ROW(ICG+1) = NROW

               ENDDO UNPACK_MATRIX_REST

          END SELECT OMESH_UNPACK_SELECT
      ENDIF OMESH_UNPACK_IF
   ENDDO OMESH_UNPACK_LOOP
ENDDO MESH_UNPACK_LOOP

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_EXCHANGE_SEND


!> ------------------------------------------------------------------------------------------------
!> Check if difference of two values is less than a given tolerance
!> ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB), INTENT(IN) :: VAL1, VAL2
REAL (EB) :: TOL
TOL = 1.0E-10_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
RETURN
END FUNCTION MATCH


!> ------------------------------------------------------------------------------------------------
!> Print out timings for ScaRC - not updated at the moment
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_TIMINGS

IF (.NOT.BVERBOSE_MUCH) RETURN

IF (MYID == 0) THEN
   IF (N_MPI_PROCESSES == 1) THEN
      WRITE(LU_OUTPUT,1002) 
      WRITE(LU_OUTPUT,1001) 
      WRITE(LU_OUTPUT,1002) 
   ELSE
      WRITE(LU_OUTPUT,2002) 
      WRITE(LU_OUTPUT,2001) 
      WRITE(LU_OUTPUT,2002) 
   ENDIF
ENDIF

IF (MYID == 0) THEN
   IF (N_MPI_PROCESSES == 1) THEN
      WRITE(LU_OUTPUT,1002) 
   ELSE
      WRITE(LU_OUTPUT,2002) 
   ENDIF
ENDIF

IF (.NOT.BDEBUG_NONE) CLOSE(MSG%LU_DEBUG)

1001 FORMAT('| Scarc routine       |    Time (s)    |')
1002 FORMAT('|---------------------|----------------|')
2001 FORMAT('| Scarc routine       |  Time_min (s)  |  Time_max (s)  | Time_mean (s)  |')
2002 FORMAT('|---------------------|----------------|----------------|----------------|')
END SUBROUTINE SCARC_TIMINGS


!> -----------------------------------------------------------------------------
!> Setup DISCRET-array
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DISCRETIZATION

INTEGER :: NM
INTEGER :: I, J, K
INTEGER, PARAMETER :: IMPADD = 1
INTEGER, PARAMETER :: SHFTM(1:3,1:6) = RESHAPE((/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_DISCRETIZATION')
!>
!> Initialize all cells as GASPHASE cells
!>
MESHES_LOOP1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NLEVEL_MIN)

   CALL SCARC_ALLOCATE_INT3(L%CELL%DOF,   0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_NONE, 'DOF')
   CALL SCARC_ALLOCATE_INT3(L%CELL%STATE, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_NONE, 'TYP')
 
   L%CELL%DOF   = NSCARC_UNDEFINED_INT
   L%CELL%STATE = NSCARC_UNDEFINED_INT
   L%CELL%STATE(1:L%NX, 1:L%NY, 1:L%NZ) = NSCARC_DISCRET_GASPHASE 

ENDDO MESHES_LOOP1

!>
!> Identify and mark OBST SOLID cells in CGSC-part of DISCRET 
!>
MESHES_LOOP2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   L => SCARC(NM)%LEVEL(NLEVEL_MIN)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (M%SOLID(M%CELL_INDEX(I, J, K))) L%CELL%STATE(I, J, K) = NSCARC_DISCRET_SOLID
         ENDDO
      ENDDO
   ENDDO
ENDDO MESHES_LOOP2

!>
!> Define local cell numbers for Poisson equation
!>
MESHES_LOOP3 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NLEVEL_MIN)

   !> consider gasphase and solid cells
   IF (PRES_ON_WHOLE_DOMAIN) THEN 
      DO K=1,L%NZ
         DO J=1,L%NY
            DO I=1,L%NX
               L%CELL%NC_LOCAL(NM) = L%CELL%NC_LOCAL(NM) + 1
               L%CELL%DOF(I,J,K)   = L%CELL%NC_LOCAL(NM)
            ENDDO
         ENDDO
      ENDDO

   !> consider only gasphase cells
   ELSE 
      DO K=1,L%NZ
         DO J=1,L%NY
            DO I=1,L%NX
               IF (L%CELL%STATE(I,J,K) == NSCARC_DISCRET_GASPHASE ) THEN
                  L%CELL%NC_LOCAL(NM) = L%CELL%NC_LOCAL(NM) + 1
                  L%CELL%DOF(I,J,K)   = L%CELL%NC_LOCAL(NM)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

   ENDIF 

   L%NCS = L%CELL%NC_LOCAL(NM)

ENDDO MESHES_LOOP3

CALL SCARC_LEAVE_ROUTINE()
END SUBROUTINE SCARC_SETUP_DISCRETIZATION


!> -----------------------------------------------------------------------------
!> MV: Setup DISCRET on coarser levels
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DISCRETIZATION_LEVEL(NL)

INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IXF, IYF, IZF, IX, IY, IZ, NSTEP
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF, LC
TYPE (SCARC_CELL_TYPE) , POINTER :: CF, CC

CALL SCARC_ENTER_ROUTINE('SCARC_SETUP_DISCRETIZATION_LEVEL')

MESHES_LOOP1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   LF => SCARC(NM)%LEVEL(NLEVEL_MIN)                  !> fine level
   LC => SCARC(NM)%LEVEL(NL)                          !> coarse level

   CF => SCARC(NM)%LEVEL(NLEVEL_MIN)%CELL             !> fine level
   CC => SCARC(NM)%LEVEL(NL)%CELL                     !> coarse level

   CALL SCARC_ALLOCATE_INT3(CC%DOF,   0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_NONE, 'DOF')
   CALL SCARC_ALLOCATE_INT3(CC%STATE, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_NONE, 'TYP')
 
   CC%DOF   = NSCARC_UNDEFINED_INT
   CC%STATE = NSCARC_UNDEFINED_INT
   CC%STATE(1:LC%NX, 1:LC%NY, 1:LC%NZ) = NSCARC_DISCRET_GASPHASE 

   NSTEP = 2**(NL - NLEVEL_MIN)

   !> In case of structured approach:
   IF (PRES_ON_WHOLE_DOMAIN) THEN
      DO IZ = 1, LC%NZ
         IZF = (IZ-1)*NSTEP + 1
         DO IY = 1, LC%NY
            IYF = (IY-1)*NSTEP + 1
            DO IX = 1, LC%NX
               IXF = (IX-1)*NSTEP + 1
               CC%STATE(IX,IY,IZ) = CF%STATE(IXF, IYF, IZF)
               CC%NC_LOCAL(NM)  = CC%NC_LOCAL(NM) + 1
               CC%DOF(IX,IY,IZ) = CC%NC_LOCAL(NM)
            ENDDO
         ENDDO
      ENDDO

   !> In case of unstructured approach:
   ELSE
      DO IZ = 1, LC%NZ
         IZF = (IZ-1)*NSTEP + 1
         DO IY = 1, LC%NY
            IYF = (IY-1)*NSTEP + 1
            DO IX = 1, LC%NX
               IXF = (IX-1)*NSTEP + 1
               CC%STATE(IX,IY,IZ) = CF%STATE(IXF, IYF, IZF)
               IF (CF%STATE(IXF, IYF, IZF) == NSCARC_DISCRET_GASPHASE) THEN
                  CC%NC_LOCAL(NM)  = CC%NC_LOCAL(NM) + 1
                  CC%DOF(IX,IY,IZ) = CC%NC_LOCAL(NM)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   LC%NCS = CC%NC_LOCAL(NM)

ENDDO MESHES_LOOP1

CALL SCARC_LEAVE_ROUTINE()

END SUBROUTINE SCARC_SETUP_DISCRETIZATION_LEVEL


!> ----------------------------------------------------------------------------------------------------
!> Filter out mean value
!> --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_FILTER_MEANVALUE(NVECTOR, NL)
INTEGER, INTENT(IN) :: NVECTOR, NL
INTEGER :: NM, IC, I, J, K
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_SCOPE_TYPE), POINTER :: SC 
REAL(EB), DIMENSION(:) , POINTER :: VC 


LOCAL_REAL = 0.0_EB
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L  => SCARC(NM)%LEVEL(NL)
   VC => POINT_TO_VECTOR(NVECTOR, NM, NL)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%CELL%STATE(I,J,K) /= NSCARC_DISCRET_GASPHASE) CYCLE
            IC = L%CELL%DOF(I,J,K)
            LOCAL_REAL(NM) = LOCAL_REAL(NM) + VC(IC)
         ENDDO
      ENDDO
   ENDDO
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,LOCAL_REAL,COUNTS,DISPLS,&
                       MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

GLOBAL_REAL = SUM(LOCAL_REAL(1:NMESHES))/REAL(N_CELLS_GLOBAL(NL))

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L  => SCARC(NM)%LEVEL(NL)
   SC => SCARC(NM)%SCOPE(1,NL)
   VC => POINT_TO_VECTOR(NVECTOR, NM, NL)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%CELL%STATE(I,J,K) /= NSCARC_DISCRET_GASPHASE) CYCLE
            IC = L%CELL%DOF(I,J,K)
            VC(IC) = VC(IC) - GLOBAL_REAL
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_FILTER_MEANVALUE


!> ------------------------------------------------------------------------------------------------
!> Preset right hand side in such a way that exact solution is known
!> temporarily not used
!> ------------------------------------------------------------------------------------------------
!SUBROUTINE SCARC_PRESET_EXACT (E, NL)
!INTEGER, INTENT(IN):: E, NL
!REAL (EB), POINTER, DIMENSION(:) :: VE
!INTEGER :: IC, NM, I, K
!REAL(EB), DIMENSION(:), POINTER :: XMID, ZMID
!TYPE (MESH_TYPE), POINTER :: M
!TYPE (SCARC_LEVEL_TYPE), POINTER :: L
!
!IF (.NOT.BVERBOSE_MUCH) RETURN
!
!DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!   M => MESHES(NM)
!   L => SCARC(NM)%LEVEL(NL)
!   VE => POINT_TO_VECTOR (E, NM, NL)
!   DO K = 1, L%NZ
!      DO I = 1, L%NX
!         IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%CELL%STATE(I,1,K) /= NSCARC_DISCRET_GASPHASE) CYCLE
!         IC = L%CELL%DOF(I,1,K)
!         IF (NL == NLEVEL_MIN) THEN
!            XMID => M%XC
!            ZMID => M%ZC
!         ELSE
!            XMID => L%COORD%XMID
!            ZMID => L%COORD%ZMID
!         ENDIF
!         VE(IC) = EXACT(XMID(I),ZMID(K))
!      ENDDO
!   ENDDO
!ENDDO
!
!END SUBROUTINE SCARC_PRESET_EXACT

!> ------------------------------------------------------------------------------------------------
!> Restore last cell of last mesh
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_LAST_CELL (X, NL)
INTEGER, INTENT(IN):: X, NL
REAL (EB), POINTER, DIMENSION(:) :: VX
TYPE (SCARC_TYPE), POINTER :: S

IF (UPPER_MESH_INDEX /= NMESHES .OR. TYPE_RELAX == NSCARC_RELAX_FFT) RETURN
S  => SCARC(UPPER_MESH_INDEX)
VX => POINT_TO_VECTOR (X, UPPER_MESH_INDEX, NL)
VX(S%N_CELLS) = S%RHS_END

END SUBROUTINE SCARC_RESTORE_LAST_CELL


!> ------------------------------------------------------------------------------------------------
!> Preset right hand side in such a way that exact solution is known
!> temporarily not used
!> ------------------------------------------------------------------------------------------------
!SUBROUTINE SCARC_PRESET_RHS (NVECTOR, NL)
!INTEGER, INTENT(IN):: NVECTOR, NL
!REAL (EB), POINTER, DIMENSION(:) :: VC
!INTEGER :: IC, NM, I, K
!REAL (EB) :: X, Z
!TYPE (MESH_TYPE), POINTER :: M
!TYPE (SCARC_LEVEL_TYPE), POINTER :: L
!
!IF (.NOT.BVERBOSE_MUCH) RETURN
!IF (ITE_TOTAL == 0) WRITE(*,*) 'ACHTUNG: PRESET_RHS AKTIV !!!'
!IF (NL > NLEVEL_MIN) CALL SCARC_SHUTDOWN('Wrong level for presetting RHS ', 'NONE', NL)
!
!DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!   M => MESHES(NM)
!   M%BXS = 0.0_EB
!   M%BXF = 0.0_EB
!   M%BZS = 0.0_EB
!   M%BZF = 0.0_EB
!   L => SCARC(NM)%LEVEL(NL)
!   VC => POINT_TO_VECTOR (NVECTOR, NM, NL)
!   DO K = 1, L%NZ
!      DO I = 1, L%NX
!         IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%CELL%STATE(I,1,K) /= NSCARC_DISCRET_GASPHASE) CYCLE
!         IC = L%CELL%DOF(I,1,K)
!         X  = M%XC(I)
!         Z  = M%ZC(K)
!         VC(IC) = RHS(X,Z)
!      ENDDO
!   ENDDO
!ENDDO
!
!END SUBROUTINE SCARC_PRESET_RHS
!

!> ------------------------------------------------------------------------------------------------
!> Set exact solution 
!> ------------------------------------------------------------------------------------------------
!DOUBLE PRECISION FUNCTION EXACT(X,Z)
!REAL (EB), INTENT(IN) :: X, Z
!EXACT = (X**2 - X**4) * (Z**4 - Z**2)                                    !> FUNCTION 1
!!EXACT = (X**2 - 1) * (Z**2 - 1)                                         !> FUNCTION 2
!!EXACT = - 625.0_EB/16.0_EB * X * (0.8_EB - X) * Z * (0.8_EB - Z)        !> FUNCTION 3
!END FUNCTION EXACT


!> ------------------------------------------------------------------------------------------------
!> Set right hand side
!> ------------------------------------------------------------------------------------------------
!DOUBLE PRECISION FUNCTION RHS(X,Z)
!REAL (EB), INTENT(IN) :: X, Z
!RHS = 2.0_EB*((1.0_EB - 6.0_EB*X**2)*Z**2*(1.0_EB-Z**2)+(1.0_EB-6.0_EB*Z**2)*X**2*(1.0_EB-X**2))
!!RHS = -X**2 - Z**2 +2
!!RHS = 625.0_EB/8.0_EB * (X * (0.8_EB - X) + Z * (0.8_EB - Z))
!END FUNCTION RHS


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT1(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR (MSG%HISTORY(MSG%NCURRENT), CTEXT, IERROR)
   SELECT CASE (NINIT) 
      CASE (NSCARC_INIT_UNDEFINED)
         WORKSPACE = NSCARC_UNDEFINED_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_SHUTDOWN('Inconsistent length for vector allocation ', CTEXT, -999)
ENDIF
END SUBROUTINE SCARC_ALLOCATE_INT1

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 2
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT2(WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR (MSG%HISTORY(MSG%NCURRENT), CTEXT, IERROR)
   SELECT CASE (NINIT) 
      CASE (NSCARC_INIT_UNDEFINED)
         WORKSPACE = NSCARC_UNDEFINED_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
       SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_SHUTDOWN('Inconsistent length for vector allocation ', CTEXT, -999)
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_INT2

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 3
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR (MSG%HISTORY(MSG%NCURRENT), CTEXT, IERROR)
   SELECT CASE (NINIT) 
      CASE (NSCARC_INIT_UNDEFINED)
         WORKSPACE = NSCARC_UNDEFINED_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
       SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
       SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_SHUTDOWN('Inconsistent length for vector allocation ', CTEXT, -999)
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_INT3

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL1(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR (MSG%HISTORY(MSG%NCURRENT), CTEXT, IERROR)
   SELECT CASE (NINIT) 
      CASE (NSCARC_INIT_UNDEFINED)
         WORKSPACE = NSCARC_UNDEFINED_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_SHUTDOWN('Inconsistent length for vector allocation ', CTEXT, -999)
ENDIF
END SUBROUTINE SCARC_ALLOCATE_REAL1


#ifdef WITH_MKL_FB
!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL1_FB(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(FB), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR (MSG%HISTORY(MSG%NCURRENT), CTEXT, IERROR)
   SELECT CASE (NINIT) 
      CASE (NSCARC_INIT_UNDEFINED)
         WORKSPACE = NSCARC_UNDEFINED_REAL_FB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_FB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_FB
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_SHUTDOWN('Inconsistent length for vector allocation ', CTEXT, -999)
ENDIF
END SUBROUTINE SCARC_ALLOCATE_REAL1_FB
#endif

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 2
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL2(WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR (MSG%HISTORY(MSG%NCURRENT), CTEXT, IERROR)
   SELECT CASE (NINIT) 
      CASE (NSCARC_INIT_UNDEFINED)
         WORKSPACE = NSCARC_UNDEFINED_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
       SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_SHUTDOWN('Inconsistent length for vector allocation ', CTEXT, -999)
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_REAL2

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 3
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR (MSG%HISTORY(MSG%NCURRENT), CTEXT, IERROR)
   SELECT CASE (NINIT) 
      CASE (NSCARC_INIT_UNDEFINED)
         WORKSPACE = NSCARC_UNDEFINED_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
       SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
       SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_SHUTDOWN('Inconsistent length for vector allocation ', CTEXT, -999)
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_REAL3

!> ------------------------------------------------------------------------------------------------
!> Allocate matrix with corresponding pointer and length structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_MATRIX(MAT, NAV, NAC, NAR, NSTENCIL, NINIT, CMAT)
TYPE (SCARC_MATRIX_TYPE), INTENT(INOUT) :: MAT
INTEGER, INTENT(IN) :: NAV, NAC, NAR, NSTENCIL, NINIT
CHARACTER(*), INTENT(IN) :: CMAT
CHARACTER(40) :: CINFO

WRITE(CINFO,'(A,A)') TRIM(CMAT),'.VAL'
CALL SCARC_ALLOCATE_REAL1(MAT%VAL, 1, NAV, NINIT, CINFO)

WRITE(CINFO,'(A,A)') TRIM(CMAT),'.COL'
CALL SCARC_ALLOCATE_INT1(MAT%COL, 1, NAC, NINIT, CINFO)

WRITE(CINFO,'(A,A)') TRIM(CMAT),'.ROW'
CALL SCARC_ALLOCATE_INT1(MAT%ROW, 1, NAR, NINIT, CINFO)

MAT%NAV = NAV
MAT%NAC = NAC
MAT%NAR = NAR
MAT%NSTENCIL = NSTENCIL

END SUBROUTINE SCARC_ALLOCATE_MATRIX

#ifdef WITH_MKL_FB
!> ------------------------------------------------------------------------------------------------
!> Allocate matrix with corresponding pointer and length structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_MATRIX_FB(MAT, NAV, NAC, NAR, NSTENCIL, NINIT, CMAT)
TYPE (SCARC_MATRIX_FB_TYPE), INTENT(INOUT) :: MAT
INTEGER, INTENT(IN) :: NAV, NAC, NAR, NSTENCIL, NINIT
CHARACTER(*), INTENT(IN) :: CMAT
CHARACTER(40) :: CINFO

WRITE(CINFO,'(A,A)') TRIM(CMAT),'.VAL'
CALL SCARC_ALLOCATE_REAL1_FB(MAT%VAL, 1, NAV, NINIT, CINFO)

WRITE(CINFO,'(A,A)') TRIM(CMAT),'.COL'
CALL SCARC_ALLOCATE_INT1(MAT%COL, 1, NAC, NINIT, CINFO)

WRITE(CINFO,'(A,A)') TRIM(CMAT),'.ROW'
CALL SCARC_ALLOCATE_INT1(MAT%ROW, 1, NAR, NINIT, CINFO)

MAT%NAV = NAV
MAT%NAC = NAC
MAT%NAR = NAR
MAT%NSTENCIL = NSTENCIL

END SUBROUTINE SCARC_ALLOCATE_MATRIX_FB
#endif

!> ------------------------------------------------------------------------------------------------
!> Allocate matrix with corresponding pointer and length structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_MATRIX(MAT)
TYPE (SCARC_MATRIX_TYPE), INTENT(INOUT) :: MAT

DEALLOCATE(MAT%VAL)
DEALLOCATE(MAT%COL)
DEALLOCATE(MAT%ROW)

MAT%NAV = 0
MAT%NAC = 0
MAT%NAR = 0
MAT%NSTENCIL = 0

END SUBROUTINE SCARC_DEALLOCATE_MATRIX

!> ------------------------------------------------------------------------------------------------
!> Reduce size of matrix to specified size
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_REDUCE_MATRIX(MAT, NAV, CMAT)
TYPE (SCARC_MATRIX_TYPE), INTENT(INOUT) :: MAT
INTEGER, INTENT(IN) :: NAV
CHARACTER(*), INTENT(IN) :: CMAT
CHARACTER(40) :: CINFO
TYPE (SCARC_MATRIX_TYPE) :: AUX

IF (SIZE(MAT%VAL) == NAV) THEN
   RETURN
ELSE IF (SIZE(MAT%VAL) > NAV) THEN
   WRITE(CINFO,'(A,A)') TRIM(CMAT),'.VAL'
   CALL SCARC_COPY_MATRIX(MAT, AUX, 'AUX') 
   CALL SCARC_DEALLOCATE_MATRIX(MAT)
   CALL SCARC_COPY_MATRIX(AUX, MAT, 'MAT-copied')
   CALL SCARC_DEALLOCATE_MATRIX(AUX)
ELSE
   CALL SCARC_SHUTDOWN('Length for resizing too big ', 'NONE', NAV)
ENDIF
   
END SUBROUTINE SCARC_REDUCE_MATRIX

!> ------------------------------------------------------------------------------------------------
!> Copy matrix with all corresponding pointer and length structures
!> If not already allocated, allocate complete matrix structure
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_COPY_MATRIX(MAT1, MAT2, CMAT2)
TYPE (SCARC_MATRIX_TYPE), INTENT(INOUT) :: MAT1, MAT2
CHARACTER(*), INTENT(IN) :: CMAT2
CHARACTER(40) :: CINFO

MAT2%NAV = MAT1%NAV
MAT2%NAC = MAT1%NAC
MAT2%NAR = MAT1%NAR
MAT2%NSTENCIL = MAT1%NSTENCIL

WRITE(CINFO,'(A,A)') TRIM(CMAT2),'.VAL'
CALL SCARC_ALLOCATE_REAL1(MAT2%VAL, 1, MAT2%NAV, NSCARC_INIT_NONE, CINFO)
CALL SCARC_COPY_REAL(MAT1%VAL, MAT2%VAL, 1.0_EB, MAT1%NAV)

WRITE(CINFO,'(A,A)') TRIM(CMAT2),'.COL'
CALL SCARC_ALLOCATE_INT1(MAT2%COL, 1, MAT2%NAC, NSCARC_INIT_NONE, CINFO)
CALL SCARC_COPY_INT(MAT1%COL, MAT2%COL, 1, MAT1%NAC)

WRITE(CINFO,'(A,A)') TRIM(CMAT2),'.ROW'
CALL SCARC_ALLOCATE_INT1(MAT2%ROW, 1, MAT2%NAR, NSCARC_INIT_NONE, CINFO)
CALL SCARC_COPY_INT(MAT1%ROW, MAT2%ROW, 1, MAT1%NAR)

END SUBROUTINE SCARC_COPY_MATRIX

!> ------------------------------------------------------------------------------------------------
!> Store name of currently called routine in history stack for verbosing
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ENTER_ROUTINE (CNAME)
CHARACTER(*), INTENT(IN) :: CNAME
MSG%NCURRENT = MSG%NCURRENT + 1
IF (MSG%NCURRENT > MSG%NHISTORY_MAX) &
   CALL SCARC_SHUTDOWN('Too many messages in calling stack ', CNAME, -999)
MSG%HISTORY(MSG%NCURRENT) = CNAME
END SUBROUTINE SCARC_ENTER_ROUTINE

!> ------------------------------------------------------------------------------------------------
!> Remove name of last called routine from history stack for verbosing
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_LEAVE_ROUTINE ()
MSG%HISTORY(MSG%NCURRENT) = 'NONE'
IF (MSG%NCURRENT == 0) &
   CALL SCARC_SHUTDOWN('Final level of calling history already reached', 'NONE', -999)
MSG%NCURRENT = MSG%NCURRENT - 1
END SUBROUTINE SCARC_LEAVE_ROUTINE

END MODULE SCRC
