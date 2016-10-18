#ifdef WITH_MKL
include "mkl_pardiso.f90"
include "mkl_cluster_sparse_solver.f90"
#endif

MODULE SCRC  

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: SECOND, GET_FILE_NUMBER, SHUTDOWN
USE MPI
USE TYPES, ONLY: MULTIPLIER_TYPE

#ifdef WITH_MKL
  USE MKL_PARDISO
  USE MKL_CLUSTER_SPARSE_SOLVER
#endif 

IMPLICIT NONE
!PRIVATE

! Public structures (needed in main, read, divg, dump)

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
PUBLIC SCARC_DISCRETIZATION                    !> Selected ScaRC discretization type (structured/unstructured)
PUBLIC SCARC_METHOD                            !> Selected ScaRC method (Krylov/Multigrid/MKL)
PUBLIC SCARC_INITIAL                           !> Initial solution for ScaRC

PUBLIC SCARC_RESIDUAL                          !> Residual of iterative solver
PUBLIC SCARC_ITERATIONS                        !> Number of iterations
PUBLIC SCARC_CAPPA                             !> Convergence rate
PUBLIC SCARC_ACCURACY                          !> Chosen accuracy type (relative/absolute)
PUBLIC SCARC_ACCURACY_DIVERGENCE               !> Divergence accuracy
PUBLIC SCARC_ACCURACY_RELATIVE                 !> Relative accuracy

PUBLIC SCARC_MKL                               !> Type of MKL method (LOCAL->Pardiso/GLOBAL->Cluster_Sparse_Solver)
PUBLIC SCARC_MKL_MTYPE                         !> Type of matrix (SYMMETRIC/NONSYMMETRIC)

PUBLIC SCARC_KRYLOV                            !> Type of Krylov method
PUBLIC SCARC_KRYLOV_ITERATIONS                 !> Maximum number of iterations for Krylov method
PUBLIC SCARC_KRYLOV_ACCURACY                   !> Requested accuracy for Krylov method

PUBLIC SCARC_MULTIGRID                         !> Type of multigrid method
PUBLIC SCARC_MULTIGRID_LEVEL                   !> Multigrid level
PUBLIC SCARC_MULTIGRID_CYCLE                   !> Type of multigrid cycle (V/W/F)
PUBLIC SCARC_MULTIGRID_COARSENING              !> Coarsening method for multigrid (AMG only)
PUBLIC SCARC_MULTIGRID_INTERPOL                !> Interpolation method for multigrid (AMG only)
PUBLIC SCARC_MULTIGRID_ITERATIONS              !> Maximum number of iterations for multigrid method
PUBLIC SCARC_MULTIGRID_ACCURACY                !> Requested accuracy for multigrid method

PUBLIC SCARC_SMOOTH                            !> Smoother for multigrid method
PUBLIC SCARC_SMOOTH_ITERATIONS                 !> Maximum number of iterations for smoothing method
PUBLIC SCARC_SMOOTH_ACCURACY                   !> Requested accuracy for smoothing method
PUBLIC SCARC_SMOOTH_OMEGA                      !> Damping parameter for smoothing method

PUBLIC SCARC_PRECON                            !> Preconditioner of defect correction method
PUBLIC SCARC_PRECON_ITERATIONS                 !> Maximum number of iterations for perconditioning method
PUBLIC SCARC_PRECON_ACCURACY                   !> Requested accuracy for preconditioning method

PUBLIC SCARC_COARSE                            !> Coarse grid solver for multigrid method
PUBLIC SCARC_COARSE_ITERATIONS                 !> Maximum number of iterations for coarse grid solver
PUBLIC SCARC_COARSE_ACCURACY                   !> Requested accuracy for coarse grid solver
PUBLIC SCARC_COARSE_MTYPE                      !> Type of direct coarse grid matrix (symmetric/nonsymmetric)

PUBLIC SCARC_CASE                              !> Predefined ScaRC case
PUBLIC SCARC_DEBUG                             !> Debugging parameter
PUBLIC SCARC_LAYER                             !> Number of layers for data exchange
PUBLIC LU_SCARC                                !> Debug file unit
<<<<<<< HEAD
PUBLIC TYPE_DEBUG                              !> Debug file unit
PUBLIC PRES_ON_WHOLE_DOMAIN
=======
PUBLIC PRES_ON_WHOLE_DOMAIN_SCRC
>>>>>>> origin/scarc

!> ------------------------------------------------------------------------------------------------
!> Corresponding declarations (with default settings)
!> ------------------------------------------------------------------------------------------------

!! General definitions
CHARACTER(40) :: SCARC_METHOD    = 'null'                   !> Requested solver method (KRYLOV/MULTIGRID)
CHARACTER(40) :: SCARC_INITIAL   = 'null'                   !> Initial solution (currently only default is used)

!! General iteration parameters
REAL (EB)     :: SCARC_RESIDUAL             = -1.0_EB       !> Residual of global selected solver
INTEGER       :: SCARC_ITERATIONS           =  0            !> Number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_CAPPA                =  1.0_EB       !> Convergence rate of selected ScarC solver
REAL (EB)     :: SCARC_ACCURACY_DIVERGENCE  =  1.E+6_EB     !> Divergence epsilon for all solvers
REAL (EB)     :: SCARC_ACCURACY_RELATIVE    =  1.E-2_EB     !> Minimum relative accuracy for all solvers
CHARACTER(40) :: SCARC_ACCURACY             = 'ABSOLUTE'    !> Accuracy type (ABSOLUTE/RELATIVE)

!! Parameter for DISCRET
CHARACTER(40) :: SCARC_DISCRETIZATION   = 'UNSTRUCTURED'    !> Discretization of whole domain or with OBSTs 

!! Parameter for MKL solver
CHARACTER(40) :: SCARC_MKL       = 'GLOBAL'                 !> Type of MKL solver (LOCAL->Pardiso/GLOBAL->Cluster_Sparse_solver)
CHARACTER(40) :: SCARC_MKL_MTYPE = 'SYMMETRIC'              !> Type of MKL solver (LOCAL->Pardiso/GLOBAL->Cluster_Sparse_solver)

!! Parameters for multigrid-type methods
CHARACTER(40) :: SCARC_MULTIGRID            = 'GEOMETRIC'   !> Type of MG-method (GEOMETRIC/ALGEBRAIC)
CHARACTER(40) :: SCARC_MULTIGRID_COARSENING = 'FALGOUT'     !> Coarsening strategy  (RS3/A1/A2/PMIS/FDS...)
CHARACTER(40) :: SCARC_MULTIGRID_INTERPOL   = 'DIRECT'      !> Interpolation strategy (DIRECT/RS/STANDARD)
CHARACTER(1)  :: SCARC_MULTIGRID_CYCLE      = 'V'           !> Cycling type  (F/V/W)
INTEGER       :: SCARC_MULTIGRID_LEVEL      = -1            !> User defined number of MG-levels (optionally)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS = 1000          !> Max number of iterations
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY   = 1.E-8_EB      !> Requested accuracy for convergence

!! Parameters for Krylov-type methods
CHARACTER(40) :: SCARC_KRYLOV            = 'CG'             !> Type of Krylov-method (CG/BICG)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 1000             !> Max number of iterations
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-12_EB        !> Requested accuracy for convergence

!! Parameters for smoothing method (used in multigrids-methods)
CHARACTER(40) :: SCARC_SMOOTH            = 'SSOR'           !> Smoother for MG (JACOBI/SSOR)
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 10               !> Max number of iterations
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-12_EB        !> Requested accuracy for convergence
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.90E+0_EB       !> Relaxation parameter

!! Parameters for preconditioning method (used in Krylov-methods)
CHARACTER(40) :: SCARC_PRECON            = 'SSOR'           !> Preconditioner for CG/BICG (JACOBI/SSOR/MG)
INTEGER       :: SCARC_PRECON_ITERATIONS = 1000             !> Max number of iterations
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-12_EB        !> Requested accuracy for convergence

!! Parameters for coarse grid method
CHARACTER(40) :: SCARC_COARSE            = 'DIRECT'         !> Coarse grid solver (iterative CG-solver/direct MKL-solver)
CHARACTER(40) :: SCARC_COARSE_MTYPE      = 'SYMMETRIC'      !> Type of coarse grid matrix (nonsymmetric/symmetric)
INTEGER       :: SCARC_COARSE_ITERATIONS = 100              !> Max number of iterations for iterative solver
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-12_EB        !> Requested accuracy for iterative solver

!! Debugging parameters
CHARACTER(40) :: SCARC_DEBUG = 'NONE'                       !> Debugging level (NONE/LESS/MEDIUM/MUCH)
CHARACTER(40) :: SCARC_FN, SCARC_FN_DUMP                    !> File name for ScaRC debug messages
INTEGER       :: LU_SCARC, LU_SCARC_DUMP                    !> Unit number for ScaRC debug file

CHARACTER(40) :: SCARC_CASE = 'NONE'                        !> Debugging level (NONE/LESS/MEDIUM/MUCH)
INTEGER       :: SCARC_LAYER = 1                            !> Unit number for ScaRC debug file

REAL (EB)     :: SCARC_WALL_TIME = -1.0_EB                  !> Time measurement parameter
REAL (EB)     :: SCARC_WALL_TIME_START = -1.0_EB            !> Initial time for time measurment

CHARACTER(100):: SCARC_MESSAGE                              !> ScaRC messages
CHARACTER(40) :: SCARC_ROUTINE                              !> name of active SCARC-routine

!! order for the treatment of the single mesh faces
INTEGER :: FACE_ORDER_STD(6) = (/-3,-1,-1,1,2,3/)           !> Standard order of mesh faces
INTEGER :: FACE_ORDER_XYZ(6) = (/1,-1,2,-2,3,-3/)           !> Coordinate direction related order of mesh faces

INTEGER :: MSOL
CHARACTER(60) :: CSOL

! Public declarations

!> ------------------------------------------------------------------------------------------------
!> Global constants
!> ------------------------------------------------------------------------------------------------
INTEGER, PARAMETER :: NSCARC_DISCRET_NONE           = -1, &
                      NSCARC_DISCRET_STRUCTURED     =  2, &
                      NSCARC_DISCRET_UNSTRUCTURED   =  3

INTEGER, PARAMETER :: NSCARC_DIMENSION_NONE         = -1, &
                      NSCARC_DIMENSION_TWO          =  2, &    !> two-dimensional problem
                      NSCARC_DIMENSION_THREE        =  3       !> three-dimensional problem

INTEGER, PARAMETER :: NSCARC_SCOPE_NONE             = -1, &
                      NSCARC_SCOPE_MAIN             =  1, &    !> method used as main solver
                      NSCARC_SCOPE_SMOOTH           =  2, &    !> method used as smoother
                      NSCARC_SCOPE_PRECON           =  3, &    !> method used as preconditiner
                      NSCARC_SCOPE_COARSE           =  4       !> method used as coarse grid solver

INTEGER, PARAMETER :: NSCARC_METHOD_NONE            = -1, &
                      NSCARC_METHOD_KRYLOV          =  1, &    !> Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID       =  2, &    !> multigrid-method as global solver
                      NSCARC_METHOD_MKL             =  3       !> multigrid-method as global solver

INTEGER, PARAMETER :: NSCARC_KRYLOV_NONE            = -1, &
                      NSCARC_KRYLOV_CG              =  1, &    !> CG   as Krylov solver
                      NSCARC_KRYLOV_BICG            =  2       !> BICG as Krylov solver

INTEGER, PARAMETER :: NSCARC_MULTIGRID_NONE         = -1, &
                      NSCARC_MULTIGRID_GEOMETRIC    =  1, &    !> geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC    =  2       !> algebraic multigrid

INTEGER, PARAMETER :: NSCARC_MKL_NONE               = -1, &
                      NSCARC_MKL_LOCAL              =  1, &    !> local use of MKL solver (Pardiso)
                      NSCARC_MKL_GLOBAL             =  2, &    !> global use of MKL solver (Cluster_Sparse_Solver)
                      NSCARC_MKL_COARSE             =  3       !> MKL solver on coarse grid level

INTEGER, PARAMETER :: NSCARC_EXCHANGE_NONE          = -1, &
                      NSCARC_EXCHANGE_BASIC         =  1, &    !> initialize wall information
                      NSCARC_EXCHANGE_WALLINFO      =  2, &    !> initialize wall information
                      NSCARC_EXCHANGE_VECTOR        =  3, &    !> matrix-vector communication
                      NSCARC_EXCHANGE_PRESSURE      =  4, &    !> vector values along internal boundaries
                      NSCARC_EXCHANGE_MEASURE       =  5, &    !> measure values along internal boundaries
                      NSCARC_EXCHANGE_MEASURE_ADD   =  6, &    !> measure values along internal boundaries
                      NSCARC_EXCHANGE_CELL_TYPE      =  7, &    !> cell types along internal boundaries
                      NSCARC_EXCHANGE_CELL_TYPE2     =  8, &    !> cell types II along internal boundaries
                      NSCARC_EXCHANGE_PROLONGATION  =  9, &    !> internal transfer weights
                      NSCARC_EXCHANGE_RESTRICTION   = 10, &    !> internal transfer weights
                      NSCARC_EXCHANGE_STENCIL       = 11, &    !> internal subdiagonal matrix values
                      NSCARC_EXCHANGE_MATRIX_SIZE   = 12, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_SUBDIAG= 13, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_STENCIL= 14, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_SYSTEM = 15, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_PROL   = 16, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_REST   = 17, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_TRANSFER_SIZE = 18, &    !> neighboring transfer matrix size
                      NSCARC_EXCHANGE_WIDTHINFO     = 19, &    !> neighboring grid resolution
                      NSCARC_EXCHANGE_MESHINFO      = 20, &    !> neighboring mesh information
                      NSCARC_EXCHANGE_GRAPH         = 21, &    !> graph along internal boundaries
                      NSCARC_EXCHANGE_CCVAR         = 22       !> exchange CCVAR along internal boundaries

INTEGER, PARAMETER :: NSCARC_SMOOTH_NONE            = -1, &
                      NSCARC_SMOOTH_JACOBI          =  1, &    !> smoothing by JACOBI-method
                      NSCARC_SMOOTH_SSOR            =  2, &    !> smoothing by SSOR-method
                      NSCARC_SMOOTH_FFT             =  3, &    !> smoothing by FFT-method
                      NSCARC_SMOOTH_PARDISO         =  4, &    !> smoothing by PARDISO-method
                      NSCARC_SMOOTH_CLUSTER         =  5       !> smoothing by PARDISO-method

INTEGER, PARAMETER :: NSCARC_PRECON_NONE            = -1, &
                      NSCARC_PRECON_JACOBI          =  1, &    !> preconditioning by JACOBI-method
                      NSCARC_PRECON_SSOR            =  2, &    !> preconditioning by SSOR-method
                      NSCARC_PRECON_FFT             =  3, &    !> preconditioning by FFT-method
                      NSCARC_PRECON_PARDISO         =  4, &    !> preconditioning by PARDISO-method
                      NSCARC_PRECON_CLUSTER         =  5, &    !> preconditioning by CLUSTER-method
                      NSCARC_PRECON_MULTIGRID       =  6       !> preconditioning by MG-method

INTEGER, PARAMETER :: NSCARC_CYCLE_NONE             = -1, &
                      NSCARC_CYCLE_F                =  0, &    !> F-cycle for mg-method
                      NSCARC_CYCLE_V                =  1, &    !> V-cycle for mg-method
                      NSCARC_CYCLE_W                =  2, &    !> W-cycle for mg-method
                      NSCARC_CYCLE_SETUP            =  3, &    !> initialize cycle counts
                      NSCARC_CYCLE_RESET            =  4, &    !> reset cycle counts
                      NSCARC_CYCLE_PROCEED          =  5, &    !> proceed cycle counts
                      NSCARC_CYCLE_PRESMOOTH        =  6, &    !> presmoothing cycle
                      NSCARC_CYCLE_POSTSMOOTH       =  7, &    !> postsmoothing cycle
                      NSCARC_CYCLE_NEXT             =  8, &    !> perform next cycling loop
                      NSCARC_CYCLE_EXIT             =  9       !> exit cycling loop

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED           =  0, &    !> proceed loop
                      NSCARC_STATE_CONV              =  1, &    !> convergence
                      NSCARC_STATE_DIVG              =  2       !> divergence

INTEGER, PARAMETER :: NSCARC_DEBUG_NONE             = -1, &    !> no debugging requested
                      NSCARC_DEBUG_INFO0            =  0, &    !> info1  level of debugging requested
                      NSCARC_DEBUG_INFO1            =  1, &    !> info1  level of debugging requested
                      NSCARC_DEBUG_INFO2            =  2, &    !> info2  level of debugging requested
                      NSCARC_DEBUG_LESS             =  3, &    !> low    level of debugging requested
                      NSCARC_DEBUG_MEDIUM           =  4, &    !> medium level of debugging requested
                      NSCARC_DEBUG_MUCH             =  5, &    !> strong level of debugging requested
                      NSCARC_DEBUG_EXTREME          = 22, &    !> extreme level of debugging requested
                      NSCARC_DEBUG_MATRIX           =  6, &    !> show matrix
                      NSCARC_DEBUG_MATRIXS          =  7, &    !> show matrix
                      NSCARC_DEBUG_MATRIXE          =  8, &    !> show matrix
                      NSCARC_DEBUG_IJKW             =  9, &    !> show IJKW
                      NSCARC_DEBUG_WALLINFO         = 10, &    !> show WALLINFO
                      NSCARC_DEBUG_FACEINFO         = 11, &    !> show FACEINFO
                      NSCARC_DEBUG_BCINDEX          = 12, &    !> show PRESSURE_BC_INDEX
                      NSCARC_DEBUG_ACELL            = 13, &    !> show WALL_CELL
                      NSCARC_DEBUG_GCELL            = 14, &    !> show GHOST_CELL
                      NSCARC_DEBUG_NCELL            = 15, &    !> show NOM_CELL
                      NSCARC_DEBUG_SUBDIVISION      = 16, &    !> show SUBDIVISION
                      NSCARC_DEBUG_MEASURE          = 17, &    !> show MEASURE
                      NSCARC_DEBUG_CELL_TYPE         = 18, &    !> show CELL_TYPE
                      NSCARC_DEBUG_GRAPH            = 19, &    !> show CELL_TYPE
                      NSCARC_DEBUG_COARSE           = 20, &    !> show coarse grid
                      NSCARC_DEBUG_PROLONGATION     = 21, &    !> show prolongation matrix
                      NSCARC_DEBUG_RESTRICTION      = 22       !> show restriction matrix

INTEGER, PARAMETER :: NSCARC_COARSENING_NONE        = -1, &
                      NSCARC_COARSENING_BASIC       =  1, &    !> basic coarsening
                      NSCARC_COARSENING_FALGOUT     =  2, &    !> parallel Falgout
                      NSCARC_COARSENING_RS3         =  3, &    !> parallel RS3
                      NSCARC_COARSENING_A1          =  4, &    !> aggressive 1 (path=1, length=2)
                      NSCARC_COARSENING_A2          =  5, &    !> aggressive 2 (path=2, length=2)
                      NSCARC_COARSENING_PMIS        =  6, &    !> PMIS
                      NSCARC_COARSENING_PMISG       =  7, &    !> PMIS
                      NSCARC_COARSENING_FDSRS3      =  8, &    !> FDSRS3 : FDS variant similar to RS3
                      NSCARC_COARSENING_FDSA1       =  9, &    !> FDSA1  : FDS variant similar to A1
                      NSCARC_COARSENING_FDSA2       = 10, &    !> FDSA2  : FDS variant similar to A2
                      NSCARC_COARSENING_BDRY        = 11, &    !> FDSA2  : FDS variant similar to A2
                      NSCARC_COARSENING_GMG         = 12, &    !> GMG    : GMG as AMG-variang
                      NSCARC_COARSENING_GMG3        = 13       !> GMG3   : for cell numbers divisable by 3

INTEGER, PARAMETER :: NSCARC_COARSE_NONE            = -1, &
                      NSCARC_COARSE_ITERATIVE       =  1, &    !> iterative solution of coarse grid problem
                      NSCARC_COARSE_DIRECT          =  2       !> direct solution of coarse grid problem

INTEGER, PARAMETER :: NSCARC_SIZE_NONE              = -1, &
                      NSCARC_SIZE_MATRIX            =  2, &    !> size of system matrix for compact system
                      NSCARC_SIZE_TRANSFER          =  3       !> size of transfer matrices for compact system

INTEGER, PARAMETER :: NSCARC_VECTOR_NONE            = -1, &
                      NSCARC_VECTOR_X               =  1, &    !> selection parameter for vector X
                      NSCARC_VECTOR_F               =  2, &    !> selection parameter for vector F
                      NSCARC_VECTOR_Y               =  3, &    !> selection parameter for vector Y
                      NSCARC_VECTOR_G               =  4, &    !> selection parameter for vector G
                      NSCARC_VECTOR_W               =  5, &    !> selection parameter for vector R
                      NSCARC_VECTOR_D               =  6, &    !> selection parameter for vector D
                      NSCARC_VECTOR_Z               =  7, &    !> selection parameter for vector Z
                      NSCARC_VECTOR_X2              =  8, &    !> selection parameter for vector X2
                      NSCARC_VECTOR_G2              =  9, &    !> selection parameter for vector D2
                      NSCARC_VECTOR_D2              = 10, &    !> selection parameter for vector D2
                      NSCARC_VECTOR_W2              = 11, &    !> selection parameter for vector R2
                      NSCARC_VECTOR_Y2              = 12, &    !> selection parameter for vector Y2
                      NSCARC_VECTOR_Z2              = 13, &    !> selection parameter for vector Y2
                      NSCARC_VECTOR_H               = 14, &    !> selection parameter for vector Y2
                      NSCARC_VECTOR_HS              = 15, &    !> selection parameter for vector Y2
                      NSCARC_VECTOR_MEASURE         = 16, &    !> selection parameter for vector MEASURE
                      NSCARC_VECTOR_CELL_TYPE        = 17, &    !> selection parameter for vector CELL_TYPE
                      NSCARC_VECTOR_GRAPH           = 18, &    !> selection parameter for vector GRAPH
                      NSCARC_VECTOR_PTR             = 19       !> selection parameter for vector CELL_TYPE

INTEGER, PARAMETER :: NSCARC_MATRIX_NONE            = -9999999, &
                      NSCARC_MATRIX_SUBDIAG         =  1, &    !> exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SUBDIAG_LOWER   =  2, &    !> exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SUBDIAG_UPPER   =  3, &    !> exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SYSTEM          =  4, &    !> exchange subdiagonal matrix entries
                      NSCARC_MATRIX_TRANSFER        =  5, &    !> exchange prolongation matrix
                      NSCARC_MATRIX_STENCIL         =  6, &    !> exchange prolongation matrix
                      NSCARC_MATRIX_PROLONGATION    =  7, &    !> exchange prolongation matrix
                      NSCARC_MATRIX_RESTRICTION     =  8, &    !> exchange restriction matrix
                      NSCARC_MATRIX_GMG             =  9       !> exchange prolongation matrix

INTEGER, PARAMETER :: NSCARC_ACCURACY_NONE          = -1, &
                      NSCARC_ACCURACY_ABSOLUTE      =  1, &    !> absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE      =  2       !> relative accuracy must be reached

REAL(EB), PARAMETER:: NSCARC_MEASURE_NONE           =  0.0_EB, &
                      NSCARC_MEASURE_ONE            =  1.0_EB, &  !> coarse-grid cell
                      NSCARC_MEASURE_COARSE         =  6.0_EB, &  !> coarse-grid cell
                      NSCARC_MEASURE_FINE           =  5.0_EB, &  !> fine-grid cell
                      NSCARC_MEASURE_SFINE          =  5.0_EB, &  !> strongly coupled fine-grid cell
                      NSCARC_MEASURE_WFINE          =  4.0_EB, &  !> weakly   coupled fine-grid cell
                      NSCARC_MEASURE_BDRY           =  1.0_EB     !> boundry weight

REAL(EB), PARAMETER:: NSCARC_GRAPH_NONE             =  0

INTEGER, PARAMETER :: NSCARC_CELL_TYPE_NONE          =  0, &
                      NSCARC_CELL_TYPE_COARSE        =  1, &   !> coarse-grid cell
                      NSCARC_CELL_TYPE_COARSE0       =  3, &   !> special coarse-grid cell
                      NSCARC_CELL_TYPE_COMMON        =  3, &   !> common cell
                      NSCARC_CELL_TYPE_FINE          = -1, &   !> fine-grid cell
                      NSCARC_CELL_TYPE_FINE0         = -3, &   !> special fine-grid cell
                      NSCARC_CELL_TYPE_SFINE         = -1, &   !> strongly coupled fine-grid cell
                      NSCARC_CELL_TYPE_WFINE         = -2, &   !> weakly   coupled fine-grid cell
                      NSCARC_CELL_TYPE_FPNT          = -1, &   !> special f-point
                      NSCARC_CELL_TYPE_ZPNT          = -2, &   !> special z-point
                      NSCARC_CELL_TYPE_SFPNT         = -3, &   !> special sf-point
                      NSCARC_CELL_TYPE_CPNT          =  2      !> special c-point

INTEGER, PARAMETER :: NSCARC_INTERPOL_NONE          =  0, &
                      NSCARC_INTERPOL_STANDARD      =  1, &    !> standard interpolation
                      NSCARC_INTERPOL_CLASSICAL     =  2, &    !> classical interpolation
                      NSCARC_INTERPOL_CLASSICAL2    =  3, &    !> classical interpolation
                      NSCARC_INTERPOL_DIRECT        =  4, &    !> direct interpolation
                      NSCARC_INTERPOL_DIRECT_BDRY   =  5, &    !> direct interpolation with special boundary
                      NSCARC_INTERPOL_MULTIPASS     =  6, &    !> multipass interpolation
                      NSCARC_INTERPOL_GMG           =  7, &    !> GMG-like interpolation
                      NSCARC_INTERPOL_GMG3          =  8       !> GMG3-like interpolation

INTEGER, PARAMETER :: NSCARC_LATEX_NONE             = -1, &    !> no latex information requested
                      NSCARC_LATEX_STAGGERED        =  1, &    !> show staggered latex information
                      NSCARC_LATEX_EQUAL            =  2, &    !> show equal latex information
                      NSCARC_LATEX_NUMBER           =  3, &    !> show number latex information
                      NSCARC_LATEX_ALL              =  4, &    !> show all latex information
                      NSCARC_LATEX_TABLE            =  5       !> show latex table

INTEGER, PARAMETER :: NSCARC_TIME_NONE              = -1, &
                      NSCARC_TIME_TOTAL             =  1, &    !> time for complete ScaRC part of FDS
                      NSCARC_TIME_SETUP             =  2, &    !> time for setup phase
                      NSCARC_TIME_PARDISO_SETUP     =  3, &    !> time for Pardiso setup
                      NSCARC_TIME_SOLVER            =  4, &    !> time for ScaRC solver
                      NSCARC_TIME_KRYLOV            =  5, &    !> time for Krylov solver
                      NSCARC_TIME_MULTIGRID         =  6, &    !> time for multigrid solver
                      NSCARC_TIME_PARDISO           =  7, &    !> time for Pardiso solver
                      NSCARC_TIME_PRECON            =  8, &    !> time for preconditioner
                      NSCARC_TIME_SMOOTH            =  9, &    !> time for smoother
                      NSCARC_TIME_COARSE            = 10, &    !> time for coarse grid solver
                      NSCARC_TIME_MATVEC            = 11, &    !> time for matrix-vector product
                      NSCARC_TIME_SCALPROD          = 12, &    !> time for scalar product
                      NSCARC_TIME_L2NORM            = 13, &    !> time for l2norm
                      NSCARC_TIME_EXCH_INIT         = 14, &    !> time for exchange initialization
                      NSCARC_TIME_EXCH_VECTOR       = 15, &    !> time for exchange of internal boundary
                      NSCARC_TIME_EXCH_MATRIX       = 16, &    !> time for exchange of internal matrix
                      NSCARC_TIME_EXCH_MEASURE      = 17       !> time for exchange of internal measure

INTEGER, PARAMETER :: NSCARC_CASE_NONE              = -1, &    !> no predefined initial solution used
                      NSCARC_CASE_CD_NSA_2D         =  1, &    !> CD_NSA_2D-case
                      NSCARC_CASE_VD_NSA_2D         =  2, &    !> VD_NSA_2D-case
                      NSCARC_CASE_CD_VA_2D          =  3, &    !> CD_VA_2D-case
                      NSCARC_CASE_ZM_GRAV_ADV_2D    =  4       !> ZM_GRAV_ADVECTED_2D

INTEGER, PARAMETER :: NSCARC_IOR_TOP_Z              = -3, &    !> top z-face
                      NSCARC_IOR_BACK_Y             = -2, &    !> back y-face
                      NSCARC_IOR_RIGHT_X            =  1, &    !> right x-face
                      NSCARC_IOR_DUMMY              =  0, &    !> dummy value
                      NSCARC_IOR_LEFT_X             =  1, &    !> left x-face
                      NSCARC_IOR_FRONT_Y            =  2, &    !> front y-face
                      NSCARC_IOR_BOTTOM_Z           =  3       !> bottom z-face

INTEGER, PARAMETER :: NSCARC_LEVEL_NONE             = -1, &    !> no predefined initial solution used
                      NSCARC_LEVEL_MIN              =  0, &    !> minimum multigrid level
                      NSCARC_LEVEL_MAX              =  15      !> maximum multigrid level

INTEGER, PARAMETER :: NSCARC_LAYER_NONE             = -1, &    !> no neighbor
                      NSCARC_LAYER_ONE              =  1, &    !> communication of one abutting layer
                      NSCARC_LAYER_TWO              =  2       !> communication of two abutting layers

INTEGER, PARAMETER :: NSCARC_DUMP_NONE              = -1, &    !> no dumping
                      NSCARC_DUMP_RHS               =  1, &    !> dump rhs
                      NSCARC_DUMP_PRES              =  2       !> dump pressure

INTEGER, PARAMETER :: NSCARC_STENCIL_NONE           = -1, &    !> no matrix stencil available
                      NSCARC_STENCIL_CENTRAL        =  1, &    !> standard 5- or 7-point stencil
                      NSCARC_STENCIL_AMG            =  2       !> arbitrary AMG-stencil

INTEGER, PARAMETER :: NSCARC_DIAG_MAIN              =  1, &    !> standard 5- or 7-point stencil
                      NSCARC_DIAG_LOWER             =  2, &    !> standard 5- or 7-point stencil
                      NSCARC_DIAG_UPPER             =  3       !> standard 5- or 7-point stencil

INTEGER, PARAMETER :: NSCARC_INITIAL_NONE           = -1       !> another initial function ?

INTEGER, PARAMETER :: NSCARC_COUPLING_MAX           = 10       !> maximum of possible couplings in stencil

INTEGER, PARAMETER :: NSCARC_DUMMY                  = -1       !> dummy variable (needed at several places)

INTEGER, PARAMETER :: NSCARC_DIRECTION_X            =  1, &    !> x-direction indicator
                      NSCARC_DIRECTION_Y            =  2, &    !> y-direction indicator
                      NSCARC_DIRECTION_Z            =  3       !> z-direction indicator

INTEGER, PARAMETER :: NSCARC_NUM_FACES              =  6       !> number of faces per mesh
INTEGER, PARAMETER :: NSCARC_MAX_FACE_NEIGHBORS     = 20       !> max number neighbors per mesh face
INTEGER, PARAMETER :: NSCARC_MAX_MESH_NEIGHBORS     =  6*NSCARC_MAX_FACE_NEIGHBORS
INTEGER, PARAMETER :: NSCARC_MAX_STENCIL            =  7

INTEGER,  PARAMETER :: NGUARD= 2 ! Two layers of guard-cells.
INTEGER,  PARAMETER :: FCELL = 1 ! Right face index.

INTEGER,  PARAMETER :: IS_GASPHASE  = -1
INTEGER,  PARAMETER :: IS_CUTCFE    =  0
INTEGER,  PARAMETER :: IS_SOLID     =  1
INTEGER,  PARAMETER :: IS_UNDEFINED =-11

INTEGER,  PARAMETER :: IS_CGSC   = 1 ! Face media type: IS_GASPHASE, IS_SOLID or IS_CUTCFE.
INTEGER,  PARAMETER :: IS_UNKH   = 2 ! H unknown number.
INTEGER,  PARAMETER :: IS_NCVARS = 2 ! Number of face variables in MESHES(NM)%CCVAR.

<<<<<<< HEAD
LOGICAL :: PRES_ON_WHOLE_DOMAIN = .TRUE.
=======
LOGICAL :: PRES_ON_WHOLE_DOMAIN_SCRC = .TRUE.
>>>>>>> origin/scarc

! Cartesian Cell centered variables, actual case initialized as CC_IBM=.FALSE.:
INTEGER :: CGSC=IS_CGSC, UNKH=IS_UNKH, NCVARS=IS_NCVARS
INTEGER :: NUNKH_LOCAL

INTEGER, SAVE :: ILO_CELL,IHI_CELL,JLO_CELL,JHI_CELL,KLO_CELL,KHI_CELL
INTEGER, SAVE :: ILO_FACE,IHI_FACE,JLO_FACE,JHI_FACE,KLO_FACE,KHI_FACE
INTEGER, SAVE :: NXB, NYB, NZB


!> --------------------------------------------------------------------------------------------
!> Global variables
!> --------------------------------------------------------------------------------------------
!> use integer types for the user defined input data (based on SCARC_TYPE_... variables)
INTEGER :: TYPE_DISCRET    = NSCARC_DISCRET_NONE          !> Type of discretization (structured/unstructured)
INTEGER :: TYPE_DIMENSION  = NSCARC_DIMENSION_NONE        !> Dimension type (2D/3D)
INTEGER :: TYPE_SCOPE      = NSCARC_SCOPE_NONE            !> Type of surrounding solver scopy
INTEGER :: TYPE_SCOPE0     = NSCARC_SCOPE_MAIN            !> Type of surrounding solver scope II
INTEGER :: TYPE_METHOD     = NSCARC_METHOD_NONE           !> Type of global ScaRC method
INTEGER :: TYPE_METHOD0    = NSCARC_METHOD_NONE           !> Type of local ScaRC method
INTEGER :: TYPE_KRYLOV     = NSCARC_KRYLOV_NONE           !> Type of Krylov method (CG/BICG)
INTEGER :: TYPE_MULTIGRID  = NSCARC_MULTIGRID_NONE        !> Type of multigrid method (GMG/AMG)
INTEGER :: TYPE_MKL        = NSCARC_MKL_NONE              !> Type of MKL method (PARDISO/CLUSTER_SPARSE_SOLVER)
INTEGER :: TYPE_MATRIX     = NSCARC_MATRIX_NONE           !> Type of matrix for data exchange
INTEGER :: TYPE_ACCURACY   = NSCARC_ACCURACY_NONE         !> Type of requested accuracy
INTEGER :: TYPE_SMOOTH     = NSCARC_SMOOTH_SSOR           !> Type of smoother for multigrid method
INTEGER :: TYPE_PRECON     = NSCARC_PRECON_SSOR           !> Type of preconditioner for global iterative solver
INTEGER :: TYPE_PRECON0    = NSCARC_PRECON_SSOR           !> Type of preconditioner for local iterative solver
INTEGER :: TYPE_CYCLE      = NSCARC_CYCLE_V               !> Type of cycling for multigrid method
INTEGER :: TYPE_COARSENING = NSCARC_COARSENING_NONE       !> Type of coarsening algorithm for AMG
INTEGER :: TYPE_INTERPOL   = NSCARC_INTERPOL_NONE         !> Type of interpolation for AMG
INTEGER :: TYPE_COARSE     = NSCARC_COARSE_NONE           !> Type of coarse grid solver for multigrid method
INTEGER :: TYPE_CASE       = NSCARC_CASE_NONE             !> Type of predefined test case
INTEGER :: TYPE_DEBUG      = NSCARC_DEBUG_NONE            !> Type of debugging level
INTEGER :: TYPE_INITIAL    = NSCARC_INITIAL_NONE          !> Type of initial solution
INTEGER :: TYPE_EXCHANGE   = NSCARC_EXCHANGE_NONE         !> Type of data exchange
INTEGER :: TYPE_VECTOR     = NSCARC_VECTOR_NONE           !> Type of vector
INTEGER :: TYPE_LATEX      = NSCARC_LATEX_ALL             !> Type of Latex output level
INTEGER :: TYPE_DUMP       = NSCARC_DUMP_RHS              !> Type of Dump level
INTEGER :: TYPE_LAYER      = NSCARC_LAYER_ONE             !> Type of layers (for overlap)

INTEGER, PARAMETER :: NSCARC_TIMER=17                     !> Number of timers for ScaRC

INTEGER :: NLEVEL, NLEVEL_MAX, NLEVEL_MIN                 !> Total, minimum and maximum number of multigrid levels
INTEGER :: NMASTER, NC_COARSE0                            !> Settings for global coarse solver
INTEGER :: NREQ_SCARC, N_EXCHANGES, TAG_SCARC             !> Variables for data exchange
INTEGER :: SNODE, RNODE                                   !> Process identifier for data exchange
INTEGER :: STATUS2_SCARC(MPI_STATUS_SIZE)                 !> Status array for data excange
REAL(EB):: SP_GLOBAL                                      !> Global scalar product

INTEGER,  ALLOCATABLE, DIMENSION (:,:):: NCS_OFFSET       !> Offsets for global numbering on different meshes
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: NCS_GLOBAL       !> Global cell number (over all meshes) over levels
INTEGER,  ALLOCATABLE, DIMENSION (:,:):: NCS_BLOCK        !> Group cell numbers (related to a process group) over levels
INTEGER,  ALLOCATABLE, DIMENSION (:,:):: NCS_LOCAL        !> Local number of cells (related to single meshes)
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: NC_COARSE        !> Global number of cells for coarse solver
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: REQ_SCARC        !> Request array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: COUNTS_SCARC     !> Counter array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: DISPLS_SCARC     !> Displacement array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: GROUP_ID         !> Displacement array for data exchange
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: SP_LOCAL         !> Local scalar procucts
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: SP_GROUP         !> Scalar procuct of process group
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TUSED_SCARC      !> Time measurement array


!! Private type declarations

!> --------------------------------------------------------------------------------------------
!> Main administration type for ScaRC-method
!> --------------------------------------------------------------------------------------------
TYPE SCARC_TYPE
!
INTEGER :: NUM_NEIGHBORS = 0                                 !> Number of adjacent neighbors of whole mesh
INTEGER, POINTER, DIMENSION(:) :: NEIGHBORS                  !> List of adjacent neighbors of whole mesh
!
INTEGER :: CYCLE_COUNT(2, NSCARC_LEVEL_MAX) = 0              !> Counter for multigrid cycling
REAL (EB), POINTER, DIMENSION (:,:) :: A_COARSE              !> System matrix for coarse grid solver
REAL (EB), POINTER, DIMENSION (:)   :: X_COARSE, X_BUF       !> Solution and RHS vecors for coarse grid solver
INTEGER  , POINTER, DIMENSION (:,:) :: PIVOT                 !> Pivot vector for coarse grid solver
!
INTEGER, POINTER, DIMENSION (:)   :: OFFSET                  !> Offset vector for global data exchange
INTEGER, POINTER, DIMENSION (:)   :: COUNTS, COUNTS1         !> Counters for data exchange
INTEGER, POINTER, DIMENSION (:)   :: DISPLS, DISPLS1         !> Displacement vectors for data exchange

INTEGER  :: NX, NY, NZ                                       !> number of cells in x-, y- and z-direction on finest level
REAL(EB) :: XS, XF, YS, YF, ZS, ZF                           !> x-, y- and z-bounds of mesh (same for all levels)
REAL(EB) :: XS_MIN, XF_MAX, YS_MIN, YF_MAX, ZS_MIN, ZF_MAX   !min and max x-, y- and z-coordinates in total


#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE)   , POINTER, DIMENSION(:) :: MKL       !> MKL-type
#endif
TYPE (SCARC_PRECON_TYPE), POINTER, DIMENSION(:) :: PRECON    !> Preconditioning type
TYPE (SCARC_LEVEL_TYPE) , POINTER, DIMENSION(:) :: LEVEL     !> Grid level type
TYPE (OSCARC_TYPE)      , POINTER, DIMENSION(:) :: OSCARC    !> ScaRC type on other mesh

END TYPE SCARC_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration other mesh data needed for the coupling of adjacent neighbors
!> --------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE
REAL (EB), POINTER, DIMENSION (:) :: SEND_REAL, RECV_REAL       !> main real send and receive buffers
INTEGER  , POINTER, DIMENSION (:) :: SEND_INTEGER, RECV_INTEGER !> main integer send and receive buffers
REAL(EB) :: SEND_REAL0(50), RECV_REAL0(50)                      !> initial real send and receive buffers
INTEGER  :: SEND_INTEGER0(50)  , RECV_INTEGER0(50)              !> initial integer send and receive buffers
INTEGER  :: NICMAX_R=0, NICMAX_S=0, NIC_R=0, NIC_S=0

TYPE (SCARC_LEVEL_TYPE), POINTER, DIMENSION(:) :: LEVEL
END TYPE OSCARC_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of face information related to wall cells and neighbors
!> --------------------------------------------------------------------------------------------
TYPE SCARC_FACE_TYPE

LOGICAL :: BXY, BXZ, BYZ                       !> type of face

INTEGER :: NFM                                 !> own mesh number
INTEGER :: NFC, NFW                            !> number of cells and wall cells along face
INTEGER :: NFX, NFY, NFZ                       !> local face dimensions
INTEGER :: NCPL = 1                            !> number of adjacent couplings
INTEGER :: NUM_NEIGHBORS = 0                   !> number of adjacent neighbors
INTEGER :: IWG_MARKER                          !> first IW number (corresponding to global numbering)
INTEGER :: IOFFSET_WALL   = 0                  !> counter for wall cells over all faces
INTEGER :: IOFF_MAT = 0                        !> matrix offset pointer

INTEGER , POINTER, DIMENSION(:) :: NEIGHBORS   !> adjacent neighbors
REAL(EB), POINTER, DIMENSION(:) :: DH          !> adjacent grid sizes

END TYPE SCARC_FACE_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of wall information related to neighbors and BC's
!> --------------------------------------------------------------------------------------------
TYPE SCARC_WALL_TYPE
!>
!> properties of the single wall cells
!>
INTEGER :: BTYPE                        !> boundary type of wall cell (1:Dirichlet/2:Neumann/3:Internal)
<<<<<<< HEAD
INTEGER :: BOUNDARY_TYPE=0              !> boundary type of wall cell (1:Dirichlet/2:Neumann/3:Internal)
=======
>>>>>>> origin/scarc
INTEGER :: IOR                          !> orientation of wall cell
INTEGER :: NOM                          !> neighbor at wall cell
INTEGER :: NCPL = 1                     !> number of couplings at wall cell (depending on resolution of neighbor)
INTEGER :: IWL                          !> corresponding local wall cell number for neighbor NOM
INTEGER :: CGSC, UNKH                   !> CCVAR-structure (for neighbors)
!>
!> Cell and coordinate numberings
!>
INTEGER:: ICW  = IS_UNDEFINED           !> internal wall cell for IW
INTEGER:: ICW2 = IS_UNDEFINED           !> internal wall cell for IW
INTEGER:: ICO, ICO2                     !> overlapping cell for IW
INTEGER:: IXG, IYG, IZG                 !> x-, y- and z-indices of ghost cells
INTEGER:: IXW, IYW, IZW                 !> x-, y- and z-indices of (internal) wall cells
INTEGER:: IXN(2), IYN(2), IZN(2)        !> x-, y- and z-indices of neighboring cells
!>
!> pointer vecors to ghost, extended and neighboring cells
!>
INTEGER, POINTER, DIMENSION(:) :: ICE, ICE2         !> extended cell for IW
INTEGER, POINTER, DIMENSION(:) :: ICG, ICG2         !> ghost cell for IW

END TYPE SCARC_WALL_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of data related to single grid level
!> --------------------------------------------------------------------------------------------
TYPE SCARC_LEVEL_TYPE

INTEGER :: SUBDIVISION(3,-3:3)=0
INTEGER :: NLEN_MATRIX_SYSTEM = 0
!>
!> Different numbers, lengths and settings
!>
INTEGER :: NX, NY, NZ           !> number of grid cells in x-, y- and z-direction
INTEGER :: NW, NWL              !> number of global and local wall cells
INTEGER :: NC, NCS              !> number of cells
INTEGER :: NCG=0, NCG0=0        !> number of ghost cells
INTEGER :: NCE=0, NCE0=0        !> number of extended cells plus auxiliary variable
INTEGER :: NCO=0, NCO0=0        !> number of overlapping cells plus auxiliary variable
INTEGER :: NA, NAE, NA0         !> number of elements in internal and extended system matrix
INTEGER :: NAS                  !> number of elements in internal and extended system matrix
INTEGER :: NP, NPE, NP0         !> number of elements in internal and extended prolongation matrix
INTEGER :: NR, NRE, NR0         !> number of elements in internal and extended restriction matrix
INTEGER :: NLAYER               !> number of layers
INTEGER :: NCPL=1, NCPL_MAX=-10 !> number of couplings
INTEGER :: NCPLS, NCPLR         !> number of couplings to send and read
INTEGER :: NCW                  !> number of coarse wall cells
INTEGER :: NCF, NCFE, NCF0      !> number of fine cells in internal and extended grid (AMG)
INTEGER :: NCC, NCCE, NCC0      !> number of coarse cells in internal and extended grid (AMG)
INTEGER :: NCCI                 !> number of internal coarse cells
INTEGER :: NOM                  !> number of adjacent neighbor
INTEGER :: NPOINTS              !> number of points in matrix stencil
INTEGER :: NOFFSET_MESH         !> offset for global mesh numbering
INTEGER :: N_EXTERNAL_WALL_CELLS
INTEGER :: N_INTERNAL_WALL_CELLS
INTEGER :: N_WALL_CELLS
INTEGER :: N_OBST
<<<<<<< HEAD
INTEGER :: CELL_COUNT
=======
>>>>>>> origin/scarc
!>
!> Different pointers
!>
INTEGER :: IWL = 0              !> local wall cell numbers
INTEGER :: ICG = 0              !> ghost cell counter
INTEGER :: ICG0 = 0             !> auxiliary ghost cell counter

INTEGER :: ICG_MARKER = 0       !> ghost cell marker
INTEGER :: ICO_MARKER = 0       !> overlapping cell marker
INTEGER :: ICE_MARKER = 0       !> extended cell marker
INTEGER :: IWL_MARKER = 0       !> local wall cell marker
INTEGER :: IWG_MARKER = 0       !> global wall cell marker
!>
!> Coordinate data
!>
REAL(EB) :: DI2(3)              !> sample vector of step sizes in x-, y- and z-direction
REAL(EB) :: DX , DY , DZ        !> step sizes in x-, y- and z-direction
REAL(EB) :: DXI, DYI, DZI       !> inversed of step sizes in x-, y- and z-direction
REAL(EB) :: DXI2, DYI2, DZI2    !> squared and inversed step sizes in x-, y- and z-direction
!>
REAL (EB), POINTER, DIMENSION (:)  :: DXL, DYL, DZL        !> step size vectors in x-, y- and z-direction
REAL (EB), POINTER, DIMENSION (:)  :: XCOR, YCOR, ZCOR     !> coordinate vectors in x-, y- and z-direction
REAL (EB), POINTER, DIMENSION (:)  :: XMID, YMID, ZMID     !> midpoint vectors in x-, y- and z-direction
!>
!> matrices with corresponding storage vectors
!>
REAL (EB), POINTER, DIMENSION (:) :: A , P , R, S     !> system, prolongation, restriction and strength matrix
INTEGER,   POINTER, DIMENSION (:) :: A_ROW , A_COL    !> row and column pointers for system matrix A
INTEGER,   POINTER, DIMENSION (:) :: S_ROW , S_COL    !> row and column pointers for strength matrix S
INTEGER,   POINTER, DIMENSION (:) :: ST_ROW, ST_COL   !> row and column pointers for transpose of strength matrix S
INTEGER,   POINTER, DIMENSION (:) :: P_ROW , P_COL    !> row and column pointers for prolongation matrix P
INTEGER,   POINTER, DIMENSION (:) :: R_ROW , R_COL    !> row and column pointers for restriction matrix A
INTEGER,   POINTER, DIMENSION (:) :: A_TAG , P_TAG    !> auxiliary arrays for mark positions in A and P
INTEGER,   POINTER, DIMENSION (:) :: A_SIZE, AUX      !> some more auxiliary arrays

INTEGER,   POINTER, DIMENSION(:,:,:,:) :: CCVAR  
INTEGER, ALLOCATABLE, DIMENSION(:) :: NUNKH_LOC, NUNKH_TOT, UNKH_IND

!>
!> Different vector definitions for global and local solvers
!>
REAL (EB), POINTER, DIMENSION (:) :: X, X2, XU        !> solution vectors
REAL (EB), POINTER, DIMENSION (:) :: F, F2, FU        !> right hand side vectors
REAL (EB), POINTER, DIMENSION (:) :: D, D2, DU        !> defect vectors
REAL (EB), POINTER, DIMENSION (:) :: Y, Y2, YU        !> auxiliary vectors
REAL (EB), POINTER, DIMENSION (:) :: G, G2, GU        !> auxiliary vectors
REAL (EB), POINTER, DIMENSION (:) :: W, W2, WU        !> auxiliary vectors
REAL (EB), POINTER, DIMENSION (:) :: Z, Z2, ZU        !> auxiliary vectors
!>
!> Different pointer vectors for wall and communication related data
!>
<<<<<<< HEAD
INTEGER, POINTER, DIMENSION(:,:)  :: WALL_INDEX         !> wall index
INTEGER, POINTER, DIMENSION(:,:,:):: CELL_INDEX         !> wall index
INTEGER, POINTER, DIMENSION(:)    :: INTERNAL_BDRY_CELL !! index of internal boundary cells
INTEGER, POINTER, DIMENSION(:)    :: EXT_PTR            !> index of external cells
=======
INTEGER, POINTER, DIMENSION(:,:):: WALL_INDEX         !> wall index
INTEGER, POINTER, DIMENSION(:)  :: INTERNAL_BDRY_CELL !! index of internal boundary cells
INTEGER, POINTER, DIMENSION(:)  :: EXT_PTR            !> index of external cells
>>>>>>> origin/scarc
!>
!> neighbor related settings for IOR, NCPL and NWL
!>
REAL(EB):: DH = 0.0_EB                                !> local step sizes
INTEGER :: IOR = 0                                    !> local orientations

!>
!> Different mappings of quantities
!>
INTEGER, POINTER, DIMENSION (:)  :: ICE_TO_IWG        !> mapping from ICE to IWG
INTEGER, POINTER, DIMENSION (:)  :: ICE_TO_IWL        !> mapping from ICE to IWL
INTEGER, POINTER, DIMENSION (:)  :: ICE_TO_ICG        !> mapping from ICE to ICG
!>
INTEGER, POINTER, DIMENSION (:)  :: ICG_TO_IWG        !> mapping from ICG to IWG
INTEGER, POINTER, DIMENSION (:)  :: ICG_TO_ICE        !> mapping from ICG to ICE
INTEGER, POINTER, DIMENSION (:)  :: ICG_TO_ICO        !> mapping from ICG to ICE
!>
INTEGER, POINTER, DIMENSION (:)  :: ICN_TO_ICE        !> mapping from ICN to ICE    ! ACHTUNG: RIESIG
!>
INTEGER, POINTER, DIMENSION (:)  :: IWL_TO_IWG        !> mapping from IWL to IWG
INTEGER, POINTER, DIMENSION (:)  :: IWL_TO_ICW        !> mapping from IWL to ICW
INTEGER, POINTER, DIMENSION (:)  :: IWL_TO_ICO        !> mapping from IWL to ICO
INTEGER, POINTER, DIMENSION (:,:):: IWL_TO_ICG        !> mapping from IWL to ICG
!>
!> special administration vectors/arrays for Algebraic multigrid (AMG)
!>
REAL(EB):: MAX_ROWSUM=0.9_EB
REAL(EB), POINTER, DIMENSION (:)  :: MEASURE          !> mapping for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:,:):: CELL_MAP         !> cell mapping for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:)  :: CELL_TYPE         !> cell types for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:)  :: CELL_MARKER      !> cell markers for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:)  :: GRAPH            !> graph vector for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:,:):: P_PTR            !> pointer vector for AMG-prolongation
!>
!> local type definitions
!>
TYPE (SCARC_WALL_TYPE), POINTER, DIMENSION(:) :: WALL !! description of single wall cells
TYPE (SCARC_FACE_TYPE), POINTER, DIMENSION(:) :: FACE !! description of single faces
TYPE (SCARC_OBST_TYPE), POINTER, DIMENSION(:) :: OBST !! description of single obstructions
!>
END TYPE SCARC_LEVEL_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of scopes for solution methods
!> --------------------------------------------------------------------------------------------
TYPE SCARC_SCOPE_TYPE
INTEGER :: X, F, Y, G, W, D, Z              !> local references of solver vectors
INTEGER :: ITE, NIT                         !> local references to iteration counters
REAL(EB) :: EPS, RES, RESIN, OMEGA          !> local references to solver parameters
CHARACTER(30) :: CROUTINE = 'null'
END TYPE SCARC_SCOPE_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of obstructions for coarser grid levels 
!> --------------------------------------------------------------------------------------------
TYPE SCARC_OBST_TYPE
INTEGER :: I1, I2, J1, J2, K1, K2
END TYPE SCARC_OBST_TYPE


!> --------------------------------------------------------------------------------------------
!> Save parent settings when starting new solution method
!> --------------------------------------------------------------------------------------------
TYPE SCARC_PARENT_TYPE
INTEGER :: TYPE_METHOD           !> type of solution method
INTEGER :: TYPE_SCOPE            !> type of related scope
INTEGER :: TYPE_PRECON           !> type of preconditioning method
INTEGER :: TYPE_SMOOTH           !> type of smoothing
INTEGER :: TYPE_ACCURACY         !> type of accuracy requirements (relative/absolute)
INTEGER :: TYPE_CYCLE            !> type of multigrid cycle
END TYPE SCARC_PARENT_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of solver information for preconditioners
!> --------------------------------------------------------------------------------------------
TYPE SCARC_PRECON_TYPE
REAL (EB), POINTER, DIMENSION (:)       :: LDX, UDX, sMUDX, MDY, MDZ
REAL (EB), POINTER, DIMENSION (:)       :: DWORK, PERIOD
REAL (EB), POINTER, DIMENSION (:, :, :) :: FFT
END TYPE SCARC_PRECON_TYPE

!> --------------------------------------------------------------------------------------------
!> Administration of MKL information
!> --------------------------------------------------------------------------------------------
#ifdef WITH_MKL
TYPE SCARC_MKL_TYPE
TYPE(MKL_PARDISO_HANDLE),               ALLOCATABLE :: PT_H(:), PT(:)
TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), ALLOCATABLE :: CT_H(:), CT(:)
INTEGER, ALLOCATABLE :: IPARM(:)
REAL (EB), POINTER, DIMENSION (:) :: ASYM                 !> symmetric and global matrix
INTEGER, POINTER, DIMENSION (:)   :: ASYM_ROW, ASYM_COL   !> row and column pointers for symmetric matrix AS
INTEGER, POINTER, DIMENSION (:)   :: AG_COL               !> column pointers for global numbering
INTEGER :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR, MSGLVL
INTEGER :: PERM(1)
END TYPE SCARC_MKL_TYPE
#endif

!> --------------------------------------------------------------------------------------------
!> Pardiso structures, temporary form
!> --------------------------------------------------------------------------------------------


!> --------------------------------------------------------------------------------------------
!> Globally used types
!> --------------------------------------------------------------------------------------------
TYPE ( SCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET ::  SCARC
TYPE (OSCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: OSCARC


CONTAINS


!> ------------------------------------------------------------------------------------------------
!> Initialize ScaRC structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP
INTEGER :: IERR
REAL(EB):: TNOW

!> ------------------------------------------------------------------------------------------------
!> Initialize time measurement and store dimension of problem
!> ------------------------------------------------------------------------------------------------
IERR = 0
TNOW = SECOND()

ALLOCATE(TUSED_SCARC(0:NSCARC_TIMER,NMESHES),STAT=IERR)
CALL ChkMemErr('SCARC_SETUP','TUSED_SCARC',IERR)
TUSED_SCARC = 0._EB

IF (TWO_D) THEN
   TYPE_DIMENSION = NSCARC_DIMENSION_TWO
ELSE
   TYPE_DIMENSION = NSCARC_DIMENSION_THREE
ENDIF

ITERATE_PRESSURE = .TRUE.  ! Although there is no need to do pressure iterations to drive down velocity error
                           ! on wall cells (i.e. the solution should give the right unique dH/dxn), leave it
                           ! .TRUE. to write out velocity error diagnostics.


!> ------------------------------------------------------------------------------------------------
!> Parse input parameters
!> ------------------------------------------------------------------------------------------------
CALL SCARC_INPUT_PARSER                  !> parse input parameters (corresponding to FDS-file)


!> ------------------------------------------------------------------------------------------------
!> Setup different components of ScaRC
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_DEBUGGING               !> open debug file if requested
CALL SCARC_SETUP_LEVELS                  !> define number of necessary grid levels
CALL SCARC_SETUP_TYPES                   !> allocate basic ScaRC-types for all necessary grid levels
CALL SCARC_SETUP_MESHES                  !> setup mesh information
CALL SCARC_SETUP_INTERFACES              !> allocate structures related to mesh interfaces
CALL SCARC_SETUP_DECOMPOSITION           !> setup information related to domain decomposition
CALL SCARC_SETUP_GLOBALS                 !> define some global variables
CALL SCARC_SETUP_EXCHANGE                !> setup information for data exchange
CALL SCARC_SETUP_SYSTEM                  !> assemble linear system of equations
CALL SCARC_SETUP_COARSENING              !> perform coarsening on different grid levels if requested (AMG only)
CALL SCARC_SETUP_SOLVER                  !> setup single solvers (allocate work and auxiliary vectors)


!> ------------------------------------------------------------------------------------------------
!> Measure time for setup routine
!> ------------------------------------------------------------------------------------------------
TUSED_SCARC(NSCARC_TIME_SETUP,MYID+1)=TUSED_SCARC(NSCARC_TIME_SETUP,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL,MYID+1)+SECOND()-TNOW

END SUBROUTINE SCARC_SETUP


!> ------------------------------------------------------------------------------------------------
!> Setup debug file if requested
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DEBUGGING
INTEGER:: NM

IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) THEN
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      WRITE (SCARC_FN, '(A,A,i3.3)') TRIM(CHID),'.scarc',NM
      WRITE(*,*) 'Opening ', SCARC_FN
      LU_SCARC = GET_FILE_NUMBER()
      OPEN (LU_SCARC, FILE=SCARC_FN)
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_DEBUGGING


!> ----------------------------------------------------------------------------------------------------
!> Determine types of input parameters
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INPUT_PARSER

SCARC_ROUTINE = 'SCARC_INPUT_PARSER'

!> ------------------------------------------------------------------------------------------------
!> set type of discretization
!> ------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_DISCRETIZATION))
   CASE ('STRUCTURED')
      TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED
   CASE ('UNSTRUCTURED')
      TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED
END SELECT
<<<<<<< HEAD
PRES_ON_WHOLE_DOMAIN = (TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED)
=======
PRES_ON_WHOLE_DOMAIN_SCRC = (TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED)
PRES_ON_WHOLE_DOMAIN = PRES_ON_WHOLE_DOMAIN_SCRC
>>>>>>> origin/scarc

!> ------------------------------------------------------------------------------------------------
!> set type of global solver
!> ------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_METHOD))

   CASE ('KRYLOV')

      TYPE_METHOD  = NSCARC_METHOD_KRYLOV
      TYPE_METHOD0 = NSCARC_METHOD_KRYLOV

      !> set type of Krylov-method (CG/BICG)
      SELECT CASE (TRIM(SCARC_KRYLOV))
         CASE ('CG')
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
         CASE ('BICG')
            TYPE_KRYLOV = NSCARC_KRYLOV_BICG
         CASE DEFAULT
            WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_KRYLOV)
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
      END SELECT

      !> set type of preconditioner (JACOBI/SSOR/MG)
      SELECT CASE (TRIM(SCARC_PRECON))
         CASE ('JACOBI')
            TYPE_PRECON = NSCARC_PRECON_JACOBI
         CASE ('SSOR')
            TYPE_PRECON = NSCARC_PRECON_SSOR
         CASE ('MULTIGRID')
            TYPE_PRECON = NSCARC_PRECON_MULTIGRID
            SELECT CASE (TRIM(SCARC_SMOOTH))
               CASE ('JACOBI')
                  TYPE_SMOOTH = NSCARC_SMOOTH_JACOBI
               CASE ('SSOR')
                  TYPE_SMOOTH = NSCARC_SMOOTH_SSOR
               CASE ('FFT')
                  TYPE_PRECON = NSCARC_SMOOTH_FFT
            END SELECT
         CASE ('FFT')
            TYPE_PRECON = NSCARC_PRECON_FFT
#ifdef WITH_MKL
         CASE ('PARDISO')
            TYPE_MKL    = NSCARC_MKL_LOCAL
            TYPE_PRECON = NSCARC_PRECON_PARDISO
         CASE ('CLUSTER')
            TYPE_MKL    = NSCARC_MKL_GLOBAL
            !IF (N_MPI_PROCESSES > 1) THEN
               TYPE_PRECON = NSCARC_PRECON_CLUSTER
            !ELSE     
            !   TYPE_PRECON = NSCARC_PRECON_PARDISO           ! in case of only 1 process, use Pardiso solver
            !ENDIF
#endif
         CASE DEFAULT
            WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_PRECON)
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
      END SELECT
      TYPE_PRECON0 = TYPE_PRECON

   CASE ('MULTIGRID')

      TYPE_METHOD  = NSCARC_METHOD_MULTIGRID
      TYPE_METHOD0 = NSCARC_METHOD_MULTIGRID

      !> set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_MULTIGRID))
         CASE ('GEOMETRIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
         CASE ('ALGEBRAIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
         CASE DEFAULT
            WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_MULTIGRID)
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
      END SELECT

      !> set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_COARSE))
         CASE ('ITERATIVE')
            TYPE_COARSE = NSCARC_COARSE_ITERATIVE
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
         CASE ('DIRECT')
            TYPE_COARSE = NSCARC_COARSE_DIRECT
            TYPE_MKL = NSCARC_MKL_GLOBAL
            IF (SCARC_COARSE_MTYPE == 'SYMMETRIC') THEN
               SCARC_MKL_MTYPE = 'SYMMETRIC'
            ELSE
               SCARC_MKL_MTYPE = 'NONSYMMETRIC'
            ENDIF
      END SELECT

      !> set type of smoother (JACOBI/SSOR)
      SELECT CASE (TRIM(SCARC_SMOOTH))                        !> use same parameters as for preconditioner
         CASE ('JACOBI')
            TYPE_PRECON = NSCARC_PRECON_JACOBI
            TYPE_SMOOTH = NSCARC_SMOOTH_JACOBI
         CASE ('SSOR')
            TYPE_PRECON = NSCARC_PRECON_SSOR
            TYPE_SMOOTH = NSCARC_SMOOTH_SSOR
         CASE ('PARDISO')
            TYPE_PRECON = NSCARC_PRECON_PARDISO
            TYPE_SMOOTH = NSCARC_SMOOTH_PARDISO
         CASE DEFAULT
            WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_SMOOTH)
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
      END SELECT

#ifdef WITH_MKL
   CASE ('MKL')

      TYPE_METHOD  = NSCARC_METHOD_MKL

      !> set type of MKL method (global/local)
      SELECT CASE (TRIM(SCARC_MKL))                      !Achtung, hier noch nacharbeiten!
         CASE ('GLOBAL')
            TYPE_MKL = NSCARC_MKL_GLOBAL
            !IF (N_MPI_PROCESSES > 1) THEN
               TYPE_MKL = NSCARC_MKL_GLOBAL
            !ELSE
            !   TYPE_MKL = NSCARC_MKL_LOCAL
            !ENDIF
         CASE ('LOCAL')
            TYPE_MKL = NSCARC_MKL_LOCAL
         CASE DEFAULT
            WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_MKL)
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
      END SELECT
#endif

   CASE DEFAULT
      WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_METHOD)
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN

END SELECT

!> ------------------------------------------------------------------------------------------------
!> if a multigrid solver is used (either as main solver or as preconditioner)
!> set types for multigrid, coarse grid solver and cycling pattern
!> ------------------------------------------------------------------------------------------------
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_PRECON_MULTIGRID) THEN

   !! set type of multigrid (GEOMETRIC/ALGEBRAIC with corresponding coarsening strategy)
   SELECT CASE (TRIM(SCARC_MULTIGRID))

      CASE ('GEOMETRIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC

      CASE ('ALGEBRAIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC

         !> set type of coarsening strategy (STANDARD/AGGRESSIVE)
         SELECT CASE (TRIM(SCARC_MULTIGRID_COARSENING))
            CASE ('BASIC')
               TYPE_COARSENING = NSCARC_COARSENING_BASIC
               TYPE_LAYER      = NSCARC_LAYER_ONE
            CASE ('FALGOUT')
               TYPE_COARSENING = NSCARC_COARSENING_FALGOUT
               TYPE_LAYER      = NSCARC_LAYER_ONE
            CASE ('RS3')
               TYPE_COARSENING = NSCARC_COARSENING_RS3
               TYPE_LAYER      = NSCARC_LAYER_ONE
            CASE ('A1')
               TYPE_COARSENING = NSCARC_COARSENING_A1
               TYPE_LAYER      = NSCARC_LAYER_ONE
            CASE ('A2')
               TYPE_COARSENING = NSCARC_COARSENING_A2
               TYPE_LAYER      = NSCARC_LAYER_ONE
            CASE ('PMIS')
               TYPE_COARSENING = NSCARC_COARSENING_PMIS
               TYPE_LAYER      = NSCARC_LAYER_TWO
            CASE ('PMISG')
               TYPE_COARSENING = NSCARC_COARSENING_PMISG
               TYPE_LAYER      = NSCARC_LAYER_TWO
            CASE ('FDSRS3')
               TYPE_COARSENING = NSCARC_COARSENING_FDSRS3
               TYPE_LAYER      = NSCARC_LAYER_ONE
            CASE ('FDSA1')
               TYPE_COARSENING = NSCARC_COARSENING_FDSA1
               TYPE_LAYER      = NSCARC_LAYER_TWO
            CASE ('FDSA2')
               TYPE_COARSENING = NSCARC_COARSENING_FDSA2
               TYPE_LAYER      = NSCARC_LAYER_TWO
            CASE ('GMG')
               TYPE_COARSENING = NSCARC_COARSENING_GMG
               TYPE_LAYER      = NSCARC_LAYER_ONE
            CASE ('GMG3')
               TYPE_COARSENING = NSCARC_COARSENING_GMG3
               TYPE_LAYER      = NSCARC_LAYER_ONE
            CASE DEFAULT
               WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_MULTIGRID_COARSENING)
               CALL SHUTDOWN(SCARC_MESSAGE); RETURN
         END SELECT

      CASE DEFAULT
         WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_MULTIGRID)
         CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   END SELECT

   !! set type of coarse grid solver (CG/GE)
   SELECT CASE (TRIM(SCARC_COARSE))
      CASE ('ITERATIVE')
         TYPE_COARSE = NSCARC_COARSE_ITERATIVE
         TYPE_KRYLOV = NSCARC_KRYLOV_CG
      CASE ('DIRECT')
         TYPE_COARSE = NSCARC_COARSE_DIRECT
#ifdef WITH_MKL
         TYPE_MKL = NSCARC_MKL_COARSE
#endif
      CASE DEFAULT
         WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_COARSE)
         CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   END SELECT

   !! set type of cycling pattern (F/V/W)
   SELECT CASE (TRIM(SCARC_MULTIGRID_CYCLE))
      CASE ('F')
         TYPE_CYCLE = NSCARC_CYCLE_F
      CASE ('V')
         TYPE_CYCLE = NSCARC_CYCLE_V
      CASE ('W')
         TYPE_CYCLE = NSCARC_CYCLE_W
      CASE DEFAULT
         WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_MULTIGRID_CYCLE)
         CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   END SELECT

   !! set type of interpolation (STANDARD/DIRECT/MULTIPASS)
   SELECT CASE (TRIM(SCARC_MULTIGRID_INTERPOL))
      CASE ('STANDARD')
         TYPE_INTERPOL = NSCARC_INTERPOL_STANDARD
      CASE ('CLASSICAL')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL
      CASE ('CLASSICAL2')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL2
      CASE ('DIRECT')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT
      CASE ('DIRECT_BDRY')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT_BDRY
      CASE ('MULTIPASS')
         TYPE_INTERPOL = NSCARC_INTERPOL_MULTIPASS
      CASE ('GMG')
         TYPE_INTERPOL = NSCARC_INTERPOL_GMG
      CASE ('GMG3')
         TYPE_INTERPOL = NSCARC_INTERPOL_GMG3
      CASE DEFAULT
         WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_MULTIGRID_INTERPOL)
         CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   END SELECT

ENDIF

!> ------------------------------------------------------------------------------------------------
!> set test case
!> ------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_CASE))
   CASE ('CD_NSA_2D')
      TYPE_CASE = NSCARC_CASE_CD_NSA_2D
   CASE ('VD_NSA_2D')
      TYPE_CASE = NSCARC_CASE_CD_NSA_2D
   CASE ('CD_VA_2D')
      TYPE_CASE = NSCARC_CASE_CD_VA_2D
   CASE ('ZM_GRAV_ADV_2D')
      TYPE_CASE = NSCARC_CASE_ZM_GRAV_ADV_2D
   CASE DEFAULT
      TYPE_CASE = NSCARC_CASE_NONE
END SELECT


!> ------------------------------------------------------------------------------------------------
!> set type of accuracy (ABSOLUTE/RELATIVE)
!> ------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_ACCURACY))
   CASE ('ABSOLUTE')
      TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
   CASE ('RELATIVE')
      TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
   CASE DEFAULT
      WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_ACCURACY)
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
END SELECT

!> ------------------------------------------------------------------------------------------------
!> set level of debugging (NONE/LESS/MEDIUM/MUCH)
!> ------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_DEBUG))
   CASE ('NONE')
      TYPE_DEBUG = NSCARC_DEBUG_NONE
   CASE ('INFO0')
      TYPE_DEBUG = NSCARC_DEBUG_INFO0
   CASE ('INFO1')
      TYPE_DEBUG = NSCARC_DEBUG_INFO1
   CASE ('INFO2')
      TYPE_DEBUG = NSCARC_DEBUG_INFO2
   CASE ('LESS')
      TYPE_DEBUG = NSCARC_DEBUG_LESS
   CASE ('MEDIUM')
      TYPE_DEBUG = NSCARC_DEBUG_MEDIUM
   CASE ('MUCH')
      TYPE_DEBUG = NSCARC_DEBUG_MUCH
   CASE ('EXTREME')
      TYPE_DEBUG = NSCARC_DEBUG_EXTREME
   CASE DEFAULT
      WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_DEBUG)
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
END SELECT

!> ------------------------------------------------------------------------------------------------
!> set type of initial solution (not yet used, may be used to define own initial vector)
!> ------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_INITIAL))
   CASE ('null')
      TYPE_INITIAL = NSCARC_INITIAL_NONE
   CASE DEFAULT
      TYPE_INITIAL = NSCARC_INITIAL_NONE
END SELECT

END SUBROUTINE SCARC_INPUT_PARSER




!> ------------------------------------------------------------------------------------------------
!> Determine number of grid levels  (1 for CG/BICG-method, NLEVEL for MG-method)
!> Note: NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LEVELS
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM

SELECT CASE (TYPE_MULTIGRID)

   !! ----------------------------------------------------------------------------------------------
   !! predefined hierarchy of levels in case of geometric multigrid-method
   !! ----------------------------------------------------------------------------------------------
   CASE (NSCARC_MULTIGRID_GEOMETRIC)

      NLEVEL = NSCARC_LEVEL_MAX
      KLEVEL = NSCARC_LEVEL_MAX

      DO NM=1,NMESHES

         KLEVEL(1)=SCARC_GET_MAXLEVEL(MESHES(NM)%IBAR,1)
         IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) KLEVEL(2)=SCARC_GET_MAXLEVEL(MESHES(NM)%JBAR,2)
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

      IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
         WRITE(LU_SCARC,*) 'KLEVEL=',KLEVEL
         WRITE(LU_SCARC,*) 'NLEVEL=',NLEVEL
         WRITE(LU_SCARC,*) 'NLEVEL_MAX=',NLEVEL_MAX
         WRITE(LU_SCARC,*) 'TYPE_MULTIGRID=',TYPE_MULTIGRID
         WRITE(LU_SCARC,*) 'TYPE_COARSE=',TYPE_COARSE
      ENDIF

   !! ----------------------------------------------------------------------------------------------
   !! first, only finest level is set, further levels are defined during coarsening process
   !! ----------------------------------------------------------------------------------------------
   CASE (NSCARC_MULTIGRID_ALGEBRAIC)

      NLEVEL_MIN = 1
      IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
         NLEVEL_MAX  = SCARC_MULTIGRID_LEVEL
      ELSE
         NLEVEL_MAX  = NSCARC_LEVEL_MAX
      ENDIF
      NLEVEL = SCARC_MULTIGRID_LEVEL


   !! ----------------------------------------------------------------------------------------------
   !! no multigrid-hierachy needed in case of a pure Krylov-method: use only one level
   !! ----------------------------------------------------------------------------------------------
   CASE DEFAULT

      NLEVEL     = 1
      NLEVEL_MIN = 1
      NLEVEL_MAX = 1

END SELECT

END SUBROUTINE SCARC_SETUP_LEVELS


!> ------------------------------------------------------------------------------------------------
!> Determine maximum number of possible levels on direction IOR0 of mesh NM
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_GET_MAXLEVEL(NC, IOR0)
INTEGER, INTENT(IN) :: NC, IOR0
INTEGER :: NC0, NL

SCARC_ROUTINE = 'SCARC_GET_MAXLEVEL'

!> In case of the GMG-method, NC must be divisable by 2 at least one time
IF (MOD(NC,2)/=0 .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
   SELECT CASE (IOR0)
      CASE (1)
         WRITE(*,1000) 'IBAR', NC
      CASE (2)
         WRITE(*,1000) 'JBAR', NC
      CASE (3)
         WRITE(*,1000) 'KBAR', NC
   END SELECT
   WRITE(SCARC_MESSAGE,'(2A)') TRIM(SCARC_ROUTINE),': Step size not divisable by 2'
   CALL SHUTDOWN(SCARC_MESSAGE); RETURN
ENDIF

!> Divide by 2 as often as possible or till user defined max-level is reached
NC0=NC
DO NL=1,NSCARC_LEVEL_MAX
   NC0=NC0/2
   IF (MOD(NC0,2)/=0) EXIT                !> NC no longer divisable by two
   IF (NL==SCARC_MULTIGRID_LEVEL) EXIT    !> max number of levels defined by user
   IF (NC0==1) EXIT                       !> NC is power of two, minimum has been reached
ENDDO

SCARC_GET_MAXLEVEL=NL
RETURN
1000 FORMAT(A,'=',I3,' must be divisable by two for ScaRC-Multigrid!')
END FUNCTION SCARC_GET_MAXLEVEL


!> ------------------------------------------------------------------------------------------------
!> Allocate ScaRC-structures for all needed levels
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TYPES
INTEGER :: IERR, NM
TYPE (SCARC_TYPE),       POINTER :: S

!> ------------------------------------------------------------------------------------------------
!> Allocate global ScaRC-structure
!> ------------------------------------------------------------------------------------------------
ALLOCATE (SCARC(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERR)


!> ------------------------------------------------------------------------------------------------
!> Allocate local OSCARC, MESHES and PRECON structures for single meshes
!> ------------------------------------------------------------------------------------------------
MESHES_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   S => SCARC(NM)

   !> allocate Scarc-types for other meshes
   ALLOCATE (S%OSCARC(NMESHES), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSCARC', IERR)

   !> allocate Scarc-types for different levels
   ALLOCATE (SCARC(NM)%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'LEVEL', IERR)

   !> allocate additional FFT-administration type in case of FFT-preconditioning
   IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
      ALLOCATE (SCARC(NM)%PRECON(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'PRECON', IERR)
   ENDIF

#ifdef WITH_MKL
   !> allocate additional MKL-administration type if MKL-Pardiso-solver is used
   IF (TYPE_MKL == NSCARC_MKL_LOCAL .OR. TYPE_MKL == NSCARC_MKL_GLOBAL) THEN
<<<<<<< HEAD
      ALLOCATE (SCARC(NM)%MKL(NLEVEL_MIN:NLEVEL_MIN), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'MKL', IERR)
   ELSE IF (TYPE_MKL == NSCARC_MKL_COARSE) THEN
=======
      IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'A: ALLOCATING MKL FOR LEVEL ', NLEVEL_MAX
      ALLOCATE (SCARC(NM)%MKL(NLEVEL_MIN:NLEVEL_MIN), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'MKL', IERR)
   ELSE IF (TYPE_MKL == NSCARC_MKL_COARSE) THEN
      IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'B: ALLOCATING MKL FOR LEVEL ', NLEVEL_MAX
>>>>>>> origin/scarc
      ALLOCATE (SCARC(NM)%MKL(NLEVEL_MAX:NLEVEL_MAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'MKL', IERR)
   ENDIF
#endif

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_TYPES


!> ----------------------------------------------------------------------------------------------------
!> Setup geometry information for mesh NM
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MESHES
INTEGER  :: IERR, NL, NM, IX, IY, IZ, IO
INTEGER  :: NX0, NY0, NZ0
TYPE (MESH_TYPE)       , POINTER :: M
TYPE (SCARC_TYPE)      , POINTER :: S
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

IERR=0

LEVEL_MESHES_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   M => MESHES(NM)
   S => SCARC(NM)

   !DO IX=1,M%N_EXTERNAL_WALL_CELLS
   !   WRITE(LU_SCARC,'(a,7i4)') '0:  IX=',IX,M%EXTERNAL_WALL(IX)%IIO_MIN,M%EXTERNAL_WALL(IX)%JJO_MIN,M%EXTERNAL_WALL(IX)%KKO_MIN,&
   !                                          M%EXTERNAL_WALL(IX)%IIO_MAX,M%EXTERNAL_WALL(IX)%JJO_MAX,M%EXTERNAL_WALL(IX)%KKO_MAX
   !ENDDO

   !> store bounds of mesh in SCARC-structure
   S%XS = M%XS                 
   S%XF = M%XF                 
   S%YS = M%YS                 
   S%YF = M%YF                 
   S%ZS = M%ZS                 
   S%ZF = M%ZF                 

   S%NX = M%IBAR
   S%NY = M%JBAR
   S%NZ = M%KBAR

   NX0 = M%IBAR
   NY0 = M%JBAR
   NZ0 = M%KBAR

   IF (N_MPI_PROCESSES>1) THEN

      CALL MPI_ALLREDUCE(S%XS,S%XS_MIN, 1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL SCARC_HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for S%XS',IERR)
      CALL MPI_ALLREDUCE(S%XF,S%XF_MAX, 1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL SCARC_HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for S%XF',IERR)

      CALL MPI_ALLREDUCE(S%YS,S%YS_MIN, 1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL SCARC_HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for S%YS',IERR)
      CALL MPI_ALLREDUCE(S%YF,S%YF_MAX, 1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL SCARC_HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for S%YF',IERR)

      CALL MPI_ALLREDUCE(S%ZS,S%ZS_MIN, 1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL SCARC_HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for S%ZS',IERR)
      CALL MPI_ALLREDUCE(S%ZF,S%ZF_MAX, 1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL SCARC_HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for S%ZF',IERR)

   ENDIF

   LEVEL_LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

      IF (NL > NLEVEL_MIN .AND. TYPE_MULTIGRID /= NSCARC_MULTIGRID_GEOMETRIC) EXIT LEVEL_LEVEL_LOOP

      !> let SM point to SCARC(NM)%MESHES(NL)
      SL => S%LEVEL(NL)

      !> numbers of cells in x-, y- and z-direction for level 'NL'
      SL%NX = NX0
      SL%NY = NY0
      SL%NZ = NZ0

      SL%NC = SL%NX * SL%NY * SL%NZ

      SL%N_EXTERNAL_WALL_CELLS = M%N_EXTERNAL_WALL_CELLS
      SL%N_INTERNAL_WALL_CELLS = M%N_INTERNAL_WALL_CELLS
      SL%NW = SL%N_EXTERNAL_WALL_CELLS + SL%N_INTERNAL_WALL_CELLS 

      IF (NL == NLEVEL_MIN) THEN
         SL%N_OBST = M%N_OBST
         ALLOCATE(SL%OBST(SL%N_OBST), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','OBST',IERR)
         DO IO = 1, SL%N_OBST
<<<<<<< HEAD
            SL%OBST(IO)%I1  = M%OBSTRUCTION(IO)%I1
            SL%OBST(IO)%I2  = M%OBSTRUCTION(IO)%I2
            SL%OBST(IO)%J1  = M%OBSTRUCTION(IO)%J1
            SL%OBST(IO)%J2  = M%OBSTRUCTION(IO)%J2
            SL%OBST(IO)%K1  = M%OBSTRUCTION(IO)%K1
            SL%OBST(IO)%K2  = M%OBSTRUCTION(IO)%K2
            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
               WRITE(LU_SCARC,*) 'SHOWING OBSTRUCTIONS'
               WRITE(LU_SCARC,*) IO, &
                                 M%OBSTRUCTION(IO)%I1, M%OBSTRUCTION(IO)%I2, &
                                 M%OBSTRUCTION(IO)%J1, M%OBSTRUCTION(IO)%J2, &
                                 M%OBSTRUCTION(IO)%K1, M%OBSTRUCTION(IO)%K2
            ENDIF
=======
            SL%OBST(IO)%I1  = M%OBSTRUCTION(IO)%I1+1
            SL%OBST(IO)%I2  = M%OBSTRUCTION(IO)%I2
            SL%OBST(IO)%J1  = M%OBSTRUCTION(IO)%J1+1
            SL%OBST(IO)%J2  = M%OBSTRUCTION(IO)%J2
            SL%OBST(IO)%K1  = M%OBSTRUCTION(IO)%K1+1
            SL%OBST(IO)%K2  = M%OBSTRUCTION(IO)%K2
         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
            WRITE(LU_SCARC,*) 'SHOWING OBSTRUCTIONS'
            WRITE(LU_SCARC,*) IO, &
                              M%OBSTRUCTION(IO)%I1, M%OBSTRUCTION(IO)%I2, &
                              M%OBSTRUCTION(IO)%J1, M%OBSTRUCTION(IO)%J2, &
                              M%OBSTRUCTION(IO)%K1, M%OBSTRUCTION(IO)%K2
         ENDIF
>>>>>>> origin/scarc
         ENDDO
      ENDIF

      !> step widths in x-, y- and z-direction for level 'NL'
      SL%DX = (S%XF-S%XS)/REAL(SL%NX,EB)
      SL%DY = (S%YF-S%YS)/REAL(SL%NY,EB)
      SL%DZ = (S%ZF-S%ZS)/REAL(SL%NZ,EB)

      SL%DXI = 1.0_EB/SL%DX
      SL%DYI = 1.0_EB/SL%DY
      SL%DZI = 1.0_EB/SL%DZ

      SL%DXI2 = SL%DXI**2
      SL%DYI2 = SL%DYI**2
      SL%DZI2 = SL%DZI**2

      SL%DI2(1) = SL%DXI2
      SL%DI2(2) = SL%DYI2
      SL%DI2(3) = SL%DZI2

      !> needed in case of GMG with multiple grid levels
      NX0=NX0/2
      IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) NY0=NY0/2
      NZ0=NZ0/2

      !> Allocate vectors for coordinate information
      IF (NL == NLEVEL_MIN) THEN

         SL%XCOR => M%X
         SL%YCOR => M%Y
         SL%ZCOR => M%Z

         SL%XMID => M%XC
         SL%YMID => M%YC
         SL%ZMID => M%ZC

      ELSE

         !> compute coordinates in x-, y- and z-direction
         ALLOCATE(SL%XCOR(0:SL%NX), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','XCOR',IERR)
         ALLOCATE(SL%YCOR(0:SL%NY), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','YCOR',IERR)
         ALLOCATE(SL%ZCOR(0:SL%NZ), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','ZCOR',IERR)

         DO IX = 0, SL%NX
            SL%XCOR(IX) = S%XS + IX*SL%DX
         ENDDO
         DO IY = 0, SL%NY
            SL%YCOR(IY) = S%YS + IY*SL%DY
         ENDDO
         DO IZ = 0, SL%NZ
            SL%ZCOR(IZ) = S%ZS + IZ*SL%DZ
         ENDDO

         !> compute midpoints in x-, y- and z-direction
         ALLOCATE(SL%XMID(0:SL%NX+1), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','XMID',IERR)
         ALLOCATE(SL%YMID(0:SL%NY+1), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','YMID',IERR)
         ALLOCATE(SL%ZMID(0:SL%NZ+1), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','ZMID',IERR)

         SL%XMID(0) = S%XS - 0.5_EB*SL%DX
         DO IX = 1, SL%NX
            SL%XMID(IX) = 0.5_EB*(SL%XCOR(IX-1) + SL%XCOR(IX))
         ENDDO
         SL%XMID(SL%NX+1) = S%XF + 0.5_EB*SL%DX

         SL%YMID(0) = S%YS - 0.5_EB*SL%DY
         DO IY = 1, SL%NY
            SL%YMID(IY) = 0.5_EB*(SL%YCOR(IY-1) + SL%YCOR(IY))
         ENDDO
         SL%YMID(SL%NY+1) = S%YF + 0.5_EB*SL%DY

         SL%ZMID(0) = S%ZS - 0.5_EB*SL%DZ
         DO IZ = 1, SL%NZ
            SL%ZMID(IZ) = 0.5_EB*(SL%ZCOR(IZ-1) + SL%ZCOR(IZ))
         ENDDO
         SL%ZMID(SL%NZ+1) = S%ZF + 0.5_EB*SL%DZ

      ENDIF

      !> Allocate vectors for step sizes in different directions
      ALLOCATE(SL%DXL(0:SL%NX), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','DXL',IERR)
      SL%DXL = 0.0_EB
      ALLOCATE(SL%DYL(0:SL%NY), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','DYL',IERR)
      SL%DYL = 0.0_EB
      ALLOCATE(SL%DZL(0:SL%NZ), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','DZL',IERR)
      SL%DZL = 0.0_EB

      !> set step sizes between cell midpoints, use interior step sizes for ghost cells as initial values
      !> correct sizes for ghost cells are exchanged later
      DO IX = 1, SL%NX-1
         SL%DXL(IX) = SL%XMID(IX+1) - SL%XMID(IX)
      ENDDO
      SL%DXL(0)     = SL%DXL(1)
      SL%DXL(SL%NX) = SL%DXL(SL%NX-1)

      DO IY = 1, SL%NY-1
         SL%DYL(IY) = SL%YMID(IY+1) - SL%YMID(IY)
      ENDDO
      SL%DYL(0)     = SL%DYL(1)
      SL%DYL(SL%NY) = SL%DYL(SL%NY-1)

      DO IZ = 1, SL%NZ-1
         SL%DZL(IZ) = SL%ZMID(IZ+1) - SL%ZMID(IZ)
      ENDDO
      SL%DZL(0)     = SL%DZL(1)
      SL%DZL(SL%NZ) = SL%DZL(SL%NZ-1)

<<<<<<< HEAD
      IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
      IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
         WRITE(LU_SCARC,*) '============ SETUP_MESHES ==============='
         WRITE(LU_SCARC,*) 'NX     =',SL%NX
         WRITE(LU_SCARC,*) 'NY     =',SL%NY
         WRITE(LU_SCARC,*) 'NZ     =',SL%NZ
         WRITE(LU_SCARC,*) 'NC     =',SL%NC
         WRITE(LU_SCARC,*) 'N_EXTERNAL_WALL_CELLS     =',SL%N_EXTERNAL_WALL_CELLS
         WRITE(LU_SCARC,*) 'N_INTERNAL_WALL_CELLS     =',SL%N_INTERNAL_WALL_CELLS
         WRITE(LU_SCARC,*) 'NW     =',SL%NW
         WRITE(LU_SCARC,*) 'XS     =',S%XS
         WRITE(LU_SCARC,*) 'XF     =',S%XF
         WRITE(LU_SCARC,*) 'YS     =',S%YS
         WRITE(LU_SCARC,*) 'YF     =',S%YF
         WRITE(LU_SCARC,*) 'ZS     =',S%ZS
         WRITE(LU_SCARC,*) 'ZF     =',S%ZF
         WRITE(LU_SCARC,*) 'DX     =',SL%DX
         WRITE(LU_SCARC,*) 'DY     =',SL%DY
         WRITE(LU_SCARC,*) 'DZ     =',SL%DZ
         WRITE(LU_SCARC,*) 'DXI    =',SL%DXI
         WRITE(LU_SCARC,*) 'DYI    =',SL%DYI
         WRITE(LU_SCARC,*) 'DZI    =',SL%DZI
         WRITE(LU_SCARC,*) 'DXI2   =',SL%DXI2
         WRITE(LU_SCARC,*) 'DYI2   =',SL%DYI2
         WRITE(LU_SCARC,*) 'DZI2   =',SL%DZI2
         WRITE(LU_SCARC,*) 'DI2(1) =',SL%DI2(1)
         WRITE(LU_SCARC,*) 'DI2(2) =',SL%DI2(2)
         WRITE(LU_SCARC,*) 'DI2(3) =',SL%DI2(3)
         WRITE(LU_SCARC,'(a,20f8.3)') 'XCOR=',SL%XCOR(0:SL%NX)
         WRITE(LU_SCARC,'(a,20f8.3)') 'YCOR=',SL%YCOR(0:SL%NY)
         WRITE(LU_SCARC,'(a,20f8.3)') 'ZCOR=',SL%ZCOR(0:SL%NZ)
         WRITE(LU_SCARC,'(a,20f8.3)') 'XMID=',SL%XMID(0:SL%NX+1)
         WRITE(LU_SCARC,'(a,20f8.3)') 'YMID=',SL%YMID(0:SL%NY+1)
         WRITE(LU_SCARC,'(a,20f8.3)') 'ZMID=',SL%ZMID(0:SL%NZ+1)
         WRITE(LU_SCARC,'(a,20f8.3)') 'DXL =',SL%DXL(0:SL%NX)
         WRITE(LU_SCARC,'(a,20f8.3)') 'DYL =',SL%DYL(0:SL%NY)
         WRITE(LU_SCARC,'(a,20f8.3)') 'DZL =',SL%DZL(0:SL%NZ)
      ENDIF

   ENDDO LEVEL_LEVEL_LOOP

ENDDO LEVEL_MESHES_LOOP


END SUBROUTINE SCARC_SETUP_MESHES


!> ------------------------------------------------------------------------------------------------
!> Setup communication structure for data exchange
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERFACES
INTEGER :: NM, NOM, NL
INTEGER :: IERR
TYPE (MESH_TYPE)          , POINTER :: M
TYPE (SCARC_TYPE)         , POINTER :: S
TYPE (OMESH_TYPE)         , POINTER :: OM
TYPE (OSCARC_TYPE)        , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLF, SLC, OSLF, OSLC

IERR = 0

!> Initialize communication counter for ScaRC, use same TAG for all communications
NREQ_SCARC  =  0
N_EXCHANGES =  0
TAG_SCARC   = 99

!> ------------------------------------------------------------------------------------------------
!> Allocate basic WALL and FACE types on mesh NM for all requested grid levels
!> ------------------------------------------------------------------------------------------------
MESHES_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   M => MESHES(NM)
   S => SCARC(NM)

   SLF => S%LEVEL(NLEVEL_MIN)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'ALLOCATING WALL IN LENGTH ', SLF%NW

   !! Allocate basic wall and face type for finest level
   ALLOCATE(SLF%WALL(SLF%NW), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_INTERFACES','WALL',IERR)

   ALLOCATE(SLF%FACE(-3:3), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_INTERFACES','FACE',IERR)

   !> Allocate basic wall and face type for coarser levels if requested
   IF (NLEVEL_MAX > NLEVEL_MIN .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN

      DO NL=NLEVEL_MIN+1,NLEVEL_MAX

         SLC => S%LEVEL(NL)

         ALLOCATE(SLC%FACE(-3:3), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_INTERFACES','WALL',IERR)
      ENDDO

   ENDIF

   !>
   !> Get communication lengths for other meshes
   !>
   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      IF (NOM == NM) CYCLE OTHER_MESHES_LOOP

      OM => M%OMESH(NOM)
      OS => S%OSCARC(NOM)

      OS%NIC_R    = OM%NIC_R
      OS%NICMAX_R = OM%NIC_R

      OS%NIC_S    = OM%NIC_S
      OS%NICMAX_S = OM%NIC_S

<<<<<<< HEAD
      IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
>>>>>>> origin/scarc
         WRITE(LU_SCARC,*) 'NIC_R=',OS%NIC_R,': NICMAX_R=',OS%NICMAX_R
         WRITE(LU_SCARC,*) 'NIC_S=',OS%NIC_R,': NICMAX_R=',OS%NICMAX_R
      ENDIF

      IF (OS%NICMAX_R==0 .AND. OS%NICMAX_S==0)  CYCLE OTHER_MESHES_LOOP

      N_EXCHANGES  = N_EXCHANGES+1

   ENDDO OTHER_MESHES_LOOP
ENDDO MESHES_LOOP


!> ------------------------------------------------------------------------------------------------
!> Initialize level structures on neighboring meshes
!> ------------------------------------------------------------------------------------------------
LEVEL_MESHES_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   M => MESHES(NM)
   S => SCARC(NM)

   LEVEL_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      IF (NOM == NM) CYCLE LEVEL_OTHER_MESHES_LOOP

      OS => S%OSCARC(NOM)
      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0)  CYCLE LEVEL_OTHER_MESHES_LOOP

      !> Allocate OSCARC level structure for mesh NM
!WRITE(*,*) NM,': ALLOCATING SCARC(',NM,')%OSCARC(',NOM,')%LEVEL(',NLEVEL_MIN,':',NLEVEL_MAX,')'

      ALLOCATE (S%OSCARC(NOM)%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OS%LEVEL', IERR)

      !> point to grid structure of OSCARC(NM) on finest level
      OSLF => S%OSCARC(NOM)%LEVEL(NLEVEL_MIN)

      OSLF%NX = MESHES(NOM)%IBAR
      OSLF%NY = MESHES(NOM)%JBAR
      OSLF%NZ = MESHES(NOM)%KBAR

      OSLF%N_EXTERNAL_WALL_CELLS = MESHES(NOM)%N_EXTERNAL_WALL_CELLS
      OSLF%N_INTERNAL_WALL_CELLS = MESHES(NOM)%N_INTERNAL_WALL_CELLS

      OSLF%NW  = OSLF%N_EXTERNAL_WALL_CELLS + OSLF%N_INTERNAL_WALL_CELLS
      OSLF%NC  = OSLF%NX*OSLF%NY*OSLF%NZ
      OSLF%NA  = 0
      OSLF%NCG = 0

      IF (S%OSCARC(NOM)%NICMAX_S == 0 .AND. S%OSCARC(NOM)%NICMAX_R == 0) CYCLE LEVEL_OTHER_MESHES_LOOP

      !> In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN

         DO NL=NLEVEL_MIN+1,NLEVEL_MAX

            OSLC => S%OSCARC(NOM)%LEVEL(NL)                          !> pointer to coarser level
            OSLF => S%OSCARC(NOM)%LEVEL(NL-1)                        !> pointer to finer level

            !> get number of internal cells and external wall cells on neighbor NOM for level NL
            OSLC%NX = OSLF%NX/2
            SELECT CASE (TYPE_DIMENSION)
               CASE (NSCARC_DIMENSION_TWO)
                  OSLC%NY = 1
               CASE (NSCARC_DIMENSION_THREE)
                  OSLC%NY = OSLF%NY/2
            END SELECT
            OSLC%NZ = OSLF%NZ/2

            OSLC%N_EXTERNAL_WALL_CELLS = 2*OSLC%NX * 2*OSLC%NZ + &
                                         2*OSLC%NX * 2*OSLC%NY + &
                                         2*OSLC%NY * 2*OSLC%NZ 

            OSLC%NC  = OSLC%NX * OSLC%NY * OSLC%NZ
            OSLC%NW  = OSLC%N_EXTERNAL_WALL_CELLS 
            OSLC%NA  = 0
            OSLC%NCG = 0

            !OSLC%N_INTERNAL_WALL_CELLS = OSLF%N_INTERNAL_WALL_CELLS/2              ! FALSCH !!!
            !OSLC%NW  = OSLC%N_EXTERNAL_WALL_CELLS + OSLC%N_INTERNAL_WALL_CELLS     ! GEBRAUCHT ??

         ENDDO
      ENDIF

   ENDDO LEVEL_OTHER_MESHES_LOOP
ENDDO LEVEL_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_INTERFACES


!> ------------------------------------------------------------------------------------------------
!> Initialize arrays for data exchange
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE
INTEGER :: NL, NM, NOM, NLEN_SEND, NLEN_RECV
INTEGER :: NLMIN, NLMAX, IX, IY, IZ
INTEGER :: IERR=0, INBR
TYPE ( SCARC_TYPE)     , POINTER :: S
TYPE (OSCARC_TYPE)     , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, OSL

!>
!> Allocate request array for data exchanges
!>
IERR = 0
IF (NMESHES>1) THEN
   ALLOCATE (REQ_SCARC(N_EXCHANGES*40), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'REQ_SCARC', IERR)
   REQ_SCARC = MPI_REQUEST_NULL
ENDIF

!>
!> Exchange basic information about wall sizes (needed for the dimensioning of !the exchange buffers)
!>
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_BASIC , NLEVEL_MIN)

!>
!> Allocate send and receive buffers (real and integer) in correct lengths
!>
DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   S => SCARC(NM)
   DO INBR = 1, S%NUM_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)

      SL  => SCARC(NM)%LEVEL(NLEVEL_MIN)
      OS  => SCARC(NM)%OSCARC(NOM)
      OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)

      !> -------------------------- real buffers
      NLEN_SEND = NSCARC_MAX_STENCIL*MAX(OSL%NWL, OSL%NCG)+10
      NLEN_RECV = NLEN_SEND

      ALLOCATE (OS%SEND_REAL(NLEN_SEND), STAT = IERR)
      CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'SEND_REAL', IERR)
      OS%SEND_REAL = 0.0_EB

      ALLOCATE (OS%RECV_REAL(NLEN_RECV), STAT = IERR)
      CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'RECV_REAL', IERR)
      OS%RECV_REAL = 0.0_EB

<<<<<<< HEAD
      IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
      IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
>>>>>>> origin/scarc
         WRITE(LU_SCARC,*) 'NOM=',NOM,': OSL%NWL=',OSL%NWL
         WRITE(LU_SCARC,*) 'NOM=',NOM,': OSL%NCG=',OSL%NCG
         WRITE(LU_SCARC,*) 'NOM=',NOM,': SL%NCPL_MAX=',SL%NCPL_MAX
         !WRITE(LU_SCARC,*) 'NOM=',NOM,': ALLOCATING REAL    SEND-VECTORS in length ',NLEN_SEND
         !WRITE(LU_SCARC,*) 'NOM=',NOM,': ALLOCATING REAL    RECV-VECTORS in length ',NLEN_RECV
      ENDIF

      !> -------------------------- integer buffers
      OSL%NCPLS = OSL%NCPL
      IF (OSL%NCG == OSL%NWL) THEN
         NLEN_SEND = 15*OSL%NWL
         NLEN_RECV = 15*OSL%NCG
         OSL%NCPLR = OSL%NCPL
      ELSE IF (OSL%NCG == 2*OSL%NWL) THEN
         NLEN_SEND = 13*OSL%NWL + 2*OSL%NCG
         NLEN_RECV = 15*OSL%NCG
         OSL%NCPLR =  1
      ELSE IF (OSL%NWL == 2*OSL%NCG) THEN
         NLEN_SEND = 15*OSL%NWL
         NLEN_RECV = 13*OSL%NCG + 2*OSL%NWL
         OSL%NCPLR =  2
      ENDIF

      ALLOCATE (OS%SEND_INTEGER(NLEN_SEND), STAT = IERR)
      CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'SEND_INTEGER', IERR)
      OS%SEND_INTEGER = 0

      ALLOCATE (OS%RECV_INTEGER(NLEN_RECV), STAT = IERR)
      CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'RECV_INTEGER', IERR)
      OS%RECV_INTEGER = 0

<<<<<<< HEAD
      IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
      IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
>>>>>>> origin/scarc
         WRITE(LU_SCARC,*) 'NOM=',NOM,': ALLOCATING INTEGER SEND-INTEGER in length ',NLEN_SEND
         WRITE(LU_SCARC,*) 'NOM=',NOM,': ALLOCATING INTEGER RECV-INTEGER in length ',NLEN_SEND
         WRITE(LU_SCARC,*) 'NOM=',NOM,': ALLOCATING WALL in length ',OSL%NCG
      ENDIF

      !> -------------------- neighboring wall (and face?) structures for common wall cells
      ALLOCATE (OSL%WALL(OSL%NCG), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OSL%WALL', IERR)

      !ALLOCATE (OSL%FACE(-3:3), STAT=IERR)
      !CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OSL%FACE', IERR)

      !> In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN
         DO NL=NLEVEL_MIN+1,NLEVEL_MAX

            OSL => S%OSCARC(NOM)%LEVEL(NL)                          !> pointer to coarser level

            ALLOCATE (OSL%WALL(OSL%NCG), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OSL%WALL', IERR)

            !ALLOCATE (OSL%FACE(-3:3), STAT=IERR)
            !CALL CHKMEMERR ('SCARC_SETUP_INTERFACES', 'OSL%FACE', IERR)

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
   IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) NLMAX = NLEVEL_MAX
   DO NL = NLMIN, NLMAX
<<<<<<< HEAD
      IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'EXCHANGE_CCVAR FOR LEVEL ', NL
=======
      IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'EXCHANGE_CCVAR FOR LEVEL ', NL
>>>>>>> origin/scarc
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CCVAR, NL)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WIDTHINFO, NL)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WALLINFO, NL)
   ENDDO
ENDIF


!>
!> Correct boundary types for cells adjacent to obstructions on ghost cells
!>
DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE
<<<<<<< HEAD
   CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NLEVEL_MIN)
   IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX
         CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
      ENDDO
   ENDIF
ENDDO
=======
   CALL SCARC_CORRECT_GHOSTCELL_TYPES(NM, NLEVEL_MIN)
   IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX
         CALL SCARC_CORRECT_GHOSTCELL_TYPES(NM, NL)
      ENDDO
   ENDIF
ENDDO

!>
!> ONLY TEMPORARILY: Debug CCVAR-info
!>
DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   NLMIN = NLEVEL_MIN
   NLMAX = NLEVEL_MIN
   IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) NLMAX = NLEVEL_MAX

   DO NL = NLMIN, NLMAX
      IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(77,*) '============================================================='
         WRITE(77,*) '     SETUP_CCVAR_EXCHANGE ', NL
         WRITE(77,*) '============================================================='
         IF (SL%NX==8) THEN
            WRITE(77,*) 'CCVAR(.,.,.,CGSC)'
            WRITE(77,'(10I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
            WRITE(77,*) 'CCVAR(.,.,.,UNKH)'
            WRITE(77,'(10I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         ELSE
            WRITE(77,*) 'CCVAR(.,.,.,CGSC)'
            WRITE(77,'(6I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
            WRITE(77,*) 'CCVAR(.,.,.,UNKH)'
            WRITE(77,'(6I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         ENDIF
         WRITE(77,*) 'NUNKH_LOC=', SL%NUNKH_LOC
         WRITE(77,*) 'NCS=', SL%NCS
         WRITE(77,*) 'NX =', SL%NX
         WRITE(77,*) 'NY =', SL%NY
         WRITE(77,*) 'NZ =', SL%NZ
       ENDIF
   ENDDO

ENDDO

END SUBROUTINE SCARC_SETUP_EXCHANGE
>>>>>>> origin/scarc

!>
!> ONLY TEMPORARILY: Debug CCVAR-info
!>
DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   NLMIN = NLEVEL_MIN
   NLMAX = NLEVEL_MIN
   IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) NLMAX = NLEVEL_MAX

   DO NL = NLMIN, NLMAX
      IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,*) '============================================================='
         WRITE(LU_SCARC,*) '     SETUP_CCVAR_EXCHANGE ', NL
         WRITE(LU_SCARC,*) '============================================================='
         IF (SL%NX==8) THEN
            WRITE(LU_SCARC,*) 'CCVAR(.,.,.,CGSC)'
            WRITE(LU_SCARC,'(10I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
            WRITE(LU_SCARC,*) 'CCVAR(.,.,.,UNKH)'
            WRITE(LU_SCARC,'(10I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         ELSE
            WRITE(LU_SCARC,*) 'CCVAR(.,.,.,CGSC)'
            WRITE(LU_SCARC,'(6I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
            WRITE(LU_SCARC,*) 'CCVAR(.,.,.,UNKH)'
            WRITE(LU_SCARC,'(6I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         ENDIF
         WRITE(LU_SCARC,*) 'NUNKH_LOC=', SL%NUNKH_LOC
         WRITE(LU_SCARC,*) 'NCS=', SL%NCS
         WRITE(LU_SCARC,*) 'NX =', SL%NX
         WRITE(LU_SCARC,*) 'NY =', SL%NY
         WRITE(LU_SCARC,*) 'NZ =', SL%NZ
       ENDIF
   ENDDO

<<<<<<< HEAD
ENDDO
=======
SCARC_CELL_INDEX = (KK-1)*SL%NX*SL%NY + (JJ-1)*SL%NX + II
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'SCARC_CELL_INDEX(',II,',',JJ,',',KK,')=', SCARC_CELL_INDEX
ENDIF
>>>>>>> origin/scarc

END SUBROUTINE SCARC_SETUP_EXCHANGE


!> ----------------------------------------------------------------------------------------------------
!> Setup neighborship structures and boundary conditions
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DECOMPOSITION
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES
INTEGER :: NM, NL
INTEGER :: IREFINE, IERR, IFACE
INTEGER :: INBR_FACE, INBR_MESH
INTEGER :: IWG, IWL, IWC
INTEGER :: III, JJJ, KKK
INTEGER :: NOM_LAST , NOM
INTEGER :: NCPL_LAST, NCPL
INTEGER :: IOR_LAST, IOR0
INTEGER :: MLATEX
INTEGER :: FACE_NEIGHBORS(NSCARC_MAX_FACE_NEIGHBORS, -3:3)
INTEGER :: MESH_NEIGHBORS(NSCARC_MAX_MESH_NEIGHBORS)
INTEGER :: NUM_FACE_NEIGHBORS(-3:3)
INTEGER :: NUM_MESH_NEIGHBORS
CHARACTER (40) :: CLATEX
LOGICAL :: BFOUND_NOM
TYPE (MESH_TYPE)        , POINTER :: M
TYPE (SCARC_TYPE)       , POINTER :: S
TYPE (OSCARC_TYPE)      , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SLF , SLC
TYPE (SCARC_LEVEL_TYPE) , POINTER :: OSLF, OSLC

IERR=0

CALL SCARC_SETUP_CCVAR(NLEVEL_MIN)

MESHES_LOOP1: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   M => MESHES(NM)
   S => SCARC(NM)

   FACE_NEIGHBORS = -1
   MESH_NEIGHBORS = -1

   NUM_FACE_NEIGHBORS = 0
   NUM_MESH_NEIGHBORS = 0

   !! For all solvers:
   !! Determine array WALL and PRESSURE_BC_INDEX on finest level SLF
   SLF => S%LEVEL(NLEVEL_MIN)

   SLF%NCE = SLF%NCS
   SLF%NCO = SLF%NCS

   SLF%N_EXTERNAL_WALL_CELLS = M%N_EXTERNAL_WALL_CELLS
   SLF%N_INTERNAL_WALL_CELLS = M%N_INTERNAL_WALL_CELLS

   SLF%NLAYER = 1

   IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
      WRITE(LU_SCARC,*) '========== DECOMPOSITION: PROCESSING LEVEL ',NLEVEL_MIN
      WRITE(LU_SCARC,*) 'SLF%NCS:',SLF%NCS
      WRITE(LU_SCARC,*) 'SLF%NCE:',SLF%NCE
      WRITE(LU_SCARC,*) 'SLF%NLAYER:',SLF%NLAYER
      WRITE(LU_SCARC,*) 'SLF%N_EXTERNAL_WALL_CELLS:',SLF%N_EXTERNAL_WALL_CELLS
      WRITE(LU_SCARC,*) 'SLF%N_INTERNAL_WALL_CELLS:',SLF%N_INTERNAL_WALL_CELLS
      WRITE(LU_SCARC,*) 'CELL_INDEX:'
      DO KKK=SLF%NZ+1,0,-1
         WRITE(LU_SCARC,'(6i6)') ((CELL_INDEX(III,JJJ,KKK),III=0,SLF%NX+1),JJJ=0,SLF%NY+1)
      ENDDO
      WRITE(LU_SCARC,*) 'SOLID:'
      DO KKK=SLF%NZ+1,0,-1
         WRITE(LU_SCARC,'(6i6)') ((SOLID(CELL_INDEX(III,JJJ,KKK)),III=0,SLF%NX+1),JJJ=0,SLF%NY+1)
      ENDDO
   ENDIF


   !> -----------------------------------------------------------------------
   !> Determine basic data for single faces (orientation, dimensions)
   !> -----------------------------------------------------------------------
   FACES_OF_MESH_LOOP1: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE FACES_OF_MESH_LOOP1
      SLF%FACE(IOR0)%NFM = NM

      !> information about face orientation and local dimensions
      SELECT CASE (ABS(IOR0))
         CASE (1)
            SLF%FACE(IOR0)%BXY  =  .FALSE.                 !> face in x-y-direction ?
            SLF%FACE(IOR0)%BXZ  =  .FALSE.                 !> face in x-z-direction ?
            SLF%FACE(IOR0)%BYZ  =  .TRUE.                  !> face in y-z-direction ?
            SLF%FACE(IOR0)%NFC  =  SLF%NX                  !> number of cells between opposite mesh faces
            SLF%FACE(IOR0)%NFX  =  1                       !> number of cells in x-direction
            SLF%FACE(IOR0)%NFY  =  SLF%NY                  !> number of cells in y-direction
            SLF%FACE(IOR0)%NFZ  =  SLF%NZ                  !> number of cells in z-direction
            SLF%FACE(IOR0)%NFW  =  SLF%NY*SLF%NZ           !> number of wall cells at that face
            SLF%FACE(IOR0)%DH   => SLF%DXL                 !> step size vector between opposite mesh faces
            SLF%FACE(IOR0)%IOFF_MAT = 1                    !> related offset for matrix diagonals
         CASE (2)
            SLF%FACE(IOR0)%BXY  =  .FALSE.                 !> same as above ...
            SLF%FACE(IOR0)%BXZ  =  .TRUE.
            SLF%FACE(IOR0)%BYZ  =  .FALSE.
            SLF%FACE(IOR0)%NFC  =  SLF%NY
            SLF%FACE(IOR0)%NFX  =  SLF%NX
            SLF%FACE(IOR0)%NFY  =  1
            SLF%FACE(IOR0)%NFZ  =  SLF%NZ
            SLF%FACE(IOR0)%NFW  =  SLF%NX*SLF%NZ
            SLF%FACE(IOR0)%DH   => SLF%DYL
            SLF%FACE(IOR0)%IOFF_MAT = SLF%NX
         CASE (3)
            SLF%FACE(IOR0)%BXY  =  .TRUE.                  !> same as above ...
            SLF%FACE(IOR0)%BXZ  =  .FALSE.
            SLF%FACE(IOR0)%BYZ  =  .FALSE.
            SLF%FACE(IOR0)%NFC  =  SLF%NZ
            SLF%FACE(IOR0)%NFX  =  SLF%NX
            SLF%FACE(IOR0)%NFY  =  SLF%NY
            SLF%FACE(IOR0)%NFZ  =  1
            SLF%FACE(IOR0)%NFW  =  SLF%NX*SLF%NY
            SLF%FACE(IOR0)%DH   => SLF%DZL
            SLF%FACE(IOR0)%IOFF_MAT = SLF%NX * SLF%NY
      END SELECT
      IF (IOR0 > 0) SLF%FACE(IOR0)%IOFF_MAT = -SLF%FACE(IOR0)%IOFF_MAT      !> take negative offset for IOR>0

   ENDDO FACES_OF_MESH_LOOP1

   !! store local IWG-number for each face
   IWG = 1
   FACE_ORDER_LOOP: DO IFACE = 1, 6
      IOR0 = FACE_ORDER_XYZ(IFACE)
      SLF%FACE(IOR0)%IWG_MARKER = IWG
      IWG = IWG + SLF%FACE(IOR0)%NFW
   ENDDO FACE_ORDER_LOOP


   !> -----------------------------------------------------------------------
   !> First loop over global IW's:
   !> store basic data and determine number of adajacent neighbors to each
   !> face with corresponding number of IW's
   !> -----------------------------------------------------------------------
   IOR_LAST =  0
   NOM_LAST  = -1
   NCPL_LAST = -1
   IWL = 0
   BFOUND_NOM = .FALSE.

   IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
      WRITE(LU_SCARC,*) 'M%N_EXTERNAL_WALL_CELLS=',M%N_EXTERNAL_WALL_CELLS
      WRITE(LU_SCARC,*) 'M%N_INTERNAL_WALL_CELLS=',M%N_INTERNAL_WALL_CELLS
      DO IWG = 1, SLF%N_EXTERNAL_WALL_CELLS
         WRITE(LU_SCARC,*) 'M%WALL(',IWG,')%ONE_D%IOR=',M%WALL(IWG)%ONE_D%IOR
      ENDDO
      DO IWG = 1, SLF%N_EXTERNAL_WALL_CELLS
         WRITE(LU_SCARC,*) 'M%EXTERNAL_WALL(',IWG,')%NOM=',M%EXTERNAL_WALL(IWG)%NOM
      ENDDO
      DO IWG = SLF%N_EXTERNAL_WALL_CELLS+1, SLF%N_EXTERNAL_WALL_CELLS+SLF%N_INTERNAL_WALL_CELLS
         WRITE(LU_SCARC,*) 'M%INTERNAL_WALL(',IWG,')=',M%WALL(IWG)%ONE_D%IOR
      ENDDO
   ENDIF

   !> First process external wall cells
<<<<<<< HEAD
   EXTERNAL_WALL_CELLS_LOOP1: DO IWG = 1, SLF%N_EXTERNAL_WALL_CELLS
=======
   EXTERNAL_WALLCELLS_LOOP1: DO IWG = 1, SLF%N_EXTERNAL_WALL_CELLS
>>>>>>> origin/scarc

      !> Determine and store neighbors, orientation and number of couplings for a single wall cell
      NOM  =  M%EXTERNAL_WALL(IWG)%NOM
      IOR0 =  M%WALL(IWG)%ONE_D%IOR
      NCPL = (M%EXTERNAL_WALL(IWG)%IIO_MAX - M%EXTERNAL_WALL(IWG)%IIO_MIN + 1) * &
             (M%EXTERNAL_WALL(IWG)%JJO_MAX - M%EXTERNAL_WALL(IWG)%JJO_MIN + 1) * &
             (M%EXTERNAL_WALL(IWG)%KKO_MAX - M%EXTERNAL_WALL(IWG)%KKO_MIN + 1)

      SLF%WALL(IWG)%NOM  = NOM                              !> store data in SCARC-structure
      SLF%WALL(IWG)%IOR  = IOR0
      SLF%WALL(IWG)%NCPL = NCPL

      IWL = IWL + 1                                         !> count local wall cells for that face

      !> for wall cells with a neighbor build extended and overlapping structures
      IF (NOM /= 0) THEN

         SLF%NCE = SLF%NCE + NCPL                           !> increase number of extended grid cells
         SLF%NCO = SLF%NCO + 1                              !> increase number of overlapping grid cells

         OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)
         OSLF%NCPL = NCPL                                   !> initialize own counter for local wall cells
         OSLF%IOR  = -IOR0                                  !> initialize own orientation variable

         !> check if this neighbor already was accounted for
         DO INBR_FACE = 1, NUM_FACE_NEIGHBORS(IOR0)
            IF (FACE_NEIGHBORS(INBR_FACE, IOR0) == NOM) BFOUND_NOM=.TRUE.
         ENDDO

         IF (.NOT.BFOUND_NOM) THEN
            OSLF%NWL = 1                                    !> initialize own counter for local wall cells
            OSLF%NCG = NCPL                                 !> initialize counter for local ghost cells
            SLF%NCPL_MAX  = MAX(SLF%NCPL_MAX, NCPL)         !> get max NCPL ever used on this mesh

            IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
               WRITE(LU_SCARC,'(11(a,i6))') & 'NM=',MYID+1,&
                     ' A: IWG=',IWG,': NOM=',NOM,': OSLF%NWL=', OSLF%NWL, ': SLF%NCO=', SLF%NCO, &
                     ': SLF%NCE=',SLF%NCE,': IWL=',IWL, ': IOR0=',IOR0,': NCPL=',OSLF%NCPL
            ENDIF

         ELSE                                               !> already known neighbor ???
            OSLF%NWL = OSLF%NWL + 1                         !> increase own counter for local wall cells
            OSLF%NCG = OSLF%NCG + NCPL                      !> increase counter for local ghost cells

            IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
               WRITE(LU_SCARC,'(11(a,i6))') & 'NM=',MYID+1,&
                     ' B: IWG=',IWG,': NOM=',NOM,': OSLF%NWL=', OSLF%NWL, ': SLF%NCO=', SLF%NCO, &
                     ': SLF%NCE=',SLF%NCE,': IWL=',IWL, ': IOR0=',IOR0,': NCPL=',OSLF%NCPL
            ENDIF

         ENDIF
      ENDIF

      IF (NOM /= 0 .AND. NOM_LAST /= NOM) THEN                            !> new neighbor at that face?

         NUM_FACE_NEIGHBORS(IOR0) = NUM_FACE_NEIGHBORS(IOR0) + 1          !> increase neighbor counter for face
         NUM_MESH_NEIGHBORS       = NUM_MESH_NEIGHBORS       + 1          !> increase neighbor counter for mesh

         FACE_NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0), IOR0) = NOM             !> store number of neighbor for face
         MESH_NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0))       = NOM             !> store number of neighbor for mesh

      ENDIF

      IOR_LAST  = IOR0                                                    !> save former values
      NCPL_LAST = NCPL
      NOM_LAST  = NOM

<<<<<<< HEAD
   ENDDO EXTERNAL_WALL_CELLS_LOOP1
 

   !> Then process internal wall cells
   INTERNAL_WALL_CELLS_LOOP1: DO IWG = SLF%N_EXTERNAL_WALL_CELLS+1, SLF%N_EXTERNAL_WALL_CELLS+SLF%N_INTERNAL_WALL_CELLS

      SLF%WALL(IWG)%IOR  = M%WALL(IWG)%ONE_D%IOR
      SLF%WALL(IWG)%NOM  = 0

      SLF%WALL(IWG)%BTYPE  = NEUMANN
      SLF%WALL(IWG)%BOUNDARY_TYPE  = M%WALL(IWG)%BOUNDARY_TYPE

      SLF%WALL(IWG)%IXG =  M%WALL(IWG)%ONE_D%II                                 !> ghost cell indices
      SLF%WALL(IWG)%IYG =  M%WALL(IWG)%ONE_D%JJ
      SLF%WALL(IWG)%IZG =  M%WALL(IWG)%ONE_D%KK

      SLF%WALL(IWG)%IXW =  M%WALL(IWG)%ONE_D%IIG                                !> (internal) wall cell indices
      SLF%WALL(IWG)%IYW =  M%WALL(IWG)%ONE_D%JJG
      SLF%WALL(IWG)%IZW =  M%WALL(IWG)%ONE_D%KKG

   ENDDO INTERNAL_WALL_CELLS_LOOP1

=======
   ENDDO EXTERNAL_WALLCELLS_LOOP1
 

   !> Then process internal wall cells
   INTERNAL_WALLCELLS_LOOP1: DO IWG = SLF%N_EXTERNAL_WALL_CELLS+1, SLF%N_EXTERNAL_WALL_CELLS+SLF%N_INTERNAL_WALL_CELLS

      SLF%WALL(IWG)%IOR  = M%WALL(IWG)%ONE_D%IOR
      SLF%WALL(IWG)%NOM  = 0

      SLF%WALL(IWG)%BTYPE  = NEUMANN

      SLF%WALL(IWG)%IXG =  M%WALL(IWG)%ONE_D%II                                 !> ghost cell indices
      SLF%WALL(IWG)%IYG =  M%WALL(IWG)%ONE_D%JJ
      SLF%WALL(IWG)%IZG =  M%WALL(IWG)%ONE_D%KK

      SLF%WALL(IWG)%IXW =  M%WALL(IWG)%ONE_D%IIG                                !> (internal) wall cell indices
      SLF%WALL(IWG)%IYW =  M%WALL(IWG)%ONE_D%JJG
      SLF%WALL(IWG)%IZW =  M%WALL(IWG)%ONE_D%KKG

   ENDDO INTERNAL_WALLCELLS_LOOP1

>>>>>>> origin/scarc

   !> allocate array which stores numbers of all neighboring meshes
   IF (NUM_MESH_NEIGHBORS /= 0) THEN
      ALLOCATE(S%NEIGHBORS(NUM_MESH_NEIGHBORS), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION','MESH_NEIGHBORS',IERR)
   ENDIF
   S%NUM_NEIGHBORS = NUM_MESH_NEIGHBORS


   !> -----------------------------------------------------------------------
   !> Store information about adjacent neighbors on different faces
   !> Allocate corresponding index arrays in OSCARC-structures
   !> -----------------------------------------------------------------------
   !! allocate administrative mapping arrays for own mesh
   CALL SCARC_SETUP_MAPPINGS(NM, NLEVEL_MIN)

   INBR_MESH = 1
   NEIGHBORS_OF_FACE_LOOP: DO IOR0 = -3, 3
      IF (IOR0 == 0) CYCLE NEIGHBORS_OF_FACE_LOOP

      !> if there are neighbors at face IOR0 store information about them
      IF (NUM_FACE_NEIGHBORS(IOR0) /= 0) THEN

         SLF%FACE(IOR0)%NUM_NEIGHBORS = NUM_FACE_NEIGHBORS(IOR0)        !> store number of neighbors on face

         IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) &
            WRITE(LU_SCARC,'(a,i4,a,i4,a,i4)') 'IOR0=',IOR0,': NUM_FACE_NEIGHBORS=',NUM_FACE_NEIGHBORS(IOR0)

         !> allocate array for storing the numbers of the single neighbors
         ALLOCATE(SLF%FACE(IOR0)%NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0)), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION','NOMS',IERR)

         !> store every neighbor and allocate corresponding administration arrays
         DO INBR_FACE = 1, NUM_FACE_NEIGHBORS(IOR0)

            NOM = FACE_NEIGHBORS(INBR_FACE,IOR0)

            SLF%FACE(IOR0)%NEIGHBORS(INBR_FACE) = NOM                   !> store NOM as a neighbor of that face
            S%NEIGHBORS(INBR_MESH)              = NOM                   !> store NOM as a neighbor of that mesh
            INBR_MESH = INBR_MESH + 1

            !! allocate administrative arrays for neighboring meshes
            CALL SCARC_SETUP_OMAPPINGS(NM, NOM, NLEVEL_MIN)

         ENDDO
      ENDIF
   ENDDO NEIGHBORS_OF_FACE_LOOP


   !> --------------------------------------------------------------------------
   !> Second loop over external wall cells:
   !> store detailed coordinate and cell data and get type of boundary condition
   !> --------------------------------------------------------------------------
   SLF%ICE_MARKER = SLF%NCS
   SLF%ICO_MARKER = SLF%NCS

   IOR_LAST = 0
<<<<<<< HEAD
   WALL_CELLS_LOOP2: DO IWG = 1, SLF%N_EXTERNAL_WALL_CELLS
=======
   WALLCELLS_LOOP2: DO IWG = 1, SLF%N_EXTERNAL_WALL_CELLS
>>>>>>> origin/scarc

      !> Determine neighbors and orientation again
      NOM  = SLF%WALL(IWG)%NOM
      IOR0 = SLF%WALL(IWG)%IOR
      NCPL = SLF%WALL(IWG)%NCPL

      !> Determine boundary type for IW
      IF (M%EXTERNAL_WALL(IWG)%NOM /= 0) THEN
         SLF%WALL(IWG)%BTYPE = INTERNAL
      ELSE IF (M%WALL(IWG)%PRESSURE_BC_INDEX == DIRICHLET) THEN
         SLF%WALL(IWG)%BTYPE = DIRICHLET
      ELSE
         SLF%WALL(IWG)%BTYPE = NEUMANN
      ENDIF

<<<<<<< HEAD
      SLF%WALL(IWG)%BOUNDARY_TYPE  = M%WALL(IWG)%BOUNDARY_TYPE

=======
>>>>>>> origin/scarc
      SLF%WALL(IWG)%IXG =  M%WALL(IWG)%ONE_D%II                                 !> ghost cell indices
      SLF%WALL(IWG)%IYG =  M%WALL(IWG)%ONE_D%JJ
      SLF%WALL(IWG)%IZG =  M%WALL(IWG)%ONE_D%KK

      SLF%WALL(IWG)%IXW =  M%WALL(IWG)%ONE_D%IIG                                !> (internal) wall cell indices
      SLF%WALL(IWG)%IYW =  M%WALL(IWG)%ONE_D%JJG
      SLF%WALL(IWG)%IZW =  M%WALL(IWG)%ONE_D%KKG


      !> If there exists a neighbor for that wall cell, setup corresponding neighborship information
      IF (NOM /= 0) THEN
         CALL SCARC_SETUP_WALLCELL_NEIGHBOR(M%EXTERNAL_WALL(IWG)%IIO_MIN, &
                                            M%EXTERNAL_WALL(IWG)%IIO_MAX, &
                                            M%EXTERNAL_WALL(IWG)%JJO_MIN, &
                                            M%EXTERNAL_WALL(IWG)%JJO_MAX, &
                                            M%EXTERNAL_WALL(IWG)%KKO_MIN, &
                                            M%EXTERNAL_WALL(IWG)%KKO_MAX, &
                                            IWG, IOR0, NM, NOM, NLEVEL_MIN)
      ENDIF

      NOM_LAST  = NOM
      IOR_LAST = IOR0

   ENDDO WALL_CELLS_LOOP2

ENDDO MESHES_LOOP1
<<<<<<< HEAD
=======

!> -----------------------------------------------------------------------------------------
!> Set CCVAR information 
!> -----------------------------------------------------------------------------------------
!> on finest level
CALL SCARC_SETUP_CCVAR_GLOBAL(NLEVEL_MIN)

!> if requested also on all coarser levels
IF (TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN
   DO NL = NLEVEL_MIN+1, NLEVEL_MAX
      CALL SCARC_SETUP_CCVAR_LEVEL(NL)
      CALL SCARC_SETUP_CCVAR_GLOBAL(NL)
   ENDDO
ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   DO NM = 1, NMESHES
      IF (PROCESS(NM) /= MYID) CYCLE
      M   => MESHES(NM)
      SLF => SCARC(NM)%LEVEL(NLEVEL_MIN)
      WRITE(LU_SCARC,*) 'M%N_OBST=',M%N_OBST
      !WRITE(LU_SCARC,*) 'M%OBST_INDEX_C(IC)'
      !WRITE(LU_SCARC,*) (M%OBST_INDEX_C(IWG),IWG=1,SLF%NW)
      WRITE(LU_SCARC,*) 'M%OBST ... I1, I2, J1, J2, K1, K2'
      DO IWG = 1, M%N_OBST
         WRITE(LU_SCARC,'(6i6)') M%OBSTRUCTION(IWG)%I1,M%OBSTRUCTION(IWG)%I2,&
                                 M%OBSTRUCTION(IWG)%J1,M%OBSTRUCTION(IWG)%J2,&
                                 M%OBSTRUCTION(IWG)%K1,M%OBSTRUCTION(IWG)%K2
      ENDDO
      WRITE(LU_SCARC,*) 'M%N_OBST=',M%N_OBST
      WRITE(LU_SCARC,*) 'M%N_WALL_CELLS=',M%N_WALL_CELLS
      WRITE(LU_SCARC,*) 'M%CELL_INDEX(.,.,.):'
      WRITE(LU_SCARC,'(6i4)') (((M%CELL_INDEX(III,JJJ,KKK), III=0,5,1), JJJ=1,1,-1),KKK=5,0,-1)
      WRITE(LU_SCARC,*) 'M%SOLID(.):'
      WRITE(LU_SCARC,'(8i4)') (M%SOLID(IWG), IWG=1, SLF%NCS)
      WRITE(LU_SCARC,*) 'M%WALL(.)%BOUNDARY_TYPE:'
      WRITE(LU_SCARC,'(8i4)') (M%WALL(IWG)%BOUNDARY_TYPE, IWG=1, SLF%NW)
      WRITE(LU_SCARC,*) 'M%WALL(.)% II, JJ, KK:'
      WRITE(LU_SCARC,'(i4,a,3i8)') ( IWG, ' : ', &
                               M%WALL(IWG)%ONE_D%II,  &
                               M%WALL(IWG)%ONE_D%JJ,  &
                               M%WALL(IWG)%ONE_D%KK,  &
                               IWG=1, SLF%NW)
      WRITE(LU_SCARC,*) 'M%WALL(.)% IIG, JJG, KKG:'
      WRITE(LU_SCARC,'(i4,a,3i8)') ( IWG, ' : ', &
                               M%WALL(IWG)%ONE_D%IIG,  &
                               M%WALL(IWG)%ONE_D%JJG,  &
                               M%WALL(IWG)%ONE_D%KKG,  &
                               IWG=1, SLF%NW)
      WRITE(LU_SCARC,*) 'M%WALL_INDEX:'
      WRITE(LU_SCARC,'(i4,a,6i8)') ( IWG, ' : ', &
                               M%WALL_INDEX(IWG, 1),  &
                               M%WALL_INDEX(IWG,-1),  &
                               M%WALL_INDEX(IWG, 2),  &
                               M%WALL_INDEX(IWG,-2),  &
                               M%WALL_INDEX(IWG, 3),  &
                               M%WALL_INDEX(IWG,-3),  &
                               IWG=1, SLF%NCS)
      WRITE(LU_SCARC,*) 'M%WALL(.)%ONE_D% IOR,II,JJ,KK:'
      DO IWG = 1, SLF%NW
         WRITE(LU_SCARC,'(5I6)') &
            IWG,M%WALL(IWG)%ONE_D%IOR,M%WALL(IWG)%ONE_D%II,M%WALL(IWG)%ONE_D%JJ,M%WALL(IWG)%ONE_D%KK
      ENDDO
      WRITE(LU_SCARC,*) 'M%WALL(.)%ONE_D% IOR,IIG,JJG,KKG:'
      DO IWG = 1, SLF%NW
         WRITE(LU_SCARC,'(5I6)') &
            IWG,M%WALL(IWG)%ONE_D%IOR,M%WALL(IWG)%ONE_D%IIG,M%WALL(IWG)%ONE_D%JJG,M%WALL(IWG)%ONE_D%KKG
      ENDDO
      WRITE(LU_SCARC,*) 'CCVAR(..., CGSC):'
      DO KKK=SLF%NZ+1,0,-1
         WRITE(LU_SCARC,'(6i6)') ((SLF%CCVAR(III,JJJ,KKK, CGSC),III=0,SLF%NX+1),JJJ=0,SLF%NY+1)
      ENDDO
      WRITE(LU_SCARC,*) 'CCVAR(..., UNKH):'
      DO KKK=SLF%NZ+1,0,-1
         WRITE(LU_SCARC,'(6i6)') ((SLF%CCVAR(III,JJJ,KKK, UNKH),III=0,SLF%NX+1),JJJ=0,SLF%NY+1)
      ENDDO
      WRITE(LU_SCARC,*) 'M%CELL_INDEX:'
      DO KKK=SLF%NZ+1,0,-1
         WRITE(LU_SCARC,'(6i6)') ((CELL_INDEX(III,JJJ,KKK),III=0,SLF%NX+1),JJJ=0,SLF%NY+1)
      ENDDO
   ENDDO
ENDIF

!> -----------------------------------------------------------------------------------------
!! Only in case of MG-method:
!! Determine WALL, FACE and OSCARC types for coarser levels
!> -----------------------------------------------------------------------------------------
LEVEL_GMG_IF: IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN

   MESHES_LOOP2: DO NM = 1, NMESHES

      IF (PROCESS(NM) /= MYID) CYCLE
      S => SCARC(NM)
>>>>>>> origin/scarc

!> -----------------------------------------------------------------------------------------
!> Set CCVAR information 
!> -----------------------------------------------------------------------------------------
!> on finest level
CALL SCARC_SETUP_CCVAR_GLOBAL(NLEVEL_MIN)

!> if requested also on all coarser levels
IF (TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN
   DO NL = NLEVEL_MIN+1, NLEVEL_MAX
      CALL SCARC_SETUP_CCVAR_LEVEL(NL)
      CALL SCARC_SETUP_CCVAR_GLOBAL(NL)
   ENDDO
ENDIF

<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   DO NM = 1, NMESHES
      IF (PROCESS(NM) /= MYID) CYCLE
      M   => MESHES(NM)
      SLF => SCARC(NM)%LEVEL(NLEVEL_MIN)
      WRITE(LU_SCARC,*) 'M%N_OBST=',M%N_OBST
      WRITE(LU_SCARC,*) 'M%OBST ... I1, I2, J1, J2, K1, K2'
      DO IWG = 1, M%N_OBST
         WRITE(LU_SCARC,'(6i6)') M%OBSTRUCTION(IWG)%I1,M%OBSTRUCTION(IWG)%I2,&
                                 M%OBSTRUCTION(IWG)%J1,M%OBSTRUCTION(IWG)%J2,&
                                 M%OBSTRUCTION(IWG)%K1,M%OBSTRUCTION(IWG)%K2
      ENDDO
      WRITE(LU_SCARC,*) 'M%N_OBST=',M%N_OBST
      WRITE(LU_SCARC,*) 'M%N_WALL_CELLS=',M%N_WALL_CELLS
      DO JJJ = SLF%NY+1,0,-1
         WRITE(LU_SCARC,*) 'M%CELL_INDEX(.,',JJJ,',.):'
         WRITE(LU_SCARC,'(18i6)') ((M%CELL_INDEX(III,JJJ,KKK), III=0,SLF%NX+1,1), KKK=SLF%NZ+1,0,-1)
      ENDDO
      WRITE(LU_SCARC,*) 'M%WALL(.)%BOUNDARY_TYPE:'
      WRITE(LU_SCARC,'(16i4)') (M%WALL(IWG)%BOUNDARY_TYPE, IWG=1, SLF%NW)
      !WRITE(LU_SCARC,*) 'M%WALL(.)% II, JJ, KK:'
      !WRITE(LU_SCARC,'(i4,a,3i8)') ( IWG, ' : ', &
      !                         M%WALL(IWG)%ONE_D%II,  &
      !                         M%WALL(IWG)%ONE_D%JJ,  &
      !                         M%WALL(IWG)%ONE_D%KK,  &
      !                         IWG=1, SLF%NW)
      !WRITE(LU_SCARC,*) 'M%WALL(.)% IIG, JJG, KKG:'
      !WRITE(LU_SCARC,'(i4,a,3i8)') ( IWG, ' : ', &
      !                         M%WALL(IWG)%ONE_D%IIG,  &
      !                         M%WALL(IWG)%ONE_D%JJG,  &
      !                         M%WALL(IWG)%ONE_D%KKG,  &
      !                         IWG=1, SLF%NW)
      WRITE(LU_SCARC,*) 'M%WALL_INDEX:'
      WRITE(LU_SCARC,'(i4,a,6i8)') ( IWG, ' : ', &
                               M%WALL_INDEX(IWG, 1),  &
                               M%WALL_INDEX(IWG,-1),  &
                               M%WALL_INDEX(IWG, 2),  &
                               M%WALL_INDEX(IWG,-2),  &
                               M%WALL_INDEX(IWG, 3),  &
                               M%WALL_INDEX(IWG,-3),  &
                               IWG=1, CELL_COUNT(NM))
      WRITE(LU_SCARC,*) 'M%WALL(.)%ONE_D% IOR,II,JJ,KK:'
      DO IWG = 1, SLF%NW
         WRITE(LU_SCARC,'(5I6)') &
            IWG,M%WALL(IWG)%ONE_D%IOR,M%WALL(IWG)%ONE_D%II,M%WALL(IWG)%ONE_D%JJ,M%WALL(IWG)%ONE_D%KK
      ENDDO
      WRITE(LU_SCARC,*) 'M%WALL(.)%ONE_D% IOR,IIG,JJG,KKG:'
      DO IWG = 1, SLF%NW
         WRITE(LU_SCARC,'(5I6)') &
            IWG,M%WALL(IWG)%ONE_D%IOR,M%WALL(IWG)%ONE_D%IIG,M%WALL(IWG)%ONE_D%JJG,M%WALL(IWG)%ONE_D%KKG
      ENDDO
      WRITE(LU_SCARC,*) 'CCVAR(..., CGSC):'
      DO JJJ=SLF%NY+1,0,-1
         IF (SLF%NX==8) THEN
            WRITE(LU_SCARC,'(10i6)') ((SLF%CCVAR(III,JJJ,KKK, CGSC),III=0,SLF%NX+1),KKK=SLF%NZ+1,0,-1)
         ELSE
            WRITE(LU_SCARC,'(6i6)') ((SLF%CCVAR(III,JJJ,KKK, CGSC),III=0,SLF%NX+1),KKK=SLF%NZ+1,0,-1)
         ENDIF
         WRITE(LU_SCARC,*)
      ENDDO
      WRITE(LU_SCARC,*) 'CCVAR(..., UNKH):'
      DO JJJ=SLF%NY+1,0,-1
         IF (SLF%NX==8) THEN
            WRITE(LU_SCARC,'(10i6)') ((SLF%CCVAR(III,JJJ,KKK, UNKH),III=0,SLF%NX+1),KKK=SLF%NZ+1,0,-1)
         ELSE
            WRITE(LU_SCARC,'(6i6)') ((SLF%CCVAR(III,JJJ,KKK, UNKH),III=0,SLF%NX+1),KKK=SLF%NZ+1,0,-1)
         ENDIF
         WRITE(LU_SCARC,*)
      ENDDO
      WRITE(LU_SCARC,*) 'M%CELL_INDEX:'
      DO JJJ=SLF%NY+1,0,-1
         IF (SLF%NX==8) THEN
            WRITE(LU_SCARC,'(10i6)') ((CELL_INDEX(III,JJJ,KKK),III=0,SLF%NX+1),KKK=SLF%NZ+1,0,-1)
         ELSE
            WRITE(LU_SCARC,'(6i6)') ((CELL_INDEX(III,JJJ,KKK),III=0,SLF%NX+1),KKK=SLF%NZ+1,0,-1)
         ENDIF
         WRITE(LU_SCARC,*)
      ENDDO
   ENDDO
ENDIF
=======
         !> point to SCARC grid structures on coarser and finer level
         SLF => S%LEVEL(NL-1)
         SLC => S%LEVEL(NL)
>>>>>>> origin/scarc

!> -----------------------------------------------------------------------------------------
!! Only in case of MG-method:
!! Determine WALL, FACE and OSCARC types for coarser levels
!> -----------------------------------------------------------------------------------------
LEVEL_GMG_IF: IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN

<<<<<<< HEAD
   MESHES_LOOP2: DO NM = 1, NMESHES

      IF (PROCESS(NM) /= MYID) CYCLE
      S => SCARC(NM)

      IREFINE=1
      LEVEL_GMG_LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX

         !> point to SCARC grid structures on coarser and finer level
         SLF => S%LEVEL(NL-1)
         SLC => S%LEVEL(NL)

         IREFINE=IREFINE*2

         CALL SCARC_CHECK_DIVISIBILITY(SLF%NCE-SLF%NCS, 'SLC%NCE', 'SCARC_SETUP_DECOMPOSITION')
         CALL SCARC_CHECK_DIVISIBILITY(SLF%NCO-SLF%NCS, 'SLC%NCO', 'SCARC_SETUP_DECOMPOSITION')

         SLC%NCE = SLC%NCS + (SLF%NCE-SLF%NCS)/2
         SLC%NCO = SLC%NCS + (SLF%NCO-SLF%NCS)/2

         SLC%ICE_MARKER = SLC%NCS
         SLC%ICO_MARKER = SLC%NCS

         SLC%N_EXTERNAL_WALL_CELLS = SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
         SLC%N_INTERNAL_WALL_CELLS = SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)
=======
         CALL SCARC_CHECK_DIVISIBILITY(SLF%NCE-SLF%NCS, 'SLC%NCE', 'SCARC_SETUP_DECOMPOSITION')
         CALL SCARC_CHECK_DIVISIBILITY(SLF%NCO-SLF%NCS, 'SLC%NCO', 'SCARC_SETUP_DECOMPOSITION')

         SLC%NCE = SLC%NCS + (SLF%NCE-SLF%NCS)/2
         SLC%NCO = SLC%NCS + (SLF%NCO-SLF%NCS)/2

         SLC%ICE_MARKER = SLC%NCS
         SLC%ICO_MARKER = SLC%NCS

         SLC%N_EXTERNAL_WALL_CELLS = 2*SLC%NX*SLC%NY + 2*SLC%NX*SLC%NZ + 2*SLC%NY*SLC%NZ
         SLC%N_INTERNAL_WALL_CELLS = SCARC_INTERNAL_WALL_CELLS(NM, NL)

         SLC%NW = SLC%N_EXTERNAL_WALL_CELLS + SLC%N_INTERNAL_WALL_CELLS

   IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) '========== DECOMPOSITION: PROCESSING LEVEL ',NL
      WRITE(LU_SCARC,*) 'SLC%NX :',SLC%NX 
      WRITE(LU_SCARC,*) 'SLC%NY :',SLC%NY
      WRITE(LU_SCARC,*) 'SLC%NZ :',SLC%NZ
      WRITE(LU_SCARC,*) 'SLC%NCS:',SLC%NCS
      WRITE(LU_SCARC,*) 'SLC%NCE:',SLC%NCE
      WRITE(LU_SCARC,*) 'SLC%NCO:',SLC%NCO
      WRITE(LU_SCARC,*) 'SLC%N_EXTERNAL_WALL_CELLS=',SLC%N_EXTERNAL_WALL_CELLS
      WRITE(LU_SCARC,*) 'SLC%N_INTERNAL_WALL_CELLS=',SLC%N_INTERNAL_WALL_CELLS
      WRITE(LU_SCARC,*) 'SLC%NW=',SLC%NW
      WRITE(LU_SCARC,*) 'SLF%N_OBST:',SLF%N_OBST
      WRITE(LU_SCARC,*) 'OBSTS:'
      DO KKK=1,SLF%N_OBST
         WRITE(LU_SCARC,'(6i6)') SLF%OBST(KKK)%I1,SLF%OBST(KKK)%I2, &
                                 SLF%OBST(KKK)%J1,SLF%OBST(KKK)%J2, &
                                 SLF%OBST(KKK)%K1,SLF%OBST(KKK)%K2
      ENDDO
      WRITE(LU_SCARC,*) 'CELL_INDEX:'
      DO KKK=SLF%NZ+1,0,-1
         WRITE(LU_SCARC,'(6i6)') ((CELL_INDEX(III,JJJ,KKK),III=0,SLF%NX+1),JJJ=0,SLF%NY+1)
      ENDDO
      WRITE(LU_SCARC,*) 'SOLID:'
      DO KKK=SLF%NZ+1,0,-1
         WRITE(LU_SCARC,'(6i6)') ((SOLID(CELL_INDEX(III,JJJ,KKK)),III=0,SLF%NX+1),JJJ=0,SLF%NY+1)
      ENDDO
ENDIF

         !! allocate WALL array for coarse level
         ALLOCATE(SLC%WALL(SLC%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_INTERFACES','WALL',IERR)
>>>>>>> origin/scarc

         SLC%NW = SLC%N_EXTERNAL_WALL_CELLS + SLC%N_INTERNAL_WALL_CELLS

         ALLOCATE(SLC%WALL(SLC%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_INTERFACES','WALL',IERR)

   IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
      WRITE(LU_SCARC,*) '========== DECOMPOSITION: PROCESSING LEVEL ',NL
      WRITE(LU_SCARC,*) 'SLC%NX :',SLC%NX 
      WRITE(LU_SCARC,*) 'SLC%NY :',SLC%NY
      WRITE(LU_SCARC,*) 'SLC%NZ :',SLC%NZ
      WRITE(LU_SCARC,*) 'SLC%NCS:',SLC%NCS
      WRITE(LU_SCARC,*) 'SLC%NCE:',SLC%NCE
      WRITE(LU_SCARC,*) 'SLC%NCO:',SLC%NCO
      WRITE(LU_SCARC,*) 'SLC%N_EXTERNAL_WALL_CELLS=',SLC%N_EXTERNAL_WALL_CELLS
      WRITE(LU_SCARC,*) 'SLC%N_INTERNAL_WALL_CELLS=',SLC%N_INTERNAL_WALL_CELLS
      WRITE(LU_SCARC,*) 'SLC%NW=',SLC%NW
      WRITE(LU_SCARC,*) 'SLF%N_OBST:',SLF%N_OBST
      WRITE(LU_SCARC,*) 'OBSTS:'
      DO KKK=1,SLF%N_OBST
         WRITE(LU_SCARC,'(6i6)') SLF%OBST(KKK)%I1,SLF%OBST(KKK)%I2, &
                                 SLF%OBST(KKK)%J1,SLF%OBST(KKK)%J2, &
                                 SLF%OBST(KKK)%K1,SLF%OBST(KKK)%K2
      ENDDO
      WRITE(LU_SCARC,*) 'CELL_INDEX:'
      DO KKK=SLF%NZ+1,0,-1
         WRITE(LU_SCARC,'(6i6)') ((CELL_INDEX(III,JJJ,KKK),III=0,SLF%NX+1),JJJ=0,SLF%NY+1)
      ENDDO
      WRITE(LU_SCARC,*) 'SOLID:'
      DO KKK=SLF%NZ+1,0,-1
         WRITE(LU_SCARC,'(6i6)') ((SOLID(CELL_INDEX(III,JJJ,KKK)),III=0,SLF%NX+1),JJJ=0,SLF%NY+1)
      ENDDO
ENDIF

IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'ALLOCATING SLC%WALL IN LENGTH ', SLC%NW

         !> ACHTUNG: HIER ÜBERARBEITEN !!!
         CALL SCARC_SETUP_CELL_INDEX(NM, NL)
         CALL SCARC_SETUP_INTERNAL_WALL_COORDS(NM, NL)
         CALL SCARC_SETUP_WALL_INDEX(NM, NL)

         !> ACHTUNG: HIER ÜBERARBEITEN !!!
         !SLC%WALL(SLC%N_EXTERNAL_WALL_CELLS+1:SLC%N_EXTERNAL_WALL_CELLS+SLC%N_INTERNAL_WALL_CELLS)%BOUNDARY_TYPE=SOLID_BOUNDARY
         !SLC%WALL(SLC%N_EXTERNAL_WALL_CELLS+1:SLC%N_EXTERNAL_WALL_CELLS+SLC%N_INTERNAL_WALL_CELLS)%BTYPE        =NEUMANN

         !! Allocate administrative mapping arrays for own level
         CALL SCARC_SETUP_MAPPINGS(NM, NL)

         !> compute FACE and WALL information for all faces of coarser level
         IWC = 1
         IWG = 1
         DO IFACE = 1, 6

            IOR0 = FACE_ORDER_XYZ(IFACE)

            !> compute mesh dimensions of coarser mesh level
            CALL SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)

            !> for every neighbor do:
            IF (SLF%FACE(IOR0)%NUM_NEIGHBORS /= 0) THEN
               DO INBR_FACE = 1, SLF%FACE(IOR0)%NUM_NEIGHBORS

                  NOM = SLF%FACE(IOR0)%NEIGHBORS(INBR_FACE)

                  OSLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL-1)

                  OSLC%NCPL = OSLF%NCPL

                  CALL SCARC_CHECK_DIVISIBILITY(OSLF%NWL, 'OSLF%NWL', 'SCARC_SETUP_DECOMPOSITION')
                  CALL SCARC_CHECK_DIVISIBILITY(OSLF%NCG, 'OSLF%NCG', 'SCARC_SETUP_DECOMPOSITION')

                  IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) THEN
                     OSLC%NWL = OSLF%NWL/4
                     OSLC%NCG = OSLF%NCG/4
                  ELSE
                     OSLC%NWL = OSLF%NWL/2
                     OSLC%NCG = OSLF%NCG/2
                  ENDIF

                  IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
                     WRITE(LU_SCARC,'(5(a,i4))') 'HIII_FINE   : NOM=',NOM,': NCPL=',OSLC%NCPL,&
                                                 ': NWL=',OSLF%NWL,': NL=',NL,': NCG=',OSLF%NCG
                     WRITE(LU_SCARC,'(5(a,i4))') 'HIII_COARSE: NOM=',NOM,': NCPL=',OSLC%NCPL,&
                                                 ': NWL=',OSLC%NWL,': NL=',NL,': NCG=',OSLC%NCG
                  ENDIF

                  !> compute exchange dimensions of coarser mesh level (analogues to NIC, I_MIN ... )
                  CALL SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)

                  !> allocate administrative mapping arrays
                  CALL SCARC_SETUP_OMAPPINGS(NM, NOM, NL)

               ENDDO
            ENDIF

            !> setup complete face information for coarser mesh
            CALL SCARC_SETUP_FACE(IOR0, IWC, IREFINE, NM, NL)

         ENDDO

<<<<<<< HEAD
         IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
           WRITE(LU_SCARC,*) 'SLC%N_INTERNAL_WALL_CELLS=', SLC%N_INTERNAL_WALL_CELLS
=======
         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
           WRITE(LU_SCARC,*) 'SLC%N_INTERNAL_WALL_CELLS=', SLC%N_INTERNAL_WALL_CELLS

         CALL SCARC_ASSIGN_INTERNAL_COORDS(NM, NL)
>>>>>>> origin/scarc

      ENDDO LEVEL_GMG_LEVEL_LOOP
   ENDDO MESHES_LOOP2
ENDIF LEVEL_GMG_IF
<<<<<<< HEAD



=======



>>>>>>> origin/scarc
!> -----------------------------------------------------------------------------------------
!> ---- ONLY TEMPORARILY FOR DEBUGGING PURPOSES !!!
!> -----------------------------------------------------------------------------------------
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
   DO NL=NLEVEL_MIN, NLEVEL_MAX
!      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL, 'SETUP_WALL_AA', 'WALL')
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACEINFO , NL, 'SETUP_WALL', 'FACE')
   ENDDO
ELSE
!   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NLEVEL_MIN, 'SETUP_WALL_BB', 'WALL')
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACEINFO , NLEVEL_MIN, 'SETUP_WALL', 'FACE')
ENDIF


IF (TYPE_DEBUG>NSCARC_DEBUG_EXTREME .AND. TYPE_LATEX  > NSCARC_LATEX_ALL) THEN
   MLATEX = 98
   DO NM = 1, NMESHES
      IF (PROCESS(NM) /= MYID) CYCLE
      CLATEX = "tables/wall  .tex"
      WRITE(CLATEX(12:13),'(i2.2)') NM
      OPEN (MLATEX, FILE=CLATEX)
      SLF => SCARC(NM)%LEVEL(NLEVEL_MIN)
      LATEX_OMESH_LOOP: DO NOM = 1, NMESHES
         IF (NOM == NM) CYCLE LATEX_OMESH_LOOP
         OS => SCARC(NM)%OSCARC(NOM)
         IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0)  CYCLE LATEX_OMESH_LOOP
         OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)
         CALL SCARC_LATEX_NOM(MLATEX,NM, NOM, NLEVEL_MIN)
      ENDDO LATEX_OMESH_LOOP
      CLOSE (MLATEX)
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_DECOMPOSITION


<<<<<<< HEAD
!> -----------------------------------------------------------------------------------------
!> --- SCARC_SETUP_CELL_INDEX
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CELL_INDEX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: I, J, K, IZERO, NO
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL
TYPE (SCARC_OBST_TYPE) , POINTER :: SO

SL => SCARC(NM)%LEVEL(NL)

ALLOCATE(SL%CELL_INDEX(0:SL%NX+1,0:SL%NY+1,0:SL%NZ+1),STAT=IZERO)
CALL ChkMemErr('SCARC_SETUP_CELL_INDEX','CELL_INDEX',IZERO)  

SL%CELL_INDEX = 0
SL%CELL_COUNT = 0

DO K=0,SL%NZ+1
   DO J=0,SL%NY+1
      DO I=0,1
         IF (SL%CELL_INDEX(I,J,K)==0) THEN
            SL%CELL_COUNT = SL%CELL_COUNT + 1
            SL%CELL_INDEX(I,J,K) = SL%CELL_COUNT
         ENDIF
      ENDDO
      DO I=SL%NX,SL%NX+1
         IF (SL%CELL_INDEX(I,J,K)==0) THEN
            SL%CELL_COUNT = SL%CELL_COUNT + 1
            SL%CELL_INDEX(I,J,K) = SL%CELL_COUNT
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO K=0,SL%NZ+1
   DO I=0,SL%NX+1
      DO J=0,1
         IF (SL%CELL_INDEX(I,J,K)==0) THEN
            SL%CELL_COUNT = SL%CELL_COUNT + 1
            SL%CELL_INDEX(I,J,K) = SL%CELL_COUNT
         ENDIF
      ENDDO
      DO J=SL%NY,SL%NY+1
         IF (SL%CELL_INDEX(I,J,K)==0) THEN
            SL%CELL_COUNT = SL%CELL_COUNT + 1
            SL%CELL_INDEX(I,J,K) = SL%CELL_COUNT
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO J=0,SL%NY+1
   DO I=0,SL%NX+1
      DO K=0,1
         IF (SL%CELL_INDEX(I,J,K)==0) THEN
            SL%CELL_COUNT = SL%CELL_COUNT + 1
            SL%CELL_INDEX(I,J,K) = SL%CELL_COUNT
         ENDIF
      ENDDO
      DO K=SL%NZ,SL%NZ+1
         IF (SL%CELL_INDEX(I,J,K)==0) THEN
            SL%CELL_COUNT = SL%CELL_COUNT + 1
            SL%CELL_INDEX(I,J,K) = SL%CELL_COUNT
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO NO=1,SL%N_OBST
   SO=>SCARC(NM)%LEVEL(NL)%OBST(NO)
   DO K=SO%K1,SO%K2+1
      DO J=SO%J1,SO%J2+1
         DO I=SO%I1,SO%I2+1
            IF (SL%CELL_INDEX(I,J,K)==0) THEN
               SL%CELL_COUNT = SL%CELL_COUNT + 1
               SL%CELL_INDEX(I,J,K) = SL%CELL_COUNT
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) 'SETTING UP CELL_INDEX FOR LEVEL :', NL
   DO J=0, SL%NY+1
      WRITE(LU_SCARC,'(6i6)') ((SL%CELL_INDEX(I,J,K),I=0,SL%NX+1),K=SL%NZ+1,0,-1)
      WRITE(LU_SCARC,*)
   ENDDO
ENDIF
END SUBROUTINE SCARC_SETUP_CELL_INDEX


!> -----------------------------------------------------------------------------------------
!> --- SCARC_SETUP_WALL_INDEX
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_INDEX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: I, J, K, ICG, IZERO, IW, IOR0
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

ALLOCATE(SL%WALL_INDEX(SL%CELL_COUNT, -3:3), STAT=IZERO) 
CALL ChkMemErr('SCARC_SETUP_WALL_INDEX','WALL_INDEX',IZERO)  
SL%WALL_INDEX = 0

DO IW = 1, SL%NW

   I = SL%WALL(IW)%IXW
   J = SL%WALL(IW)%IYW
   K = SL%WALL(IW)%IZW

   IOR0 = SL%WALL(IW)%IOR
   ICG = SL%CELL_INDEX(I,J,K)
  
   IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,'(A,i4,a,i4,a,3i4,a,i4)')  &
      'SETTING WALL_INDEX FOR IW=',IW,': ICG=',ICG,': I,J,K=',I,J,K,': IOR=',IOR0

   SL%WALL_INDEX(ICG,-IOR0) = IW

ENDDO

IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) 'SETTING UP WALL_INDEX FOR LEVEL :', NL, ': NW=', SL%NW
   DO IW=1, SL%NW
      I = SL%WALL(IW)%IXW
      J = SL%WALL(IW)%IYW
      K = SL%WALL(IW)%IZW
      ICG = SL%CELL_INDEX(I,J,K)
      WRITE(LU_SCARC,'(i6,A,7i6)') IW, ':', (SL%WALL_INDEX(ICG,IOR0), IOR0=-3,3)
   ENDDO
ENDIF
END SUBROUTINE SCARC_SETUP_WALL_INDEX


=======
>>>>>>> origin/scarc
!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> ACHTUNG: FUNKTIONIERT NUR FÜR SPEZIALFÄLLE, DIE AUCH FÜR GMG LAUFEN !!!
!> ACHTUNG: MUSS DRINGEND NOCH ERWEITERT WERDEN
!> -------------------------------------------------------------------------------------------------
<<<<<<< HEAD
SUBROUTINE SCARC_SETUP_INTERNAL_WALL_COORDS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, IO, IWC
INTEGER :: I, J, K
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLC, SLF
=======
SUBROUTINE SCARC_ASSIGN_INTERNAL_COORDS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWC, NW_INT
INTEGER :: IC, IO, IERR
INTEGER :: I, J, K
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLF, SLC
>>>>>>> origin/scarc

SLF => SCARC(NM)%LEVEL(NL-1)
SLC => SCARC(NM)%LEVEL(NL)

<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) 'SLF%N_EXTERNAL_WALL_CELLS=',SLF%N_EXTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) 'SLF%N_INTERNAL_WALL_CELLS=',SLF%N_INTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) 'SLF%NW=',SLF%NW
   WRITE(LU_SCARC,*) 'SLF%N_OBST=',SLF%N_OBST
=======
SLC%N_OBST = SLF%N_OBST
ALLOCATE(SLC%OBST(SLC%N_OBST), STAT=IERR)
CALL ChkMemErr('SCARC_ASSIGN_INTERNAL_COORDS','OBST',IERR)

IF (SLC%N_INTERNAL_WALL_CELLS == 0) RETURN

NW_INT = 0
IWC = SLC%N_EXTERNAL_WALL_CELLS + 1

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'ASSIGNING: BTYPE FINE'
   WRITE(LU_SCARC,'(4i4)') (SLF%WALL(I)%BTYPE, I=1, SLF%NW)
   WRITE(LU_SCARC,*) 'ASSIGNING: NOM FINE'
   WRITE(LU_SCARC,'(4i4)') (SLF%WALL(I)%NOM, I=1, SLF%NW)
   WRITE(LU_SCARC,*) 'ASSIGNING: BTYPE COARSE'
   WRITE(LU_SCARC,'(4i4)') (SLC%WALL(I)%BTYPE, I=1, SLC%NW)
   WRITE(LU_SCARC,*) 'ASSIGNING: NOM COARSE'
   WRITE(LU_SCARC,'(4i4)') (SLC%WALL(I)%NOM, I=1, SLC%NW)
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'SLC%N_EXTERNAL_WALL_CELLS=',SLC%N_EXTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) 'SLC%N_INTERNAL_WALL_CELLS=',SLC%N_INTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) 'SLC%NW=',SLC%NW
   WRITE(LU_SCARC,*) 'SLC%N_OBST=',SLC%N_OBST
<<<<<<< HEAD
   WRITE(LU_SCARC,*) 'SETUP_INTERNAL_WALL_COORDS: BOUNDARY_TYPE FINE'
   WRITE(LU_SCARC,'(16i4)') (SLF%WALL(I)%BOUNDARY_TYPE, I=1, SLF%NW)
   WRITE(LU_SCARC,*) 'SETUP_INTERNAL_WALL_COORDS: BTYPE FINE'
   WRITE(LU_SCARC,'(16i4)') (SLF%WALL(I)%BTYPE, I=1, SLF%NW)
   !WRITE(LU_SCARC,*) 'SETUP_INTERNAL_WALL_COORDS: NOM FINE'
   !WRITE(LU_SCARC,'(16i4)') (SLF%WALL(I)%NOM, I=1, SLF%NW)
   WRITE(LU_SCARC,*) 'SETUP_INTERNAL_WALL_COORDS: IOR FINE'
   WRITE(LU_SCARC,'(16i4)') (SLF%WALL(I)%IOR, I=1, SLF%NW)
   WRITE(LU_SCARC,*) 'SETUP_INTERNAL_WALL_COORDS: BOUNDARY_TYPE COARSE'
   WRITE(LU_SCARC,'(8i4)') (SLC%WALL(I)%BOUNDARY_TYPE, I=1, SLC%NW)
   WRITE(LU_SCARC,*) 'SETUP_INTERNAL_WALL_COORDS: BTYPE COARSE'
   WRITE(LU_SCARC,'(8i4)') (SLC%WALL(I)%BTYPE, I=1, SLC%NW)
   WRITE(LU_SCARC,*) 'SETUP_INTERNAL_WALL_COORDS: NOM COARSE'
   !WRITE(LU_SCARC,'(8i4)') (SLC%WALL(I)%NOM, I=1, SLC%NW)
   !WRITE(LU_SCARC,*) 'SETUP_INTERNAL_WALL_COORDS: IOR COARSE'
   WRITE(LU_SCARC,'(8i4)') (SLC%WALL(I)%IOR, I=1, SLC%NW)
ENDIF
IWC = SLC%N_EXTERNAL_WALL_CELLS + 1

!> Number of obstructions on coarse level is the same as on fine level
DO IO = 1, SLC%N_OBST

   !> Analyze IOR = 1
   I = SLC%OBST(IO)%I1
   DO K = SLC%OBST(IO)%K1+1, SLC%OBST(IO)%K2
      DO J = SLC%OBST(IO)%J1+1, SLC%OBST(IO)%J2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'A: Processing ', IC, I,J,K, 1, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         SLC%WALL(IWC)%IXW = I+1
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IXG = I
         SLC%WALL(IWC)%IYG = J
         SLC%WALL(IWC)%IZG = K
         SLC%WALL(IWC)%IOR = 1
         SLC%WALL(IWC)%BTYPE = NEUMANN
         SLC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
=======
   WRITE(LU_SCARC,*) 'IWC=',IWC
ENDIF

!> Number of obstructions on coarse level is the same as on fine level
DO IO = 1, SLF%N_OBST

   SLC%OBST(IO)%I1 = (SLF%OBST(IO)%I1+1)/2
   SLC%OBST(IO)%I2 = SLF%OBST(IO)%I2/2

   IF (TYPE_DIMENSION == NSCARC_DIMENSION_TWO) THEN
      SLC%OBST(IO)%J1 = 1
      SLC%OBST(IO)%J2 = 1
   ELSE
      SLC%OBST(IO)%J1 = (SLF%OBST(IO)%J1+1)/2
      SLC%OBST(IO)%J2 = SLF%OBST(IO)%J2/2
   ENDIF

   SLC%OBST(IO)%K1 = (SLF%OBST(IO)%K1+1)/2
   SLC%OBST(IO)%K2 = SLF%OBST(IO)%K2/2
   
   IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) '=========================================================='
      WRITE(LU_SCARC,*) 'FINE: OBSTRUCTIONS (',NM,',',NL,')', IWC, IO
      WRITE(LU_SCARC,'(7i6)') IO, &
                        SLF%OBST(IO)%I1, SLF%OBST(IO)%I2, &
                        SLF%OBST(IO)%J1, SLF%OBST(IO)%J2, &
                        SLF%OBST(IO)%K1, SLF%OBST(IO)%K2
      WRITE(LU_SCARC,*) 'COARSE: OBSTRUCTIONS (',NM,',',NL,')', IO
      WRITE(LU_SCARC,'(7i6)') IO, &
                        SLC%OBST(IO)%I1, SLC%OBST(IO)%I2, &
                        SLC%OBST(IO)%J1, SLC%OBST(IO)%J2, &
                        SLC%OBST(IO)%K1, SLC%OBST(IO)%K2
   ENDIF

   !> Analyze IOR = 1
   I = SLC%OBST(IO)%I1-1
   DO K = SLC%OBST(IO)%K1, SLC%OBST(IO)%K2
      DO J = SLC%OBST(IO)%J1, SLC%OBST(IO)%J2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
            WRITE(LU_SCARC,*) 'A: Processing ', IC, I,J,K, 1, IWC, NW_INT
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IOR = 1
         SLC%WALL(IWC)%BTYPE = NEUMANN
         IWC = IWC + 1
         NW_INT = NW_INT + 1
>>>>>>> origin/scarc
      ENDDO
   ENDDO

   !> Analyze IOR = -1
<<<<<<< HEAD
   I = SLC%OBST(IO)%I2
   DO K = SLC%OBST(IO)%K1+1, SLC%OBST(IO)%K2
      DO J = SLC%OBST(IO)%J1+1, SLC%OBST(IO)%J2
         IC = SLC%CCVAR(I+1, J, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'B: Processing ', IC, I,J,K,-1, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IXG = I+1
         SLC%WALL(IWC)%IYG = J
         SLC%WALL(IWC)%IZG = K
         SLC%WALL(IWC)%IOR =-1
         SLC%WALL(IWC)%BTYPE = NEUMANN
         SLC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
=======
   I = SLC%OBST(IO)%I2+1
   DO K = SLC%OBST(IO)%K1, SLC%OBST(IO)%K2
      DO J = SLC%OBST(IO)%J1, SLC%OBST(IO)%J2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
            WRITE(LU_SCARC,*) 'B: Processing ', IC, I,J,K,-1, IWC, NW_INT
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IOR =-1
         SLC%WALL(IWC)%BTYPE = NEUMANN
         IWC = IWC + 1
         NW_INT = NW_INT + 1
>>>>>>> origin/scarc
      ENDDO
   ENDDO

   !> Analyze IOR = 2
<<<<<<< HEAD
   J = SLC%OBST(IO)%J1
   DO K = SLC%OBST(IO)%K1+1, SLC%OBST(IO)%K2
      DO I = SLC%OBST(IO)%I1+1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'C: Processing ', IC, I,J,K,2, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J+1
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IXG = I
         SLC%WALL(IWC)%IYG = J
         SLC%WALL(IWC)%IZG = K
         SLC%WALL(IWC)%IOR = 2
         SLC%WALL(IWC)%BTYPE = NEUMANN
         SLC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
=======
   J = SLC%OBST(IO)%J1-1
   DO K = SLC%OBST(IO)%K1, SLC%OBST(IO)%K2
      DO I = SLC%OBST(IO)%I1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
            WRITE(LU_SCARC,*) 'C: Processing ', IC, I,J,K,2, IWC, NW_INT
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IOR = 2
         SLC%WALL(IWC)%BTYPE = NEUMANN
         IWC = IWC + 1
         NW_INT = NW_INT + 1
>>>>>>> origin/scarc
      ENDDO
   ENDDO

   !> Analyze IOR = -2
<<<<<<< HEAD
   J = SLC%OBST(IO)%J2
   DO K = SLC%OBST(IO)%K1+1, SLC%OBST(IO)%K2
      DO I = SLC%OBST(IO)%I1+1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J+1, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'D: Processing ', IC, I,J,K,-2, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IXG = I
         SLC%WALL(IWC)%IYG = J+1
         SLC%WALL(IWC)%IZG = K
         SLC%WALL(IWC)%IOR =-2
         SLC%WALL(IWC)%BTYPE = NEUMANN
         SLC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
=======
   J = SLC%OBST(IO)%J2+1
   DO K = SLC%OBST(IO)%K1, SLC%OBST(IO)%K2
      DO I = SLC%OBST(IO)%I1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
            WRITE(LU_SCARC,*) 'D: Processing ', IC, I,J,K,-2, IWC, NW_INT
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IOR =-2
         SLC%WALL(IWC)%BTYPE = NEUMANN
         IWC = IWC + 1
         NW_INT = NW_INT + 1
>>>>>>> origin/scarc
      ENDDO
   ENDDO

   !> Analyze IOR = 3
<<<<<<< HEAD
   K = SLC%OBST(IO)%K1
   DO J = SLC%OBST(IO)%J1+1, SLC%OBST(IO)%J2
      DO I = SLC%OBST(IO)%I1+1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'E: Processing ', IC, I,J,K,3, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K+1
         SLC%WALL(IWC)%IXG = I
         SLC%WALL(IWC)%IYG = J
         SLC%WALL(IWC)%IZG = K
         SLC%WALL(IWC)%IOR = 3
         SLC%WALL(IWC)%BTYPE = NEUMANN
         SLC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
=======
   K = SLC%OBST(IO)%K1-1
   DO J = SLC%OBST(IO)%J1, SLC%OBST(IO)%J2
      DO I = SLC%OBST(IO)%I1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
            WRITE(LU_SCARC,*) 'E: Processing ', IC, I,J,K,3, IWC, NW_INT
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IOR = 3
         SLC%WALL(IWC)%BTYPE = NEUMANN
         IWC = IWC + 1
         NW_INT = NW_INT + 1
>>>>>>> origin/scarc
      ENDDO
   ENDDO

   !> Analyze IOR = -3
<<<<<<< HEAD
   K = SLC%OBST(IO)%K2
   DO J = SLC%OBST(IO)%J1+1, SLC%OBST(IO)%J2
      DO I = SLC%OBST(IO)%I1+1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K+1, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'F: Processing ', IC, I,J,K,-3, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IXG = I
         SLC%WALL(IWC)%IYG = J
         SLC%WALL(IWC)%IZG = K+1
         SLC%WALL(IWC)%IOR =-3
         SLC%WALL(IWC)%BTYPE = NEUMANN
         SLC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
=======
   K = SLC%OBST(IO)%K2+1
   DO J = SLC%OBST(IO)%J1, SLC%OBST(IO)%J2
      DO I = SLC%OBST(IO)%I1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
            WRITE(LU_SCARC,*) 'F: Processing ', IC, I,J,K,-3, IWC, NW_INT
         SLC%WALL(IWC)%IXW = I
         SLC%WALL(IWC)%IYW = J
         SLC%WALL(IWC)%IZW = K
         SLC%WALL(IWC)%IOR =-3
         SLC%WALL(IWC)%BTYPE = NEUMANN
         IWC = IWC + 1
         NW_INT = NW_INT + 1
>>>>>>> origin/scarc
      ENDDO
   ENDDO
         
ENDDO

<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) 'IWC=',IWC
   WRITE(LU_SCARC,*) 'SLC%N_INTERNAL_WALL_CELLS=',SLC%N_INTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) 'SLC%NCS=',SLC%NCS
   DO IWC = 1, SLC%NW
      WRITE(LU_SCARC,'(9I4)') SLC%WALL(IWC)%IOR, &
                        SLC%WALL(IWC)%IXW, SLC%WALL(IWC)%IYW, SLC%WALL(IWC)%IZW, &
                        SLC%WALL(IWC)%IXG, SLC%WALL(IWC)%IYG, SLC%WALL(IWC)%IZG, &
                        SLC%WALL(IWC)%BTYPE, SLC%WALL(IWC)%BOUNDARY_TYPE
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_INTERNAL_WALL_COORDS
=======
IF (NW_INT /= SLC%N_INTERNAL_WALL_CELLS) THEN
   WRITE(*,*) 'ERROR with NW_INT in ASSIGN_INTERNAL_COORDS'
   STOP
ENDIF

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'NW_INT=',NW_INT
   WRITE(LU_SCARC,*) 'SLC%N_INTERNAL_WALL_CELLS=',SLC%N_INTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) 'SLC%NCS=',SLC%NCS
   DO IWC = 1, SLC%NW
      WRITE(LU_SCARC,*) SLC%WALL(IWC)%IOR, SLC%WALL(IWC)%IXW, SLC%WALL(IWC)%IYW, SLC%WALL(IWC)%IZW, &
                        SLC%WALL(IWC)%BTYPE
   ENDDO
ENDIF

END SUBROUTINE SCARC_ASSIGN_INTERNAL_COORDS
>>>>>>> origin/scarc


!> -------------------------------------------------------------------------------------------------
!> Correct BTYPE related to internal obstructions on ghost cells
!> -------------------------------------------------------------------------------------------------
<<<<<<< HEAD
SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWG
!INTEGER, POINTER, DIMENSION(:,:,:):: CELL_INDEX
!INTEGER, POINTER, DIMENSION(:,:)  :: WALL_INDEX
INTEGER :: IX, IY, IZ, IOR0, BTYPE0
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL
TYPE (MESH_TYPE), POINTER :: M

SL => SCARC(NM)%LEVEL(NL)
M  => MESHES(NM)

!IF (NL == NLEVEL_MIN) THEN
!   CELL_INDEX => M%CELL_INDEX
!   WALL_INDEX => M%WALL_INDEX
!ELSE
!   CELL_INDEX => SL%CELL_INDEX
!   WALL_INDEX => SL%WALL_INDEX
!ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) 'SL%WALL(.)%NOM:'
   WRITE(LU_SCARC,'(4i4)') (SL%WALL(IWG)%NOM, IWG=1, SL%NW)
   WRITE(LU_SCARC,*) 'SL%WALL(.)%BTYPE:'
   WRITE(LU_SCARC,'(4i4)') (SL%WALL(IWG)%BTYPE, IWG=1, SL%NW)
   WRITE(LU_SCARC,*) 'SL%WALL(.)%BOUNDARY_TYPE:'
   WRITE(LU_SCARC,'(4i4)') (SL%WALL(IWG)%BOUNDARY_TYPE, IWG=1, SL%NW)
ENDIF

DO IWG = 1, SL%N_EXTERNAL_WALL_CELLS

   IF (SL%WALL(IWG)%NOM /= 0) CYCLE

   IX = SL%WALL(IWG)%IXW
   IY = SL%WALL(IWG)%IYW
   IZ = SL%WALL(IWG)%IZW

   IOR0   = SL%WALL(IWG)%IOR
   BTYPE0 = SL%WALL(IWG)%BTYPE

  ! ICG = CELL_INDEX(IX, IY, IZ)
  ! IWG = WALL_INDEX(ICG, IOR0)

  ! IF (SL%WALL(IWG)%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY) SL%WALL(IWG)%BTYPE=NEUMANN
  ! IF (SL%CCVAR(IX, IY, IZ, CGSC) == IS_SOLID) SL%WALL(IWG)%BTYPE=NEUMANN

   IF (SL%WALL(IWG)%BOUNDARY_TYPE == SOLID_BOUNDARY) SL%WALL(IWG)%BTYPE=NEUMANN

   IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
     WRITE(LU_SCARC,'(A, 6I6)') 'SETTING TO NEUMANN: ', IWG, IX, IY, IZ, BTYPE0, SL%WALL(IWG)%BTYPE

ENDDO

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) 'AFTER CORRECTING:'
   WRITE(LU_SCARC,*) 'SL%WALL(.)%BTYPE:'
   WRITE(LU_SCARC,'(4i4)') (SL%WALL(IWG)%BTYPE, IWG=1, SL%NW)
ENDIF

END SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS
=======
SUBROUTINE SCARC_CORRECT_GHOSTCELL_TYPES(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWG
INTEGER :: IX, IY, IZ
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'M%WALL(.)%NOM:'
   WRITE(LU_SCARC,'(4i4)') (SL%WALL(IWG)%NOM, IWG=1, SL%NW)
   WRITE(LU_SCARC,*) 'M%WALL(.)%BTYPE:'
   WRITE(LU_SCARC,'(4i4)') (SL%WALL(IWG)%BTYPE, IWG=1, SL%NW)
ENDIF

DO IWG = 1, SL%N_EXTERNAL_WALL_CELLS
   IF (SL%WALL(IWG)%NOM == 0) CYCLE
   IX = SL%WALL(IWG)%IXG
   IY = SL%WALL(IWG)%IYG
   IZ = SL%WALL(IWG)%IZG
   IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
     WRITE(LU_SCARC,*) '1: CORRECTING ', IWG, IX, IY, IZ, SL%CCVAR(IX, IY, IZ, CGSC), SL%WALL(IWG)%BTYPE
   IF (SL%CCVAR(IX, IY, IZ, CGSC) == IS_SOLID) SL%WALL(IWG)%BTYPE=NEUMANN
   IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
     WRITE(LU_SCARC,*) '2: CORRECTING ', IWG, IX, IY, IZ, SL%CCVAR(IX, IY, IZ, CGSC), SL%WALL(IWG)%BTYPE
ENDDO

END SUBROUTINE SCARC_CORRECT_GHOSTCELL_TYPES
>>>>>>> origin/scarc


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
<<<<<<< HEAD
INTEGER FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IXC, IYC, IZC
INTEGER :: IXF, IYF, IZF
INTEGER :: IWC, ICF(4), IWF(4), IOR0
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
INTEGER, POINTER, DIMENSION(:,:)   :: WALL_INDEX
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLF, SLC
TYPE (MESH_TYPE), POINTER :: M

SLC => SCARC(NM)%LEVEL(NL)
SLF => SCARC(NM)%LEVEL(NL-1)
M   => MESHES(NM)

IF (NL == NLEVEL_MIN+1) THEN
   CELL_INDEX => M%CELL_INDEX
   WALL_INDEX => M%WALL_INDEX
ELSE
   CELL_INDEX => SLF%CELL_INDEX
   WALL_INDEX => SLF%WALL_INDEX
ENDIF

IWC = 0
ICF = 0
IWF = 0

SELECT CASE (TYPE_DIMENSION)
   CASE (NSCARC_DIMENSION_TWO)

      IYC = 1
      IYF = 1

      !> IOR = 1
      IOR0 = 1
      DO IZC = 1, SLC%NZ 
         IZF = 2*IZC - 1
         ICF(1) = CELL_INDEX(1  , IYF  , IZF  )
         ICF(2) = CELL_INDEX(1  , IYF  , IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, 1, IYC, IZC, 1, IYF, IZF, ICF
      ENDDO

      !> IOR = -1
      IOR0 = -1
      DO IZC = 1, SLC%NZ 
         IZF = 2*IZC - 1
         ICF(1) = CELL_INDEX(SLF%NX, IYF  , IZF  )
         ICF(2) = CELL_INDEX(SLF%NX, IYF  , IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, SLC%NX, IYC, IZC, SLF%NX, IYF, IZF, ICF
      ENDDO

      !> IOR = 2
      IOR0 = 2
      DO IZC = 1, SLC%NZ 
         IZF = 2*IZC - 1
         DO IXC = 1, SLC%NX 
            IXF = 2*IXC - 1
            ICF(1) = CELL_INDEX(IXF  , IYF, IZF  )
            ICF(2) = CELL_INDEX(IXF+1, IYF, IZF  )
            ICF(3) = CELL_INDEX(IXF  , IYF, IZF+1)
            ICF(4) = CELL_INDEX(IXF+1, IYF, IZF+1)
            IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, IXC, 1, IZC, IXF, 1, IZF, ICF
         ENDDO
      ENDDO

      !> IOR = -2
      IOR0 = -2
      IXF  = 1
      DO IZC = 1, SLC%NZ 
         IZF = 2*IZC - 1
         DO IXC = 1, SLC%NX 
            IXF = 2*IXC - 1
            ICF(1) = CELL_INDEX(IXF    , IYF, IZF  )
            ICF(2) = CELL_INDEX(IXF+1  , IYF, IZF  )
            ICF(3) = CELL_INDEX(IXF    , IYF, IZF+1)
            ICF(4) = CELL_INDEX(IXF+1  , IYF, IZF+1)
            IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, IXC, SLC%NY, IZC, IXF, SLF%NY, IZF, ICF
         ENDDO
      ENDDO

      !> IOR = 3
      IOR0 = 3
      DO IXC = 1, SLC%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF    , IYF  , 1)
         ICF(2) = CELL_INDEX(IXF+1  , IYF  , 1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, IXC, IYC, 1, IXF, IYF, 1, ICF
      ENDDO

      !> IOR = -3
      IOR0 = -3
      DO IXC = 1, SLC%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF  , IYF  , SLF%NZ)
         ICF(2) = CELL_INDEX(IXF+1, IYF  , SLF%NZ)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, IXC, IYC, SLC%NZ, IXF, IYF, SLF%NZ, ICF
      ENDDO

   CASE (NSCARC_DIMENSION_THREE)

      !> IOR = 1
      IOR0 = 1
      DO IZC = 1, SLC%NZ 
         IZF = 2*IZC - 1
         DO IYC = 1, SLC%NY 
            IYF = 2*IYC - 1
            ICF(1) = CELL_INDEX(1  , IYF  , IZF  )
            ICF(2) = CELL_INDEX(1  , IYF+1, IZF  )
            ICF(3) = CELL_INDEX(1  , IYF  , IZF+1)
            ICF(4) = CELL_INDEX(1  , IYF+1, IZF+1)
            IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, 1, IYC, IZC, 1, IYF, IZF, ICF
         ENDDO
      ENDDO

      !> IOR = -1
      IOR0 = -1
      DO IZC = 1, SLC%NZ 
         IZF = 2*IZC - 1
         DO IYC = 1, SLC%NY 
            IYF = 2*IYC - 1
            ICF(1) = CELL_INDEX(SLF%NX, IYF  , IZF  )
            ICF(2) = CELL_INDEX(SLF%NX, IYF+1, IZF  )
            ICF(3) = CELL_INDEX(SLF%NX, IYF  , IZF+1)
            ICF(4) = CELL_INDEX(SLF%NX, IYF+1, IZF+1)
            IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, SLC%NX, IYC, IZC, SLF%NX, IYF, IZF, ICF
         ENDDO
      ENDDO

      !> IOR = 2
      IOR0 = 2
      DO IZC = 1, SLC%NZ 
         IZF = 2*IZC - 1
         DO IXC = 1, SLC%NX 
            IXF = 2*IXC - 1
            ICF(1) = CELL_INDEX(IXF  , 1, IZF  )
            ICF(2) = CELL_INDEX(IXF+1, 1, IZF  )
            ICF(3) = CELL_INDEX(IXF  , 1, IZF+1)
            ICF(4) = CELL_INDEX(IXF+1, 1, IZF+1)
            IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, IXC, 1, IZC, IXF, 1, IZF, ICF
         ENDDO
      ENDDO

      !> IOR = -2
      IOR0 = -2
      IXF  = 1
      DO IZC = 1, SLC%NZ 
         IZF = 2*IZC - 1
         DO IXC = 1, SLC%NX 
            IXF = 2*IXC - 1
            ICF(1) = CELL_INDEX(IXF    , SLF%NY, IZF  )
            ICF(2) = CELL_INDEX(IXF+1  , SLF%NY, IZF  )
            ICF(3) = CELL_INDEX(IXF    , SLF%NY, IZF+1)
            ICF(4) = CELL_INDEX(IXF+1  , SLF%NY, IZF+1)
            IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, IXC, SLC%NY, IZC, IXF, SLF%NY, IZF, ICF
         ENDDO
      ENDDO

      !> IOR = 3
      IOR0 = 3
      DO IYC = 1, SLC%NY 
         IYF = 2*IYC - 1
         DO IXC = 1, SLC%NX 
            IXF = 2*IXC - 1
            ICF(1) = CELL_INDEX(IXF    , IYF  , 1)
            ICF(2) = CELL_INDEX(IXF+1  , IYF  , 1)
            ICF(3) = CELL_INDEX(IXF    , IYF+1, 1)
            ICF(4) = CELL_INDEX(IXF+1  , IYF+1, 1)
            IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, IXC, IYC, 1, IXF, IYF, 1, ICF
         ENDDO
      ENDDO

      !> IOR = -3
      IOR0 = -3
      DO IYC = 1, SLC%NY 
         IYF = 2*IYC - 1
         DO IXC = 1, SLC%NX 
            IXF = 2*IXC - 1
            ICF(1) = CELL_INDEX(IXF  , IYF  , SLF%NZ)
            ICF(2) = CELL_INDEX(IXF+1, IYF  , SLF%NZ)
            ICF(3) = CELL_INDEX(IXF  , IYF+1, SLF%NZ)
            ICF(4) = CELL_INDEX(IXF+1, IYF+1, SLF%NZ)
            IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1000) NL, IOR0,IWC, IXC, IYC, SLC%NZ, IXF, IYF, SLF%NZ, ICF
         ENDDO
      ENDDO

END SELECT

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) NL, ': FINAL  NUMBER OF EXTERNAL WALL CELLS ', IWC
SCARC_COUNT_EXTERNAL_WALL_CELLS = IWC
RETURN

1000 FORMAT('LEVEL=', I3,' : IOR=', I3,' : IWC=',I6,' : IXC, IYC, IZC=',3I4,' : IXF, IYF, IZF=',3I4,': ICF=',8I6)
END FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS

!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWC, NW_INT
INTEGER :: IC, IO, IERR
INTEGER :: I, J, K
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLF, SLC

SLF => SCARC(NM)%LEVEL(NL-1)
SLC => SCARC(NM)%LEVEL(NL)

SLC%N_OBST = SLF%N_OBST
ALLOCATE(SLC%OBST(SLC%N_OBST), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_INTERNAL_WALL_COORDS','OBST',IERR)

NW_INT = 0
IWC = SLC%N_EXTERNAL_WALL_CELLS + 1


!> Number of obstructions on coarse level is the same as on fine level
DO IO = 1, SLF%N_OBST

   SLC%OBST(IO)%I1 = (SLF%OBST(IO)%I1+1)/2
   SLC%OBST(IO)%I2 = SLF%OBST(IO)%I2/2

   IF (TYPE_DIMENSION == NSCARC_DIMENSION_TWO) THEN
      SLC%OBST(IO)%J1 = 0
      SLC%OBST(IO)%J2 = 1
   ELSE
      SLC%OBST(IO)%J1 = (SLF%OBST(IO)%J1+1)/2
      SLC%OBST(IO)%J2 = SLF%OBST(IO)%J2/2
   ENDIF

   SLC%OBST(IO)%K1 = (SLF%OBST(IO)%K1+1)/2
   SLC%OBST(IO)%K2 = SLF%OBST(IO)%K2/2
   
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
      WRITE(LU_SCARC,*) '=========================================================='
      WRITE(LU_SCARC,*) 'FINE: OBSTRUCTIONS (',NM,',',NL,')', IWC, IO
      WRITE(LU_SCARC,'(7i6)') IO, &
                        SLF%OBST(IO)%I1, SLF%OBST(IO)%I2, &
                        SLF%OBST(IO)%J1, SLF%OBST(IO)%J2, &
                        SLF%OBST(IO)%K1, SLF%OBST(IO)%K2
      WRITE(LU_SCARC,*) 'COARSE: OBSTRUCTIONS (',NM,',',NL,')', IO
      WRITE(LU_SCARC,'(7i6)') IO, &
                        SLC%OBST(IO)%I1, SLC%OBST(IO)%I2, &
                        SLC%OBST(IO)%J1, SLC%OBST(IO)%J2, &
                        SLC%OBST(IO)%K1, SLC%OBST(IO)%K2
   ENDIF

   !> Analyze IOR = 1
   I = SLC%OBST(IO)%I1
   DO K = SLC%OBST(IO)%K1+1, SLC%OBST(IO)%K2
      DO J = SLC%OBST(IO)%J1+1, SLC%OBST(IO)%J2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'A: Processing ', IC, I,J,K, 1, IWC, NW_INT
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -1
   I = SLC%OBST(IO)%I2
   DO K = SLC%OBST(IO)%K1+1, SLC%OBST(IO)%K2
      DO J = SLC%OBST(IO)%J1+1, SLC%OBST(IO)%J2
         IC = SLC%CCVAR(I+1, J, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'B: Processing ', IC, I,J,K,-1, IWC, NW_INT
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 2
   J = SLC%OBST(IO)%J1
   DO K = SLC%OBST(IO)%K1+1, SLC%OBST(IO)%K2
      DO I = SLC%OBST(IO)%I1+1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'C: Processing ', IC, I,J,K,2, IWC, NW_INT
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -2
   J = SLC%OBST(IO)%J2
   DO K = SLC%OBST(IO)%K1+1, SLC%OBST(IO)%K2
      DO I = SLC%OBST(IO)%I1+1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J+1, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'D: Processing ', IC, I,J,K,-2, IWC, NW_INT
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 3
   K = SLC%OBST(IO)%K1
   DO J = SLC%OBST(IO)%J1+1, SLC%OBST(IO)%J2
      DO I = SLC%OBST(IO)%I1+1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'E: Processing ', IC, I,J,K,3, IWC, NW_INT
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -3
   K = SLC%OBST(IO)%K2
   DO J = SLC%OBST(IO)%J1+1, SLC%OBST(IO)%J2
      DO I = SLC%OBST(IO)%I1+1, SLC%OBST(IO)%I2
         IC = SLC%CCVAR(I, J, K+1, UNKH)
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'F: Processing ', IC, I,J,K,-3, IWC, NW_INT
         IF (IC == IS_UNDEFINED) CYCLE
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'YES'
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO
         
ENDDO

SCARC_COUNT_INTERNAL_WALL_CELLS = NW_INT

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) 'NW_INT=',NW_INT
   WRITE(LU_SCARC,*) 'SLC%N_INTERNAL_WALL_CELLS=',SLC%N_INTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) 'SLC%NCS=',SLC%NCS
ENDIF

END FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS



!> -------------------------------------------------------------------------------------------------
!> Count external wall cells on face IOR
!> -------------------------------------------------------------------------------------------------
LOGICAL FUNCTION IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, NCNT, NM, NL)
INTEGER, INTENT(IN) :: IOR0, NM, NL, NCNT
INTEGER, DIMENSION(:), INTENT(IN) :: ICF
INTEGER, POINTER, DIMENSION(:,:), INTENT(IN) :: WALL_INDEX
INTEGER:: I, IWF_LAST, IWF(4)
REAL(EB) :: BSUM
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLF

IS_EXTERNAL_WALLCELL = .FALSE.

SLF  => SCARC(NM)%LEVEL(NL)

DO I = 1, NCNT
   IWF(I) = WALL_INDEX(ICF(I), -IOR0)
ENDDO

BSUM = 0.0_EB
IWF_LAST = 0

DO I = 1, NCNT
   IF (IWF(I)>0) THEN
      BSUM = BSUM + REAL(SLF%WALL(IWF(I))%BOUNDARY_TYPE,EB)
      IWF_LAST = IWF(I)
      IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'BOUNDARY_TYPE(',IWF(I),')=',SLF%WALL(IWF(I))%BOUNDARY_TYPE
   ENDIF
ENDDO

IF (IWF_LAST == 0) RETURN

IF (ABS(BSUM/REAL(NCNT,EB) - REAL(SLF%WALL(IWF_LAST)%BOUNDARY_TYPE,EB)) < 1E-12) THEN
   IS_EXTERNAL_WALLCELL = .TRUE.
   RETURN
ELSE
   WRITE(*,*) 'Wrong sum of BOUNDARY_TYPE for IOR =  1'
   WRITE(LU_SCARC,*) 'Wrong sum of BOUNDARY_TYPE for IOR =  1'
   STOP
ENDIF

END FUNCTION IS_EXTERNAL_WALLCELL
=======
INTEGER FUNCTION SCARC_INTERNAL_WALL_CELLS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IX, IY, IZ, NW_INT
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

NW_INT = 0
DO IZ = 1, SL%NZ 
   DO IY = 1, SL%NY 
      DO IX = 1, SL%NX 
         IF (SL%CCVAR(IX  , IY  , IZ  , UNKH) /= IS_UNDEFINED) CYCLE
         IF (SL%CCVAR(IX-1, IY  , IZ  , UNKH) /= IS_UNDEFINED) NW_INT = NW_INT+1
         IF (SL%CCVAR(IX+1, IY  , IZ  , UNKH) /= IS_UNDEFINED) NW_INT = NW_INT+1
         IF (SL%CCVAR(IX  , IY-1, IZ  , UNKH) /= IS_UNDEFINED) NW_INT = NW_INT+1
         IF (SL%CCVAR(IX  , IY+1, IZ  , UNKH) /= IS_UNDEFINED) NW_INT = NW_INT+1
         IF (SL%CCVAR(IX  , IY  , IZ-1, UNKH) /= IS_UNDEFINED) NW_INT = NW_INT+1
         IF (SL%CCVAR(IX  , IY  , IZ+1, UNKH) /= IS_UNDEFINED) NW_INT = NW_INT+1
      ENDDO
   ENDDO
ENDDO

SCARC_INTERNAL_WALL_CELLS = NW_INT
RETURN

END FUNCTION SCARC_INTERNAL_WALL_CELLS
>>>>>>> origin/scarc


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLCELL_NEIGHBOR(NX1, NX2, NY1, NY2, NZ1, NZ2, IWG, IOR0, NM, NOM, NL)
INTEGER, INTENT(IN) :: NX1, NX2, NY1, NY2, NZ1, NZ2
INTEGER, INTENT(IN) :: IWG, IOR0, NM, NOM, NL
INTEGER :: NOMX, NOMY, NOMZ
INTEGER :: ICG, ICO, ICE, IWL, ICPL, IERR=0, IX, IY, IZ, JL
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, OSL

SL  => SCARC(NM)%LEVEL(NL)
OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

!> store basic information about neighbor and orientation
OSL%IOR  = IOR0
OSL%NOM  = NOM

ICO  = SL%ICO_MARKER
ICE  = SL%ICE_MARKER

ICG  = OSL%ICG_MARKER
IWL  = OSL%IWL_MARKER

OSL%NCPLS = OSL%NCPL
IF (OSL%NCG == OSL%NWL) THEN
   OSL%NCPLR = OSL%NCPL
ELSE IF (OSL%NCG == 2*OSL%NWL) THEN
   OSL%NCPLR =  1
ELSE IF (OSL%NWL == 2*OSL%NCG) THEN
   OSL%NCPLR =  2
ENDIF

<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
   WRITE(LU_SCARC,'(i4,a,2i4)') NM, ' ========== CALLING WALLCELL_NEIGHBOR FOR ', IOR0, OSL%NCG
   WRITE(LU_SCARC,*) 'XXX: OSL%NCPL=',OSL%NCPL
   WRITE(LU_SCARC,*) 'XXX: OSL%IOR =',OSL%IOR
   WRITE(LU_SCARC,*) 'XXX: OSL%NCG =',OSL%NCG
   WRITE(LU_SCARC,*) 'XXX: OSL%NWL =',OSL%NWL
   WRITE(LU_SCARC,*) 'XXX: OSL%NCPLR =',OSL%NCPLR
   WRITE(LU_SCARC,*) 'XXX: OSL%NCPLS =',OSL%NCPLS
   WRITE(LU_SCARC,*) 'XXX: ICE     =', ICE
   WRITE(LU_SCARC,*) 'XXX: ICO     =', ICO
   WRITE(LU_SCARC,*) 'XXX: ICG     =', ICG
   WRITE(LU_SCARC,*) 'XXX: IWL     =', IWL
   WRITE(LU_SCARC,*) 'XXX: IWG     =', IWG
   WRITE(LU_SCARC,*) 'XXX: NX1, NX2     =', NX1, NX2
   WRITE(LU_SCARC,*) 'XXX: NY1, NY2     =', NY1, NY2
   WRITE(LU_SCARC,*) 'XXX: NZ1, NZ2     =', NZ1, NZ2
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'1: ICE_MARKER=',SL%ICE_MARKER
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'1: ICO_MARKER=',SL%ICO_MARKER
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'1: ICG_MARKER=',OSL%ICG_MARKER
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'1: IWL_MARKER=',OSL%IWL_MARKER
ENDIF

!> set neighboring coordinates
SL%WALL(IWG)%IXN(1) = NX1
SL%WALL(IWG)%IXN(2) = NX2
SL%WALL(IWG)%IYN(1) = NY1
SL%WALL(IWG)%IYN(2) = NY2
SL%WALL(IWG)%IZN(1) = NZ1
SL%WALL(IWG)%IZN(2) = NZ2

!> allocate pointer arrays for extended, ghost and neighboring cells
ALLOCATE(SL%WALL(IWG)%ICE(OSL%NCPL), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION','ICE',IERR)
SL%WALL(IWG)%ICE = NSCARC_DUMMY

ALLOCATE(SL%WALL(IWG)%ICG(OSL%NCPL), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION','ICG',IERR)
SL%WALL(IWG)%ICG = NSCARC_DUMMY

IWL = IWL + 1
ICO = ICO + 1

NOMX = MESHES(NOM)%IBAR
NOMY = MESHES(NOM)%JBAR
NOMZ = MESHES(NOM)%KBAR

IF (NL > 1) THEN
   DO JL = 2, NL
      NOMX = MESHES(NOM)%IBAR/NL
      IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) NOMY = MESHES(NOM)%JBAR/NL
      NOMZ = MESHES(NOM)%KBAR/NL
   ENDDO
ENDIF

<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'MESHES(',NM,')%IBAR=',MESHES(NM)%IBAR
   WRITE(LU_SCARC,*) 'MESHES(',NM,')%JBAR=',MESHES(NM)%JBAR
   WRITE(LU_SCARC,*) 'MESHES(',NM,')%KBAR=',MESHES(NM)%KBAR
   WRITE(LU_SCARC,*) 'MESHES(',NOM,')%IBAR=',MESHES(NOM)%IBAR
   WRITE(LU_SCARC,*) 'MESHES(',NOM,')%JBAR=',MESHES(NOM)%JBAR
   WRITE(LU_SCARC,*) 'MESHES(',NOM,')%KBAR=',MESHES(NOM)%KBAR
   WRITE(LU_SCARC,*) 'NOMX=',NOMX
   WRITE(LU_SCARC,*) 'NOMY=',NOMY
   WRITE(LU_SCARC,*) 'NOMZ=',NOMZ
ENDIF

!> store information about overlapped cells and set mapping arrays
ICPL = 0
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ICPL = ICPL + 1
         ICG  = ICG  + 1
         ICE  = ICE  + 1

<<<<<<< HEAD
         IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
         IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
            WRITE(LU_SCARC,'(8(a,i4))') 'IX=',IX, ' :IY=',IY, ' :IZ=',IZ, ' :IWL=',IWL,&
                                        ' :ICE=', ICE, ' :ICG=',ICG
         ENDIF

         SL%WALL(IWG)%ICE(ICPL)    = ICE                       !> number of extended grid cell
         SL%WALL(IWG)%ICG(ICPL)    = ICG                       !> number of ghost grid cell
   
         SL%ICE_TO_IWG(ICE)        = IWG                       !> map extended cell to global wall cell
         SL%ICE_TO_IWL(ICE)        = IWL                       !> map extended cell to local wall cell
         SL%ICE_TO_ICG(ICE)        = ICG                       !> map extended cell to ghost cell

         OSL%ICG_TO_IWG(ICG)       = IWG                       !> map ghost cell to global wall cell
         OSL%ICG_TO_ICO(ICG)       = ICO                       !> map ghost cell to extended grid cell
         OSL%ICG_TO_ICE(ICG)       = SL%WALL(IWG)%ICE(ICPL)    !> map ghost cell to extended grid cell

         OSL%IWL_TO_ICG(IWL, ICPL) = ICG                       !> map extended cell to ghost cell
      ENDDO
   ENDDO
ENDDO

SL%WALL(IWG)%ICO = ICO                                   !> number of overlapping cell
SL%WALL(IWG)%IWL = IWL                                   !> number of local wall cell

OSL%IWL_TO_IWG(IWL) = IWG                                !> map local wall cell to global wall cell
OSL%IWL_TO_ICO(IWL) = ICO                                !> map local wall cell to internal grid cell
OSL%IWL_TO_ICW(IWL) = SL%WALL(IWG)%ICW                   !> map local wall cell to internal grid cell

OSL%IWL_MARKER = IWL                                     !> store local wall cell marker
OSL%ICG_MARKER = ICG                                     !> store ghost cell counter

SL%ICO_MARKER  = ICO                                     !> store overlapping cell counter
SL%ICE_MARKER  = ICE                                     !> store extended cell counter

<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'NM=',NM,': NOM=',NOM,': IWL          =',IWL
   WRITE(LU_SCARC,*) 'NM=',NM,': NOM=',NOM,': ICE          =',ICE
   WRITE(LU_SCARC,*) 'NM=',NM,': NOM=',NOM,': ICO          =',ICO
   WRITE(LU_SCARC,*) 'NM=',NM,': NOM=',NOM,': OSL%NWL     =',OSL%NWL
   WRITE(LU_SCARC,*) 'NM=',NM,': NOM=',NOM,': OSL%IOR     =',OSL%IOR
   WRITE(LU_SCARC,*) 'NM=',NM,': NOM=',NOM,': OSL%NCG     =',OSL%NCG
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,' ICO=',SL%WALL(IWG)%ICO
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,' ICE=',SL%WALL(IWG)%ICE
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,' ICE=',SL%WALL(IWG)%ICG
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,',2: ICE_MARKER=',SL%ICE_MARKER
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,',2: ICO_MARKER=',SL%ICO_MARKER
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,',2: ICG_MARKER=',OSL%ICG_MARKER
   WRITE(LU_SCARC,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,',2: IWL_MARKER=',OSL%IWL_MARKER
ENDIF

END SUBROUTINE SCARC_SETUP_WALLCELL_NEIGHBOR


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MAPPINGS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: NCS, NCO, NCE, NW, IERR
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL  => SCARC(NM)%LEVEL(NL)
NCS = SL%NCS
NCE = SL%NCE
NCO = SL%NCO
NW  = SL%NW

<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
   WRITE(LU_SCARC,'(8(a,i4))') 'SETUP MAPPINGS FOR NCS=',NCS, ': NCE=',NCE,': NCO=',NCO,': NW=',NW,': NL=',NL
   WRITE(LU_SCARC,'(a,2i4)') 'ALLOCATING ICE_TO_IWG in length ', NCS+1, NCE
ENDIF

!> Allocate mappings for ICE
ALLOCATE(SL%ICE_TO_IWG(NCS+1:NCE), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','ICE_TO_IWG',IERR)
SL%ICE_TO_IWG=0

ALLOCATE(SL%ICE_TO_IWL(NCS+1:NCE), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','ICE_TO_IWL',IERR)
SL%ICE_TO_IWL=0

ALLOCATE(SL%ICE_TO_ICG(NCS+1:NCE), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','ICE_TO_ICG',IERR)
SL%ICE_TO_ICG=0


!> Allocate remaining mappings (outdated)
ALLOCATE(SL%INTERNAL_BDRY_CELL(NCS), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','INTERNAL_BDRY_CELL',IERR)
SL%INTERNAL_BDRY_CELL=0

<<<<<<< HEAD
=======
ALLOCATE(SL%WALL_INDEX(NCS, -3:3), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','WALL_INDEX',IERR)
SL%WALL_INDEX = 0

>>>>>>> origin/scarc

END SUBROUTINE SCARC_SETUP_MAPPINGS


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OMAPPINGS(NM, NOM, NL)
INTEGER, INTENT(IN) :: NM, NOM, NL
INTEGER:: NWL, NCG, NCPL, IERR
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSL

!> Get corresponding dimensions
OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

NWL  = OSL%NWL
NCG  = OSL%NCG
NCPL = OSL%NCPL

!> initialize some counters
OSL%IWL = 0
OSL%ICG = 0

!WRITE(LU_SCARC,'(8(a,i4))') &
!  'SETUP OMAPPINGS FOR NWL=',NWL, ' :NCG=',NCG,': NCPL=',NCPL,': NOM=',NOM,': NL=',NL

!> Allocate mapping array for local to global wall cells
ALLOCATE(OSL%IWL_TO_IWG(NWL), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','IWL_TO_IWG',IERR)
OSL%IWL_TO_IWG = 0

ALLOCATE(OSL%IWL_TO_ICW(NWL), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','IWL_TO_ICW',IERR)
OSL%IWL_TO_ICW = 0

ALLOCATE(OSL%IWL_TO_ICO(NWL), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','IWL_TO_ICO',IERR)
OSL%IWL_TO_ICO = 0

<<<<<<< HEAD
IF (PRES_ON_WHOLE_DOMAIN) THEN
=======
IF (PRES_ON_WHOLE_DOMAIN_SCRC) THEN
>>>>>>> origin/scarc
   ALLOCATE(OSL%ICN_TO_ICE(1:OSL%NCS), STAT=IERR)             !> ACHTUNG: kann verbessert werden??? riesig!
   CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICN_TO_ICE',IERR)
   OSL%ICN_TO_ICE = 0
ENDIF

ALLOCATE(OSL%IWL_TO_ICG(NWL, NCPL), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','SL%IWL_TO_ICG',IERR)
OSL%IWL_TO_ICG = 0


!> Allocate mapping arrays for ghost cells
ALLOCATE(OSL%ICG_TO_IWG(NCG), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_IWG',IERR)
OSL%ICG_TO_IWG = 0

ALLOCATE(OSL%ICG_TO_ICO(NCG), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_ICO',IERR)
OSL%ICG_TO_ICO = 0

ALLOCATE(OSL%ICG_TO_ICE(NCG), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_ICE',IERR)
OSL%ICG_TO_ICE = 0

END SUBROUTINE SCARC_SETUP_OMAPPINGS


!> ----------------------------------------------------------------------------------------------------
!> Check divisibility by 2 of a given number of elements (in one grid direction)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_DIVISIBILITY(NN, CDIR, CROUTINE)
INTEGER, INTENT(IN) :: NN
CHARACTER(*) , INTENT(IN) :: CDIR, CROUTINE
IF (MOD(NN,2) /= 0) THEN
   WRITE(SCARC_MESSAGE,'(4A,I8,A)') TRIM(CROUTINE),': Parameter ',CDIR,' = ',NN,' not divisable by 2'
   CALL SHUTDOWN(SCARC_MESSAGE); RETURN
ENDIF
END SUBROUTINE SCARC_CHECK_DIVISIBILITY


!> ----------------------------------------------------------------------------------------------------
!> Set wall cell information on coarse level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)
INTEGER, INTENT(IN)    :: IOR0, NM, NL
INTEGER, INTENT(INOUT) :: IWG
INTEGER:: INBR_FACE, IERR
TYPE (SCARC_FACE_TYPE) , POINTER               :: FC, FF              !> coarse and fine FACE types
TYPE (SCARC_LEVEL_TYPE), POINTER               :: SLC, SLF            !> coarse and fine LEVEL types

!> reference coarse and fine LEVEL type
SLC => SCARC(NM)%LEVEL(NL)
SLF => SCARC(NM)%LEVEL(NL-1)

!> reference coarse and fine FACE type
FC => SCARC(NM)%LEVEL(NL)%FACE(IOR0)
FF => SCARC(NM)%LEVEL(NL-1)%FACE(IOR0)

!> initialize FACE type for coarser mesh
FC%NFM = FF%NFM

FC%BXY = FF%BXY
FC%BXZ = FF%BXZ
FC%BYZ = FF%BYZ

FC%IWG_MARKER = IWG

FC%NUM_NEIGHBORS = FF%NUM_NEIGHBORS
IF (FC%NUM_NEIGHBORS /= 0) THEN
   ALLOCATE(FC%NEIGHBORS(FC%NUM_NEIGHBORS), STAT=IERR)
   CALL ChkMemErr('SCARC_FACE_DIMENSIONS','FC%NEIGHBORS',IERR)
ENDIF
DO INBR_FACE= 1, FC%NUM_NEIGHBORS
   FC%NEIGHBORS(INBR_FACE) = FF%NEIGHBORS(INBR_FACE)
ENDDO

SELECT CASE (ABS(IOR0))

   CASE (1)

      FC%DH => SLC%DXL
      FC%IOFF_MAT = 1

      FC%NFX = 1                                                         !> no extension in x-direction
      FC%NFC = SLC%NX                                                    !> number of cells between opposite faces

      IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) THEN                 !> only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y', 'SCARC_SETUP_FACE')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF

      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z', 'SCARC_SETUP_FACE')     !> number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2

   CASE (2)

      FC%DH => SLC%DYL
      FC%IOFF_MAT = SLC%NX

      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X', 'SCARC_SETUP_FACE')     !> number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2

      FC%NFY = 1                                                         !> no extension in y-direction
      FC%NFC = SLC%NY                                                    !> number of cells between opposite faces

      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z', 'SCARC_SETUP_FACE')     !> number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2

   CASE (3)

      FC%DH => SLC%DZL
      FC%IOFF_MAT = SLC%NX * SLC%NY

      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X', 'SCARC_SETUP_FACE')     !> number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2

      IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) THEN                 !> only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y', 'SCARC_SETUP_FACE')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF

      FC%NFZ = 1                                                         !> no extension in y-direction
      FC%NFC = SLC%NZ                                                    !> number of cells between opposite faces
END SELECT
IF (IOR0 > 0) FC%IOFF_MAT = -FC%IOFF_MAT

FC%NFW = FC%NFX * FC%NFY * FC%NFZ                                        !> get number of wall cells for that face

IWG = IWG + FC%NFW                                                       !> increase global wall cell counter

IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
   WRITE(LU_SCARC,*) '========================== SETUP_FACE_DIMENSIONS:  IOR0=',IOR0
   WRITE(LU_SCARC,*) 'FC%NFW=',FC%NFW
   WRITE(LU_SCARC,*) 'FC%NFX=',FC%NFX
   WRITE(LU_SCARC,*) 'FC%NFY=',FC%NFY
   WRITE(LU_SCARC,*) 'FC%NFZ=',FC%NFZ
   WRITE(LU_SCARC,*) 'FF%IWG_MARKER=',FF%IWG_MARKER
   WRITE(LU_SCARC,*) 'IWG=',IWG
ENDIF

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
TYPE (SCARC_WALL_TYPE) , POINTER, DIMENSION(:) :: WC, WF              !> coarse and fine WALL types
TYPE (SCARC_FACE_TYPE) , POINTER               :: FF, FC              !> coarse and fine FACE types
TYPE (SCARC_LEVEL_TYPE), POINTER               :: SLC, SLF, OSLC      !> coarse and fine LEVEL types

SCARC_ROUTINE = 'SCARC_SETUP_FACE'

!> reference coarse and fine LEVEL type
SLC => SCARC(NM)%LEVEL(NL)
SLF => SCARC(NM)%LEVEL(NL-1)

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
         NX1 = SLC%NX+1
         NX2 = SLC%NX+1
      ENDIF
      NY1 = 1
      NY2 = SLC%NY
      NZ1 = 1
      NZ2 = SLC%NZ
   CASE (2)
      NX1 = 1                                                      !> set dimensions for wall cell counting
      NX2 = SLC%NX
      IF (IOR0 > 0) THEN
         NY1 = 0
         NY2 = 0
      ELSE
         NY1 = SLC%NY+1
         NY2 = SLC%NY+1
      ENDIF
      NZ1 = 1
      NZ2 = SLC%NZ
   CASE (3)
      NX1 = 1                                                      !> set dimensions for wall cell counting
      NX2 = SLC%NX
      NY1 = 1
      NY2 = SLC%NY
      IF (IOR0 > 0) THEN
         NZ1 = 0
         NZ2 = 0
      ELSE
         NZ1 =SLC%NZ+1
         NZ2 =SLC%NZ+1
      ENDIF
END SELECT

IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) '========================== SETUP_FACE:  IOR0=',IOR0
   WRITE(LU_SCARC,*) 'NX1=',NX1
   WRITE(LU_SCARC,*) 'NX2=',NX2
   WRITE(LU_SCARC,*) 'NY1=',NY1
   WRITE(LU_SCARC,*) 'NY2=',NY2
   WRITE(LU_SCARC,*) 'NZ1=',NZ1
   WRITE(LU_SCARC,*) 'NZ2=',NZ2
ENDIF

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
               WC(IWC)%ICW = (IZ-1)*SLC%NX*SLC%NY + (IY-1)*SLC%NX + IX + 1
            CASE (-1)
               WC(IWC)%ICW = (IZ-1)*SLC%NX*SLC%NY + (IY-1)*SLC%NX + IX - 1
            CASE (2)
               WC(IWC)%ICW = (IZ-1)*SLC%NX*SLC%NY +  IY   *SLC%NX + IX
            CASE (-2)
               WC(IWC)%ICW = (IZ-1)*SLC%NX*SLC%NY + (IY-2)*SLC%NX + IX
            CASE (3)
               WC(IWC)%ICW =  IZ   *SLC%NX*SLC%NY + (IY-1)*SLC%NX + IX
            CASE (-3)
               WC(IWC)%ICW = (IZ-2)*SLC%NX*SLC%NY + (IY-1)*SLC%NX + IX
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

         IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
            WRITE(LU_SCARC,*) 'WC(',IWC,')%IOR0=',WC(IWC)%IOR
            WRITE(LU_SCARC,*) 'WC(',IWC,')%ICW =',WC(IWC)%ICW
            WRITE(LU_SCARC,*) 'WC(',IWC,')%IXG =',WC(IWC)%IXG
            WRITE(LU_SCARC,*) 'WC(',IWC,')%IYG =',WC(IWC)%IYG
            WRITE(LU_SCARC,*) 'WC(',IWC,')%IZG =',WC(IWC)%IZG
            WRITE(LU_SCARC,*) 'WC(',IWC,')%IXW =',WC(IWC)%IXW
            WRITE(LU_SCARC,*) 'WC(',IWC,')%IYW =',WC(IWC)%IYW
            WRITE(LU_SCARC,*) 'WC(',IWC,')%IZW =',WC(IWC)%IZW
            WRITE(LU_SCARC,*) 'FF%IWG_MARKER=',FF%IWG_MARKER
         ENDIF

         SELECT CASE (TYPE_DIMENSION)

            !> ------------------------------------------------------------------------------------
            !> 2D-version
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_DIMENSION_TWO)

               !> determine fine IW's, which must be merged to one coarse IW
               SELECT CASE (ABS(IOR0))
                  CASE ( 1)
                     IWF(1) = FF%IWG_MARKER + 2*(IZ-1)
                  CASE ( 2)
                     IWF(1) = FF%IWG_MARKER + 2*(IZ-1)*SLF%NX + 2*(IX - 1)
                  CASE ( 3)
                     IWF(1) = FF%IWG_MARKER + 2*(IX-1)
               END SELECT
               IWF(2) = IWF(1)+1

               IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
                  WRITE(LU_SCARC,*) 'IWF(1)=',IWF(1)
                  WRITE(LU_SCARC,*) 'IWF(2)=',IWF(2)
               ENDIF

               !> set fine cell neighbors (they must be the same for all fine IW's)
               NOMF(1) = WF(IWF(1))%NOM
               NOMF(2) = WF(IWF(2))%NOM
               IF (NOMF(1) /= NOMF(2)) THEN
                  WRITE(SCARC_MESSAGE,'(2A,2I4)') TRIM(SCARC_ROUTINE),': Inconsistent neighbors ',NOMF(1),NOMF(2)
                  CALL SHUTDOWN(SCARC_MESSAGE); RETURN
               ENDIF

               IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
                  WRITE(LU_SCARC,*) 'NOMF(1)=',NOMF(1)
                  WRITE(LU_SCARC,*) 'NOMF(2)=',NOMF(2)
               ENDIF

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

<<<<<<< HEAD
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
                  WRITE(*,*) 'WRONG BOUNDARY TYPE'
               ENDIF

=======
>>>>>>> origin/scarc
               !> in case of an internal boundary set neighboring WALL cells
               IF (NOMF(1) > 0) THEN

                  OSLC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)

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
                           WRITE(SCARC_MESSAGE,'(A,i4)') 'SCARC_SETUP_FACE: Wrong resolutions for IOR0=',IOR0
                           CALL SHUTDOWN(SCARC_MESSAGE); RETURN
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
                           WRITE(SCARC_MESSAGE,'(A,i4)') 'SCARC_SETUP_FACE: Wrong resolutions for IOR0=',IOR0
                           CALL SHUTDOWN(SCARC_MESSAGE); RETURN
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

                  !!!
                  !> Allocate and specify ICN and ICE arrays for OSLC
                  !!!
                  NCPL = (IZ2-IZ1+1)*(IX2-IX1+1)
                  OSLC%NCPL = NCPL

                  CALL SCARC_SETUP_WALLCELL_NEIGHBOR(IX1, IX2, 1, 1, IZ1, IZ2, IWC, IOR0, NM, NOMF(1), NL)

               ENDIF


            !> ------------------------------------------------------------------------------------
            !> 3D-version
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_DIMENSION_THREE)

               !> determine fine IW's, which must be merged to one coarse IW
               SELECT CASE (ABS(IOR0))
                  CASE (1)
                     IWF(1) = FF%IWG_MARKER + (2*IZ-2)*SLF%NY + 2*IY - 2
                     IWF(3) = FF%IWG_MARKER + (2*IZ-1)*SLF%NY + 2*IY - 2
                  CASE (2)
                     IWF(1) = FF%IWG_MARKER + (2*IZ-2)*SLF%NX + 2*IX - 2
                     IWF(3) = FF%IWG_MARKER + (2*IZ-1)*SLF%NX + 2*IX - 2
                  CASE (3)
                     IWF(1) = FF%IWG_MARKER + (2*IY-2)*SLF%NX + 2*IX - 2
                     IWF(3) = FF%IWG_MARKER + (2*IY-1)*SLF%NX + 2*IX - 2
               END SELECT
               IWF(2) = IWF(1)+1
               IWF(4) = IWF(3)+1

               !> set fine cell neighbors (they must be the same for all fine IW's)
               DO I=1,4
                  NOMF(I) = WF(IWF(I))%NOM
               ENDDO

               IF (NOMF(1)/=NOMF(2) .OR. NOMF(1)/=NOMF(3) .OR. NOMF(1)/=NOMF(4)) THEN
                  WRITE(SCARC_MESSAGE,'(2A,I4,A)') &
                     TRIM(SCARC_ROUTINE),': Inconsistent neighbors on IOR=', IOR0,' not allowed!'
                  CALL SHUTDOWN(SCARC_MESSAGE); RETURN
               ENDIF
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

<<<<<<< HEAD
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
                  WRITE(*,*) 'WRONG BOUNDARY TYPE'
               ENDIF

=======
>>>>>>> origin/scarc
               !> in case of an internal boundary set WALL(10:15,IWC)
               IF (NOMF(1) > 0) THEN

                  OSLC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)

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
                           WRITE(SCARC_MESSAGE,'(2A,I4)') TRIM(SCARC_ROUTINE),': Wrong resolutions for IOR=', IOR0
                           CALL SHUTDOWN(SCARC_MESSAGE); RETURN
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
                           WRITE(SCARC_MESSAGE,'(2A,I4)') TRIM(SCARC_ROUTINE),': Wrong resolutions for IOR=', IOR0
                           CALL SHUTDOWN(SCARC_MESSAGE); RETURN
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
                           WRITE(SCARC_MESSAGE,'(2A,I4)') TRIM(SCARC_ROUTINE),': Wrong resolutions for IOR=', IOR0
                           CALL SHUTDOWN(SCARC_MESSAGE); RETURN
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

                  !> Allocate and specify ICN and ICE arrays for OSLC
                  NCPL = (IZ2-IZ1+1)*(IY2-IY1+1)*(IX2-IX1+1)
                  OSLC%NCPL = NCPL

                  CALL SCARC_SETUP_WALLCELL_NEIGHBOR(IX1, IX2, IY1, IY2, IZ1, IZ2, IWC, IOR0, NM, NOMF(1), NL)

               ENDIF

         END SELECT

         IWC = IWC + 1

      ENDDO
   ENDDO
ENDDO

FF%IOFFSET_WALL = FF%IOFFSET_WALL + FF%NFW

END SUBROUTINE SCARC_SETUP_FACE


!> -------------------------------------------------------------------------------------------------
!> Set analogues to I_MIN, IMAX, J_MIN, J_MAX, K_MIN, K_MAX on coarse level (only GMG!)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)
INTEGER, INTENT(IN) :: IREFINE, NM, NOM, NL
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IW
LOGICAL :: FOUND
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, OSL

SL  => SCARC(NM)%LEVEL(NL)
OS  => SCARC(NM)%OSCARC(NOM)
OSL => OS%LEVEL(NL)

IMIN=0
IMAX=MESHES(NOM)%IBAR/IREFINE+1

SELECT CASE (TYPE_DIMENSION)
   CASE (NSCARC_DIMENSION_TWO)
      JMIN=0
      JMAX=2
   CASE (NSCARC_DIMENSION_THREE)
      JMIN=0
      JMAX=MESHES(NOM)%JBAR/IREFINE+1
END SELECT

KMIN=0
KMAX=MESHES(NOM)%KBAR/IREFINE+1

OS%NIC_R = 0
FOUND = .FALSE.

SEARCH_LOOP: DO IW=1, OSL%NCG

   !> neighborship structure already known from finest level
   IF (SL%WALL(IW)%NOM/=NOM) CYCLE SEARCH_LOOP
   OS%NIC_R = OS%NIC_R + 1
   FOUND = .TRUE.

   SELECT CASE (SL%WALL(IW)%IOR)
      CASE ( 1)
         IMIN=MAX(IMIN,SL%WALL(NM)%IXN(1)-1)
      CASE (-1)
         IMAX=MIN(IMAX,SL%WALL(NM)%IXN(2)+1)
      CASE ( 2)
         JMIN=MAX(JMIN,SL%WALL(NM)%IYN(1)-1)
      CASE (-2)
         JMAX=MIN(JMAX,SL%WALL(NM)%IYN(2)+1)
      CASE ( 3)
         KMIN=MAX(KMIN,SL%WALL(NM)%IZN(1)-1)
      CASE (-3)
         KMAX=MIN(KMAX,SL%WALL(NM)%IZN(2)+1)
   END SELECT
ENDDO SEARCH_LOOP

N_EXCHANGES = N_EXCHANGES+1

END SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS


!> ----------------------------------------------------------------------------------------------------
!> Allocate several global structures for data exchange
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NL, NP, IERR=0
INTEGER :: IREFINE, IOFFSET, NC_SUM
INTEGER, DIMENSION(:), ALLOCATABLE :: NC_AUX

!>
!> Allocate counter and displacement vector for global data exchanges
!>
ALLOCATE(COUNTS_SCARC(0:N_MPI_PROCESSES-1), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'COUNTS_SCARC', IERR)

ALLOCATE(DISPLS_SCARC(0:N_MPI_PROCESSES-1), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'COUNTS_SCARC', IERR)

ALLOCATE(NC_AUX(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'NC_AUX', IERR)

COUNTS_SCARC = 0
DO NP = 0, N_MPI_PROCESSES-1
   DO NM = 1, NMESHES
      IF (PROCESS(NM)==NP) COUNTS_SCARC(NP) = COUNTS_SCARC(NP) + 1
   ENDDO
ENDDO

DISPLS_SCARC(0) = 0
DO NP = 1, N_MPI_PROCESSES-1
   DISPLS_SCARC(NP) = COUNTS_SCARC(NP-1) + DISPLS_SCARC(NP-1)
ENDDO

!> Allocate arrays which are used for (global) communciations
ALLOCATE(NCS_OFFSET(NMESHES, NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBALS', 'NCS_OFFSET', IERR)
NCS_OFFSET = 0
<<<<<<< HEAD

ALLOCATE(NCS_GLOBAL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBALS', 'NCS_GLOBAL', IERR)
NCS_GLOBAL = 0

ALLOCATE(NCS_BLOCK(NMESHES, NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBALS', 'NCS_BLOCK', IERR)
NCS_BLOCK = 0

=======

ALLOCATE(NCS_GLOBAL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBALS', 'NCS_GLOBAL', IERR)
NCS_GLOBAL = 0

ALLOCATE(NCS_BLOCK(NMESHES, NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBALS', 'NCS_BLOCK', IERR)
NCS_BLOCK = 0

>>>>>>> origin/scarc
ALLOCATE(NCS_LOCAL(NMESHES, NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBALS', 'NCS_LOCAL', IERR)
NCS_LOCAL = 0

ALLOCATE(SP_LOCAL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBALS', 'SP_LOCAL', IERR)
SP_LOCAL = 0.0_EB

ALLOCATE(SP_GROUP(N_MPI_PROCESSES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBALS', 'SP_GROUP', IERR)
SP_GROUP = 0.0_EB


!> Compute global number of cells for all levels
IREFINE = 0
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

   NCS_GLOBAL(NL)  = 0
   NCS_BLOCK(:,NL) = 0

   IOFFSET = 0
   DO NM = 1, NMESHES

      IF (PROCESS(NM) /= MYID) CYCLE

      NCS_LOCAL(NM, NL) = SCARC(NM)%LEVEL(NL)%NCS
      NCS_BLOCK(NM, NL) = NCS_BLOCK(NM, NL) + SCARC(NM)%LEVEL(NL)%NCS   
      
<<<<<<< HEAD
      IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
      IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
         WRITE(LU_SCARC,*) 'NCS_BLOCK(',NM,',',NL,')==', NCS_BLOCK(NM, NL)
         WRITE(LU_SCARC,*) 'NCS_LOCAL(',NM,',',NL,')==', NCS_LOCAL(NM, NL)
      ENDIF
   ENDDO

   !! Determine global number of cells for all levels
   IF (N_MPI_PROCESSES > 1) THEN

      CALL MPI_ALLREDUCE(NCS_LOCAL(MYID+1,NL), NCS_GLOBAL(NL), 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)

      !WRITE(LU_SCARC,*) 'COUNTS(',MYID,')=',COUNTS_SCARC(MYID)
      !WRITE(LU_SCARC,*) 'COUNTS=',COUNTS_SCARC
      !WRITE(LU_SCARC,*) 'DISPLS=',DISPLS_SCARC

      !CALL MPI_ALLGATHERV(NCS_LOCAL(MYID+1,NL),COUNTS_SCARC(MYID),MPI_INTEGER,NCS_GLOBAL(NL),COUNTS_SCARC,DISPLS_SCARC,&
      !                    MPI_INTEGER,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHERV(NCS_LOCAL(MYID+1, NL),1,MPI_INTEGER,NC_AUX,COUNTS_SCARC,DISPLS_SCARC,&
                          MPI_INTEGER,MPI_COMM_WORLD,IERR)
      NCS_LOCAL(1:NMESHES, NL) = NC_AUX(1:NMESHES)

      NC_SUM  = 0
      DO NM = 2, NMESHES
         NC_SUM  = NC_SUM  + NCS_LOCAL(NM-1, NL)
         NCS_OFFSET(NM, NL)  = NCS_OFFSET(NM, NL)  + NC_SUM
      ENDDO
   ELSE
      NCS_GLOBAL(NL)  = NCS_BLOCK(1,NL)
   ENDIF

   IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
      WRITE(LU_SCARC,*) 'NCS_LOCAL  =', NCS_LOCAL
      WRITE(LU_SCARC,*) 'NCS_BLOCK  =', NCS_BLOCK
      WRITE(LU_SCARC,*) 'NCS_GLOBAL =', NCS_GLOBAL
      WRITE(LU_SCARC,*) 'NCS_OFFSET =', NCS_OFFSET
   ENDIF

ENDDO LEVEL_LOOP

DEALLOCATE(NC_AUX, STAT=IERR)

END SUBROUTINE SCARC_SETUP_GLOBALS



!> ----------------------------------------------------------------------------------------------------
!> Setup system of equation:
!> Define matrix stencils and initialize matrices and boundary conditions on all needed levels
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM
INTEGER :: NM, NL

SCARC_ROUTINE = 'SCARC_SETUP_SYSTEM'

!> In case of AMG: setup sizes for transfer matrices
SELECT CASE (TYPE_MULTIGRID)
   CASE (NSCARC_MULTIGRID_GEOMETRIC)                                !> set sizes for all levels
      DO NL=NLEVEL_MIN, NLEVEL_MAX
         CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NL)
      ENDDO
   CASE DEFAULT                                                     !> set sizes only for finest level
      CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NLEVEL_MIN)
END SELECT

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'SETTING UP SYSTEM: HALLO'

MESHES_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   SELECT_SOLVER: SELECT CASE (TYPE_METHOD)

      !> ---------------------------------------------------------------------------------------------
      !> Krylov method (CG/BICG) as main solver, different preconditioners possible
      !> ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_KRYLOV)

         SELECT_PRECON: SELECT CASE (TYPE_PRECON)

            !> ---------------------------------------------------------------------------------------
            !> in case of multigrid as preconditioner:
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_PRECON_MULTIGRID)

               SELECT_PRECON_MG: SELECT CASE (TYPE_MULTIGRID)

                  !> Geometric multigrid:
                  !>    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
                  CASE (NSCARC_MULTIGRID_GEOMETRIC)

                     DO NL = NLEVEL_MIN, NLEVEL_MAX
                        CALL SCARC_SETUP_MATRIX  (NM, NL)
                        CALL SCARC_SETUP_BOUNDARY(NM, NL)
                     ENDDO

                     IF (TYPE_MKL == NSCARC_MKL_COARSE) THEN
                        CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MAX)
                     ENDIF

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
            !> in case of Pardiso/Cluster_Sparse_Solver as preconditioner: 
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_PRECON_PARDISO, NSCARC_PRECON_CLUSTER)

               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MIN)
#endif


            !> ---------------------------------------------------------------------------------------
            !> in case of one-level preconditioners (JACOBI/SSOR/FFT)
            !>    -  assemble standard n-point-matrix on finest level
            !> ---------------------------------------------------------------------------------------
            CASE DEFAULT

               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

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

               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
                  CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MAX)
               ENDIF

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
         !> ---------------------------------------------------------------------------------------
         !> in case of Pardiso as smoother: setup symmetric matrix information
         !> ---------------------------------------------------------------------------------------
         IF (TYPE_SMOOTH == NSCARC_SMOOTH_PARDISO .OR. TYPE_SMOOTH == NSCARC_SMOOTH_CLUSTER) THEN

            DO NL = NLEVEL_MIN, NLEVEL_MAX
               CALL SCARC_SETUP_MATRIX_MKL(NM, NL)
            ENDDO

         ENDIF
#endif


#ifdef WITH_MKL
      !> ---------------------------------------------------------------------------------------------
      !> MKL-Pardiso method 
      !> ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_MKL)

<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'SETTING UP MKL MATRIX: BEFORE SETUP_MATRIX'

         CALL SCARC_SETUP_MATRIX(NM, NLEVEL_MIN)
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'SETTING UP MKL MATRIX: BEFORE SETUP_BOUNDARY'
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'SETTING UP MKL MATRIX: BEFORE SETUP_MATRIX_MKL'
         CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MIN)
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'SETTING UP MKL MATRIX: AFTER SETUP_MATRIX_MKL'
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'SETTING UP MKL MATRIX'

         CALL SCARC_SETUP_MATRIX(NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MIN)
>>>>>>> origin/scarc
#endif

   END SELECT SELECT_SOLVER

ENDDO MESHES_LOOP


!> ---------------------------------------------------------------------------------------------
!> Only temporarily: Show different matrices in debug file 
!> ---------------------------------------------------------------------------------------------
<<<<<<< HEAD
!CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NLEVEL_MIN, 'SETUP_SYSTEM1', 'MATRIX')
=======
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NLEVEL_MIN, 'SETUP_SYSTEM1', 'MATRIX')
>>>>>>> origin/scarc
!#ifdef WITH_MKL
!IF (TYPE_MKL == NSCARC_MKL_LOCAL .OR. TYPE_MKL == NSCARC_MKL_GLOBAL) THEN
!   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIXS, NLEVEL_MIN, 'SETUP_MKL_SYSTEM', 'MKL MATRIX')
!ENDIF
!#endif


!> ------------------------------------------------------------------------------------------------
!> Exchange matrix stencil information in case of AMG
!> ------------------------------------------------------------------------------------------------
IF (NMESHES>1 .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_ALGEBRAIC) THEN

   !   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NLEVEL_MIN, 'SETUP_SYSTEM', 'MATRIX000')
   !   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NLEVEL_MIN, 'SETUP_SYSTEM', 'WALLINFO')

   !> set sizes for matrix stencils on overlapping parts
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_STENCIL, NLEVEL_MIN)
   CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_STENCIL, NLEVEL_MIN)

   !> exchange matrix entries on overlapping parts
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SYSTEM, NLEVEL_MIN)
   CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_SYSTEM, NLEVEL_MIN)

ENDIF


!> ------------------------------------------------------------------------------------------------
!> Only temporarily: Show all relevant arrays and structures in debug file
!> ------------------------------------------------------------------------------------------------
SELECT CASE (TYPE_MULTIGRID)
   CASE (NSCARC_MULTIGRID_GEOMETRIC)
      DO NL=NLEVEL_MIN, NLEVEL_MAX
         CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NL, 'SETUP_SYSTEM', 'MATRIX000')
      !   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL, 'SETUP_SYSTEM', 'WALLINFO')
      ENDDO
   CASE (NSCARC_MULTIGRID_ALGEBRAIC)
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NLEVEL_MIN, 'SETUP_SYSTEM', 'MATRIX000')
      !CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NLEVEL_MIN, 'SETUP_SYSTEM', 'WALLINFO')
END SELECT
!CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NLEVEL_MIN, 'SETUP_SYSTEM', 'MATRIX000')
!CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_BCINDEX, NLEVEL_MIN, 'SETUP_SYSTEM', 'PRESSURE_BC_INDEX')
!CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_ACELL  , NLEVEL_MIN, 'SETUP_SYSTEM', 'WALL_CELL')
!CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_GCELL  , NLEVEL_MIN, 'SETUP_SYSTEM', 'GHOST_CELL')
!CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_NCELL  , NLEVEL_MIN, 'SETUP_SYSTEM', 'NOM_CELL')

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
<<<<<<< HEAD
INTEGER :: IX, IY, IZ, IW
INTEGER :: ICNUM, ICTYPE, IWNUM, IWTYPE
INTEGER :: IC, ICI, IP, IERR, IOR0
!INTEGER :: IWM1, IWP1, IWM2, IWP2, IWM3, IWP3
!INTEGER :: IBM1, IBP1, IBM2, IBP2, IBM3, IBP3
INTEGER :: NCELLS
INTEGER, POINTER, DIMENSION(:,:,:)   :: CELL_INDEX
INTEGER, POINTER, DIMENSION(:,:)     :: WALL_INDEX
=======
INTEGER :: IX, IY, IZ
INTEGER :: ICNUM, ICTYPE, IWNUM, IWTYPE
INTEGER :: IC, ICI, IP, IERR, IOR0
INTEGER :: IWM1, IWP1, IWM2, IWP2, IWM3, IWP3
INTEGER :: IBM1, IBP1, IBM2, IBP2, IBM3, IBP3
>>>>>>> origin/scarc
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: SM
#endif

SL => SCARC(NM)%LEVEL(NL)
M => MESHES(NM)

ALLOCATE (SL%A_ROW(SL%NCS+10), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%A_ROW', IERR)
SL%A_ROW = 0

ALLOCATE (SL%A(SL%NA+10), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%A', IERR)
SL%A = 0.0_EB

ALLOCATE (SL%A_COL(SL%NA+10), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%A_COL', IERR)
SL%A_COL = 0

IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
   WRITE(LU_SCARC,*) 'TYPE_MKL=',TYPE_MKL
   WRITE(LU_SCARC,*) 'TYPE_MKL_GLOBAL=',NSCARC_MKL_GLOBAL
   WRITE(LU_SCARC,*) 'TYPE_MKL_COARSE=',NSCARC_MKL_COARSE
   WRITE(LU_SCARC,*) 'NL=',NL
ENDIF


#ifdef WITH_MKL
IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX)) THEN

   SM => SCARC(NM)%MKL(NL)

<<<<<<< HEAD
=======
   IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'ALLOCATING AG_COL IN LENGTH ', SL%NA+10
>>>>>>> origin/scarc
   WRITE(*,*) 'ACHTUNG: WIEDER DEALLOKIEREN !!'
   ALLOCATE (SM%AG_COL(SL%NA+10), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SM%AG_COL', IERR)
   SM%AG_COL = 0

ENDIF
#endif

<<<<<<< HEAD
IF (NL == NLEVEL_MIN) THEN
   CELL_INDEX => M%CELL_INDEX
   WALL_INDEX => M%WALL_INDEX
   NCELLS     =  CELL_COUNT(NM)
ELSE
   CELL_INDEX => SL%CELL_INDEX
   WALL_INDEX => SL%WALL_INDEX
   NCELLS     =  SL%CELL_COUNT
ENDIF

IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) '------------------ CALLING SETUP_MATRIX FOR MESH ',NM
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SL%NA=',SL%NA
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SL%NCS=',SL%NCS
=======
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) '------------------ CALLING SETUP_MATRIX FOR MESH ',NM
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SL%NA=',SL%NA
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SL%NC=',SL%NC
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'NL   =',NL
   WRITE(LU_SCARC,*) 'TYPE_MKL   =',TYPE_MKL
   WRITE(LU_SCARC,*) 'TYPE_PRECON=',TYPE_PRECON
   WRITE(LU_SCARC,*) 'SIZE(CCVAR)=',SIZE(SL%CCVAR)
<<<<<<< HEAD
   WRITE(LU_SCARC,*) 'SL%NW=',SL%NW
   WRITE(LU_SCARC,*) 'BOUNDARY_TYPE:'
   IF (SL%NX==8) THEN
      WRITE(LU_SCARC,'(8I6)') (SL%WALL(IW)%BOUNDARY_TYPE, IW=1,SL%NW)
   ELSE
      WRITE(LU_SCARC,'(4I6)') (SL%WALL(IW)%BOUNDARY_TYPE, IW=1,SL%NW)
   ENDIF
   WRITE(LU_SCARC,*) 'WALL_INDEX:'
   WRITE(LU_SCARC,'(I6,A,7I6)') (IC,':', (WALL_INDEX(IC,IOR0), IOR0=-3,3),IC=1,NCELLS)
   WRITE(LU_SCARC,*) 'CELL_INDEX:'
   DO IY=0,SL%NY+1
   IF (SL%NX==8) THEN
   WRITE(LU_SCARC,'(I6,A,10I6)') (IY,':',(CELL_INDEX(IX,IY,IZ),IX=0,SL%NX+1),IZ=SL%NZ+1,0,-1)
   ELSE
   WRITE(LU_SCARC,'(I6,A,6I6)') (IY,':',(CELL_INDEX(IX,IY,IZ),IX=0,SL%NX+1),IZ=SL%NZ+1,0,-1)
   ENDIF
   ENDDO
ENDIF

IWTYPE=-9999
=======
ENDIF
>>>>>>> origin/scarc

!> Compute single matrix entries and corresponding row and column pointers
!> Along internal boundaries use placeholders for the neighboring matrix entries
!> which will be communicated in a following step
SELECT CASE (TYPE_DIMENSION)

   !! --------------------- 2D -----------------
   CASE (NSCARC_DIMENSION_TWO)

      IP  = 1
      IY  = 1
      DO IZ = 1, SL%NZ
         DO IX = 1, SL%NX

<<<<<<< HEAD
            IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. SL%CCVAR(IX, IY, IZ, CGSC)/=IS_GASPHASE) CYCLE

            IWTYPE = 0
            IWNUM  = 0

            IC  = SL%CCVAR(IX, IY, IZ, UNKH)
            ICI = CELL_INDEX(IX, IY, IZ)

            !> ACHTUNG, warum sind hier die IOR-Richtungen vertauscht ????
            !IWP1 = WALL_INDEX(ICI,  1)
            !IWM1 = WALL_INDEX(ICI, -1)

            !IWP3 = WALL_INDEX(ICI,  3)
            !IWM3 = WALL_INDEX(ICI, -3)

            !IF (IWP1 /= 0) THEN
            !   IBP1 = SL%WALL(IWP1)%BOUNDARY_TYPE
            !ELSE
            !   IBP1 = 0
            !ENDIF
            !IF (IWM1 /= 0) THEN
            !   IBM1 = SL%WALL(IWM1)%BOUNDARY_TYPE
            !ELSE
            !   IBM1 = 0
            !ENDIF

            !IF (IWP3 /= 0) THEN
            !   IBP3 = SL%WALL(IWP3)%BOUNDARY_TYPE
            !ELSE
            !   IBP3 = 0
            !ENDIF
            !IF (IWM3 /= 0) THEN
            !   IBM3 = SL%WALL(IWM3)%BOUNDARY_TYPE
            !ELSE
            !   IBM3 = 0
            !ENDIF

            !IWP2 = -99
            !IWM2 = -99

            !IBP2 = -99
            !IBM2 = -99

            !IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
            !  WRITE(LU_SCARC,1001) IX,IY,IZ,IC,ICI,IWP1,IWM1,IWP2,IWM2,IWP3,IWM3,IBP1,IBM1,IBP2,IBM2,IBP3,IBM3
=======
            IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC .AND. SL%CCVAR(IX, IY, IZ, CGSC)/=IS_GASPHASE) CYCLE

            IC = SL%CCVAR(IX, IY, IZ, UNKH)
            ICI = M%CELL_INDEX(IX, IY, IZ)

            !> ACHTUNG, warum sind hier die IOR-Richtungen vertauscht ????
            IWP1 = M%WALL_INDEX(ICI,  1)
            IWM1 = M%WALL_INDEX(ICI, -1)

            IWP3 = M%WALL_INDEX(ICI,  3)
            IWM3 = M%WALL_INDEX(ICI, -3)

            IBP1 = M%WALL(IWP1)%BOUNDARY_TYPE
            IBM1 = M%WALL(IWM1)%BOUNDARY_TYPE

            IBP3 = M%WALL(IWP3)%BOUNDARY_TYPE
            IBM3 = M%WALL(IWM3)%BOUNDARY_TYPE

            IWP2 = -99
            IWM2 = -99

            IBP2 = -99
            IBM2 = -99

            IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
              WRITE(LU_SCARC,1001) IX,IY,IZ,IC,ICI,IWP1,IWM1,IWP2,IWM2,IWP3,IWM3,IBP1,IBM1,IBP2,IBM2,IBP3,IBM3
>>>>>>> origin/scarc
            CALL SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)

            !> IOR = 1
            IOR0 = 1
<<<<<<< HEAD
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
            IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = SL%CCVAR(IX-1, IY, IZ, UNKH)
               ICTYPE = SL%CCVAR(IX-1, IY, IZ, CGSC)
               IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX-1,IY,IZ,IOR0,IC,IWNUM,ICNUM,ICTYPE
=======
            IWNUM  = M%WALL_INDEX(ICI, -IOR0)
            IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
            IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = SL%CCVAR(IX-1, IY, IZ, UNKH)
               ICTYPE = SL%CCVAR(IX-1, IY, IZ, CGSC)
               IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX-1,IY,IZ,IOR0,IC,IWNUM,ICNUM,ICTYPE
>>>>>>> origin/scarc
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
            ENDIF

            !> IOR = -1
            IOR0 = -1
<<<<<<< HEAD
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
            IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = SL%CCVAR(IX+1, IY, IZ, UNKH)
               ICTYPE = SL%CCVAR(IX+1, IY, IZ, CGSC)
               IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX+1,IY,IZ,IOR0,IC,IWNUM,ICNUM,ICTYPE
=======
            IWNUM  = M%WALL_INDEX(ICI, -IOR0)
            IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
            IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = SL%CCVAR(IX+1, IY, IZ, UNKH)
               ICTYPE = SL%CCVAR(IX+1, IY, IZ, CGSC)
               IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX+1,IY,IZ,IOR0,IC,IWNUM,ICNUM,ICTYPE
>>>>>>> origin/scarc
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
            ENDIF

            !> IOR = 3
            IOR0 = 3
<<<<<<< HEAD
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
            IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = SL%CCVAR(IX, IY, IZ-1, UNKH)
               ICTYPE = SL%CCVAR(IX, IY, IZ-1, CGSC)
               IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX,IY,IZ-1,IOR0,IC,IWNUM,ICNUM,ICTYPE
=======
            IWNUM  = M%WALL_INDEX(ICI, -IOR0)
            IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
            IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = SL%CCVAR(IX, IY, IZ-1, UNKH)
               ICTYPE = SL%CCVAR(IX, IY, IZ-1, CGSC)
               IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX,IY,IZ-1,IOR0,IC,IWNUM,ICNUM,ICTYPE
>>>>>>> origin/scarc
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
            ENDIF

            !> IOR = -3
            IOR0 = -3
<<<<<<< HEAD
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
            IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = SL%CCVAR(IX, IY, IZ+1, UNKH)
               ICTYPE = SL%CCVAR(IX, IY, IZ+1, CGSC)
               IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX,IY,IZ+1,IOR0,IC,IWNUM,ICNUM,ICTYPE
               CALL SCARC_SETUP_MATRIX_SUBDIAG  (IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
            ENDIF
=======
            IWNUM  = M%WALL_INDEX(ICI, -IOR0)
            IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
            IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = SL%CCVAR(IX, IY, IZ+1, UNKH)
               ICTYPE = SL%CCVAR(IX, IY, IZ+1, CGSC)
               IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX,IY,IZ+1,IOR0,IC,IWNUM,ICNUM,ICTYPE
               CALL SCARC_SETUP_MATRIX_SUBDIAG  (IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
            ENDIF

>>>>>>> origin/scarc

         ENDDO
      ENDDO

   !! --------------------- 3D -----------------
   CASE (NSCARC_DIMENSION_THREE)

IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'SETTING UP MATRIX FOR UNSTRUCTURED'
<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) 'CCVAR:'
   IF (SL%NX==16) THEN
      DO IZ=SL%NZ,0,-1
         WRITE(LU_SCARC,'(18i5)') ((SL%CCVAR(IX,IY,IZ,UNKH),IX=0,SL%NX+1),IY=SL%NY+1,0,-1)
         WRITE(LU_SCARC,*)
      ENDDO
   ELSE
      DO IZ=SL%NZ,0,-1
         WRITE(LU_SCARC,'(10i5)') ((SL%CCVAR(IX,IY,IZ,UNKH),IX=0,SL%NX+1),IY=SL%NY+1,0,-1)
         WRITE(LU_SCARC,*)
      ENDDO
   ENDIF
   WRITE(LU_SCARC,*) 'CELL_INDEX:'
   IF (SL%NX==16) THEN
      DO IZ=SL%NZ,0,-1
         WRITE(LU_SCARC,'(18i5)') ((CELL_INDEX(IX,IY,IZ),IX=0,SL%NX+1),IY=SL%NY+1,0,-1)
         WRITE(LU_SCARC,*)
      ENDDO
   ELSE
      DO IZ=SL%NZ,0,-1
         WRITE(LU_SCARC,'(10i5)') ((CELL_INDEX(IX,IY,IZ),IX=0,SL%NX+1),IY=SL%NY+1,0,-1)
         WRITE(LU_SCARC,*)
      ENDDO
   ENDIF
ENDIF
=======
>>>>>>> origin/scarc

      IP  = 1
      DO IZ = 1, SL%NZ
         DO IY = 1, SL%NY
            DO IX = 1, SL%NX

<<<<<<< HEAD
               IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. SL%CCVAR(IX, IY, IZ, CGSC)/=IS_GASPHASE) CYCLE

               IWTYPE = 0
               IWNUM  = 0

               IC  = SL%CCVAR(IX, IY, IZ, UNKH)
               ICI = CELL_INDEX(IX, IY, IZ)

               !IWP1 = WALL_INDEX(ICI,  1)
               !IWM1 = WALL_INDEX(ICI, -1)

               !IWP2 = WALL_INDEX(ICI,  2)
               !IWM2 = WALL_INDEX(ICI, -2)

               !IWP3 = WALL_INDEX(ICI,  3)
               !IWM3 = WALL_INDEX(ICI, -3)

               !IBP1 = SL%WALL(IWP1)%BOUNDARY_TYPE
               !IBM1 = SL%WALL(IWM1)%BOUNDARY_TYPE

               !IBP2 = SL%WALL(IWP2)%BOUNDARY_TYPE
               !IBM2 = SL%WALL(IWM2)%BOUNDARY_TYPE

               !IBP3 = SL%WALL(IWP3)%BOUNDARY_TYPE
               !IBM3 = SL%WALL(IWM3)%BOUNDARY_TYPE

               !IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
               !  WRITE(LU_SCARC,1001) IX,IY,IZ,IC,ICI,IWP1,IWM1,IWP2,IWM2,IWP3,IWM3,IBP1,IBM1,IBP2,IBM2,IBP3,IBM3
=======
               IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC .AND. SL%CCVAR(IX, IY, IZ, CGSC)/=IS_GASPHASE) CYCLE

               IC  = SL%CCVAR(IX, IY, IZ, UNKH)
               ICI = M%CELL_INDEX(IX, IY, IZ)
>>>>>>> origin/scarc

               IWP1 = M%WALL_INDEX(ICI,  1)
               IWM1 = M%WALL_INDEX(ICI, -1)

               IWP2 = M%WALL_INDEX(ICI,  2)
               IWM2 = M%WALL_INDEX(ICI, -2)

               IWP3 = M%WALL_INDEX(ICI,  3)
               IWM3 = M%WALL_INDEX(ICI, -3)

               IBP1 = M%WALL(IWP1)%BOUNDARY_TYPE
               IBM1 = M%WALL(IWM1)%BOUNDARY_TYPE

               IBP2 = M%WALL(IWP2)%BOUNDARY_TYPE
               IBM2 = M%WALL(IWM2)%BOUNDARY_TYPE

               IBP3 = M%WALL(IWP3)%BOUNDARY_TYPE
               IBM3 = M%WALL(IWM3)%BOUNDARY_TYPE

               IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) &
                 WRITE(LU_SCARC,1001) IX,IY,IZ,IC,ICI,IWP1,IWM1,IWP2,IWM2,IWP3,IWM3,IBP1,IBM1,IBP2,IBM2,IBP3,IBM3
               CALL SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)

<<<<<<< HEAD
               !> IOR = 1
               IOR0 = 1
               IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
               IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
=======
!! ACHTUNG, WARUM IST BEI WALL_INDEX DIE IOR-RICHTUNG VERTAUSCHT ???

               !> IOR = 1
               IOR0 = 1
               IWNUM  = M%WALL_INDEX(ICI, -IOR0)
               IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
>>>>>>> origin/scarc
                  ICNUM  = SL%CCVAR(IX-1, IY, IZ, UNKH)
                  ICTYPE = SL%CCVAR(IX-1, IY, IZ, CGSC)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX-1,IY,IZ,IOR0,IC,IWNUM,ICNUM,ICTYPE
                  CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
               ENDIF

               !> IOR = -1
               IOR0 = -1
<<<<<<< HEAD
               IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
               IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
=======
               IWNUM  = M%WALL_INDEX(ICI, -IOR0)
               IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
>>>>>>> origin/scarc
                  ICNUM  = SL%CCVAR(IX+1, IY, IZ, UNKH)
                  ICTYPE = SL%CCVAR(IX+1, IY, IZ, CGSC)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX+1,IY,IZ,IOR0,IC,IWNUM,ICNUM,ICTYPE
                  CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
               ENDIF

               !> IOR = 2
               IOR0 = 2
<<<<<<< HEAD
               IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
               IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
=======
               IWNUM  = M%WALL_INDEX(ICI, -IOR0)
               IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
>>>>>>> origin/scarc
                  ICNUM  = SL%CCVAR(IX, IY-1, IZ, UNKH)
                  ICTYPE = SL%CCVAR(IX, IY-1, IZ, CGSC)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX,IY-1,IZ,IOR0,IC,IWNUM,ICNUM,ICTYPE
                  CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IY, ICTYPE, ICNUM, NM, NL)
               ENDIF

               !> IOR = -2
               IOR0 = -2
<<<<<<< HEAD
               IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
               IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
=======
               IWNUM  = M%WALL_INDEX(ICI, -IOR0)
               IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
>>>>>>> origin/scarc
                  ICNUM  = SL%CCVAR(IX, IY+1, IZ, UNKH)
                  ICTYPE = SL%CCVAR(IX, IY+1, IZ, CGSC)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX,IY+1,IZ,IOR0,IC,IWNUM,ICNUM,ICTYPE
                  CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IY, ICTYPE, ICNUM, NM, NL)
               ENDIF

               !> IOR = 3
               IOR0 = 3
<<<<<<< HEAD
               IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
               IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
=======
               IWNUM  = M%WALL_INDEX(ICI, -IOR0)
               IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
>>>>>>> origin/scarc
                  ICNUM  = SL%CCVAR(IX, IY, IZ-1, UNKH)
                  ICTYPE = SL%CCVAR(IX, IY, IZ-1, CGSC)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX,IY,IZ-1,IOR0,IC,IWNUM,ICNUM,ICTYPE
                  CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
               ENDIF

               !> IOR = -3
               IOR0 = -3
<<<<<<< HEAD
               IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
               IF (IWNUM /= 0) IWTYPE = SL%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
=======
               IWNUM  = M%WALL_INDEX(ICI, -IOR0)
               IWTYPE = M%WALL(IWNUM)%BOUNDARY_TYPE
               IF (IWNUM==0.OR.IWTYPE==INTERPOLATED_BOUNDARY) THEN
>>>>>>> origin/scarc
                  ICNUM  = SL%CCVAR(IX, IY, IZ+1, UNKH)
                  ICTYPE = SL%CCVAR(IX, IY, IZ+1, CGSC)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,1002) IX,IY,IZ+1,IOR0,IC,IWNUM,ICNUM,ICTYPE
                  CALL SCARC_SETUP_MATRIX_SUBDIAG  (IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
               ENDIF

            ENDDO
         ENDDO
      ENDDO

END SELECT
SL%A_ROW(SL%NCS+1) = IP
SL%NA              = IP -1                                   !> set correct number of matrix entries
SL%NAS             = IP -1                                   !> set correct number of matrix entries


IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
   WRITE(LU_SCARC,*) 'IN SETUP_MATRIX with FULL ENTRIES'
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SL%NC =',SL%NC
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SL%NCS=',SL%NCS
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SL%NA=',SL%NA
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SL%NAS=',SL%NAS
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SIZE(SL%A)=',SIZE(SL%A)
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SIZE(SL%A_ROW)=',SIZE(SL%A_ROW)
   WRITE(LU_SCARC,*) 'SETUP_MATRIX: SIZE(SL%A_COL)=',SIZE(SL%A_COL)
ENDIF

1001 FORMAT('SETTING UP ENTRY FOR DIAG   : IX,IY,IZ=',3I4,': IC  =',I4,': ICI=',I4,': IWs=',6I4,': IBs=',6I4)
1002 FORMAT('SETTING UP ENTRY FOR SUBDIAG: IX,IY,IZ=',3I4,': IOR = ',I4,' : IC=',I4,&
            ': IWNUM=',I4,': ICNUM=',I4,': ICTYPE=',I4)
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
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: SM
#endif

SL => SCARC(NM)%LEVEL(NL)

SL%A(IP) = SL%A(IP) - 2.0_EB/(SL%DXL(IX-1)*SL%DXL(IX))

IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) &
   SL%A(IP) = SL%A(IP) - 2.0_EB/(SL%DYL(IY-1)*SL%DYL(IY))

SL%A(IP) = SL%A(IP) - 2.0_EB/(SL%DZL(IZ-1)*SL%DZL(IZ))

SL%A_ROW(IC) = IP
SL%A_COL(IP) = IC
      
#ifdef WITH_MKL
IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX)) THEN
<<<<<<< HEAD
   SM => SCARC(NM)%MKL(NL)
   SM%AG_COL(IP) = SL%A_COL(IP) + NCS_OFFSET(NM, NL)
   IF (TYPE_DEBUG > NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(A,5I4,E14.6)') &
       'MAINDIAG:', NL,IC, IP, IC, SL%A_COL(IP)+SL%NUNKH_LOC(NM), SL%A(IP)
ENDIF
=======
  SM => SCARC(NM)%MKL(NL)
  SM%AG_COL(IP) = SL%A_COL(IP) + NCS_OFFSET(NM, NL)
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,'(A,5I4,E14.6)') &
    'MAINDIAG:', NL,IC, IP, IC, SL%A_COL(IP)+SL%NUNKH_LOC(NM), SL%A(IP)
      ENDIF
>>>>>>> origin/scarc
#endif
      
IP = IP + 1

END SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG


SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG (IOR0, IC, IP, IPTR, NBR_TYPE, NBR_NUM, NM, NL)
INTEGER, INTENT(IN)  :: IC, IOR0, IPTR, NBR_TYPE, NBR_NUM
INTEGER, INTENT(IN)  :: NM, NL
INTEGER, INTENT(INOUT) :: IP
INTEGER  :: IW, IX, IY, IZ
REAL(EB) :: DSCAL, DH1, DH2
LOGICAL  :: BINTERNAL
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL
TYPE (MESH_TYPE)       , POINTER :: M

#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: SM
#endif

SL => SCARC(NM)%LEVEL(NL)
M  => MESHES(NM)

DH1 = SL%FACE(IOR0)%DH(IPTR-1)
DH2 = SL%FACE(IOR0)%DH(IPTR)
   
<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_EXTREME) THEN
   WRITE(LU_SCARC,*) '------------- SETUP_MATRIX_SUBDIAG:', IOR0
   WRITE(LU_SCARC,*) 'IOR0=',IOR0,': IC=',IC,': IP=',IP,': IPTR=',IPTR,': NBR_TYPE=',NBR_TYPE
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,'(a,i3)') '------------- SETUP_MATRIX_SUBDIAG:', IOR0
   WRITE(LU_SCARC,'(5(a,i4))') 'IOR0=',IOR0,': IC=',IC,': IP=',IP,': IPTR=',IPTR,': NBR_TYPE=',NBR_TYPE
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'DH1 = ', DH1
   WRITE(LU_SCARC,*) 'DH2 = ', DH2
   WRITE(LU_SCARC,*) 'NBR_TYPE = ', NBR_TYPE
   WRITE(LU_SCARC,*) 'IPTR     = ', IPTR    
   WRITE(LU_SCARC,*) 'SL%FACE(',IOR0,')%NUM_NEIGHBORS=',SL%FACE(IOR0)%NUM_NEIGHBORS
   WRITE(LU_SCARC,*) 'SL%FACE(',IOR0,')%NFC          =',SL%FACE(IOR0)%NFC
<<<<<<< HEAD
   WRITE(LU_SCARC,*) 'PRES_ON_WHOLE_DOMAIN = ', PRES_ON_WHOLE_DOMAIN
=======
   WRITE(LU_SCARC,*) 'PRES_ON_WHOLE_DOMAIN_SCRC = ', PRES_ON_WHOLE_DOMAIN_SCRC
>>>>>>> origin/scarc
ENDIF
   
!> set bounds and step sizes depending on orientation of face (lower or upper subdiagonal)
IF (IOR0 > 0) THEN
   BINTERNAL = IPTR > 1
   DSCAL = 2.0_EB/(DH1*(DH1+DH2))
ELSE
   BINTERNAL = IPTR < SL%FACE(IOR0)%NFC
   DSCAL = 2.0_EB/(DH2*(DH1+DH2))
ENDIF

<<<<<<< HEAD
!> if IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
IF (BINTERNAL) THEN

   IF (PRES_ON_WHOLE_DOMAIN .OR. NBR_TYPE == IS_GASPHASE) THEN
   
      SL%A(IP)     = SL%A(IP) + DSCAL
      SL%A_COL(IP) = NBR_NUM
IF (TYPE_DEBUG > NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,*) 'INTERNAL: IC=',IC,': A_COL(',IP,') =',SL%A_COL(IP)
         
#ifdef WITH_MKL
      IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
          (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX)) THEN
         SM => SCARC(NM)%MKL(NL)
         SM%AG_COL(IP)= SL%A_COL(IP) + NCS_OFFSET(NM, NL)
IF (TYPE_DEBUG > NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,*) 'INTERNAL: IC=',IC,': AG_COL(',IP,')=',SM%AG_COL(IP)
      ENDIF
      IP = IP + 1

=======
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'BINTERNAL: ',BINTERNAL

!> if IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
IF (BINTERNAL) THEN

   IF (PRES_ON_WHOLE_DOMAIN_SCRC .OR. NBR_TYPE == IS_GASPHASE) THEN
   
      SL%A(IP)     = SL%A(IP) + DSCAL
      SL%A_COL(IP) = NBR_NUM
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'INTERNAL: IC=',IC,': A_COL(',IP,') =',SL%A_COL(IP)
         
#ifdef WITH_MKL
      IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
          (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX)) THEN
         SM => SCARC(NM)%MKL(NL)
         SM%AG_COL(IP)= SL%A_COL(IP) + NCS_OFFSET(NM, NL)
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'INTERNAL: IC=',IC,': AG_COL(',IP,')=',SM%AG_COL(IP)
      ENDIF
      IP = IP + 1

>>>>>>> origin/scarc
   ENDIF
#endif
         
!> if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
ELSE IF (SL%FACE(IOR0)%NUM_NEIGHBORS /= 0) THEN
   
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'EXTERNAL: HAS NEIGHBOR??? '
   IF (HAS_NEIGHBOR(SL%CCVAR, SL%WALL, IC, IW, SL%FACE(IOR0)%IWG_MARKER, SL%FACE(IOR0)%NFW)) THEN
      SL%A(IP)     = SL%A(IP) + DSCAL
      SL%A_COL(IP) = SL%WALL(IW)%ICE(1)
      !SL%A_COL(IP) = SL%WALL(IW)%ICO
<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,*) 'EXTERNAL: IC=',IC,': A_COL(',IP,') =',SL%A_COL(IP)
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'EXTERNAL: IC=',IC,': A_COL(',IP,') =',SL%A_COL(IP)
>>>>>>> origin/scarc
         
#ifdef WITH_MKL
      IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
          (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX)) THEN
         SM => SCARC(NM)%MKL(NL)
         IX = SL%WALL(IW)%IXG
         IY = SL%WALL(IW)%IYG
         IZ = SL%WALL(IW)%IZG
         SM%AG_COL(IP) = SL%CCVAR(IX, IY, IZ, UNKH) + NCS_OFFSET(SL%WALL(IW)%NOM, NL)
<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,*) 'EXTERNAL: IC=',IC,': AG_COL(',IP,')=',SM%AG_COL(IP)
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'EXTERNAL: IC=',IC,': AG_COL(',IP,')=',SM%AG_COL(IP)
>>>>>>> origin/scarc
      ENDIF
#endif
         
      IP = IP + 1
   ENDIF
   
ENDIF

<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_EXTREME) THEN
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
   WRITE(LU_SCARC,'(a, 10f12.6)') 'DH=',SL%FACE(IOR0)%DH
   WRITE(LU_SCARC,*) 'DSCAL = ', DSCAL
   WRITE(LU_SCARC,*) 'IOR0 = ', IOR0
   WRITE(LU_SCARC,*) 'SL%FACE(IOR0)%NUM_NEIGHBORS=',SL%FACE(IOR0)%NUM_NEIGHBORS
   WRITE(LU_SCARC,*) 'SL%FACE(IOR0)%IWG_MARKER=',SL%FACE(IOR0)%IWG_MARKER
   WRITE(LU_SCARC,*) 'SL%FACE(IOR0)%NFW=',SL%FACE(IOR0)%NFW
   WRITE(LU_SCARC,*) 'SL%FACE(IOR0)%NFC=',SL%FACE(IOR0)%NFC
<<<<<<< HEAD
=======
   WRITE(LU_SCARC,'(a,i4,a,f12.6)') 'A(',IP-1,')    =',SL%A(IP-1)
   WRITE(LU_SCARC,'(a,i4,a,i4)') 'A_ROW(',IC,')=',SL%A_ROW(IC)
   WRITE(LU_SCARC,'(a,i4,a,i4)') 'A_COL(',IP-1,')=',SL%A_COL(IP-1)
>>>>>>> origin/scarc
#ifdef WITH_MKL
   IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
       (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX)) THEN
      SM => SCARC(NM)%MKL(NL)
      WRITE(LU_SCARC,'(a,i4,a,i4)') 'AG_COL(',IP-1,')=',SM%AG_COL(IP-1)
   ENDIF
#endif
ENDIF
END SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG




!> ------------------------------------------------------------------------------------------------
!> Determine if cell IC has a neighbor and, if yes, return corresponding IW-value
!> ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION HAS_NEIGHBOR(CCVAR0, WALL0, IC, IW, IW0, IL0)
TYPE (SCARC_WALL_TYPE), DIMENSION(:), INTENT(IN):: WALL0
INTEGER, DIMENSION(0:,0:,0:,:), INTENT(IN):: CCVAR0
INTEGER, INTENT(IN)    :: IC, IW0, IL0
INTEGER :: IXW0, IYW0, IZW0
INTEGER :: IXG0, IYG0, IZG0
INTEGER :: IC0
INTEGER, INTENT(INOUT) :: IW

HAS_NEIGHBOR = .FALSE.
<<<<<<< HEAD
SEARCH_WALL_CELLS_LOOP: DO IW = IW0, IW0+IL0-1
=======
SEARCH_WALLCELLS_LOOP: DO IW = IW0, IW0+IL0-1
>>>>>>> origin/scarc

  IXW0 = WALL0(IW)%IXW
  IYW0 = WALL0(IW)%IYW
  IZW0 = WALL0(IW)%IZW

  IXG0 = WALL0(IW)%IXG
  IYG0 = WALL0(IW)%IYG
  IZG0 = WALL0(IW)%IZG

  IC0 = CCVAR0(IXW0, IYW0, IZW0, UNKH) 
IF (TYPE_DEBUG>NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'HAS_NEIGHBOR? IC, IW, IX0, IY0, IZ0, IXG0, IYG0, IZG0, IC0:', &
    IC, IW, IXW0, IYW0, IZW0, IXG0, IYG0, IZG0, IC0

  IF (IC == IC0 .AND. (WALL0(IW)%NOM /= 0.AND.CCVAR0(IXG0, IYG0, IZG0, CGSC)/=IS_SOLID)) THEN
     HAS_NEIGHBOR = .TRUE.
IF (TYPE_DEBUG>NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'HAS_NEIGHBOR      ==============> YES'
<<<<<<< HEAD
     EXIT SEARCH_WALL_CELLS_LOOP
=======
     EXIT SEARCH_WALLCELLS_LOOP
>>>>>>> origin/scarc
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
INTEGER :: ICOL, JCOL, IAS, IERR, IP
INTEGER :: ISYM, JSYM, NSYM
REAL(EB) :: VAL, VALS, DIFF
LOGICAL  :: BSYM
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL
TYPE (SCARC_MKL_TYPE) , POINTER :: SM
INTEGER, DIMENSION(:), ALLOCATABLE :: ICOL_AUX, JC_AUX

SL => SCARC(NM)%LEVEL(NL)
SM => SCARC(NM)%MKL(NL)

!> ------------------------------------------------------------------------------------------------
!> Only temporarily : Debug matrix
!> ------------------------------------------------------------------------------------------------
<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'Checking SYMMETRY'
   WRITE(LU_SCARC,*) '---------------------- A_ROW:', SL%NCS
   WRITE(LU_SCARC,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NCS+1)
   WRITE(LU_SCARC,*) '---------------------- A_COL:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%A_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- A_GCOL:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SM%AG_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- A:', SL%A_ROW(IC+1)- SL%A_ROW(IC)
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SL%A(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Store only symmetric parts of matrix (diagonal and upper part)
!> ------------------------------------------------------------------------------------------------
IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN

   !> First check whether symmetry of system matrix is guaranteed
   DO IC = 1, SL%NCS
   
      COLUMN_LOOP: DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
         ICS = SL%A_COL(ICOL)
         VAL = SL%A(ICOL)
         IF (ICS > IC .AND. ICS <= SL%NCS) THEN
            BSYM = .FALSE.
            DO JCOL = SL%A_ROW(ICS)+1, SL%A_ROW(ICS+1)-1
               JCS = SL%A_COL(JCOL)
               IF (JCS == IC) THEN
                  VALS = SL%A(JCOL)
                  DIFF = ABS(VAL-VALS)
                  !WRITE(LU_SCARC,*) IC, ICS,': DIFF=',DIFF
                  IF (ABS(VAL - VALS) < 1E-10) THEN
                     BSYM=.TRUE.
<<<<<<< HEAD
                     IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
                     IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
>>>>>>> origin/scarc
                        WRITE(LU_SCARC,'(a,i4,a,i4,a,i4,a,i4,a,f12.6)') 'A(',IC,',',ICS,')= A(',ICS,',',JCS,')=',VAL
                     ENDIF
                     CYCLE COLUMN_LOOP
                  ENDIF
               ENDIF
            ENDDO
            IF (.NOT.BSYM) THEN
<<<<<<< HEAD
               IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
               IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
>>>>>>> origin/scarc
                  WRITE(LU_SCARC,'(a,i4,a,i4,a,f25.16,a,i4,a,i4,a,f25.16)')  &
                                 'ERROR: A(',IC,',',ICS,')=',VAL,'/= A(',ICS,',',JCS,')=',VALS
               ENDIF
               WRITE(SCARC_MESSAGE,'(3A)') TRIM('SCARC_CHECK_SYMMETRY'),': System matrix not symmetric'
               CALL SHUTDOWN(SCARC_MESSAGE); RETURN
            ENDIF
         ENDIF
      ENDDO COLUMN_LOOP
   
   ENDDO

   !> Compute number of entries in symmetric matrix
   SL%NAS = 0
   DO IC = 1, SL%NCS
      DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
         IF (TYPE_MKL == NSCARC_MKL_LOCAL) THEN
            JC = SL%A_COL(ICOL)
            IF (JC >= IC .AND. JC <= SL%NCS) SL%NAS = SL%NAS+1   !really neccesary ???
         ELSE
            JC = SM%AG_COL(ICOL)
            IF (JC >= IC + NCS_OFFSET(NM,NL)) SL%NAS = SL%NAS+1
<<<<<<< HEAD
            IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
=======
            IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
>>>>>>> origin/scarc
              WRITE(LU_SCARC,*) 'COUNTING NAS: IC=',IC,', ICOL=',ICOL,', JC=',JC,', NAS=',SL%NAS
         ENDIF
      ENDDO
   ENDDO

!> Non-symmetric case: use complete matrix (not only symmetric part) of matrix
ELSE
   SL%NAS = SL%NA
ENDIF

IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
   WRITE(LU_SCARC,*) 'Setting up MKL-matrix for level ', NL
   WRITE(LU_SCARC,*) 'TYPE_MKL=',TYPE_MKL
   WRITE(LU_SCARC,*) 'SL%NCS=',SL%NCS
   WRITE(LU_SCARC,*) 'SL%NA=',SL%NA
   WRITE(LU_SCARC,*) 'SL%NAS=',SL%NAS
ENDIF


!> ------------------------------------------------------------------------------------------------
!> allocate storage for symmetric matrix and its column and row pointers
!> ------------------------------------------------------------------------------------------------

ALLOCATE (SM%ASYM(SL%NAS), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SM%ASYM', IERR)
SM%ASYM = 0.0_EB

ALLOCATE (SM%ASYM_COL(SL%NAS), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SM%ASYM_COL', IERR)
SM%ASYM_COL = 0

ALLOCATE (SM%ASYM_ROW(SL%NCS+1), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SM%ASYM_ROW', IERR)
SM%ASYM_ROW = 0


IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX)) THEN
   ALLOCATE (ICOL_AUX(SL%NPOINTS), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'ICOL_AUX', IERR)
   ALLOCATE (JC_AUX(SL%NPOINTS), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'JC_AUX', IERR)
ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
   WRITE(LU_SCARC,*) 'SYMMETRIC MATRIX:'
   WRITE(LU_SCARC,*) 'NCS_LOCAL=',NCS_LOCAL
ENDIF

!> ------------------------------------------------------------------------------------------------
!> extract symmetric matrix part from usual system matrix
!> ------------------------------------------------------------------------------------------------
IAS = 1
DO IC = 1, SL%NCS
   SM%ASYM_ROW(IC) = IAS

   !> blockwise use of MKL solver
   IF (TYPE_MKL == NSCARC_MKL_LOCAL) THEN    !really neccesary ??

      DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
         JC = SL%A_COL(ICOL)

            IF (JC >= IC .AND. JC <= SL%NCS) THEN
               SM%ASYM_COL(IAS) = SL%A_COL(ICOL) 
               SM%ASYM(IAS)     = SL%A(ICOL)
<<<<<<< HEAD
               IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
=======
               IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
>>>>>>> origin/scarc
                  WRITE(LU_SCARC,'(a,i6,a,i6,a,i6,a,i6,a,i6,a,f12.6)') 'IC:',IC,' JC=',JC,' IAS:',IAS,' ROW:',&
                        SM%ASYM_ROW(IC),' COL:',SM%ASYM_COL(IAS),' AS:',SM%ASYM(IAS)
               ENDIF
               IAS = IAS + 1
            ENDIF
      ENDDO

   !> global use of MKL solver
   ELSE

      !> store indices of all diagonal and upper-diagonal entries
      ICOL_AUX = 0
      JC_AUX   = 99999999
      ISYM = 1
      JC0 = SM%AG_COL(SL%A_ROW(IC))
      DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
         !JC = SL%A_COL(ICOL)
         JC = SM%AG_COL(ICOL)
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
      IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
        WRITE(LU_SCARC,*) 'NSYM = ', NSYM
        WRITE(LU_SCARC,*) 'ICOL_AUX=', ICOL_AUX
        WRITE(LU_SCARC,*) 'JC_AUX=', JC_AUX
      ENDIF
   
      !> sort them in increasing order (for the use of DSS-PARDISO functionality)
      JSYM = 1
      SORT_LOOP: DO WHILE (JSYM <= NSYM)
         DO ISYM = 1, NSYM
            JC = JC_AUX(ISYM)
            IF (JC == 99999999) CYCLE
            IF (JC <= MINVAL(JC_AUX)) THEN
               ICOL = ICOL_AUX(ISYM)
               SM%ASYM(IAS) = SL%A(ICOL)
               !SM%ASYM_COL(IAS) = SL%A_COL(ICOL)
               SM%ASYM_COL(IAS) = SM%AG_COL(ICOL)
<<<<<<< HEAD
               IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
               IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
                  WRITE(LU_SCARC,'(a,i6,a,i6,a,i6,a,i6,a,i6,a,f12.6)') 'IC:',IC,' JC=',JC,' IAS:',IAS,' ROW:',&
                        SM%ASYM_ROW(IC),' COL:',SM%ASYM_COL(IAS),' AS:',SM%ASYM(IAS)
               ENDIF
               JC_AUX(ISYM) = 99999999            ! mark entry as already used
               IAS  = IAS  + 1
            ENDIF
         ENDDO
         JSYM = JSYM + 1
      ENDDO SORT_LOOP
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
        WRITE(LU_SCARC,*) 'NSYM = ', NSYM
        WRITE(LU_SCARC,*) 'ICOL_AUX=', ICOL_AUX
      ENDIF

   ENDIF

ENDDO

SM%ASYM_ROW(SL%NCS+1) = IAS
!WRITE(LU_SCARC,*) 'SIZE(ASYM_ROW)=',SIZE(SM%ASYM_ROW), SM%ASYM_ROW(SL%NCS)
!WRITE(LU_SCARC,*) 'SIZE(ASYM_ROW)=',SIZE(SM%ASYM_ROW), SM%ASYM_ROW(SL%NCS+1)
!WRITE(LU_SCARC,'(A,17I4)') 'ASYM_ROW=',SM%ASYM_ROW
SM%ASYM_COL = SM%ASYM_COL 

IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX)) THEN
   DEALLOCATE (ICOL_AUX, STAT=IERR)
   DEALLOCATE (JC_AUX, STAT=IERR)
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Only temporarily - debug information
!> ------------------------------------------------------------------------------------------------
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
   WRITE(LU_SCARC,*) 'SYMMETRIC MATRIX:'
   WRITE(LU_SCARC,*) 'NCS_LOCAL=',NCS_LOCAL
   WRITE(LU_SCARC,*) '---------------------- A_ROW:', SL%NCS
   WRITE(LU_SCARC,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NCS+1)
   WRITE(LU_SCARC,*) '---------------------- A_COL:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%A_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- AG_COL:'
   IF (TYPE_MKL == NSCARC_MKL_GLOBAL) THEN
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SM%AG_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
   ENDIF
   WRITE(LU_SCARC,*) '---------------------- A:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SL%A(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- ASYM_ROW:'
   WRITE(LU_SCARC,'(4i9)') (SM%ASYM_ROW(IC), IC=1,SL%NCS+1)
   WRITE(LU_SCARC,*) '---------------------- ASYM_COL:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SM%ASYM_COL(IP),IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- AS:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SM%ASYM(IP),IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1)
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_MKL
#endif


!> ----------------------------------------------------------------------------------------------------
!> Extract overlapping matrix data
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OVERLAPS (NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM, LAST_NOM
INTEGER :: IERR, IW, IWG, IWC, IG, IC, JC, JCA, ICW, ICN, ICE, ICE0, ICE2, ICG, ICG2
INTEGER :: IX, IY, IZ
INTEGER :: ICOL, ICCE, ICPL, NCOL
INTEGER :: ICC, ICOLG, ICOLE, IOR0, NCPL, NCE0
TYPE (SCARC_TYPE)         , POINTER :: S
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SLF, SLC
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSLF, OSLC

SCARC_ROUTINE = 'SCARC_SETUP_OVERLAPS'

IERR = 0
SELECT_TYPE: SELECT CASE (NTYPE)

   !! --------------------------------------------------------------------------------------------
   !! set stencil sizes on overlapping parts
   !! --------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_STENCIL)

      !WRITE(LU_SCARC,*) 'MATRIX_STENCIL:'

      STENCIL_MESHES_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         S => SCARC(NM)
         SLF => S%LEVEL(NL)

         STENCIL_EXTERNAL_CELL_LOOP: DO ICE = SLF%NCS+1, SLF%NCE

           IW = ABS(SLF%ICE_TO_IWG(ICE))

           IF (IW > 0) THEN
              ICG = SLF%WALL(IW)%ICG(1)
              !WRITE(LU_SCARC,*) 'A: ICE=',ICE,': IW=',IW,': ICG=',ICG
           ELSE
              IW  = ABS(IW)
              ICG = SLF%WALL(ABS(IW))%ICG2(1)
              !WRITE(LU_SCARC,*) 'B: ICE=',ICE,': IW=',IW,': ICG=',ICG
           ENDIF

           NOM = SLF%WALL(IW)%NOM
           OSLF => S%OSCARC(NOM)%LEVEL(NL)

           SLF%A_ROW(ICE+1) = SLF%A_ROW(ICE) + OSLF%A_SIZE(ICG)
           !WRITE(LU_SCARC,*) 'NOM=',NOM,': SLF%A_ROW(',ICE+1,')=',SLF%A_ROW(ICE+1)

         ENDDO STENCIL_EXTERNAL_CELL_LOOP

      ENDDO STENCIL_MESHES_LOOP

   !! --------------------------------------------------------------------------------------------
   !! after exchange of system matrix
   !! --------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SYSTEM)

      SYSTEM_MESHES_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         S => SCARC(NM)
         SLF => S%LEVEL(NL)

         NCE0 = SLF%NCE+1
         DO ICE = SLF%NCS+1, SLF%NCE

            IC = ICE
            IWG = SLF%ICE_TO_IWG(ICE)
            IX  = SLF%WALL(IWG)%IXG
            IY  = SLF%WALL(IWG)%IYG
            IZ  = SLF%WALL(IWG)%IZG
            ICN = SLF%CCVAR(IX, IY, IZ, UNKH)

            NOM = SLF%WALL(IWG)%NOM
            OSLF => S%OSCARC(NOM)%LEVEL(NL)

            IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
               WRITE(LU_SCARC,'(10(a,i4))') 'ICE=',ICE,': IWG=',IWG,': ICN=',ICN,': NOM=',NOM
               WRITE(LU_SCARC,'(a,10i4)') '============= OVERLAPS 0: A_ROW:',&
               (SLF%A_COL(ICOL), ICOL = SLF%A_ROW(ICE), SLF%A_ROW(ICE+1)-1)
               WRITE(LU_SCARC,'(a,10i4)') '============= OVERLAPS 0: A_COL:',&
               (SLF%A_COL(ICOL), ICOL = SLF%A_ROW(ICE), SLF%A_ROW(ICE+1)-1)
            ENDIF

            ICPL = 0
            DO ICOL = SLF%A_ROW(ICE),SLF%A_ROW(ICE+1)-1
               JC  = SLF%A_COL(ICOL)
               JCA = ABS(SLF%A_COL(ICOL))
               !WRITE(LU_SCARC,'(10(a,i4))') 'ICOL=',ICOL,': JC=',JC,': JCA=',JCA,': NOM=',NOM
WRITE(*,*) 'ACHTUNG HIER WEGEN ICN SCHAUEN!'
               IF (JC < 0 .AND. JCA <= OSLF%NC) THEN
                  !IF (OSLF%ICN_TO_ICE(JCA)==0) THEN
                  IF (JC==-1) THEN
                     OSLF%ICN_TO_ICE(JCA) = NCE0
                     SLF%CELL_MAP(2, ICE) = NCE0
                     NCE0 = NCE0 + 1
                  ELSE
                     SLF%CELL_MAP(2, ICE) = OSLF%ICN_TO_ICE(JCA)
                  ENDIF
                  SLF%CELL_MAP(3, ICE) = JC
                  SLF%A_COL(ICOL) = SLF%CELL_MAP(2,ICE)

                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,3i3,a,i3)') &
                  !'OVERLAPS A: ICOL=',ICOL,': ICE=',ICE,': NOM =',NOM,': ICN_TO_ICE(',JCA,')=',&
                  !OSLF%ICN_TO_ICE(JCA),': CELL_MAP=',&
                  !SLF%CELL_MAP(1,IC),SLF%CELL_MAP(2,IC),SLF%CELL_MAP(2,IC), ': NCE0=',NCE0

               ELSE
                  SLF%A_COL(ICOL) = JC

                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,2i3)') &
                  !'OVERLAPS B: ICOL=',ICOL,': ICE=',ICE,': NOM =',NOM,&
                  !': A_COL(',ICOL,') =',SLF%A_COL(ICOL)

               ENDIF
            ENDDO

            !WRITE(LU_SCARC,'(a,10i4)') '============= OVERLAPS C: A_COL:',&
            ! (SLF%A_COL(ICOL), ICOL = SLF%A_ROW(IC), SLF%A_ROW(IC+1)-1)
         ENDDO

         !SLF%NCE2 = NCE0-1              !> needed where ????? OUTDATED ???

      ENDDO SYSTEM_MESHES_LOOP

   !! --------------------------------------------------------------------------------------------
   !! after exchange of prolongation matrix
   !! --------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_PROLONGATION)

      PROLONGATION_MESHES_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         S => SCARC(NM)
         SLF => S%LEVEL(NL)
         SLC => S%LEVEL(NL+1)

         IG   = 0
         IWC  = 0
         ICCE = SLF%NCC
         PROLONGATION_WALL_LOOP1: DO IW = 1, SLF%NW

            NOM  = SLF%WALL(IW)%NOM
            IF (NOM == 0) CYCLE PROLONGATION_WALL_LOOP1

            OSLF => S%OSCARC(NOM)%LEVEL(NL)
            OSLC => S%OSCARC(NOM)%LEVEL(NL+1)

            IOR0 = SLF%WALL(IW)%IOR
            ICW  = SLF%WALL(IW)%ICW
            ICC  = SLF%CELL_TYPE(ICW)

            ICE   = SLF%WALL(IW)%ICE(1)
            ICG   = SLF%WALL(IW)%ICG(1)
            NCOL  = OSLF%P_ROW(ICG +1)-OSLF%P_ROW(ICG)         !> wozu wird das nochmal gebraucht ??

            IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
               WRITE(LU_SCARC,*) '============= SETUP_OVERLAPS: IW=',IW
               WRITE(LU_SCARC,*) 'ICW  =',ICW
               WRITE(LU_SCARC,*) 'ICC  =',ICC
               WRITE(LU_SCARC,*) 'ICE  =',ICE
               WRITE(LU_SCARC,*) 'ICG  =',ICG
               WRITE(LU_SCARC,*) 'NCOL =',NCOL
            ENDIF

            IF (ICC >= NSCARC_CELL_TYPE_COARSE) THEN
               IWC = IWC + 1
               !SLC%INTERNAL_BDRY_CELL(ICC) = NSCARC_LAYER_ONE
               SLC%WALL_INDEX(ICC, IOR0)   = IWC
               SLC%WALL(IWC)%NOM = NOM
               SLC%WALL(IWC)%ICW = ICC
               SLC%WALL(IWC)%IOR = IOR0

               IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
                  WRITE(LU_SCARC,*) 'A: SLC%INTERNAL_BDRY_CELL(',ICC,')=',SLC%INTERNAL_BDRY_CELL(ICC)
                  WRITE(LU_SCARC,*) 'A: SLC%WALL_INDEX(',ICC,',', IOR0,')=',SLC%WALL_INDEX(ICC, IOR0)
                  WRITE(LU_SCARC,*) 'A: SLC%WALL(',IWC,')%NOM=',SLC%WALL(IWC)%NOM
                  WRITE(LU_SCARC,*) 'A: SLC%WALL(',IWC,')%ICW=',SLC%WALL(IWC)%ICW
                  WRITE(LU_SCARC,*) 'A: SLC%WALL(',IWC,')%IOR=',SLC%WALL(IWC)%IOR
               ENDIF

            ENDIF

            IF (SLF%CELL_TYPE(ICE) >= NSCARC_CELL_TYPE_COARSE) THEN
               ICCE = ICCE + 1
               SLF%CELL_TYPE(ICE) = ICCE
               !SLC%EXT_PTR(ICCE) = OSLF%CELL_TYPE(ICG)
               ICOL = OSLF%P_ROW(ICG)
               IG  = IG + 1
               !OSLF%ICG0 = OSLF%ICG0 + 1
               OSLF%ICG0 = IG
               OSLC%ICN_TO_ICE(OSLF%P_COL(ICOL))   = ICCE

               IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
                  WRITE(LU_SCARC,*) 'B: ICCE=',ICCE
                  WRITE(LU_SCARC,*) 'B: ICOL=',ICOL
                  WRITE(LU_SCARC,*) 'B: OSLF%ICG0=',OSLF%ICG0
                  WRITE(LU_SCARC,*) 'B: SLF%CELL_TYPE(',ICE,')=',SLF%CELL_TYPE(ICE)
                  WRITE(LU_SCARC,*) 'B: SLC%EXT_PTR(',ICCE,')=',OSLF%CELL_TYPE(ICG)
                  WRITE(LU_SCARC,*) 'B: OSLC%ICN_TO_ICE(',OSLF%P_COL(ICOL),')=',OSLC%ICN_TO_ICE(OSLF%P_COL(ICOL))
                  WRITE(LU_SCARC,*) 'OSLC%ICG_TO_ICE(',OSLF%ICG0,')=',OSLC%ICG_TO_ICE(OSLF%ICG0)
               ENDIF

               OSLC%ICG_TO_ICE(OSLF%ICG0) = ICCE

            ENDIF

            IF (NL/=NLEVEL_MIN .OR. TYPE_LAYER/=NSCARC_LAYER_TWO) CYCLE PROLONGATION_WALL_LOOP1

            IF (SLF%CELL_TYPE(ICE2) >= NSCARC_CELL_TYPE_COARSE) THEN
               SLF%CELL_TYPE(ICE2) = ICCE
               SLC%EXT_PTR(ICCE)  = OSLF%CELL_TYPE(ICG2)
               ICOL = OSLF%P_ROW(ICG2)
               IG  = IG + 1
               OSLF%ICG0 = OSLF%ICG0 + 1
               OSLC%ICN_TO_ICE(OSLF%P_COL(ICOL)) = ICCE
               ICCE = ICCE + 1
            ENDIF

         ENDDO PROLONGATION_WALL_LOOP1

         IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
            DO IC = 1, SLC%NC
               WRITE(LU_SCARC,*) 'SLC%INTERNAL_BDRY_CELL(',IC,')=',SLC%INTERNAL_BDRY_CELL(IC)
            ENDDO
            DO IC = 1, SLC%NC
               WRITE(LU_SCARC,'(a,i3,a,8i4)') 'SLC%WALL_INDEX(',IC,', -3:3)=',SLC%WALL_INDEX(IC,-3:3)
            ENDDO
         ENDIF

         !> Replace negative cell numbers of neighboring cells for internal wall cells
         !> Note that this only holds true for fine grid cells
         !> (weights of coarse cells don't reach in neighboring mesh)
         IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) THEN

         PROLONGATION_WALL_LOOP2: DO IW = 1, SLF%NW

            NOM  = SLF%WALL(IW)%NOM
            IF (NOM == 0) CYCLE PROLONGATION_WALL_LOOP2

            OSLF => S%OSCARC(NOM)%LEVEL(NL)

            ICW  = SLF%WALL(IW)%ICW
            ICC  = SLF%CELL_TYPE(ICW)

            !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3)') &
            !     'B0: ============= : IW=',IW,': NOM=',NOM,': ICW=',ICW,': ICC=',ICC

            IF (ICC < NSCARC_CELL_TYPE_COARSE) THEN

               DO ICOL = SLF%P_ROW(ICW), SLF%P_ROW(ICW+1)-1
                  JC   = SLF%P_COL(ICOL)

                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3)') &
                  !   'B1: ============= : ICC=',ICC,': JC=',JC

                  !> Additionally identify coarse cells from second layer
                  !> adjacent to internal boundary
WRITE(*,*) 'ACHTUNG HIER WEGEN ICN SCHAUEN'
                  IF (JC < 0) THEN
                     IC = OSLF%ICN_TO_ICE(ABS(JC))
                     SLF%P_COL(ICOL) = SLF%CELL_TYPE(IC)

                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
                  !   'B2: ============= : IC=',IC,': SLF%P_COL(',ICOL,')=',SLF%P_COL(ICOL)

                  ELSE IF (JC <= SLF%NCC) THEN
                     IF (SLC%INTERNAL_BDRY_CELL(JC) > NSCARC_LAYER_ONE) THEN
                        IWC = IWC + 1
                        SLC%INTERNAL_BDRY_CELL(JC) = NSCARC_LAYER_TWO
                        SLC%WALL_INDEX(ICC, IOR0)  = -IWC
                        SLC%WALL(IWC)%NOM = NOM
                        SLC%WALL(IWC)%ICW = JC
                        SLC%WALL(IWC)%IOR = SLF%WALL(IW)%IOR

                        !WRITE(LU_SCARC,*) 'B3: IWC=',IWC
                        !WRITE(LU_SCARC,*) 'B3: SLC%INTERNAL_BDRY_CELL(',JC,')=',SLC%INTERNAL_BDRY_CELL(JC)
                        !WRITE(LU_SCARC,*) 'B3: SLC%WALL(',IWC,')%NOM=',NOM
                        !WRITE(LU_SCARC,*) 'B3: SLC%WALL(',IWC,')%ICW=',JC
                        !WRITE(LU_SCARC,*) 'B3: SLC%WALL(',IWC,')%IOR=',SLF%WALL(IWC)%IOR

                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO PROLONGATION_WALL_LOOP2

   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL 001')
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION', 'FINAL CELL_TYPE 001')

         PROLONGATION_ECELL_LOOP: DO ICE0 = SLF%NCS+1, SLF%NCE

            IW  = SLF%ICE_TO_IWG(ICE0)
            NOM = SLF%WALL(ABS(IW))%NOM
            OSLF => S%OSCARC(NOM)%LEVEL(NL)
            OSLC => S%OSCARC(NOM)%LEVEL(NL+1)

            !WRITE(LU_SCARC,*) 'C: =======================  IW=',IW, ': NOM=',NOM
            !WRITE(LU_SCARC,'(a,5i4)') 'OSLF%P_ROW=',(OSLF%P_ROW(ICG), ICG=1,5)
            !WRITE(LU_SCARC,'(a,10i4)') 'OSLF%P_COL=',(OSLF%P_COL(ICG), ICG=1,10)

            !> Replace positive cell numbers of neighboring cells in first layer

            !WRITE(*,*) 'G: HIER WALL%ICN �ndern auf OSLF!'
            !WRITE(*,*) 'G: MUSS NOCH UEBERARBEITET WERDEN !!!'

            IF (IW > 0) THEN
               ICG = SLF%WALL(IW)%ICG(1)
               ICE = SLF%WALL(IW)%ICE(1)
               ICOLE = SLF%P_ROW(ICE)

               !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4)') 'C00: ICG=',ICG,': ICE=',ICE,': ICOLE=',ICOLE

               DO ICOLG = OSLF%P_ROW(ICG), OSLF%P_ROW(ICG+1)-1
                  SLF%P_COL(ICOLE) = OSLF%P_COL(ICOLG)
                  SLF%P_COL(ICOLE) = OSLF%ICN_TO_ICE(OSLF%P_COL(ICOLG))
                  SLF%P(ICOLE) = OSLF%P(ICOLG)
                  ICOLE = ICOLE + 1

                  !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4,a,i4)') 'C01: ICOLG=',ICOLG,': JC=',JC,': OSLF%NCC=',OSLF%NCC

!>                 IF (JC > 0 .AND. JC <= OSLF%NCC) THEN
!>                    IC = OSLC%ICN_TO_ICE(JC)
!>                    !OSLF%P_COL(ICOLG) = IC
!>                    SLF%P_COL(ICOLE)  = IC
!>                    !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4,a,i4)') &
!>                    !'C10: OSLF%P_COL(',ICOLG,')=',IC,': SLF%P_COL(',ICOLE,')=',SLF%P_COL(ICOLE)
!>                 ELSE
!>                    !OSLF%P_COL(ICOLG) = SLF%CELL_TYPE(ABS(JC))
!>                    SLF%P_COL(ICOLE) = SLF%CELL_TYPE(ABS(JC))
!>                    !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4,a,i4)') &
!>                    !'C20: OSLF%P_COL(',ICOLG,')=',IC,': SLF%P_COL(',ICOLE,')=',SLF%P_COL(ICOLE)
!>                 ENDIF
!>                 SLF%P(ICOLE) = OSLF%P(ICOLG)
!>                 ICOLE = ICOLE + 1
               ENDDO
               SLF%P_ROW(ICE+1) = ICOLE

            !> Replace positive cell numbers of neighboring cells in second layer (if requested)
            ELSE
               IW = ABS(IW)

               WRITE(*,*) 'D: HIER WALL%ICN �ndern auf OSLF!'
               WRITE(*,*) 'D: MUSS NOCH UEBERARBEITET WERDEN !!'

               IF (NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) THEN
                  !ICG2 = SLF%WALL(IW)%ICG2(1)
                  !ICE2 = SLF%WALL(IW)%ICE2(1)
                  ICOLE = SLF%P_ROW(ICE2)
                  DO ICOLG = OSLF%P_ROW(ICG2), OSLF%P_ROW(ICG2+1)-1
                     JC   = OSLF%P_COL(ICOLG)
                     IF (JC > 0 .AND. JC <= OSLF%NCC) THEN
                        IC = OSLC%ICN_TO_ICE(JC)
                        IF (IC > 0) THEN
                           OSLF%P_COL(ICOLG) = IC
                           SLF%P_COL(ICOLE) = IC
                        ELSE
                           OSLF%P_COL(ICOLG) = -OSLF%P_COL(ICOLG)
                           SLF%P_COL(ICOLE)  =  OSLF%P_COL(ICOLG)
                        ENDIF
                     ELSE
                        OSLF%P_COL(ICOLG) = SLF%CELL_TYPE(ABS(JC))
                        SLF%P_COL(ICOLE) = SLF%CELL_TYPE(ABS(JC))
                     ENDIF
                     SLF%P(ICOLE) = OSLF%P(ICOLG)
                     ICOLE = ICOLE + 1
                  ENDDO
                  SLF%P_ROW(ICE2+1) = ICOLE
               ENDIF
            ENDIF

         ENDDO PROLONGATION_ECELL_LOOP

         SLC%NW = IWC
         !WRITE(LU_SCARC,*) 'SLC%NW=',SLC%NW

         ENDIF     !GMG-ENDIF

      ENDDO PROLONGATION_MESHES_LOOP

 !! --------------------------------------------------------------------------------------------
 !! after exchange of GMG-matrix
 !! --------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_GMG)

      GMG_MESHES_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         S => SCARC(NM)
         SLF => S%LEVEL(NL)
         SLC => S%LEVEL(NL+1)

         IWC  = 0
         ICCE = SLF%NCC
         LAST_NOM = 0
         GMG_WALL_LOOP1: DO IW = 1, SLF%NW

            NOM  = SLF%WALL(IW)%NOM
            IF (NOM == 0) CYCLE GMG_WALL_LOOP1
            IF (NOM /= LAST_NOM) IG = 0
            LAST_NOM = NOM

            !WRITE(LU_SCARC,*) '============ IW=',IW,': NOM=',NOM, '; IG=',IG

            OSLF => S%OSCARC(NOM)%LEVEL(NL)
            OSLC => S%OSCARC(NOM)%LEVEL(NL+1)

            IOR0 = SLF%WALL(IW)%IOR
            ICW  = SLF%WALL(IW)%ICW
            ICC  = SLF%CELL_TYPE(ICW)

            WRITE(*,*) 'I: HIER WALL%ICN �ndern auf OSLF!'
            WRITE(*,*) 'I: MUSS NOCH UEBERARBEITET WERDEN !!'

            !ICE   = SLF%WALL(IW)%ICE(1)
            !ICG   = SLF%WALL(IW)%ICG(1)
            NCOL  = OSLF%P_ROW(ICG +1)-OSLF%P_ROW(ICG)

            !WRITE(LU_SCARC,*) 'SETUP_OVERLAPS:'
            !WRITE(LU_SCARC,*) 'ICW  =',ICW
            !WRITE(LU_SCARC,*) 'ICC  =',ICC
            !WRITE(LU_SCARC,*) 'ICE  =',ICE
            !WRITE(LU_SCARC,*) 'ICG  =',ICG
            !WRITE(LU_SCARC,*) 'NCOL =',NCOL

            IF (ICC >= NSCARC_CELL_TYPE_COARSE) THEN

               IWC  = IWC + 1
               ICCE = ICCE + 1
               NCPL = 1

               SLC%INTERNAL_BDRY_CELL(ICC) = NSCARC_LAYER_ONE
               SLC%WALL_INDEX(ICC, IOR0)   = IWC

               WRITE(*,*) 'J: HIER WALL%ICN �ndern auf OSLF!'
               WRITE(*,*) 'J: MUSS NOCH UEBERARBEITET WERDEN !!'

               SLC%WALL(IWC)%NOM  = NOM
               !SLC%WALL(IWC)%ICW  = ICC
               SLC%WALL(IWC)%IOR  = IOR0
               !SLC%WALL(IWC)%NCPL = NCPL

               !WRITE(LU_SCARC,*) 'SLC%WALL(',IWC,')%NOM=',SLC%WALL(IWC)%NOM
               !WRITE(LU_SCARC,*) 'SLC%WALL(',IWC,')%ICW=',SLC%WALL(IWC)%ICW
               !WRITE(LU_SCARC,*) 'SLC%WALL(',IWC,')%IOR=',SLC%WALL(IWC)%IOR , SLC%INTERNAL_BDRY_CELL(ICC)

               !IF (SLF%CELL_TYPE(ICE) < NSCARC_CELL_TYPE_COARSE) THEN
               !   WRITE(SCARC_MESSAGE,'(2A)') TRIM(SCARC_ROUTINE),': Wrong cominbation of coarse cells in OVERLAP_GMG, stop!'
               !   CALL SHUTDOWN(SCARC_MESSAGE); RETURN
               !NDIF

               WRITE(*,*) 'H: HIER WALL%ICN �ndern auf OSLF!'
               WRITE(*,*) 'H: MUSS NOCH UEBERARBEITET WERDEN !!'

               !ALLOCATE(SLC%WALL(IWC)%ICG(NCPL), STAT=IERR)
               !CALL ChkMemErr('SCARC_SETUP_OVERLAP_GMG','ICG',IERR)

               !ALLOCATE(SLC%WALL(IWC)%ICE(NCPL), STAT=IERR)
               !CALL ChkMemErr('SCARC_SETUP_SETUP_OVERLAP_GMG','ICE',IERR)

               !ALLOCATE(SLC%WALL(IWC)%ICN(NCPL), STAT=IERR)
               !CALL ChkMemErr('SCARC_SETUP_SETUP_OVERLAP_GMG','ICN',IERR)

               SLF%CELL_TYPE(ICE) = ICCE
               SLC%EXT_PTR(ICCE) = OSLF%CELL_TYPE(ICG)
               ICOL = OSLF%P_ROW(ICG)

               IG  = IG + 1
               OSLC%ICN_TO_ICE(OSLF%P_COL(ICOL)) = ICCE
               OSLC%ICG_TO_ICE(IG)             = ICCE

               !SLC%WALL(IWC)%ICE(1) = ICCE
               !SLC%WALL(IWC)%ICG(1) = IG
               !SLC%WALL(IWC)%ICN(1) = OSLF%CELL_TYPE(ICG)

               !WRITE(LU_SCARC,*) 'ICCE=',ICCE
               !WRITE(LU_SCARC,*) 'SLF%CELL_TYPE(',ICE,')=',SLF%CELL_TYPE(ICE)
               !WRITE(LU_SCARC,*) 'SLC%EXT_PTR(',ICCE,')=',SLC%EXT_PTR(ICCE)
               !WRITE(LU_SCARC,*) 'OSLC%ICN_TO_ICE(',OSLF%P_COL(ICOL),')=',OSLC%ICN_TO_ICE(OSLF%P_COL(ICOL))
               !!WRITE(LU_SCARC,*) 'SLC%WALL(',IWC,')%ICE(1)=',SLC%WALL(IWC)%ICE(1)
               !!WRITE(LU_SCARC,*) 'SLC%WALL(',IWC,')%ICG(1)=',SLC%WALL(IWC)%ICG(1)
               !!WRITE(LU_SCARC,*) 'SLC%WALL(',IWC,')%ICN(1)=',SLC%WALL(IWC)%ICN(1)
            ENDIF

         ENDDO GMG_WALL_LOOP1

         !DO IC = 1, SLC%NC
         !   WRITE(LU_SCARC,*) 'SLC%INTERNAL_BDRY_CELL(',IC,')=',SLC%INTERNAL_BDRY_CELL(IC)
         !ENDDO
         !DO IC = 1, SLC%NC
         !   WRITE(LU_SCARC,'(a,i3,a,8i4)') 'SLC%WALL_INDEX(',IC,', -3:3)=',SLC%WALL_INDEX(IC,-3:3)
         !ENDDO

      ENDDO GMG_MESHES_LOOP

END SELECT SELECT_TYPE

END SUBROUTINE SCARC_SETUP_OVERLAPS


!> ------------------------------------------------------------------------------------------------
!> Set pointer for different structures on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP, NW, ITYPE
REAL(EB) :: DBC, A_OLD
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL


SL => SCARC(NM)%LEVEL(NL)

IF (PRES_ON_WHOLE_DOMAIN) THEN
   NW = SL%N_EXTERNAL_WALL_CELLS
ELSE
   NW = SL%N_EXTERNAL_WALL_CELLS + SL%N_INTERNAL_WALL_CELLS
ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) '==============================================================='
   WRitE(LU_SCARC,*) ' SETUP BOUNDARY '
   WRITE(LU_SCARC,*) '==============================================================='
   WRITE(LU_SCARC,*) 'NM = ', NM
   WRITE(LU_SCARC,*) 'NL = ', NL
   WRITE(LU_SCARC,*) 'NW = ', NW
   WRITE(LU_SCARC,*) 'NCS= ', SL%NCS
   WRITE(LU_SCARC,*) 'NA = ', SL%NA
   WRITE(LU_SCARC,*) 'N_EXTERNAL_WALL_CELLS= ', SL%N_EXTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) 'N_INTERNAL_WALL_CELLS= ', SL%N_INTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) '---------------------- WALL(.)%BTYPE:'
   IF (SL%NX==8) THEN
   WRITE(LU_SCARC,'(8i9)') (SL%WALL(IW)%BTYPE, IW=1,SL%NW)
   ELSE
   WRITE(LU_SCARC,'(4i9)') (SL%WALL(IW)%BTYPE, IW=1,SL%NW)
   ENDIF
   WRITE(LU_SCARC,*) '---------------------- WALL(.)%BOUNDARY_TYPE:'
   IF (SL%NX==8) THEN
   WRITE(LU_SCARC,'(8i9)') (SL%WALL(IW)%BOUNDARY_TYPE, IW=1,SL%NW)
   ELSE
   WRITE(LU_SCARC,'(4i9)') (SL%WALL(IW)%BOUNDARY_TYPE, IW=1,SL%NW)
   ENDIF
   WRITE(LU_SCARC,*) '---------------------- CCVAR(...,CGSC):'
   IF (SL%NX==8) THEN
   WRITE(LU_SCARC,'(10i9)') ((SL%CCVAR(I,1,K,CGSC), I=0,SL%NX+1), K=0,SL%NZ+1)
   ELSE
   WRITE(LU_SCARC,'(6i9)') ((SL%CCVAR(I,1,K,CGSC), I=0,SL%NX+1), K=0,SL%NZ+1)
   ENDIF
   WRITE(LU_SCARC,*) '---------------------- CCVAR(...,UNKH):'
   IF (SL%NX==8) THEN
   WRITE(LU_SCARC,'(10i9)') ((SL%CCVAR(I,1,K,UNKH), I=0,SL%NX+1), K=0,SL%NZ+1)
   ELSE
   WRITE(LU_SCARC,'(6i9)') ((SL%CCVAR(I,1,K,CGSC), I=0,SL%NX+1), K=0,SL%NZ+1)
   ENDIF
   WRITE(LU_SCARC,*) '---------------------- A_ROW:', SL%NCS, SIZE(SL%A_ROW)
   WRITE(LU_SCARC,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NCS+1)
   WRITE(LU_SCARC,*) '---------------------- A_COL:', SIZE(SL%A_COL)
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%A_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- A:', SL%A_ROW(IC+1)- SL%A_ROW(IC), SIZE(SL%A)
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SL%A(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
ENDIF

IF (PRES_ON_WHOLE_DOMAIN_SCRC) THEN
   NW = SL%N_EXTERNAL_WALL_CELLS
ELSE
   NW = SL%N_EXTERNAL_WALL_CELLS + SL%N_INTERNAL_WALL_CELLS
ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) '==============================================================='
   WRitE(LU_SCARC,*) ' SETUP BOUNDARY '
   WRITE(LU_SCARC,*) '==============================================================='
   WRITE(LU_SCARC,*) 'NM = ', NM
   WRITE(LU_SCARC,*) 'NL = ', NL
   WRITE(LU_SCARC,*) 'NW = ', NW
   WRITE(LU_SCARC,*) 'NCS= ', SL%NCS
   WRITE(LU_SCARC,*) 'NA = ', SL%NA
   WRITE(LU_SCARC,*) 'N_EXTERNAL_WALL_CELLS= ', SL%N_EXTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) 'N_INTERNAL_WALL_CELLS= ', SL%N_INTERNAL_WALL_CELLS
   WRITE(LU_SCARC,*) '---------------------- WALL(.)%BTYPE:'
   WRITE(LU_SCARC,'(6i9)') (SL%WALL(IW)%BTYPE, IW=1,SL%NW)
   WRITE(LU_SCARC,*) '---------------------- CCVAR:'
   WRITE(LU_SCARC,'(6i9)') ((SL%CCVAR(I,1,K,CGSC), I=0,SL%NX+1), K=0,SL%NZ+1)
   WRITE(LU_SCARC,*) '---------------------- A_ROW:', SL%NCS, SIZE(SL%A_ROW)
   WRITE(LU_SCARC,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NCS+1)
   WRITE(LU_SCARC,*) '---------------------- A_COL:', SIZE(SL%A_COL)
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%A_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- A:', SL%A_ROW(IC+1)- SL%A_ROW(IC), SIZE(SL%A)
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SL%A(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
ENDIF

SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)

   !! ------------------ 2D -----------------------
   CASE (NSCARC_DIMENSION_TWO)

<<<<<<< HEAD
      WALL_CELLS_LOOP2D: DO IW = 1, NW
=======
      WALLCELLS_LOOP2D: DO IW = 1, NW
>>>>>>> origin/scarc

         IOR0 = SL%WALL(IW)%IOR
         IF (ABS(IOR0) == 2) CYCLE          !> 2D: cycle boundaries in y-direction

         I    = SL%WALL(IW)%IXW
         J    = 1
         K    = SL%WALL(IW)%IZW

         NOM  = SL%WALL(IW)%NOM

         ITYPE = SL%CCVAR(I,J,K,CGSC)

<<<<<<< HEAD
         IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.ITYPE/=IS_GASPHASE) CYCLE
=======
         IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC.AND.ITYPE/=IS_GASPHASE) CYCLE
>>>>>>> origin/scarc

         SL%WALL(IW)%ICW = SL%CCVAR(I, J, K, UNKH)
         IC = SL%CCVAR(I, J, K, UNKH)

<<<<<<< HEAD
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS)  WRITE(LU_SCARC,'(A,8i4)') &
   'BOUNDARY, PROCESSING: IW, I, J, K, NOM, IOR: ', IW, I, J, K, NOM, IOR0, ITYPE, IC

>>>>>>> origin/scarc
         SELECT CASE (ABS(IOR0))
            CASE (1)
               DBC= SL%DXI2
            CASE (3)
               DBC= SL%DZI2
         END SELECT

         IP = SL%A_ROW(IC)
         A_OLD = SL%A(IP)
         SELECT CASE (SL%WALL(IW)%BTYPE)
            CASE (DIRICHLET)                      !> set Dirichlet BC's along open boundary cells
               SL%A(IP) = SL%A(IP) - DBC
<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
    WRITE(LU_SCARC,'(a,i5,a,i5,a,3f15.6,i5)') 'IW=',IW,': DIRICHLET: SL%A(',IP,')=',SL%A(IP), A_OLD, DBC, IC
            !CASE (INTERNAL)                      !> do nothing along internal boundaries (only debugging)
            CASE (NEUMANN)                        !> set Neumann BC's at all other nodes
               SL%A(IP) = SL%A(IP) + DBC
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
    WRITE(LU_SCARC,'(a,i5,a,i5,a,3f15.6,i5)') 'IW=',IW,': NEUMANN  : SL%A(',IP,')=',SL%A(IP), A_OLD, DBC, IC
=======
               !WRITE(LU_SCARC,'(a,i5,a,i5,a,3f15.6,i5)') 'IW=',IW,': DIRICHLET: SL%A(',IP,')=',SL%A(IP), A_OLD, DBC, IC
            !CASE (INTERNAL)                      !> do nothing along internal boundaries (only debugging)
            CASE (NEUMANN)                        !> set Neumann BC's at all other nodes
               SL%A(IP) = SL%A(IP) + DBC
               !WRITE(LU_SCARC,'(a,i5,a,i5,a,3f15.6,i5)') 'IW=',IW,': NEUMANN  : SL%A(',IP,')=',SL%A(IP), A_OLD, DBC, IC
>>>>>>> origin/scarc
         END SELECT

      ENDDO WALL_CELLS_LOOP2D


   !! ------------------ 2D -----------------------
   CASE (NSCARC_DIMENSION_THREE)

<<<<<<< HEAD
      WALL_CELLS_LOOP3D: DO IW = 1, NW
=======
      WALLCELLS_LOOP3D: DO IW = 1, NW
>>>>>>> origin/scarc

         IOR0 = SL%WALL(IW)%IOR

         I    = SL%WALL(IW)%IXW
         J    = SL%WALL(IW)%IYW
         K    = SL%WALL(IW)%IZW

         NOM  = SL%WALL(IW)%NOM
         ITYPE = SL%CCVAR(I,J,K,CGSC)

<<<<<<< HEAD
         IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.ITYPE/=IS_GASPHASE) CYCLE
=======
         IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC.AND.ITYPE/=IS_GASPHASE) CYCLE
>>>>>>> origin/scarc

         SL%WALL(IW)%ICW = SL%CCVAR(I, J, K, UNKH)
         IC = SL%CCVAR(I, J, K, UNKH)

<<<<<<< HEAD
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
   WRITE(LU_SCARC,*) 'BOUNDARY, PROCESSING: IW, I, J, K, NOM, IOR: ', IW, I, J, K, NOM, IOR0, ITYPE

>>>>>>> origin/scarc
         SELECT CASE (ABS(IOR0))
            CASE (1)
               DBC= SL%DXI2           ! Achtung: Wirklich richtig oder Mittelwert?
            CASE (2)
               DBC= SL%DYI2
            CASE (3)
               DBC= SL%DZI2
         END SELECT
         IP = SL%A_ROW(IC)
         SELECT CASE (SL%WALL(IW)%BTYPE)
            CASE (DIRICHLET)                      !> set Dirichlet BC's at open and null boundary cells
               SL%A(IP) = SL%A(IP) - DBC
<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
    WRITE(LU_SCARC,*) 'DIRICHLET  IW=',IW,' IOR0=',IOR0,': IC=',IC,': DBC=',DBC,': A(',IP,')=',SL%A(IP)
            !CASE (INTERNAL)                      !> do nothing along internal boundaries (only debugging)
            CASE (NEUMANN)                        !> set Neumann BC's at all other cells
               SL%A(IP) = SL%A(IP) + DBC
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
    WRITE(LU_SCARC,*) 'NEUMANN   IW=',IW,' IOR0=',IOR0,': IC=',IC,': DBC=',DBC,': A(',IP,')=',SL%A(IP)
=======
!WRITE(LU_SCARC,*) 'DIRICHLET  IW=',IW,' IOR0=',IOR0,': IC=',IC,': DBC=',DBC,': A(',IP,')=',SL%A(IP)
            !CASE (INTERNAL)                      !> do nothing along internal boundaries (only debugging)
            CASE (NEUMANN)                        !> set Neumann BC's at all other cells
               SL%A(IP) = SL%A(IP) + DBC
!WRITE(LU_SCARC,*) 'NEUMANN   IW=',IW,' IOR0=',IOR0,': IC=',IC,': DBC=',DBC,': A(',IP,')=',SL%A(IP)
>>>>>>> origin/scarc
         END SELECT

      ENDDO WALL_CELLS_LOOP3D

END SELECT SELECT_DIMENSION

!> Only temporarily
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) 'SETUP BOUNDARY '
   WRITE(LU_SCARC,*) '---------------------- A_ROW:', SL%NCS
   WRITE(LU_SCARC,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NCS+1)
   WRITE(LU_SCARC,*) '---------------------- A_COL:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%A_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- A:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SL%A(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_BOUNDARY


!> ----------------------------------------------------------------------------------------------------
!> Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SOLVER
INTEGER :: IERR, NC
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: SL
TYPE (SCARC_PRECON_TYPE) , POINTER :: SP
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE)    , POINTER :: SM
#endif
INTEGER :: NM, NL

IERR = 0

MESHES_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

      SL => SCARC(NM)%LEVEL(NL)

 ! ACHTUNG: HIER CHECKEN !!
<<<<<<< HEAD
 !     IF (PRES_ON_WHOLE_DOMAIN) THEN
=======
 !     IF (PRES_ON_WHOLE_DOMAIN_SCRC) THEN
>>>>>>> origin/scarc
 !        NC = SL%NCE
 !     ELSE
 !        NC = SL%NCS
 !     ENDIF
      NC = SL%NCE

<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'PRES_ON_WHOLE_DOMAIN=',PRES_ON_WHOLE_DOMAIN
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'ALLOCATING VECTORS IN LENGTH ', NC, SL%NCE
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'PRES_ON_WHOLE_DOMAIN_SCRC=',PRES_ON_WHOLE_DOMAIN_SCRC
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'ALLOCATING VECTORS IN LENGTH ', NC, SL%NCE
>>>>>>> origin/scarc

      SELECT_METHOD: SELECT CASE (TRIM(SCARC_METHOD))

         !> -------------------------------------------------------------
         !> working and auxiliary vectors for global CG/BICG-method
         !> -------------------------------------------------------------
         CASE ('KRYLOV')

            ALLOCATE (SL%X(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'X', IERR)
            SL%X = 0.0_EB

            ALLOCATE (SL%F(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'B', IERR)
            SL%F = 0.0_EB

            ALLOCATE (SL%D(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'D', IERR)
            SL%D = 0.0_EB

            ALLOCATE (SL%W(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'W', IERR)
            SL%W = 0.0_EB

            ALLOCATE (SL%G(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'G', IERR)
            SL%G = 0.0_EB

            ALLOCATE (SL%Y(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'Y', IERR)
            SL%Y = 0.0_EB

            ALLOCATE (SL%Z(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'Z', IERR)
            SL%Z = 0.0_EB

            IF (TYPE_PRECON == NSCARC_PRECON_MULTIGRID) THEN

               ALLOCATE (SL%X2(NC), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'X2', IERR)
               SL%X2 = 0.0_EB

               ALLOCATE (SL%D2(NC), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'D2', IERR)
               SL%D2 = 0.0_EB

               IF (NL==NLEVEL_MAX) THEN

                  ALLOCATE (SL%W2(NC), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'W2', IERR)
                  SL%W2 = 0.0_EB

                  ALLOCATE (SL%G2(NC), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'G2', IERR)
                  SL%G2 = 0.0_EB

                  ALLOCATE (SL%Y2(NC), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Y2', IERR)
                  SL%Y2 = 0.0_EB

               ENDIF

#ifdef WITH_MKL
            ELSE IF (TYPE_PRECON == NSCARC_PRECON_PARDISO) THEN

               CALL SCARC_SETUP_PARDISO(NM, NL)

            ELSE IF (TYPE_PRECON == NSCARC_PRECON_CLUSTER) THEN

               CALL SCARC_SETUP_CLUSTER(NM, NL)
#endif

            ELSE IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN

               SP => SCARC(NM)%PRECON(NL)
               ALLOCATE (SP%FFT(1:SL%NX+1, 1:SL%NY+1, 1:SL%NZ+1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'FFT', IERR)
               SP%FFT = 0.0_EB

            ENDIF


         !> -------------------------------------------------------------
         !> working and auxiliary vectors for global GMG/AMG-method
         !> -------------------------------------------------------------
         CASE ('MULTIGRID')

            ALLOCATE (SL%X(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'X', IERR)
            SL%X = 0.0_EB

            ALLOCATE (SL%F(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'F', IERR)
            SL%F = 0.0_EB

            ALLOCATE (SL%D(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'D', IERR)
            SL%D = 0.0_EB

            IF (NL==NLEVEL_MAX) THEN

               ALLOCATE (SL%W(NC), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'W', IERR)
               SL%W = 0.0_EB

               ALLOCATE (SL%G(NC), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'G', IERR)
               SL%G = 0.0_EB

               ALLOCATE (SL%Y(NC), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'Y', IERR)
               SL%Y = 0.0_EB

            ENDIF

            IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
               SP => SCARC(NM)%PRECON(NL)
               ALLOCATE (SP%FFT(1:SL%NX+1, 1:SL%NY+1, 1:SL%NZ+1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'FFT', IERR)
               SP%FFT = 0.0_EB
            ENDIF

#ifdef WITH_MKL
            IF (TYPE_SMOOTH == NSCARC_SMOOTH_PARDISO) CALL SCARC_SETUP_PARDISO(NM, NL)

            IF (TYPE_COARSE == NSCARC_COARSE_DIRECT .AND. NL== NLEVEL_MAX) THEN
               IF (N_MPI_PROCESSES > 1) THEN
                  CALL SCARC_SETUP_CLUSTER(NM, NLEVEL_MAX)
               ELSE
                  CALL SCARC_SETUP_PARDISO(NM, NLEVEL_MAX)
               ENDIF
            ENDIF
#endif

#ifdef WITH_MKL
         !> -------------------------------------------------------------
         !> working and auxiliary vectors for global GMG/AMG-method
         !> -------------------------------------------------------------
         CASE ('MKL')

            SM => SCARC(NM)%MKL(NL)

            ALLOCATE (SL%X(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'X', IERR)
            SL%X = 0.0_EB

            ALLOCATE (SL%F(NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'B', IERR)
            SL%F = 0.0_EB

            SELECT_MKL: SELECT CASE (TYPE_MKL)
               CASE (NSCARC_MKL_GLOBAL)
                  CALL SCARC_SETUP_CLUSTER(NM, NL)
               CASE (NSCARC_MKL_LOCAL)
                  CALL SCARC_SETUP_PARDISO(NM, NL)
            END SELECT SELECT_MKL
#endif

      END SELECT SELECT_METHOD


   ENDDO LEVEL_LOOP
ENDDO  MESHES_LOOP

END SUBROUTINE SCARC_SETUP_SOLVER


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CLUSTER(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, IP, IERR, IC
REAL (EB) :: TNOW, DDUMMY(1)=0.0_EB
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: SL
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE)    , POINTER :: SM
#endif

TNOW = SECOND()

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'TRYING TO SETUP CLUSTER FOR LEVEL ',NL

SL => SCARC(NM)%LEVEL(NL)
#ifdef WITH_MKL
SM => SCARC(NM)%MKL(NL)
#endif

ALLOCATE(SM%IPARM(64), STAT=IERR) 
CALL CHKMEMERR ('SCARC', 'IPARM', IERR)
SM%IPARM = 0

ALLOCATE(SM%CT(64), STAT=IERR) 
CALL CHKMEMERR ('SCARC', 'CT', IERR)
DO I=1,64
   SM%CT(I)%DUMMY = 0
ENDDO

SM%NRHS   =  1         ! one right hand side
SM%MAXFCT =  1         ! one matrix
SM%MNUM   =  1         ! number of matrix to be factorized
SM%ERROR  =  0         ! initialize error flag
SM%MSGLVL =  0         ! print statistical information

IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
   SM%MTYPE  = -2         ! Matrix type real and symmetric indefinite
ELSE
   SM%MTYPE  = 11         ! Matrix type real and non-symmetric
ENDIF

! Define control parameter vector iparm:
!SM%IPARM(1) = 1   ! no solver default
!SM%IPARM(2) = 3   ! Parallel fill-in reordering from METIS
!SM%IPARM(4) = 0   ! no iterative-direct algorithm
!SM%IPARM(5) = 0   ! no user fill-in reducing permutation
!SM%IPARM(6) = 2   ! =0 solution on the first n components of x
!SM%IPARM(8) = 2   ! numbers of iterative refinement steps
!SM%IPARM(10) = 13 ! perturb the pivot elements with 1E-13
!SM%IPARM(11) = 1  ! use nonsymmetric permutation and scaling MPS  !!!!! was 1
!SM%IPARM(13) = 1  ! maximum weighted matching algorithm is switched-off
!               !(default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
!SM%IPARM(14) = 0  ! Output: number of perturbed pivots
!SM%IPARM(18) = 0 !-1 ! Output: number of nonzeros in the factor LU
!SM%IPARM(19) = 0 !-1 ! Output: Mflops for LU factorization
!SM%IPARM(20) = 0  ! Output: Numbers of CG Iterations
!SM%IPARM(21) = 1  ! 1x1 diagonal pivoting for symmetric indefinite matrices.
!SM%IPARM(24) = 0
!SM%IPARM(27) = 1 ! Check matrix
!SM%IPARM(37) = 0 ! CSR-format
!SM%IPARM(40) = 2 ! Matrix, solution and rhs provided in distributed assembled matrix input format.


!WRITE(*,*) 'ACHTUNG HIER NOCHMAL CHECKEN'
SM%IPARM(1)  =  1      ! supply own parameters
<<<<<<< HEAD
SM%IPARM(2)  =  3      ! supply own parameters
SM%IPARM(4)  =  0      ! supply own parameters
SM%IPARM(5)  =  0      ! supply own parameters
SM%IPARM(6)  =  0      ! write solution to x
SM%IPARM(8)  =  2      ! automatic setting of iterative refinement steps
SM%IPARM(10) = 13      ! pivoting perturbation
SM%IPARM(11) =  1      ! pivoting perturbation
SM%IPARM(13) =  1      ! pivoting perturbation
SM%IPARM(14) =  0      ! pivoting perturbation
SM%IPARM(18) =  0      ! pivoting perturbation
SM%IPARM(19) =  0      ! pivoting perturbation
SM%IPARM(20) =  0      ! pivoting perturbation
SM%IPARM(21) =  1      ! Bunch-Kaufman pivoting which is default in case of IPARM(0)=0
SM%IPARM(24) =  0      ! Bunch-Kaufman pivoting which is default in case of IPARM(0)=0
SM%IPARM(27) =  1      ! use matrix checker
SM%IPARM(40) = 2       ! provide matrix in distributed format
SM%IPARM(41) = NCS_OFFSET(NM, NL) + 1                  ! first global cell number for mesh NM
SM%IPARM(42) = NCS_OFFSET(NM, NL) + NCS_LOCAL(NM, NL)  ! last global cell number for mesh NM
=======
SM%IPARM(2)  =  2      ! supply own parameters
SM%IPARM(6)  =  0      ! write solution to x
SM%IPARM(8)  =  2      ! automatic setting of iterative refinement steps
SM%IPARM(10) = 13      ! pivoting perturbation
SM%IPARM(21) =  1      ! Bunch-Kaufman pivoting which is default in case of IPARM(0)=0
SM%IPARM(27) =  1      ! use matrix checker
SM%IPARM(37) =  0      ! CSR format for matrix storage
SM%IPARM(40) = 2       ! provide matrix in distributed format
SM%IPARM(41) = NCS_OFFSET(NM, NL)+1                  ! first global cell number for mesh NM
SM%IPARM(42) = NCS_OFFSET(NM, NL)+NCS_LOCAL(NM, NL)  ! last global cell number for mesh NM
>>>>>>> origin/scarc
!SM%IPARM(39) = 2       ! provide matrix in distributed format
!SM%IPARM(40) = NCS_OFFSET(NM, NL)+1                  ! first global cell number for mesh NM
!SM%IPARM(41) = NCS_OFFSET(NM, NL)+NCS_LOCAL(NM, NL)  ! last global cell number for mesh NM


! Only reordering and symbolic factorization
   
<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'BEFORE FACTORIZATION'
   WRITE(LU_SCARC,*) '-------------- IPARM:A: -----------------'
   WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
   WRITE(LU_SCARC,*) 'PERM=',SM%PERM
   WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
   WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
   WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
   WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
   WRITE(LU_SCARC,*) 'ERROR=',SM%ERROR
<<<<<<< HEAD
   WRITE(LU_SCARC,*) 'IPARM(40:42)=',SM%IPARM(40:42)
=======
   WRITE(LU_SCARC,*) 'IPARM(39:42)=',SM%IPARM(39:42)
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'NCS_GLOBAL(NL)=',NCS_GLOBAL(NL)
   WRITE(LU_SCARC,*) 'SIZE(SM%ASYM)=',SIZE(SM%ASYM)
   WRITE(LU_SCARC,*) 'SIZE(SM%ASYM_ROW)=',SIZE(SM%ASYM_ROW)
   WRITE(LU_SCARC,*) 'SIZE(SM%ASYM_COL)=',SIZE(SM%ASYM_COL)
   WRITE(LU_SCARC,*) 'SIZE(X)=',SIZE(SL%X)
   WRITE(LU_SCARC,*) 'SIZE(F)=',SIZE(SL%F)
   WRITE(LU_SCARC,*) '---------------------- ASYM_ROW:', SL%NCS
   WRITE(LU_SCARC,'(4i9)') (SM%ASYM_ROW(IC), IC=1,SL%NCS+1)
   WRITE(LU_SCARC,*) '---------------------- ASYM_COL:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SM%ASYM_COL(IP),IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- AS:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SM%ASYM(IP),IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1)
   ENDDO
ENDIF

SM%PHASE = 11
CALL CLUSTER_SPARSE_SOLVER(SM%CT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, NCS_GLOBAL(NL), &
                           SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, SM%PERM, SM%NRHS, SM%IPARM, &
                           SM%MSGLVL, DDUMMY, DDUMMY, MPI_COMM_WORLD, SM%ERROR)
               
IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'Reordering completed ... ', SM%ERROR
   IF (SM%ERROR /= 0) THEN
      WRITE(LU_SCARC,*) 'Sym Factor: The following ERROR was detected: ', SM%ERROR
   END IF
   WRITE(LU_SCARC,*) 'Number of nonzeros in factors= ',SM%IPARM(18)
   WRITE(LU_SCARC,*) 'Number of factorization MFLOPS= ',SM%IPARM(19)
ENDIF
<<<<<<< HEAD

=======
>>>>>>> origin/scarc
               
! only Factorization.
SM%PHASE = 22 
CALL CLUSTER_SPARSE_SOLVER(SM%CT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, NCS_GLOBAL(NL), &
                           SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, SM%PERM, SM%NRHS, SM%IPARM, &
                           SM%MSGLVL, DDUMMY, DDUMMY, MPI_COMM_WORLD, SM%ERROR)
         
IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'Factorization completed ... ', SM%ERROR
   IF (SM%ERROR /= 0) THEN
      WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR, SM%IPARM(30)
   ENDIF
ENDIF
<<<<<<< HEAD


=======
WRITE(*,*) 'SETUP CLUSTER:', SECOND()-TNOW
>>>>>>> origin/scarc
TUSED_SCARC(NSCARC_TIME_PARDISO_SETUP,MYID+1)=TUSED_SCARC(NSCARC_TIME_PARDISO_SETUP,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL  ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL  ,MYID+1)+SECOND()-TNOW
END SUBROUTINE SCARC_SETUP_CLUSTER
#endif


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PARDISO(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, IERR, IDUMMY(1)=0, IC, IP
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: SL
TYPE (SCARC_MKL_TYPE)    , POINTER :: SM
REAL (EB) :: TNOW, DDUMMY(1)=0.0_EB

TNOW = SECOND()

SL => SCARC(NM)%LEVEL(NL)
SM => SCARC(NM)%MKL(NL)

ALLOCATE(SM%IPARM(64), STAT=IERR) 
CALL CHKMEMERR ('SCARC', 'IPARM', IERR)
SM%IPARM = 0

ALLOCATE(SM%PT(64), STAT=IERR) 
CALL CHKMEMERR ('SCARC', 'PT', IERR)
DO I=1,64
   SM%PT(I)%DUMMY = 0
ENDDO

SM%NRHS   = 1
SM%MAXFCT = 1
SM%MNUM   = 1

!SM%IPARM(1) = 1   ! no solver default
!SM%IPARM(2) = 2   ! fill-in reordering from METIS
!SM%IPARM(4) = 0   ! no iterative-direct algorithm
!SM%IPARM(5) = 0   ! no user fill-in reducing permutation
!SM%IPARM(6) = 2   ! =0 solution on the first n components of x
!SM%IPARM(8) = 2   ! numbers of iterative refinement steps
!SM%IPARM(10) = 13 ! perturb the pivot elements with 1E-13
!SM%IPARM(11) = 1  ! use nonsymmetric permutation and scaling MPS  !!!!! was 1
!SM%IPARM(13) = 1  ! maximum weighted matching algorithm is switched-off
!                  !(default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
!SM%IPARM(14) = 0  ! Output: number of perturbed pivots
!SM%IPARM(18) = 0  !-1 ! Output: number of nonzeros in the factor LU
!SM%IPARM(19) = 0  !-1 ! Output: Mflops for LU factorization
!SM%IPARM(20) = 0  ! Output: Numbers of CG Iterations
!SM%IPARM(21) = 1  ! 1x1 diagonal pivoting for symmetric indefinite matrices.
!SM%IPARM(24) = 0
!SM%IPARM(27) = 1  ! Check matrix
!SM%IPARM(37) = 0  ! Matrix, solution and rhs provided in distributed assembled matrix input format.  ???
!SM%IPARM(40) = 2  ! Matrix, solution and rhs provided in distributed assembled matrix input format.   ???


SM%IPARM(1)  =  1      ! no solver default
SM%IPARM(2)  =  2      ! nested dissection algorithm
SM%IPARM(4)  =  0      ! factorization computed as required by phase
SM%IPARM(5)  =  0      ! user permutation ignored
SM%IPARM(6)  =  0      ! write solution on x
SM%IPARM(8)  =  2      ! numbers of iterative refinement steps
SM%IPARM(10) = 13      ! perturb the pivot elements with 1E-13
SM%IPARM(11) =  0      ! disable scaling (default for SPD)
SM%IPARM(13) =  0      ! disable matching
SM%IPARM(18) = -1      ! Output: number of nonzeros in the factor LU
SM%IPARM(19) = -1      ! Output: number of floating points operations
SM%IPARM(20) =  1      ! Output: Numbers of CG Iterations
SM%IPARM(27) =  1      ! use matrix checker
SM%IPARM(37) =  0      ! matrix storage in CSR-format

SM%ERROR  =  0       ! initialize error flag
SM%MSGLVL =  0       ! print statistical information
SM%MTYPE  = -2       ! Matrix type real non-symmetric

! Only reordering and symbolic factorization
SM%PHASE = 11
   
IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
WRITE(LU_SCARC,*) '-------------- PARDISO-INIT -----'
WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
WRITE(LU_SCARC,*) 'PERM=',SM%PERM
WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
WRITE(LU_SCARC,*) 'PHASE=',SM%PHASE
WRITE(LU_SCARC,*) 'IDUM=',IDUMMY
WRITE(LU_SCARC,*) 'DDUM=',DDUMMY
WRITE(LU_SCARC,*) 'ERROR=',SM%ERROR
   WRITE(LU_SCARC,*) '---------------------- ASYM_ROW:', SL%NCS
   WRITE(LU_SCARC,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NCS+1)
   WRITE(LU_SCARC,*) '---------------------- ASYM_COL:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SM%ASYM_COL(IP),IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1)
   ENDDO
   WRITE(LU_SCARC,*) '---------------------- AS:'
   DO IC = 1, SL%NCS
      WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SM%ASYM(IP),IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1)
   ENDDO
ENDIF


CALL PARDISO_D(SM%PT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, SL%NCS, SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, &
               IDUMMY, SM%NRHS, SM%IPARM, SM%MSGLVL, DDUMMY, DDUMMY, SM%ERROR)
               
IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'Reordering ZZ completed ... ', SM%ERROR
   IF (SM%ERROR /= 0) THEN
      WRITE(LU_SCARC,*) 'Sym Factor: The following ERROR was detected: ', SM%ERROR
   END IF
   WRITE(LU_SCARC,*) 'Number of nonzeros in factors ZZ= ',SM%IPARM(18)
   WRITE(LU_SCARC,*) 'Number of factorization MFLOPS ZZ= ',SM%IPARM(19)
   WRITE(LU_SCARC,*) 'ERROR:1: ',SM%ERROR
ENDIF
               
! only Factorization.
SM%PHASE = 22 
CALL PARDISO_D(SM%PT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, SL%NCS, SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, &
               IDUMMY, SM%NRHS, SM%IPARM, SM%MSGLVL, DDUMMY, DDUMMY, SM%ERROR)

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'Factorization completed ... ', SM%ERROR
   IF (SM%ERROR /= 0) THEN
      WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR, SM%IPARM(30)
   ENDIF
   WRITE(LU_SCARC,*) 'ERROR:2: ',SM%ERROR
ENDIF

TUSED_SCARC(NSCARC_TIME_PARDISO_SETUP,MYID+1)=TUSED_SCARC(NSCARC_TIME_PARDISO_SETUP,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL  ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL  ,MYID+1)+SECOND()-TNOW
END SUBROUTINE SCARC_SETUP_PARDISO
#endif


!> ------------------------------------------------------------------------------------------------
!> Set sizes for transfer matrices
!> ------------------------------------------------------------------------------------------------
SUBROUTINE  SCARC_SETUP_SIZES(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: IW, IC, IG, ICP, IERR, ICCE
INTEGER :: NM, NOM
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SLF, SLC
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSLF, OSLC

SCARC_ROUTINE = 'SCARC_SETUP_SIZES'

SELECT CASE (NTYPE)

   !! --------------------------------------------------------------------------------------------------
   !! Define sizes for system matrix A (including extended regions related to overlaps)
   !! --------------------------------------------------------------------------------------------------
   CASE (NSCARC_SIZE_MATRIX)

      LEVEL_SYSTEM_MESHES_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SLF => SCARC(NM)%LEVEL(NL)

         !> Determine number of couplings in matrix stencil and matrix size in own part of mesh
         SELECT CASE (TYPE_DIMENSION)
            CASE (NSCARC_DIMENSION_TWO)
               SLF%NPOINTS = 5
            CASE (NSCARC_DIMENSION_THREE)
               SLF%NPOINTS = 7
         END SELECT
         SLF%NA  = SLF%NCS * SLF%NPOINTS
         SLF%NAS = SLF%NCS * SLF%NPOINTS
<<<<<<< HEAD


         IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
            WRITE(LU_SCARC,*) 'NPOINTS=',SLF%NPOINTS
            WRITE(LU_SCARC,*) 'NC     =',SLF%NC
            WRITE(LU_SCARC,*) 'NC     =',SLF%NC
            WRITE(LU_SCARC,*) 'NCS    =',SLF%NCS
            WRITE(LU_SCARC,*) 'NA     =',SLF%NA
            WRITE(LU_SCARC,*) 'NAS    =',SLF%NAS
         ENDIF

=======


         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
            WRITE(LU_SCARC,*) 'NPOINTS=',SLF%NPOINTS
            WRITE(LU_SCARC,*) 'NC     =',SLF%NC
            WRITE(LU_SCARC,*) 'NC     =',SLF%NC
            WRITE(LU_SCARC,*) 'NCS    =',SLF%NCS
            WRITE(LU_SCARC,*) 'NA     =',SLF%NA
            WRITE(LU_SCARC,*) 'NAS    =',SLF%NAS
         ENDIF

>>>>>>> origin/scarc
         !> Determine sizes of overlapped parts for later communication with corresponding neighbors
         DO IW = 1, SLF%NW
            NOM = SLF%WALL(IW)%NOM
            IF (NOM /= 0) THEN

               OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

               IC = SLF%WALL(IW)%ICW
               OSLF%NA0 = OSLF%NA0 + SLF%NPOINTS

               !WRITE(LU_SCARC,'(4(a,i4))') 'IW=',IW,': NOM=',NOM,': IC=',IC,': NA0=',OSLF%NA0

               !IF (TYPE_LAYER == NSCARC_LAYER_TWO) THEN
               !   IC = SLF%WALL(IW)%ICW2
               !   OSLF%NA0 = OSLF%NA0 + SLF%NPOINTS
               !ENDIF

            ENDIF
         ENDDO

      ENDDO LEVEL_SYSTEM_MESHES_LOOP

      !> Exchange sizes for AMG matrices and vectors
      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SIZE, NL)


      !> Determine extended sizes for extended prolongation and restriction matrix
      LEVEL_SYSTEM_MESHES_LOOP2: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SLF => SCARC(NM)%LEVEL(NL)
         SLF%NAE = SLF%NA

         ALLOCATE (SLF%CELL_MAP(3,SLF%NCS+1:2*SLF%NCE), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'SLF%CELL_MAP', IERR)
         SLF%CELL_MAP = 0

         LEVEL_OTHER_MESHES_LOOP2: DO NOM = 1, NMESHES

            IF (NOM == NM) CYCLE LEVEL_OTHER_MESHES_LOOP2
            OS => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE LEVEL_OTHER_MESHES_LOOP2

            OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            SLF%NAE  = SLF%NAE  + OSLF%NA

            ALLOCATE (OSLF%A_SIZE(1:OSLF%NCG), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'OSLF%A_SIZE', IERR)
            OSLF%A_SIZE = 0

         ENDDO LEVEL_OTHER_MESHES_LOOP2

      ENDDO LEVEL_SYSTEM_MESHES_LOOP2

   !! -------------------------------------------------------------------------------------------
   !! Define sizes for transfer matrices P and R (including extended regions related to overlaps)
   !! -------------------------------------------------------------------------------------------
   CASE (NSCARC_SIZE_TRANSFER)


      TRANSFER_MESHES_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SLF => SCARC(NM)%LEVEL(NL)
         SLC => SCARC(NM)%LEVEL(NL+1)

         !> Determine dimensions of restriction and prolongation matrices in own mesh
         SLF%NCF  = 0
         SLF%NCC  = 0
         SLF%NP   = 0
         SLF%NR   = 0
         SLF%NCCI = 0

         DO IC = 1, SLF%NC
            IF (SLF%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
               SLF%NCC = SLF%NCC + 1
               SLF%NP  = SLF%NP  + 1
               SLF%NR  = SLF%NP
    !>         SLF%CELL_TYPE(IC)  = SLF%NCC
            ELSE
               SLF%NCF = SLF%NCF + 1
               SLF%NP  = SLF%NP  + SLF%A_ROW(IC+1)-SLF%A_ROW(IC) - 1
               SLF%NR  = SLF%NP
            ENDIF
         ENDDO

    !>   SLF%NCCE = SLF%NCC
    !>   SLF%NCFE = SLF%NCF
    !>   DO IC = SLF%NC+1, SLF%NCE
    !>      IF (SLF%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
    !>         SLF%NCCE = SLF%NCCE + 1
    !>         SLF%CELL_TYPE(IC)  = SLF%NCCE
    !>      ELSE
    !>         SLF%NCFE = SLF%NCFE + 1
    !>      ENDIF
    !>   ENDDO

         SLF%NCCI = SLF%NCC
         SLF%NCCE = SLF%NCC
         !SLF%NCE0 = SLF%NC+1      !> really needed? OUTDATED ?

         !WRITE(LU_SCARC,*) 'NCC =',SLF%NCC
         !WRITE(LU_SCARC,*) 'NCCE=',SLF%NCCE
         !WRITE(LU_SCARC,*) 'NP  =',SLF%NP
         !WRITE(LU_SCARC,*) 'NR  =',SLF%NR
         !WRITE(LU_SCARC,*) 'CELL_TYPE  ='
         !WRITE(LU_SCARC,'(4i4)') (SLF%CELL_TYPE(IC), IC=1, SLF%NCE)

         !> Determine number of coarse and fine cells and check correctness of computation
         IF (SLF%NCC + SLF%NCF /= SLF%NC) THEN
            WRITE(SCARC_MESSAGE,'(2A,I8,A,I8,A,I4)') TRIM(SCARC_ROUTINE),&
                 ': N_CELLS_COARSE + N_CELLS_FINE = ', SLF%NCC + SLF%NCF, &
                 '  differs from N_CELLS = ', SLF%NC, ' on level ', NL
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
         ENDIF


         !> define variables for overlapping parts
         TRANSFER_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

            IF (NOM == NM) CYCLE TRANSFER_OTHER_MESHES_LOOP
            OS => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE TRANSFER_OTHER_MESHES_LOOP

            OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            OSLF%NP   = 0
            OSLF%NR   = 0
            OSLF%NP0  = 0
            OSLF%NR0  = 0
            OSLF%NCC  = 0
            OSLF%NCF  = 0
            OSLF%NCC0 = 0
            OSLF%NCF0 = 0
            OSLF%ICG0  = 0

         ENDDO TRANSFER_OTHER_MESHES_LOOP

         !> Determine sizes of overlapped parts for later communication with corresponding neighbors
         ICCE = SLF%NCC
         IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) SLF%NCW = 0

         DO IW = 1, SLF%NW
            NOM = SLF%WALL(IW)%NOM
            IF (NOM /= 0) THEN

               OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

               !WRITE(*,*) 'K: HIER WALL%ICN �ndern auf OSLF!'
               !WRITE(*,*) 'K: MUSS NOCH UEBERARBEITET WERDEN !!'

               IC  = SLF%WALL(IW)%ICW
               !ICE = SLF%WALL(IW)%ICE(1)

               IF (TYPE_COARSENING<NSCARC_COARSENING_GMG.AND.SLF%CELL_TYPE(IC)>0) SLF%NCW = SLF%NCW + 1

               !WRITE(LU_SCARC,*) 'CELL_TYPE(',IC,')=',SLF%CELL_TYPE(IC)

               !IF (SLF%CELL_TYPE(ICE) >= NSCARC_CELL_TYPE_COARSE) THEN
               IF (SLF%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
                  OSLF%NP0  = OSLF%NP0  + 1
                  OSLF%ICG0  = OSLF%ICG0  + 1
                  OSLF%NR0  = OSLF%NP0
                  OSLF%NCC0 = OSLF%NCC0 + 1

                  !WRITE(LU_SCARC,*) 'OSLF(',NOM,')%ICG0 =',OSLF%ICG0
                  !WRITE(LU_SCARC,*) 'OSLF(',NOM,')%NCC0=',OSLF%NCC0

               ELSE
                  OSLF%NP0  = OSLF%NP0  + SLF%A_ROW(IC+1)-SLF%A_ROW(IC) - 1
                  OSLF%NR0  = OSLF%NP0
                  OSLF%NCF0 = OSLF%NCF0 + 1
               ENDIF

               OSLF%NCC=OSLF%NCC0
               OSLF%NCF=OSLF%NCF0

            ENDIF
         ENDDO

         !!!
         !> Determine new numbering for coarse cells in interior of mesh
         !!!
         IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) THEN
            ICP   = 0
            DO IC = 1, SLF%NCE
               IF (SLF%CELL_TYPE(IC) == NSCARC_CELL_TYPE_COARSE) THEN
                  ICP = ICP + 1
                  SLF%CELL_TYPE(IC) = ICP
               ENDIF
            ENDDO
         ENDIF

         !WRITE(LU_SCARC,*) 'NCC =',SLF%NCC
         !WRITE(LU_SCARC,*) 'NCCE=',SLF%NCCE
         !WRITE(LU_SCARC,*) 'NP  =',SLF%NP
         !WRITE(LU_SCARC,*) 'NR  =',SLF%NR
         !WRITE(LU_SCARC,*) 'CELL_TYPE  ='
         !WRITE(LU_SCARC,'(4i4)') (SLF%CELL_TYPE(IC), IC=1, SLF%NCE)

      ENDDO TRANSFER_MESHES_LOOP

      !> Exchange sizes for AMG matrices and vectors
      IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_TRANSFER_SIZE, NL)

      !> Determine extended sizes for extended prolongation and restriction matrix
      TRANSFER_MESHES_LOOP2: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SLF => SCARC(NM)%LEVEL(NL)

         SLF%NPE  = SLF%NP
         SLF%NRE  = SLF%NR
         SLF%NCCE = SLF%NCC
         SLF%NCFE = SLF%NCF

         TRANSFER_OTHER_MESHES_LOOP2: DO NOM = 1, NMESHES

            IF (NOM == NM) CYCLE TRANSFER_OTHER_MESHES_LOOP2

            OS => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE TRANSFER_OTHER_MESHES_LOOP2

            OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            OSLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

            OSLC%NC  = OSLF%NCC0
            OSLF%NCC = OSLF%NCC0
            OSLC%NCG = OSLF%ICG0

            ALLOCATE (OSLC%ICG_TO_ICE(1:OSLC%NC), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'OSLC%ICG_TO_ICE', IERR)
            OSLC%ICG_TO_ICE = 0

WRITE(*,*) 'ACHTUNG HIER WEGEN ICN SCHAUEN'
            ALLOCATE (OSLC%ICN_TO_ICE(1:OSLF%NCCI), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'OSLC%ICN_TO_ICE', IERR)
            OSLC%ICN_TO_ICE = 0

            DO IG = 1, OSLC%NC
               OSLC%ICG_TO_ICE(IG) = SLF%NCCE + IG
            ENDDO

            SLF%NPE  = SLF%NPE  + OSLF%NP
            SLF%NRE  = SLF%NRE  + OSLF%NR
            SLF%NCCE = SLF%NCCE + OSLF%NCC
            SLF%NCFE = SLF%NCFE + OSLF%NCF

            !WRITE(LU_SCARC,*) '==============  NM=',NM,': NOM=',NOM
            !WRITE(LU_SCARC,*) 'OSLF%NCC =',OSLF%NCC,': OSLF%NC=',OSLF%NC
            !WRITE(LU_SCARC,*) 'SLF%NPE =',SLF%NPE,': SLF%NP=',SLF%NP
            !WRITE(LU_SCARC,*) 'SLF%NRE =',SLF%NRE,': SLF%NR=',SLF%NR
            !WRITE(LU_SCARC,*) 'SLF%NCCE=',SLF%NCCE,': SLF%NCCE=',SLF%NCCE
            !WRITE(LU_SCARC,*) 'SLF%NCFE=',SLF%NCFE,': SLF%NCFE=',SLF%NCFE

         ENDDO TRANSFER_OTHER_MESHES_LOOP2

      ENDDO TRANSFER_MESHES_LOOP2

      !> Determine extended sizes for extended prolongation and restriction matrix
      TRANSFER_MESHES_LOOP3: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SLF => SCARC(NM)%LEVEL(NL)
         SLC => SCARC(NM)%LEVEL(NL+1)

         SLC%NC  = SLF%NCC
         SLC%NCE = SLF%NCCE
         SLC%NW  = SLF%NCW

         !WRITE(LU_SCARC,*) 'SLC%NC =',SLC%NC
         !WRITE(LU_SCARC,*) 'SLC%NCE=',SLC%NCE
         !WRITE(LU_SCARC,*) 'SLC%NW =',SLC%NW
         !WRITE(LU_SCARC,*) 'ALLOCATING SLC%WALL in length ', SLC%NW, ' to ', SLC%NW

         ALLOCATE (SLC%WALL(1:SLC%NW), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'SLC%WALL', IERR)

         ALLOCATE (SLC%CELL_TYPE(1:SLC%NCE), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'SLC%CELL_TYPE', IERR)
         SLC%CELL_TYPE = NSCARC_CELL_TYPE_NONE

         ALLOCATE (SLC%MEASURE(1:SLC%NCE), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'SLC%MEASURE', IERR)
         SLC%MEASURE = NSCARC_MEASURE_NONE

         ALLOCATE(SLC%INTERNAL_BDRY_CELL(SLC%NC), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION','INTERNAL_BDRY_CELL',IERR)
         SLC%INTERNAL_BDRY_CELL = NSCARC_LAYER_NONE

         ALLOCATE(SLC%WALL_INDEX(SLC%NC, -3:3), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION','WALL_INDEX',IERR)
         SLC%WALL_INDEX = 0

         IF (SLC%NCE > SLC%NC) THEN
            ALLOCATE(SLC%ICE_TO_IWG(SLC%NC+1:SLC%NCE), STAT=IERR)
            CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION','ICE_TO_IWG',IERR)

            ALLOCATE(SLC%EXT_PTR(SLC%NC+1:SLC%NCE), STAT=IERR)
            CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION','EXT_PTR',IERR)
         ENDIF

      ENDDO TRANSFER_MESHES_LOOP3

END SELECT

END SUBROUTINE SCARC_SETUP_SIZES



!> ------------------------------------------------------------------------------------------------
!> Initialize global 3D-solver methods (cg/mg)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSENING
INTEGER :: IERR, NL, NM, NOM
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLF, SLC
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSLF
TYPE (OSCARC_TYPE), POINTER :: OSF

IERR = 0
IF (TYPE_MULTIGRID /= NSCARC_MULTIGRID_ALGEBRAIC) RETURN

!> Determine number of multigrid levels
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX-1
!LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MIN

   !! ---------------------------------------------------------------------------------------------
   !! Determine coarser meshes corresponding to requested coarsening strategy
   !!  --- allocate necessary arrays
   !!  --- setup measures of single cells
   !!  --- setup CELL_TYPEs of single cells
   !!  --- setup sizes of transformation matrices (prolongation/restriction)
   !! ---------------------------------------------------------------------------------------------
   !IF (NL == NLEVEL_MIN) THEN

   DO NM = 1, NMESHES

      IF (PROCESS(NM) /= MYID) CYCLE

      SLF => SCARC(NM)%LEVEL(NL)

      ALLOCATE (SLF%MEASURE(1:SLF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLF%MEASURE', IERR)
      SLF%MEASURE = NSCARC_MEASURE_NONE

      ALLOCATE (SLF%CELL_TYPE(1:SLF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLF%CELL_TYPE', IERR)
      SLF%CELL_TYPE = NSCARC_CELL_TYPE_NONE

      ALLOCATE (SLF%CELL_MARKER(1:SLF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLF%CELL_MARKER', IERR)
      SLF%CELL_MARKER = -1

      ALLOCATE (SLF%S(SLF%NA), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLF%S', IERR)
      CALL SCARC_COPY_REAL(SLF%A, SLF%S, 1.0_EB, SLF%NA)

      ALLOCATE (SLF%S_COL(SLF%NA), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLF%S_COL', IERR)
      CALL SCARC_COPY_INTEGER(SLF%A_COL, SLF%S_COL, 1, SLF%NA)

      ALLOCATE (SLF%S_ROW(SLF%NC+1), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLF%S_ROW', IERR)
      CALL SCARC_COPY_INTEGER(SLF%A_ROW, SLF%S_ROW, 1, SLF%NC+1)

      ALLOCATE (SLF%ST_COL(SLF%NAE+1), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLF%ST_COL', IERR)
      SLF%ST_COL = 0

      ALLOCATE (SLF%ST_ROW(SLF%NCE+1), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLF%ST_ROW', IERR)
      SLF%ST_ROW = 0

      ALLOCATE (SLF%GRAPH(1:SLF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLF%GRAPH', IERR)
      SLF%GRAPH = NSCARC_GRAPH_NONE

      OTHER_MESH_LOOP1: DO NOM = 1, NMESHES

         IF (NOM == NM) CYCLE OTHER_MESH_LOOP1

         OSF  => SCARC(NM)%OSCARC(NOM)
         IF (OSF%NICMAX_S==0 .AND. OSF%NICMAX_R==0) CYCLE OTHER_MESH_LOOP1

         OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

         ALLOCATE (OSLF%MEASURE(1:OSLF%NCG), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'OSLF%MEASURE', IERR)
         OSLF%MEASURE = NSCARC_MEASURE_NONE

         ALLOCATE (OSLF%CELL_TYPE(1:OSLF%NCG), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'OSLF%CELL_TYPE', IERR)
         OSLF%CELL_TYPE = NSCARC_CELL_TYPE_NONE

         ALLOCATE (OSLF%GRAPH(1:OSLF%NCG), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'OSLF%GRAPH', IERR)
         OSLF%GRAPH = NSCARC_GRAPH_NONE

      ENDDO OTHER_MESH_LOOP1

   ENDDO

   !ENDIF

   CALL SCARC_SETUP_STRENGTH_MATRIX  (NSCARC_COARSENING_RS3, NL)

   !! Then set measures and CELL_TYPEs on internal cells due to chosen coarsening strategy
   SELECT CASE (TYPE_COARSENING)
      CASE (NSCARC_COARSENING_GMG)
         CALL SCARC_SETUP_CELL_TYPES (NSCARC_COARSENING_GMG, NL)
      CASE (NSCARC_COARSENING_RS3)
         CALL SCARC_SETUP_COLORING  (NSCARC_COARSENING_RS3, NL)
      CASE (NSCARC_COARSENING_FALGOUT)
         CALL SCARC_SETUP_COLORING  (NSCARC_COARSENING_RS3, NL)
         CALL SCARC_SETUP_COLORING  (NSCARC_COARSENING_FALGOUT, NL)
   END SELECT

   !! Set dimensions for coarser mesh and define sizes of prolongation and restriction matrices
   CALL SCARC_SETUP_SIZES (NSCARC_SIZE_TRANSFER, NL)

   !CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE,  NL)
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPES', 'CELL_TYPE BBB')


   !! --------------------------------------------------------------------------------------------
   !! Allocate and define grid transfer matrices
   !!  NCE  : number of extended cells in fine grid
   !!  NCCE : number of extended cells in coarse grid
   !!  NCW  : number of wall cells in coarse grid
   !! --------------------------------------------------------------------------------------------
   DO NM = 1, NMESHES

      IF (PROCESS(NM) /= MYID) CYCLE

      SLF => SCARC(NM)%LEVEL(NL)

      !WRITE(LU_SCARC,*) 'ALLOCATING P in size ', SLF%NPE
      !WRITE(LU_SCARC,*) 'ALLOCATING R in size ', SLF%NRE
      !WRITE(LU_SCARC,*) 'CURRENT FINE   LEVEL ', NL
      !WRITE(LU_SCARC,*) 'CURRENT COARSE LEVEL ', NL+1

      !> ---------------- allocate prolongation matrix including row and column pointers
      ALLOCATE (SLF%P(SLF%NPE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P', IERR)
      SLF%P = 0.0_EB

      ALLOCATE (SLF%P_COL(SLF%NPE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_COL', IERR)
      SLF%P_COL = 0.0_EB

      ALLOCATE (SLF%P_ROW(SLF%NCE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_ROW', IERR)
      SLF%P_ROW = 0.0_EB

      !> ---------------- allocate auxiliary arrays to mark relevant positions in A and P
      ALLOCATE (SLF%A_TAG(SLF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'A_TAG', IERR)
      SLF%A_TAG = 0

      ALLOCATE (SLF%P_TAG(SLF%NCCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_TAG', IERR)
      SLF%P_TAG = 0


      !> ---------------- allocate restriction matrix including row and column pointers
      ALLOCATE (SLF%R(SLF%NRE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R', IERR)
      SLF%R = 0.0_EB

      ALLOCATE (SLF%R_COL(SLF%NRE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R_COL', IERR)
      SLF%R_COL = 0.0_EB

      ALLOCATE (SLF%R_ROW(SLF%NCCE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R_ROW', IERR)
      SLF%R_ROW = 0.0_EB


      OTHER_MESH_LOOP2: DO NOM = 1, NMESHES

         IF (NOM == NM) CYCLE OTHER_MESH_LOOP2
         OSF  => SCARC(NM)%OSCARC(NOM)

         IF (OSF%NICMAX_S==0 .AND. OSF%NICMAX_R==0) CYCLE OTHER_MESH_LOOP2
         OSLF  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

         ALLOCATE (OSLF%P(OSLF%NP+10), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P', IERR)
         OSLF%P = 0.0_EB

         ALLOCATE (OSLF%P_COL(OSLF%NP+10), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_COL', IERR)
         OSLF%P_COL = 0.0_EB

         ALLOCATE (OSLF%P_ROW(OSLF%NC+10), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_ROW', IERR)
         OSLF%P_ROW = 0.0_EB

         ALLOCATE (OSLF%R(OSLF%NP+10), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R', IERR)
         OSLF%P = 0.0_EB

         ALLOCATE (OSLF%R_COL(OSLF%NP+10), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_AMG_CnARSENING', 'R_COL', IERR)
         OSLF%P_COL = 0.0_EB

         ALLOCATE (OSLF%R_ROW(OSLF%NC+10), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R_ROW', IERR)
         OSLF%P_ROW = 0.0_EB

      ENDDO OTHER_MESH_LOOP2
   ENDDO

   !!
   !! determine prolongation and restriction matrix
   !! set corresponding overlap information between neighboring meshes for coarser grid level
   !!
   IF (TYPE_COARSENING >= NSCARC_COARSENING_GMG) THEN
      CALL SCARC_SETUP_RESTRICTION(NL)
      CALL SCARC_SETUP_PROLONGATION(NL)
      CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_GMG, NL)
   ELSE
      CALL SCARC_SETUP_PROLONGATION(NL)
      CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_PROLONGATION, NL)
      CALL SCARC_MATRIX_TRANSPOSE(SLF%P, SLF%P_ROW, SLF%P_COL, SLF%R, SLF%R_ROW, SLF%R_COL, &
                                  SLF%NCE, SLF%NCCE )
   ENDIF

   !CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL 002')
   !CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_RESTRICTION, NL, 'SETUP_RESTRICTION', 'FINAL REST 002')
   !CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL+1, 'SETUP_SYSTEM', 'WALLINFO 002')
   !CALL SCARC_LATEX_GRID(NSCARC_LATEX_STAGGERED,NL)
   !CALL SCARC_LATEX_GRID(NSCARC_LATEX_EQUAL,NL)
   !CALL SCARC_LATEX_GRID(NSCARC_LATEX_NUMBER,NL)
   !DO NM = 1, NMESHES
   !   IF (PROCESS(NM) /= MYID) CYCLE
   !   SLF => SCARC(NM)%LEVEL(NL)
   !   CALL SCARC_PRINT_MATRIX2(SLF%A, SLF%A_ROW, SLF%A_COL, SLF%NC  , SLF%NC  , NM, NL, 'A')
   !   CALL SCARC_PRINT_MATRIX2(SLF%R, SLF%R_ROW, SLF%R_COL, SLF%NCCE, SLF%NCE , NM, NL, 'R')
   !   CALL SCARC_PRINT_MATRIX2(SLF%P, SLF%P_ROW, SLF%P_COL, SLF%NCE , SLF%NCCE, NM, NL, 'P')
   !   CALL SCARC_MATLAB_MATRIX(SLF%A, SLF%A_ROW, SLF%A_COL, SLF%NC  , SLF%NC  , NM, NL, 'A')
   !   CALL SCARC_MATLAB_MATRIX(SLF%R, SLF%R_ROW, SLF%R_COL, SLF%NCCE, SLF%NCE , NM, NL, 'R')
   !   CALL SCARC_MATLAB_MATRIX(SLF%P, SLF%P_ROW, SLF%P_COL, SLF%NCE , SLF%NCCE, NM, NL, 'P')
   !ENDDO


   !! -----------------------------------------------------------------------------------------
   !! Allocate coarse grid matrix including pointer arrays
   !! Note: number of cells on coarse level corresponds to number of c-points on fine level
   !! Compute coarse grid matrix by multiplication with restriction and prolongation matrix:
   !!  A_coarse := R * A_fine * P
   !! -----------------------------------------------------------------------------------------
   DO NM = 1, NMESHES

      IF (PROCESS(NM) /= MYID) CYCLE

      SLF => SCARC(NM)%LEVEL(NL)               !> Pointer to fine level
      SLC => SCARC(NM)%LEVEL(NL+1)             !> Pointer to coarse level

      SLC%NC  = SLF%NCC
      SLC%NCE = SLF%NCCE
      SLC%NW  = SLF%NCW

      !> Allocate WALL-information for coarser grids (own and neighboring grids)
      ALLOCATE(SLC%ICG_TO_ICE(SLC%NC+1:SLC%NCE), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSENING','ICG_TO_ICE',IERR)
      SLC%ICG_TO_ICE = 0

      IF (TYPE_COARSENING >= NSCARC_COARSENING_GMG) THEN
         ALLOCATE(SLC%XCOR(0:SLC%NX), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_CELL_TYPE','SLC%XCOR',IERR)

         ALLOCATE(SLC%YCOR(0:SLC%NY), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_CELL_TYPE','SLC%YCOR',IERR)

         ALLOCATE(SLC%ZCOR(0:SLC%NZ), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_CELL_TYPE','SLC%ZCOR',IERR)

         CALL SCARC_SETUP_COORDINATES_AMG(NM, NL)
      ENDIF

      !DO IC = SLF%NC+1,SLF%NCE
      !>   ICC = SLF%CELL_TYPE(IC)
      !ENDDO


   ENDDO
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL A')

!>  IF (NMESHES>1) CALL SCARC_SETUP_DECOMPOSITION_AMG (NL)          !HIER NOCHMAL CHECKEN ----

   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL B')

   CALL SCARC_SETUP_SUBDIVISION_AMG(NL)                       !HIER NOCHMAL CHECKEN ----

   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL C')
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_SUBDIVISION  , NL+1, 'SETUP_SYSTEM', 'SUBDIVISION new level')
   !CALL SCARC_TRANSFER_MATRIX (NL)
   CALL SCARC_SETUP_SYSTEM_AMG (NL)

   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SUBDIAG, NL+1)

   !IF (NMESHES > 1) THEN
   !!   IF (TYPE_COARSENING == NSCARC_COARSENING_GMG)  THEN
   !!      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SYSTEM, NL+1)
   !!   ENDIF
   !ENDIF

   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NL+1, 'SETUP_COARSENING', 'AMG_MATRIX_FINAL')

!>! IMPORTANT : Deallocate auxiliary arrays ----!!!
!>! IMPORTANT : Deallocate auxiliary arrays ----!!!
!>! IMPORTANT : Deallocate auxiliary arrays ----!!!
!>! IMPORTANT : Deallocate auxiliary arrays ----!!!
!>! IMPORTANT : Deallocate auxiliary arrays ----!!!
!>! IMPORTANT : Deallocate auxiliary arrays ----!!!
!>! IMPORTANT : Deallocate auxiliary arrays ----!!!

   IF (SLC%NC <= 4) THEN
      NLEVEL_MAX = NL + 1
      EXIT LEVEL_LOOP
   ENDIF

ENDDO LEVEL_LOOP

!DO NL = NLEVEL_MIN, NLEVEL_MAX
!>  DO NM = 1, NMESHES
!>     IF (PROCESS(NM) /= MYID) CYCLE
!>     SLF => SCARC(NM)%LEVEL(NL)
!>     DEALLOCATE(SLF%MEASURE)
!>     DEALLOCATE(SLF%CELL_TYPE)
!>  ENDDO
!ENDDO

END SUBROUTINE SCARC_SETUP_COARSENING



!> ------------------------------------------------------------------------------------------------
!> Resize an existing array whose final size wasn't clear at the point of first allocation
!> ------------------------------------------------------------------------------------------------
!SUBROUTINE RESIZE_ARRAY(A, newSize)
!IMPLICIT NONE
!INTEGER, DIMENSION(:), INTENT(INOUT) :: A
!INTEGER, INTENT(IN) :: newSize
!INTEGER, DIMENSION(:), ALLOCATABLE :: B
!ALLOCATE(B(LBOUND(A):UBOUND(A))
!B = A
!DEALLOCATE(A)
!ALLOCATE(A(newSize))
!A(LBOUND(B):UBOUND(B)) = B
!DEALLOCATE(B)
!END SUBROUTINE RESIZE_ARRAY


!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRINT_MATRIX(A, A_ROW, A_COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , DIMENSION(:), INTENT(IN) :: A_COL
INTEGER , INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(*), INTENT(IN):: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60):: MATRIX
REAL(EB):: MATRIX_LINE(1000)

IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN

WRITE (MATRIX, '(A,A,A,A,A,i2.2,A,i2.2,A)') 'matrix/',TRIM(CHID),'_',TRIM(CNAME),'_M',NM,'_L',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=MATRIX)

DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= A_ROW(IC), A_ROW(IC+1)-1
         IF (A_COL(ICOL)==JC) MATRIX_LINE(JC)=A(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX,1001) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO

CLOSE(MMATRIX)

1001 FORMAT(40(F6.0,','))
1000 FORMAT(F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',&
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,';')
END SUBROUTINE SCARC_PRINT_MATRIX

!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRINT_MATRIX2(A, A_ROW, A_COL, NC1, NC2, NX, NY, NZ, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , DIMENSION(:), INTENT(IN) :: A_COL
INTEGER , INTENT(IN) :: NM, NL, NC1, NC2, NX, NY, NZ
CHARACTER(*), INTENT(IN):: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX, MAX_NC, IDUMMY
CHARACTER(60):: MATRIX
CHARACTER(2), DIMENSION(64):: LINE
REAL(EB):: MATRIX_LINE(1000)

IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN

WRITE (MATRIX, '(A,A,A,A,A,i2.2,A,i2.2,A)') 'matrix2/',TRIM(CHID),'_',TRIM(CNAME),'_M',NM,'_L',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=MATRIX)

MAX_NC = 64
IDUMMY = NY
IDUMMY = NZ

DO IC = 1, MIN(NC1,MAX_NC)
   MATRIX_LINE=0.0_EB
   LINE = ' .'
   DO JC = 1, MIN(NC2,MAX_NC)
      DO  ICOL= A_ROW(IC), A_ROW(IC+1)-1
         IF (A_COL(ICOL)==JC) THEN
            MATRIX_LINE(JC)=A(ICOL)/100.0
            WRITE(LINE(JC),'(I2)') NINT(MATRIX_LINE(JC))
         ENDIF
      ENDDO
   ENDDO
   !IF (NX==8) THEN
   !   WRITE(MMATRIX,1008) (NINT(MATRIX_LINE(JC)),JC=1,MIN(NC2,MAX_NC))
   !ELSE IF (NX==4) THEN
   !   WRITE(MMATRIX,1004) (NINT(MATRIX_LINE(JC)),JC=1,MIN(NC2,MAX_NC))
   !ELSE
   !   WRITE(MMATRIX,1000) (NINT(MATRIX_LINE(JC)),JC=1,MIN(NC2,MAX_NC))
   !ENDIF
   IF (NX==8) THEN
      WRITE(MMATRIX,2008) (LINE(JC),JC=1,MIN(NC2,MAX_NC))
   ELSE IF (NX==4) THEN
      WRITE(MMATRIX,2004) (LINE(JC),JC=1,MIN(NC2,MAX_NC))
   ELSE IF (NX==6) THEN
      WRITE(MMATRIX,2006) (LINE(JC),JC=1,MIN(NC2,MAX_NC))
   ELSE IF (NX==3) THEN
      WRITE(MMATRIX,2003) (LINE(JC),JC=1,MIN(NC2,MAX_NC))
   ELSE
      WRITE(MMATRIX,2000) (LINE(JC),JC=1,MIN(NC2,MAX_NC))
   ENDIF
   IF (MOD(IC, NX) == 0) WRITE(MMATRIX,*)
ENDDO

CLOSE(MMATRIX)

1000 FORMAT(40i3)
1004 FORMAT(4i3, 9(i5, 3i3))
1008 FORMAT(8i3, 9(i5, 7i3))
2000 FORMAT(40a3)
2003 FORMAT(3a3, 9(a5, 2a3))
2004 FORMAT(4a3, 9(a5, 3a3))
2006 FORMAT(6a3, 9(a5, 5a3))
2008 FORMAT(8a3, 9(a5, 7a3))
END SUBROUTINE SCARC_PRINT_MATRIX2


!> ----------------------------------------------------------------------------------------------------
!> Store subdivision information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COORDINATES_AMG(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER IXC, IYC, IZC
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLC, SLF

SLF => SCARC(NM)%LEVEL(NL)
SLC => SCARC(NM)%LEVEL(NL+1)

SELECT CASE (TYPE_DIMENSION)
   CASE (NSCARC_DIMENSION_TWO)

      DO IXC = 0, SLC%NX-1
         SLC%XCOR(IXC)  = SLF%XCOR(2*IXC)
      ENDDO
      SLC%XCOR(IXC)  = SLF%XCOR(SLF%NX)

      DO IZC = 0, SLC%NZ-1
         SLC%ZCOR(IZC)  = SLF%ZCOR(2*IZC)
      ENDDO
      SLC%ZCOR(IZC)  = SLF%ZCOR(SLF%NZ)

      !WRITE(LU_SCARC,'(a,8f8.2)') 'SLC%XCOR:', (SLC%XCOR(IXC), IXC=0,SLC%NX)
      !WRITE(LU_SCARC,'(a,8f8.2)') 'SLC%ZCOR:', (SLC%ZCOR(IZC), IZC=0,SLC%NZ)

   CASE (NSCARC_DIMENSION_THREE)

      DO IZC = 1, SLC%NZ
         DO IYC = 1, SLC%NY
            DO IXC = 1, SLC%NX
            ENDDO
         ENDDO
      ENDDO

END SELECT

END SUBROUTINE SCARC_SETUP_COORDINATES_AMG


!> ----------------------------------------------------------------------------------------------------
!> Store subdivision information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIVISION_AMG(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IOR_LAST, NOM, INBR
INTEGER :: NEIGHBORS(20,-3:3)
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLC

SCARC_ROUTINE = 'SCARC_SETUP_SUBDIVISION'

IOR_LAST    = 0
NEIGHBORS   = 0

MESHES_LOOP1: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   SLC => SCARC(NM)%LEVEL(NL+1)

   SLC%SUBDIVISION = 0

   WALL_CELLS_LOOP: DO IW = 1, SLC%NW

      IOR0 = SLC%WALL(IW)%IOR

      IF (IOR_LAST /= IOR0) SLC%SUBDIVISION(1,IOR0) = IW
      SLC%SUBDIVISION(2,IOR0) = SLC%SUBDIVISION(2,IOR0) + 1

      NOM= SLC%WALL(IW)%NOM

      IF (NOM /= 0) THEN
         NEIGHBOR_LOOP: DO INBR = 1, 20
            IF (NOM == NEIGHBORS(INBR, IOR0)) THEN
               EXIT NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) /= 0) THEN
               CYCLE NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) == 0) THEN
               NEIGHBORS(INBR, IOR0) = NOM
               SLC%SUBDIVISION(3,IOR0) = SLC%SUBDIVISION(3,IOR0) + 1
               EXIT NEIGHBOR_LOOP
            ELSE
               WRITE(SCARC_MESSAGE,'(2A)') TRIM(SCARC_ROUTINE),': More than 20 neighbors at one face not allowed yet!'
               CALL SHUTDOWN(SCARC_MESSAGE); RETURN
            ENDIF
         ENDDO NEIGHBOR_LOOP
      ENDIF

      IOR_LAST = IOR0

   ENDDO WALL_CELLS_LOOP

   !WRITE(LU_SCARC,*) 'SETTING SUBDIVISION for mesh ', NM, ' on level ', NL
   !WRITE(LU_SCARC,'(a,7i8)') 'SUBDIVISION (1,.): ', (SLC%SUBDIVISION(1, IOR0), IOR0=-3,3)
   !WRITE(LU_SCARC,'(a,7i8)') 'SUBDIVISION (2,.): ', (SLC%SUBDIVISION(2, IOR0), IOR0=-3,3)
   !WRITE(LU_SCARC,'(a,7i8)') 'SUBDIVISION (3,.): ', (SLC%SUBDIVISION(3, IOR0), IOR0=-3,3)

ENDDO MESHES_LOOP1

END SUBROUTINE SCARC_SETUP_SUBDIVISION_AMG



!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DECOMPOSITION_AMG(NL)
INTEGER , INTENT(IN) :: NL
INTEGER :: NM, NOM, ICOL, IW, ICW, NCPL, JC, IERR
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SLC
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSLC


!> -------------------------------------------------------------------------
!> Loop over all boundary cells IW of fine grid
!> Get corresponding adjacent and ghost cell
!> -------------------------------------------------------------------------
IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) THEN
MESHES_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   SLC => SCARC(NM)%LEVEL(NL+1)

   DO IW = 1, SLC%NW

      NOM = SLC%WALL(IW)%NOM
      OSLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

      ICW = SLC%WALL(IW)%ICW

      IF (TYPE_COARSENING == NSCARC_COARSENING_GMG) NCPL = 1

      WRITE(*,*) 'L: HIER WALL%ICN �ndern auf OSLF!'
      WRITE(*,*) 'L: MUSS NOCH UEBERARBEITET WERDEN !!'

!>     SLC%WALL(IW)%NCPL = NCPL

!>     ALLOCATE(SLC%WALL(IW)%ICN(NCPL), STAT=IERR)
!>     CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION_AMG','ICN',IERR)

!>     ALLOCATE(SLC%WALL(IW)%ICE(NCPL), STAT=IERR)
!>     CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION_AMG','ICE',IERR)

!>     ALLOCATE(SLC%WALL(IW)%ICG(NCPL), STAT=IERR)
!>     CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION_AMG','ICG',IERR)

      NCPL = 0
      DO ICOL = SLC%A_ROW(ICW)+1, SLC%A_ROW(ICW+1)-1
         JC = SLC%A_COL(ICOL)
         IF (JC > SLC%NC) THEN
            NCPL = NCPL + 1
!>           SLC%WALL(IW)%ICE(NCPL) = JC
!>           SLC%WALL(IW)%ICN(NCPL) = SLC%EXT_PTR(JC)
         ENDIF
      ENDDO

    ENDDO

ENDDO MESHES_LOOP

ENDIF

!CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL+1, 'SETUP_INTERFACES_AMG', 'WALL pos 0 ')

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MESHINFO, NL+1)

MESHES_LOOP2: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
      OS  => SCARC(NM)%OSCARC(NOM)

      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE OTHER_MESHES_LOOP
      OSLC  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

      ALLOCATE(OSLC%WALL(1:OSLC%NW), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_DECOMPOSITION_AMG','OSLC%WALL',IERR)

   ENDDO OTHER_MESHES_LOOP
ENDDO MESHES_LOOP2

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL+1, 'SETUP_INTERFACES_AMG', 'WALL pos 1 ')

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WALLINFO, NL+1)

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL+1, 'SETUP_INTERFACES_AMG', 'WALL pos 2')

END SUBROUTINE SCARC_SETUP_DECOMPOSITION_AMG


!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATLAB_MATRIX(A, A_ROW, A_COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , DIMENSION(:), INTENT(IN) :: A_COL
INTEGER , INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(1), INTENT(IN):: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60):: MATRIX
REAL(EB):: MATRIX_LINE(1000)

IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN

SCARC_ROUTINE = 'SCARC_MATLAB_MATRIX'

WRITE (MATRIX, '(A,A1,A,i2.2,A,i2.2,A)') 'matlab/',CNAME,'_mesh',NM,'_level',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=MATRIX)

!WRITE(LU_SCARC,*) 'PRINTING MATLAB INFORMATION FOR LEVEL ', NL
!WRITE(LU_SCARC,*) 'NC1  =',NC1
!WRITE(LU_SCARC,*) 'NC2  =',NC2
!WRITE(LU_SCARC,*) 'CNAME=',CNAME

WRITE(MMATRIX, 1000) CNAME
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= A_ROW(IC), A_ROW(IC+1)-1
         IF (A_COL(ICOL)==JC) MATRIX_LINE(JC)=A(ICOL)
      ENDDO
   ENDDO
   SELECT CASE (NC2)
      CASE (4)
         WRITE(MMATRIX,1004) (MATRIX_LINE(JC),JC=1,NC2)
      CASE (8)
         WRITE(MMATRIX,1008) (MATRIX_LINE(JC),JC=1,NC2)
      CASE (12)
         WRITE(MMATRIX,1012) (MATRIX_LINE(JC),JC=1,NC2)
      CASE (16)
         WRITE(MMATRIX,1016) (MATRIX_LINE(JC),JC=1,NC2)
      CASE (24)
         WRITE(MMATRIX,1024) (MATRIX_LINE(JC),JC=1,NC2)
      CASE (32)
         WRITE(MMATRIX,1032) (MATRIX_LINE(JC),JC=1,NC2)
      CASE (64)
         WRITE(MMATRIX,1064) (MATRIX_LINE(JC),JC=1,NC2)
      CASE DEFAULT
         WRITE(SCARC_MESSAGE,'(2A,I8)') TRIM(SCARC_ROUTINE),': Wrong value for ', NC2
         CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   END SELECT
ENDDO
WRITE(MMATRIX, 1100)

CLOSE(MMATRIX)

1000 FORMAT(a,' = [ ')
1004 FORMAT( 3(f24.16,','),f24.16,';')
1008 FORMAT( 7(f24.16,','),f24.16,';')
1012 FORMAT(11(f24.16,','),f24.16,';')
1016 FORMAT(15(f24.16,','),f24.16,';')
1024 FORMAT(23(f24.16,','),f24.16,';')
1032 FORMAT(31(f24.16,','),f24.16,';')
1064 FORMAT(63(f24.16,','),f24.16,';')
1104 FORMAT(i4,' :', 4(f8.2))
1108 FORMAT(i4,' :', 8(f8.2))
1112 FORMAT(i4,' :',12(f8.2))
1116 FORMAT(i4,' :',16(f8.2))
1124 FORMAT(i4,' :',24(f8.2))
1132 FORMAT(i4,' :',32(f8.2))
1164 FORMAT(i4,' :',64(f8.2))
1100 FORMAT(' ] ')
END SUBROUTINE SCARC_MATLAB_MATRIX



!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_LATEX_NOM(MTABLE, NM, NOM, NL)
INTEGER, INTENT(IN):: NM, NOM, NL, MTABLE
INTEGER :: IG
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSL

IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN

OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

!WRITE (CTABLE, '(A,A,A,i2.2,A,I2.2,A,i2.2,A)') 'tables/',CNAME,'_mesh',NM,'_nom',nom,'_level',NL,'.tex'
!MTABLE=GET_FILE_NUMBER()
!OPEN(MTABLE,FILE=CTABLE)

!WRITE(MTABLE,100) NM, NOM

!> NOM: ICE_TO_IWG
WRITE(MTABLE, 500) NM,':WALL\_PTR(',NOM,'), GHOST\_PTR(',NOM,'), NOM\_PTR(',NOM,')'
WRITE(LU_SCARC, 500) NM,':WALL\_PTR(',NOM,'), GHOST\_PTR(',NOM,'), NOM\_PTR(',NOM,')'
WRITE(MTABLE,1000)
WRITE(MTABLE,1500) NM,NOM
WRITE(MTABLE,2000)
WRITE(MTABLE,3000) 'IG', 'IW','ICE','IW'
WRITE(LU_SCARC,3000) 'IG', 'IW','ICE','IW'
WRITE(MTABLE,2000)
DO IG = 1, OSL%NCG
   WRITE(MTABLE,3001) IG, OSL%ICG_TO_IWG(IG), OSL%ICG_TO_ICE(IG), OSL%ICN_TO_ICE(IG)
   WRITE(LU_SCARC,3001) IG, OSL%ICG_TO_IWG(IG), OSL%ICG_TO_ICE(IG), OSL%ICN_TO_ICE(IG)
   !WRITE(MTABLE,2000)
ENDDO
WRITE(MTABLE,2000)
WRITE(MTABLE,4000)
WRITE(MTABLE,5000)



100  FORMAT( /,'\vspace{5mm}',/,'{\bf MESH',i3,': NOM ',i3,' }\\', /)
500  FORMAT( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/, &
             '%%%    ', i3, a,i3,a,i3,a,i3,a/,&
             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
1000 FORMAT( '\begin{minipage}{4cm}')
1500 FORMAT( '\flushleft{\bf {\footnotesize Mesh ',i3,' : NOM',i3,':}}\\[2ex]',/,&
             '\begin{tabular}{|c|c|c|c|}')
2000 FORMAT( '\hline')
3000 FORMAT ('{\footnotesize ',a,'} & {\footnotesize ',a, '} & {\footnotesize ',a, '} & {\footnotesize ',a,'} \\')
3001 FORMAT (i4,' &',i4, ' &',i4,' &',i4,'\\')
4000 FORMAT( '\end{tabular}')
5000 FORMAT( '\end{minipage}')
END SUBROUTINE SCARC_LATEX_NOM



!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
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

!> ------------------------------------------------------------------------------------------------
!> Determine measure of cells corresponding to requested coarsening type
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STRENGTH_MATRIX(NTYPE, NL)
INTEGER, PARAMETER  :: NCOUPLINGS = 20
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER  :: NM, IC, JC, ICOL, IDIAG, IS , IPOS
REAL(EB) :: ADIAG, ACOL, ROW_SCALE, ROW_SUM, MAX_ROW_SUM = 0.9_EB, THRESHOLD = 0.25_EB
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

!> Only dummy (NTYPE really used ?)
IC = NTYPE

STRENGTH_MESHES_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   !WRITE(LU_SCARC,*) 'SETUP_STRENGTH_MATRIX: ============== '

   SL => SCARC(NM)%LEVEL(NL)

   STRENGTH_CELL_LOOP: DO IC = 1, SL%NC

      !WRITE(LU_SCARC,*) '==============   IC =', IC, '==== NTYPE=',NTYPE

      IDIAG = SL%A_ROW(IC)
      ADIAG = SL%A(IDIAG)

      !> get row sum and scaling factor
      ROW_SCALE = 0.0_EB
      ROW_SUM   = ADIAG

      !WRITE(LU_SCARC,'(a,f12.6)') 'ADIAG    =',ADIAG
      !WRITE(LU_SCARC,'(a,f12.6)') 'ROW_SCALE=',ROW_SCALE
      !WRITE(LU_SCARC,'(a,f12.6)') 'ROW_SUM  =',ROW_SUM

      DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
         JC = SL%A_COL(ICOL)
         IF (JC <= SL%NC) THEN
            ACOL = SL%A(ICOL)
            ROW_SCALE = MAX(ROW_SCALE, ACOL)
            ROW_SUM   = ROW_SUM + ACOL
            !WRITE(LU_SCARC,'(a,i4,a,f12.6,a,f12.6, a, f12.6)') &
            !     'ICOL=',ICOL,': ACOL=',ACOL,': ROW_SCALE=',ROW_SCALE, ': ROW_SUM=',ROW_SUM
         ENDIF
      ENDDO

      !> get row entries of strength matrix S
      ROW_SUM = ABS(ROW_SUM/ADIAG)
      SL%S_COL(IDIAG) = -1

      !WRITE(LU_SCARC,'(a,f12.6)') 'ROW_SUM=',ROW_SUM
      !WRITE(LU_SCARC,'(a,i4,a,i4)') 'S_COL(',IDIAG,')=',SL%S_COL(IDIAG)

      IF ((ROW_SUM > MAX_ROW_SUM) .AND. (MAX_ROW_SUM < 1.0_EB)) THEN
         !> set all dependencies to be weak
         DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
            JC = SL%A_COL(ICOL)
            IF (JC <= SL%NC) THEN
               SL%S_COL(ICOL) = -1
               !WRITE(LU_SCARC,'(a,i4,a,i4)') 'A: S_COL(',ICOL,')=',SL%S_COL(ICOL)
            ENDIF
         ENDDO
      ELSE
         !> set dependencies to be weak related to threshold
         DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
            JC = SL%A_COL(ICOL)
            IF (JC <= SL%NC) THEN
               IF (SL%A(ICOL) <= THRESHOLD * ROW_SCALE) SL%S_COL(ICOL) = -1
               !WRITE(LU_SCARC,'(a,i4,a,f12.6,a,f12.6,a,i3,a,i4)') &
               !  'ZZZZ2: A(',ICOL,')=',SL%A(ICOL), ': T*R_S=',THRESHOLD*ROW_SCALE,' S_COL(',ICOL,')=',SL%S_COL(ICOL)
            !ELSE
            !>   SL%S_COL(ICOL) = -1
            ENDIF
          ENDDO
      ENDIF

   ENDDO STRENGTH_CELL_LOOP

   !WRITE(LU_SCARC,*) 'A: ---------------------- S_ROW:', SL%NC
   !WRITE(LU_SCARC,'(4i9)') (SL%S_ROW(IC), IC=1,SL%NC+1)
   !WRITE(LU_SCARC,*) '---------------------- S_COL:'
   !DO IC = 1, SL%NC
   !   WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%S_COL(IP),IP=SL%S_ROW(IC),SL%S_ROW(IC+1)-1)
   !ENDDO
   !WRITE(LU_SCARC,*) '---------------------- S:'
   !DO IC = 1, SL%NC
   !   WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SL%S(IP),IP=SL%S_ROW(IC),SL%S_ROW(IC+1)-1)
   !ENDDO

   !! Compress strength matrix
   IS = 1
   STRENGTH_CELL_LOOP2: DO IC = 1, SL%NC
      SL%S_ROW(IC) = IS
      DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
         !IF (SL%S_COL(ICOL) > -1.AND.SL%S_COL(ICOL)<=SL%NC) THEN
         IF (SL%S_COL(ICOL) > -1) THEN
            SL%S_COL(IS) = SL%S_COL(ICOL)
            !WRITE(LU_SCARC,*) 'S_COL(',IS,')=',SL%S_COL(ICOL)
            IS = IS + 1
         ENDIF
      ENDDO
   ENDDO STRENGTH_CELL_LOOP2
   SL%S_ROW(SL%NC+1) = IS

   !WRITE(LU_SCARC,*) 'B: ---------------------- S_ROW:', SL%NC
   !WRITE(LU_SCARC,'(4i9)') (SL%S_ROW(IC), IC=1,SL%NC+1)
   !WRITE(LU_SCARC,*) '---------------------- S_COL:'
   !DO IC = 1, SL%NC
   !   WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%S_COL(IP),IP=SL%S_ROW(IC),SL%S_ROW(IC+1)-1)
   !ENDDO
   !WRITE(LU_SCARC,*) 'SIZE(ST_ROW)=',SIZE(SL%ST_ROW)
   !WRITE(LU_SCARC,*) 'SIZE(ST_COL)=',SIZE(SL%ST_COL)

   DO IC = 1, SL%NCE+1
      SL%ST_ROW(IC) = 0
   ENDDO

   IS = SL%S_ROW(SL%NC+1)-1
   DO ICOL = 1, IS
      SL%ST_ROW(SL%S_COL(ICOL)+1) = SL%ST_ROW(SL%S_COL(ICOL)+1) + 1
      !WRITE(LU_SCARC,*) 'A: ST_ROW(',SL%S_COL(ICOL)+1,')=',SL%ST_ROW(SL%S_COL(ICOL)+1)
   ENDDO
   SL%ST_ROW(1) = 1

   DO IC = 1, SL%NCE
      SL%ST_ROW(IC+1)= SL%ST_ROW(IC+1) + SL%ST_ROW(IC)
      !WRITE(LU_SCARC,*) 'B: ST_ROW(',IC+1,')=',SL%ST_ROW(IC+1)
   ENDDO
   DO IC = 1, SL%NC
      DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
         IPOS = SL%S_COL(ICOL)
         SL%ST_COL(SL%ST_ROW(IPOS)) = IC
         !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
         !  'C: IC=',IC,': ICOL=',ICOL,': IPOS=', IPOS,': ST_COL(',SL%ST_ROW(IPOS),')=',IC
         SL%ST_ROW(IPOS) = SL%ST_ROW(IPOS) + 1
      ENDDO
   ENDDO
   DO IC = SL%NCE+1, 2, -1
      SL%ST_ROW(IC) = SL%ST_ROW(IC-1)
      !WRITE(LU_SCARC,*) 'D: ST_ROW(',IC,')=',SL%ST_ROW(IC)
   ENDDO
   SL%ST_ROW(1) = 1

   !WRITE(LU_SCARC,*) '---------------------- ST_ROW:', SL%NC
   !WRITE(LU_SCARC,'(4i9)') (SL%ST_ROW(IC), IC=1,SL%NCE+1)
   !WRITE(LU_SCARC,*) '---------------------- ST_COL:'
   !DO IC = 1, SL%NCE
   !   WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%ST_COL(IP),IP=SL%ST_ROW(IC),SL%ST_ROW(IC+1)-1)
   !ENDDO

ENDDO STRENGTH_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_STRENGTH_MATRIX


!> ------------------------------------------------------------------------------------------------
!> Determine measure of cells corresponding to requested coarsening type
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COLORING(NTYPE, NL)
INTEGER, PARAMETER  :: NCOUPLINGS = 20
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, IC, IC0, JC, KC, ICOL, JCOL, IRAND, REMAINING_CELLS, ICG, IW
INTEGER :: IGRAPH, IGRAPHE, IGRAPH_GLOBAL, ICYCLE, IERR
INTEGER :: FCELL = -1, ZCELL =-2, ICT2 = -1, ICT
LOGICAL :: BEMPTY=.FALSE., BNONEMPTY=.FALSE., BNEIGHBOR, BREAK = .TRUE., BCYCLE = .TRUE.
REAL(EB) :: RAND_NUM, MEASURE, NEW_MEASURE
REAL(EB) :: MEASURE_MAX, EPS
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

EPS = 1.0E-12
MEASURE_MAX = 0.0_EB


!> Select coarsening strategy
SELECT CASE (NTYPE)

   !! ---------------------------------------------------------------------------------------------
   !! RS3-coarsening:
   !!      - Original Ruge-Stuben method with parallel postprocessing
   !!      - Produces good C/F splittings but is inherently serial.
   !!      - May produce AMG hierarchies with relatively high operator complexities.
   !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_RS3)

      RS3_MESH_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SL => SCARC(NM)%LEVEL(NL)

         REMAINING_CELLS = 0

         !WRITE(LU_SCARC,*) 'SCARC_SETUP_COLORING: INIT: NL=',NL
         !WRITE(LU_SCARC,*) 'SCARC_SETUP_COLORING: CELL_TYPE: ', SL%NC
         !WRITE(LU_SCARC,*) (SL%CELL_TYPE(IC), IC=1, SL%NC)

         !> Currently the measures are computed as row sums of ST (number of influences for IC)
         RS3_MEASURE_LOOP0: DO IC = 1, SL%NC
            SL%MEASURE(IC) = SL%ST_ROW(IC+1)-SL%ST_ROW(IC)
         ENDDO RS3_MEASURE_LOOP0


         !> Subdivide in fine and coarse grid cells
         RS3_MEASURE_LOOP1: DO IC = 1, SL%NC

            IF (SL%S_ROW(IC+1)-SL%S_ROW(IC) == 0) THEN
               SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE0
               SL%MEASURE(IC) = 0.0_EB

               !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
               !     'A: CELL_TYPE(',IC,') =',SL%CELL_TYPE(IC),': MEASURE(',IC,')=',&
               !      SL%MEASURE(IC),': Remaining cells=',REMAINING_CELLS

             !> IF (AGGRESSIVE2) SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE0
            ELSE
               SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
               REMAINING_CELLS   = REMAINING_CELLS + 1

               !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
               !    'B: CELL_TYPE(',IC,') =',SL%CELL_TYPE(IC),': MEASURE(',IC,')=',&
               !     SL%MEASURE(IC),': Remaining cells=',REMAINING_CELLS

            ENDIF
         ENDDO RS3_MEASURE_LOOP1

         RS3_MEASURE_LOOP2: DO IC = 1, SL%NC
            MEASURE = SL%MEASURE(IC)
            IF (SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_FINE0 .AND. SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_COARSE0) THEN
               IF (SL%MEASURE(IC) > 0.0_EB) THEN

                  !SL%MEASURE(IC) = SL%MEASURE(IC)
                  !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
                  !    'C: CELL_TYPE(',IC,') =',SL%CELL_TYPE(IC),': MEASURE(',IC,')=',&
                  !     SL%MEASURE(IC),': Remaining cells=',REMAINING_CELLS

               ELSE
                  IF (SL%MEASURE(IC) < 0.0_EB) WRITE(*,*) 'SCARC_SETUP_MEASURE: Negative measure !!'
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                     JC = SL%S_COL(ICOL)
                     IF (SL%CELL_TYPE(JC) /= NSCARC_CELL_TYPE_COARSE0 .AND. SL%CELL_TYPE(JC) /= NSCARC_CELL_TYPE_FINE0) THEN
                        IF (JC < IC) THEN
                           NEW_MEASURE = SL%MEASURE(JC)
                           IF (NEW_MEASURE > 0.0_EB) SL%MEASURE(JC) = 0.0_EB
                           NEW_MEASURE = SL%MEASURE(JC)+1
                           SL%MEASURE(JC) = NEW_MEASURE

                           !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
                           !    'D: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',&
                           !     SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS

                        ELSE
                           NEW_MEASURE = SL%MEASURE(JC)+1

                           !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
                           !    'E: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',&
                           !     SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS

                        ENDIF
                     ENDIF
                  ENDDO
                  REMAINING_CELLS = REMAINING_CELLS - 1
               ENDIF
            ENDIF
         ENDDO RS3_MEASURE_LOOP2


CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_COLORING', 'MEASURE AFTER FIRST RS3 LOOP ')

         RS3_CYCLE_LOOP: DO WHILE (REMAINING_CELLS > 0)


            !> get maximum (remaining) measure for all cells
            MEASURE_MAX = MAXVAL(SL%MEASURE(1:SL%NC))
            IF (MEASURE_MAX <= EPS) EXIT RS3_CYCLE_LOOP

            !WRITE(LU_SCARC,*) 'MEASURE_MAX=',MEASURE_MAX

            RS3_CELL_LOOP: DO IC = 1, SL%NC

             !> Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, SL%MEASURE(IC))) THEN

                  SL%MEASURE(IC)  = NSCARC_MEASURE_NONE
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
                  REMAINING_CELLS = REMAINING_CELLS - 1

                  !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
                  !'F: CELL_TYPE(',IC,') =',SL%CELL_TYPE(IC),': MEASURE(',IC,')=',&
                  !SL%MEASURE(IC),': Remaining cells=',REMAINING_CELLS

                  !> Determine set of fine cells
                  DO ICOL = SL%ST_ROW(IC), SL%ST_ROW(IC+1)-1

                   !> IF JC hasn't been marked yet, set it to be a fine cell which is no longer measured
                     JC = SL%ST_COL(ICOL)

                     IF (JC > SL%NC) CYCLE
                     IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) THEN

                        SL%MEASURE(JC)  = NSCARC_MEASURE_NONE
                        SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_FINE
                        REMAINING_CELLS = REMAINING_CELLS - 1

                        !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
                        !'G: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',&
                        !SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS

                        !>  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                        DO JCOL = SL%S_ROW(JC), SL%S_ROW(JC+1)-1
                           KC = SL%S_COL(JCOL)
                           IF (SL%CELL_TYPE(KC)==NSCARC_CELL_TYPE_NONE) THEN
                              SL%MEASURE(KC) = SL%MEASURE(KC) + 1.0_EB

                              !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
                              !'H: CELL_TYPE(',KC,') =',SL%CELL_TYPE(KC),': MEASURE(',KC,')=',&
                              !SL%MEASURE(KC),': Remaining cells=',REMAINING_CELLS

                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO

!WRITE(LU_SCARC,*) '====================== REMAINING_CELLS = ', REMAINING_CELLS
!CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELL_TYPE', 'RS3_MEASURE ')
!CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPE', 'RS3_CELL_TYPE ')
!WRITE(LU_SCARC,*) '====================== REMAINING_CELLS = ', REMAINING_CELLS
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_COLORING', 'MEASURE AFTER SECOND RS3 LOOP ')
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE AFTER SECOND RS3 LOOP ')

                  EXIT RS3_CELL_LOOP
               ENDIF
            ENDDO RS3_CELL_LOOP

            DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
               JC = SL%S_COL(ICOL)
               !WRITE(LU_SCARC,*) 'ICOL=',ICOL, ': JC=',JC
               IF (JC > SL%NC) CYCLE
               !WRITE(LU_SCARC,*) 'noch da ...', SL%CELL_TYPE(JC)
               IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) THEN
               !WRITE(LU_SCARC,*) 'immer noch da ...'
                  MEASURE = SL%MEASURE(JC) - 1
                  SL%MEASURE(JC) = MEASURE
                  IF (MEASURE > 0.0_EB) THEN
                     SL%MEASURE(JC) = MEASURE

                     !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
                     !'I: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',&
                     !SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS

                  ELSE
                     SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_FINE
                     REMAINING_CELLS = REMAINING_CELLS - 1

                     !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
                     !'J: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',&
                     !SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS

                     DO JCOL = SL%S_ROW(JC), SL%S_ROW(JC+1)-1
                        KC = SL%S_COL(JCOL)
                        IF (SL%CELL_TYPE(KC)==NSCARC_CELL_TYPE_NONE) THEN
                           SL%MEASURE(KC) = SL%MEASURE(KC) + 1.0_EB
                           MEASURE_MAX = MAX(MEASURE_MAX, SL%MEASURE(KC))

                           !WRITE(LU_SCARC,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
                           !'K: CELL_TYPE(',KC,') =',SL%CELL_TYPE(KC),': MEASURE(',KC,')=',&
                           !SL%MEASURE(KC),': Remaining cells=',REMAINING_CELLS

                        ENDIF
                     ENDDO
                  ENDIF

               ENDIF
            ENDDO

CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_COLORING', 'MEASURE AFTER THIRD RS3 LOOP ')
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE AFTER THIRD RS3 LOOP ')

         ENDDO RS3_CYCLE_LOOP
         SL%NCW = 0

         DO IC = 1, SL%NC
            IF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_COARSE0) SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
         ENDDO

         !> exchange information with neighboring meshes
         IF (NMESHES > 1) THEN

            DO IC = 1, SL%NCE
               SL%GRAPH(IC) = -1
            ENDDO

            IC0 = 1
            DO IC = 1, SL%NC

               !WRITE(LU_SCARC,*) 'GRAPH::: =============== IC =',IC
               !WRITE(LU_SCARC,*) 'GRAPH:::0  CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)

               IF (ICT2 /= IC) ICT = -1

               !WRITE(LU_SCARC,*) 'GRAPH:::0  ICT=',ICT

               IF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE) THEN

                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                     JC = SL%S_COL(ICOL)
                     IF (JC <= SL%NC) THEN
                        IF (SL%CELL_TYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                           SL%GRAPH(JC) = IC
                           !WRITE(LU_SCARC,*) 'GRAPH::: GRAPH(',JC,')=',SL%GRAPH(JC)
                        ENDIF
                     ENDIF
                  ENDDO

                  !> Hier fehlt noch die Abfrage nach Nachbarmesh !!!

                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                     JC = SL%S_COL(ICOL)
                     IF (JC <= SL%NC) THEN
                     IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_FINE) THEN

                     !WRITE(LU_SCARC,*) 'GRAPH:::1  CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)

                        !> DIESER PART WIRD ERST FUER ANDERE VERFEINERUNGEN AKTIV
                        !> ACHTUNG: DANN NOCHMAL UEBERPRUEFEN!
                        BEMPTY = .TRUE.
                        DO JCOL = SL%S_ROW(JC), SL%S_ROW(JC+1)-1
                           KC = SL%S_COL(JCOL)
                           IF (SL%GRAPH(KC) == IC) THEN
                              BEMPTY = .FALSE.
                              EXIT
                           ENDIF
                           IF (BEMPTY) THEN
                              IF (BNONEMPTY) THEN
                                 SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
                                 !WRITE(LU_SCARC,*) 'GRAPH:::2  CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
                                 IF (ICT > -1) THEN
                                    SL%CELL_TYPE(ICT) = NSCARC_CELL_TYPE_FINE
                                    !WRITE(LU_SCARC,*) 'GRAPH:::3  CELL_TYPE(',ICT,')=',SL%CELL_TYPE(ICT)
                                    ICT = -1
                                 ENDIF
                                 !> Hier fehlt noch Nachbaranteil
                                 BNONEMPTY = .FALSE.
                                 EXIT
                              ELSE
                                 ICT  = JC
                                 ICT2 = IC
                                 SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_COARSE
                                 !WRITE(LU_SCARC,*) 'GRAPH::: CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
                                 BNONEMPTY = .FALSE.
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


         !DO IC = 1, SL%NC
         !   WRITE(LU_SCARC,*) 'GRAPH(',IC,')=',SL%GRAPH(IC)
         !ENDDO
         !DO IC = 1, SL%NC
         !   WRITE(LU_SCARC,*) 'CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
         !ENDDO
      ENDDO RS3_MESH_LOOP

    !> Exchange CELL_TYPE-data along internal boundaries
      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)


      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE BEFORE LAST EXCHANGE')
      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE AFTER LAST EXCHANGE')

    !> Third pass: identify fine cells along boundary and get their coarse neighbors
      RS3_MESH_LOOP2: DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         DO IC = 1, SL%NC
            SL%GRAPH(IC) = -1
         ENDDO
      ENDDO RS3_MESH_LOOP2


!> ----------------------------------------------------------------------------------------
!> Falgout coarsening
!> ----------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_FALGOUT)

      FCELL = ZCELL
      FALGOUT_MESHES_LOOP1: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SL => SCARC(NM)%LEVEL(NL)

         !> Reset the measures as row sums of ST (number of influences for IC)
         !> plus random number to make them unique
         FALGOUT_CELL_LOOP1: DO IC = 1, SL%NC
            RAND_NUM = 0.01_EB
            DO IRAND = 1, 5
               CALL RANDOM_NUMBER(RAND_NUM)
               RAND_NUM = RAND_NUM + RAND_NUM/10**(IRAND-1)
               SL%MEASURE(IC) = SL%MEASURE(IC) + RAND_NUM
            ENDDO
            SL%MEASURE(IC) = SL%S_ROW(IC+1)-SL%S_ROW(IC) + RAND_NUM
         ENDDO FALGOUT_CELL_LOOP1

      ENDDO FALGOUT_MESHES_LOOP1

      !> Initial exchange of measure array
      IF (NMESHES > 1) THEN
         TYPE_VECTOR = NSCARC_VECTOR_MEASURE
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR , NL)
      ENDIF

      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_FALGOUT', 'FIRST MEASURE PRINTOUT ')
      FALGOUT_INIT_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE
         !> reset CELL_TYPE for cells with neighbors in other meshes
         FALGOUT_INTERNAL_CELL_LOOP1: DO IC = 1, SL%NC
            BNEIGHBOR = .FALSE.
            DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
               JC = SL%S_COL(ICOL)
               IF (JC > SL%NC) THEN
                  BNEIGHBOR = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (BNEIGHBOR .OR. SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE) SL%CELL_TYPE(IC)=NSCARC_CELL_TYPE_NONE
         ENDDO FALGOUT_INTERNAL_CELL_LOOP1

         !> initialize GRAPH and reset CELL_TYPE on ghost cells
         FALGOUT_EXTENDED_CELL_LOOP: DO IC = SL%NC+1, SL%NCE
            SL%GRAPH(IC) = IC
            SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
         ENDDO FALGOUT_EXTENDED_CELL_LOOP

         !> reset CELL_TYPE on internal wall cells
         DO IW = 1, SL%NW
            IF (SL%WALL(IW)%NOM /= 0) THEN
               IC = SL%WALL(IW)%ICW
               SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
               !WRITE(LU_SCARC,*) 'INTERNAL WALL CELLS: CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
            ENDIF
         ENDDO

         !> reset CELL_TYPE on internal fine cells
         ICG = 1
         FALGOUT_INTERNAL_CELL_LOOP2: DO IC = 1, SL%NC
            IF (SL%CELL_TYPE(IC)<NSCARC_CELL_TYPE_NONE) THEN
               SL%CELL_TYPE(IC)=NSCARC_CELL_TYPE_NONE
               !WRITE(LU_SCARC,*) 'FINE CELLS: CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
            ENDIF
            IF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_ZPNT) THEN
               IF (SL%MEASURE(IC) >= 1.0_EB .OR. (SL%S_ROW(IC+1)-SL%S_ROW(IC)) > 0) THEN
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
                  SL%GRAPH(ICG) = IC
                  !WRITE(LU_SCARC,*) 'ZPNTS : CELL_TYPE(',&
                  !  IC,')=',SL%CELL_TYPE(IC),': GRAPH(',ICG,')=',SL%GRAPH(ICG)
                  ICG = ICG + 1
               ELSE
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
               ENDIF
            ELSE IF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_SFPNT) THEN
               SL%MEASURE(IC) = NSCARC_MEASURE_NONE
            ELSE
               SL%GRAPH(ICG) = IC
               !WRITE(LU_SCARC,*) 'GRAPH(',IC,')=',SL%GRAPH(ICG)
               ICG = ICG + 1
            ENDIF
         ENDDO FALGOUT_INTERNAL_CELL_LOOP2

         IGRAPH  = ICG-1
         IGRAPHE = SL%NCE-SL%NC

      ENDDO FALGOUT_INIT_LOOP

      FALGOUT_EXTERNAL_LOOP: DO IC = SL%NC+1, SL%NCE
         SL%MEASURE(IC) = NSCARC_MEASURE_NONE
      ENDDO FALGOUT_EXTERNAL_LOOP

      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_FALGOUT', 'CELL_TYPE AFTER FALGOUT-INIT ')
      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_FALGOUT', 'MEASURE AFTER FALGOUT-INIT ')

      !> Coloring loop until all cells are coarse or fine
      FALGOUT_COLORING_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         ICYCLE = 0
         FALGOUT_GRAPH_LOOP: DO WHILE (BCYCLE)

            IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MEASURE_ADD, NL)
            !IF (NMESHES > 1) THEN
            !   TYPE_VECTOR = NSCARC_VECTOR_MEASURE
            !   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)
            !ENDIF

            CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_COLORING', 'FIRST COLORING')

            !WRITE(LU_SCARC,*) 'ICYCLE=',ICYCLE
            IF (ICYCLE > 0) THEN

               ICG = 1
               !WRITE(LU_SCARC,*) '====================================================='
               !WRITE(LU_SCARC,*) '==== ICYCLE=',ICYCLE, ': ICG=',ICG,': IGRAPH=',IGRAPH
               !WRITE(LU_SCARC,*) '====================================================='
               FALGOUT_GRAPH_LOOP1: DO WHILE (ICG <= IGRAPH)

                  IC = SL%GRAPH(ICG)
                  !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4,a,i4,a,f12.6)') &
                  !   'FALGOUT_LOOP: ICG=',ICG,': IGRAPH=',IGRAPH,&
                  !   ': IC=',IC,': CELL_TYPE=',SL%CELL_TYPE(IC),': MEASURE=',SL%MEASURE(IC)

                  !> if cell isn't marked as coarse cell yet and has measure
                  !> less than 1, mark it as fine cell
                  !> take care that all dependencies have been taken into account
                  IF (SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_COARSE .AND. &
                      SL%MEASURE(IC)  <  NSCARC_MEASURE_ONE) THEN
                     SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
                     DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                        JC = SL%S_COL(ICOL)
                        IF (JC < 0) CYCLE
                        SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
                     ENDDO
                     !WRITE(LU_SCARC,'(a,i4,a,f12.6,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4)') &
                     !   'CASE 1: IC=',IC,': Measure=',SL%MEASURE(IC),&
                     !   ':CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC), &
                     !   ':GRAPH(',ICG,')=',SL%GRAPH(ICG), &
                     !   ':GRAPH(',IGRAPH,')=',SL%GRAPH(IGRAPH)
                  ENDIF

                  !> if cell is already marked as fine or coarse, set its measure to zero
                  !> and extract it from the graph (put it at the end of the graph
                  !> array and decrease number of relevant graph entries)
                  IF (SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_NONE) THEN
                     SL%MEASURE(IC) = NSCARC_MEASURE_NONE
                     SL%GRAPH(ICG) = SL%GRAPH(IGRAPH)
                     SL%GRAPH(IGRAPH) = IC
                     !WRITE(LU_SCARC,'(a,i4,a,f12.6,a,i4,a,i4,a,i4,a,i4,a,i4,a,3i4)') &
                     !   'CASE 2: IC=',IC,': Measure=',SL%MEASURE(IC),&
                     !   ':CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC), &
                     !   ':GRAPH(',ICG,')=',SL%GRAPH(ICG), &
                     !   ':GRAPH(',IGRAPH,')=',SL%GRAPH(IGRAPH), IGRAPH-1, ICG-1
                     IGRAPH = IGRAPH - 1
                     ICG = ICG - 1
                  ENDIF

                  ICG = ICG + 1
               ENDDO FALGOUT_GRAPH_LOOP1

               CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_COLORING', 'END OF FIRST COLORING')
               !WRITE(LU_SCARC,*) 'GRAPH:'
               !WRITE(LU_SCARC,'(16i4)') (SL%GRAPH(III),III=1,16)

            ENDIF

            !> exchanges measure and CELL_TYPEs of neighboring cells
            IF (NMESHES > 1) THEN
               TYPE_VECTOR = NSCARC_VECTOR_MEASURE
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)
            ENDIF

      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_FALGOUT', 'BEFORE MPI_ALLREDUCE ')
      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_FALGOUT', 'BEFORE MPI_ALLREDUCE ')

            !> Get global number of relevant graph entries (was ist im seriellen Fall??)
            IF (NMESHES>1 .AND. N_MPI_PROCESSES>1) THEN
               CALL MPI_ALLREDUCE (IGRAPH, IGRAPH_GLOBAL, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
            ENDIF
            !WRITE(LU_SCARC,*) 'AFTER MPI_ALLREDUCE: IGRAPH=',IGRAPH,': IGRAPH_GLOBAL=',IGRAPH_GLOBAL

            !> If all cells have been encountered for, leave loop
            IF (IGRAPH_GLOBAL == 0) EXIT FALGOUT_GRAPH_LOOP

            !> search for an independent set of points with maximal measure
            IF (ICYCLE > 0) CALL SCARC_SETUP_GRAPHSET(IGRAPH, NM, NL)
            ICYCLE = ICYCLE + 1

            !> exchanges CELL_TYPEs of neighboring cells (with new information from graphset)
            IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)

            ICG = SL%NC+1
            !WRITE(LU_SCARC,*) 'ICG=',ICG,': IGRAPHE=',IGRAPHE

      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_FALGOUT', 'BEFORE GRAPH_LOOP2 ')
      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_FALGOUT', 'BEFORE GRAPH_LOOP2 ')

            FALGOUT_GRAPH_LOOP2: DO WHILE (ICG <= SL%NC+IGRAPHE)
               IC = SL%GRAPH(ICG)
               !WRITE(LU_SCARC,*) 'IN SECOND FALGOUT_WHILE, ICG=',ICG, ': IC=',IC

               IF (IC < 0) CYCLE
               IF (SL%CELL_TYPE(IC) < NSCARC_CELL_TYPE_NONE) THEN
                  SL%GRAPH(ICG) = SL%GRAPH(SL%NC+IGRAPHE)
                  SL%GRAPH(SL%NC+IGRAPHE) = IC
                  !WRITE(LU_SCARC,'(a,i4,a,f12.6,a,i4,a,i4,a,i4,a,i4,a,i4,a,3i4)') &
                  !   'CASE 3: IC=',IC,': Measure=',SL%MEASURE(IC),&
                  !   ':CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC), &
                  !   ':GRAPH(',ICG,')=',SL%GRAPH(ICG), &
                  !   ':GRAPH(',SL%NC+IGRAPHE,')=',SL%GRAPH(SL%NC+IGRAPHE), IGRAPHE-1, ICG
                  IGRAPHE = IGRAPHE - 1
               ENDIF
               ICG = ICG + 1
            ENDDO FALGOUT_GRAPH_LOOP2

            SL%MEASURE(SL%NC+1:SL%NCE) = NSCARC_MEASURE_NONE

      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_FALGOUT', 'BEFORE GRAPH_LOOP3 ')
      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_FALGOUT', 'BEFORE GRAPH_LOOP3 ')
            FALGOUT_GRAPH_LOOP3: DO ICG = 1, IGRAPH

               IC = SL%GRAPH(ICG)

   !WRITE(LU_SCARC,*) '--------------------------------------------------'
   !WRITE(LU_SCARC,*) '---- PROCESSING NODE ',IC
   !WRITE(LU_SCARC,*) '--------------------------------------------------'

               !> Coarse cells don't interpolate from influencing neighbors
               IF (SL%CELL_TYPE(IC) > NSCARC_CELL_TYPE_NONE) THEN

                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE

                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                     JC = SL%S_COL(ICOL)

                     !> remove edge from S and decrement measures of unmarked neighbors
                     IF (JC > -1) THEN
                        SL%S_COL(ICOL) = - SL%S_COL(ICOL) - 1
                        IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) SL%MEASURE(JC) = SL%MEASURE(JC) - 1
   !WRITE(LU_SCARC,'(A,i4,A,I4,A,i4,a,i4,A,I4,A,I4,a,i4,a,f12.6)') &
   !' 1. IF    S_COL(',ICOL,')=',SL%S_COL(ICOL),&
   !': CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC),&
   !': CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC),&
   !': MEASURE(',JC,')=',SL%MEASURE(JC)
                     ENDIF
                  ENDDO

               ELSE

                  !> dependencies which are already marked
                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1

                     JC = SL%S_COL(ICOL)
                     IF (JC < 0) JC = -JC - 1

                     !> remove edge from S and temporarily reset CELL_TYPE
                     IF (SL%CELL_TYPE(JC) > NSCARC_CELL_TYPE_NONE) THEN

                        IF (SL%S_COL(ICOL) > -1) SL%S_COL(ICOL) = -SL%S_COL(ICOL) - 1
                        SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_CPNT
                        !WRITE(LU_SCARC,*) 'SETTING TO CPNT : CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)

                     ELSE IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_SFPNT) THEN    !> necessary ??
                        IF (SL%S_COL(ICOL) > -1) SL%S_COL(ICOL) = -SL%S_COL(ICOL) - 1
                     ENDIF

   !WRITE(LU_SCARC,'(A,i4,A,I4,A,i4,a,i4,A,I4,A,I4,a,i4,a,f12.6)') &
   !' 2. IF    S_COL(',ICOL,')=',SL%S_COL(ICOL),&
   !': CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC),&
   !': CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC),&
   !': MEASURE(',JC,')=',SL%MEASURE(JC)

                  ENDDO

                  !> dependencies which aren't marked yet
                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1

                     JC = SL%S_COL(ICOL)
                     !WRITE(LU_SCARC,*) 'JC=',JC,': BREAK=',BREAK
                     IF (JC > -1 .AND. JC<=SL%NC) THEN

                        BREAK = .TRUE.

                        !> check if there are common C-points
                        DO JCOL = SL%S_ROW(JC), SL%S_ROW(JC+1)-1
                           KC = SL%S_COL(JCOL)
                           IF (KC < 0) KC = -KC - 1
                     !WRITE(LU_SCARC,*) 'JCOL=',JCOL,': KC=',KC

                           IF (KC <= SL%NC) THEN
                              !> remove edge from S and update measure
                              IF (SL%CELL_TYPE(KC) == NSCARC_CELL_TYPE_CPNT) THEN
                                 SL%S_COL(ICOL) = - SL%S_COL(ICOL) - 1
                                 SL%MEASURE(JC) = SL%MEASURE(JC) - 1
                                 BREAK = .FALSE.
                                 EXIT
                              ENDIF
   !WRITE(LU_SCARC,'(A,i4,A,I4,A,i4,a,i4,A,I4,A,I4,a,i4,a,f12.6)') &
   !' 3. IF    S_COL(',JCOL,')=',SL%S_COL(JCOL),&
   !': CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC),&
   !': CELL_TYPE(',KC,')=',SL%CELL_TYPE(KC),&
   !': MEASURE(',JC,')=',SL%MEASURE(JC)
                           ENDIF
                        ENDDO

                        IF (BREAK) THEN
                           DO JCOL = SL%S_ROW(JC), SL%S_ROW(JC+1)-1
                              KC = SL%S_COL(JCOL)
                              IF (KC < 0) KC = -KC - 1
                              IF (KC > SL%NC) THEN
                                 !> remove edge from S and update measure
                                 IF (SL%CELL_TYPE(KC) == NSCARC_CELL_TYPE_CPNT) THEN
                                    SL%S_COL(ICOL) = - SL%S_COL(ICOL) - 1
                                    SL%MEASURE(JC) = SL%MEASURE(JC) - 1
                                    BREAK = .FALSE.
                                    EXIT
                                 ENDIF
                              ENDIF
   !WRITE(LU_SCARC,'(A,i4,A,I4,A,i4,a,i4,A,I4,A,I4,a,i4,a,f12.6)') &
   !' 4. IF    S_COL(',JCOL,')=',SL%S_COL(JCOL),&
   !': CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC),&
   !': CELL_TYPE(',KC,')=',SL%CELL_TYPE(KC),&
   !': MEASURE(',JC,')=',SL%MEASURE(JC)
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF

             !> reset CELL_TYPEs
               DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                  JC = SL%S_COL(ICOL)
                  IF (JC < 1) JC = -JC - 1
                  IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_CPNT) SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_COARSE
               ENDDO

            ENDDO FALGOUT_GRAPH_LOOP3
            CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_COLORING', 'MEASURE1 IN FALGOUT WHILE')
            CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE1 IN FALGOUT WHILE')
            CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_GRAPH   , NL, 'SETUP_COLORING', 'GRAPH1 IN FALGOUT WHILE')

         ENDDO FALGOUT_GRAPH_LOOP

       !> reset S-matrix
         DO ICOL = 1, SL%S_ROW(SL%NC+1)
            IF (SL%S_COL(ICOL) < 0) SL%S_COL(ICOL) = -SL%S_COL(ICOL)-1
         ENDDO

      ENDDO FALGOUT_COLORING_LOOP

      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_COLORING', 'MEASURE2 IN FALGOUT WHILE')
      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE2 IN FALGOUT WHILE')
      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_GRAPH   , NL, 'SETUP_COLORING', 'GRAPH2 IN FALGOUT WHILE')

END SELECT


END SUBROUTINE SCARC_SETUP_COLORING


!> ----------------------------------------------------------------------------------------
!> Select an independent set from the graph
!> ----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRAPHSET(IGRAPH, NM, NL)
INTEGER, INTENT(IN):: IGRAPH, NM, NL
INTEGER :: IG, IC, ICOL, JC
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

!> First mark every cell from the 'internal' graphset with measure bigger
!> than one as a coarse cell
DO IG = 1, IGRAPH
   IC = SL%GRAPH(IG)
   IF (SL%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
      !WRITE(LU_SCARC,*) 'SETUP_GRAPHSET1: IGRAPH=',IGRAPH,': CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
   ENDIF
ENDDO

!> Do the same with cells on the overlapping areas
DO IG = SL%NC+1, SL%NCE
   IC = SL%GRAPH(IG)
   IF (IC < 0) CYCLE
   IF (SL%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
      !WRITE(LU_SCARC,*) 'SETUP_GRAPHSET2: IGRAPH=',IGRAPH,': CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
   ENDIF
ENDDO

!> remove nodes from the initial independent set depending on their measure
!> For each cell consider every connection in stencil and set cell with
!> 1G/smaller measure to Zero-Type
DO IG = 1, IGRAPH
   IC = SL%GRAPH(IG)

   IF (IC < 0) CYCLE

   IF (SL%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
         JC = SL%A_COL(ICOL)
         IF (JC < 0) CYCLE
         IF (SL%MEASURE(JC) > NSCARC_MEASURE_ONE) THEN
            IF (SL%MEASURE(IC) > SL%MEASURE(JC)) THEN
               SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_NONE
      !WRITE(LU_SCARC,*) 'SETUP_GRAPHSET3: IGRAPH=',IGRAPH,': CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
            ELSE IF (SL%MEASURE(JC) > SL%MEASURE(IC)) THEN
               SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
      !WRITE(LU_SCARC,*) 'SETUP_GRAPHSET4: IGRAPH=',IGRAPH,': CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
            ENDIF
         ENDIF
      ENDDO
   ENDIF
ENDDO


END SUBROUTINE SCARC_SETUP_GRAPHSET

!> ------------------------------------------------------------------------------------------------
!> Numbering of single patches
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PATCH_CELL_TYPES(CELL_TYPE, NX, NZ, NX1, NX2, NY1, NY2, IZ, INCRX, INCRY, INCRZ)
INTEGER, INTENT(IN):: NX, NZ, NX1, NX2, NY1, NY2, INCRX, INCRY, INCRZ, IZ
INTEGER, DIMENSION(:), INTENT(OUT) :: CELL_TYPE
INTEGER :: IX, INX, INY, INZ, PX(3), PY(3), PZ(3), PC(3,3), INX0, INZ0

!> Only dummy (check, if variables are still needed!)
IX = NY1
IX = NY2
IX = INCRY

INX0 = 1
INZ0 = 1
INY  = 1
PY   = 1
DO IX = NX1, NX2, INCRX

   DO INZ = 1, INCRZ

      PZ(INZ) = IZ + INZ - 1
      IF (NZ>3.AND.(PZ(INZ)==1 .OR. PZ(INZ)== NZ)) INZ0 = INZ

      DO INX = 1, INCRX

         PX(INX) = IX + INX - 1
         IF (NX>3.AND.(PX(INX)==1 .OR. PX(INX)== NX)) INX0 = INX

         PC(INX,INZ) = (PZ(INZ)-1)*NX + PX(INX)
         !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4)') 'PC(',INX,',',INZ,')=',PC(INX,INZ)
         CELL_TYPE(PC(INX,INZ)) = NSCARC_CELL_TYPE_SFINE
      ENDDO
   ENDDO
   !WRITE(LU_SCARC,*) 'INX0=',INX0,': INZ0=',INZ0, NY1, NY2, INCRY
   CELL_TYPE(PC(INX0,INZ0)) = NSCARC_CELL_TYPE_COARSE

ENDDO

END SUBROUTINE SCARC_PATCH_CELL_TYPES


!> ------------------------------------------------------------------------------------------------
!> Setup CELL_TYPEs of mesh corresponding to requested coarsening strategy
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CELL_TYPES(NTYPE, NL)
INTEGER, PARAMETER  :: NCOUPLINGS = 20, NCYC_MAX=1000
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER  :: NM, REMAINING_CELLS, MEASURE_TYPE
INTEGER  :: IA, IC, JC, KC, IG, ICOL, JCOL, ICASE, IOR0, IZ, ILOOP, IGRAPH, ICOUNT, IROW, JROW
INTEGER  :: IW, INCRX, INCRZ
REAL(EB) :: MEASURE_MAX, EPS
LOGICAL  :: BEVENX, BEVENZ, BTWO_X, BTHREE_X, BTWO_Z, BTHREE_Z, BREMOVE=.FALSE.
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, SLF, SLC


SCARC_ROUTINE = 'SCARC_SETUP_CELL_TYPES'

EPS = 1.0E-12
MEASURE_MAX = 0.0_EB


!> Define CELL_TYPEs for corresponding coarsening strategy
SELECT CASE (NTYPE)

   !! ---------------------------------------------------------------------------------------------
   !! standard Rue-Stueben coarsening
   !!  - identify first point IC with maximum measure and set it to be a C-point
   !!  - identify surrounding strongly connected points JC and set them to be F-points
   !!  - increase measures of strongly connected neighbours KC of JC
   !!  - decrease measures of strongly connected neighbours of IC
   !! ---------------------------------------------------------------------------------------------

   CASE (NSCARC_COARSENING_RS3)

      MEASURE_TYPE = 1
      RS3_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SL => SCARC(NM)%LEVEL(NL)
         REMAINING_CELLS = SL%NC

         RS3_CYCLE_LOOP: DO WHILE (REMAINING_CELLS > 0)

            !> get maximum (remaining) measure for all cells
            MEASURE_MAX = MAXVAL(SL%MEASURE(1:SL%NC))
            IF (MEASURE_MAX <= EPS) EXIT RS3_CYCLE_LOOP

            RS3_CELL_LOOP: DO IC = 1, SL%NC

               !> Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, SL%MEASURE(IC))) THEN

                  SL%MEASURE(IC)  = NSCARC_MEASURE_NONE
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
                  REMAINING_CELLS = REMAINING_CELLS - 1

                  !> Determine set of fine cells
                  DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1

                     IF (STRONGLY_COUPLED(SL%A, SL%A_ROW, IC, ICOL)) THEN

                        !> IF JC hasn't been marked yet, set it to be a fine cell which is no longer measured
                        JC = SL%A_COL(ICOL)
                        IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) THEN

                           SL%MEASURE(JC)  = NSCARC_MEASURE_NONE
                           SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_FINE
                           REMAINING_CELLS = REMAINING_CELLS - 1

                           !>  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                           DO JCOL = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                              IF (STRONGLY_COUPLED(SL%A, SL%A_ROW, JC, JCOL)) THEN
                                 KC = SL%A_COL(JCOL)
                                 IF (KC /= 0) THEN
                                    IF (KC /= IC .AND. (SL%CELL_TYPE(KC)==NSCARC_CELL_TYPE_NONE)) THEN
                                       SL%MEASURE(KC) = SL%MEASURE(KC) + 1.0_EB
                                       MEASURE_MAX = MAX(MEASURE_MAX, SL%MEASURE(KC))
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDDO

                        ENDIF
                     ENDIF
                  ENDDO
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELL_TYPE', 'RS3_MEASURE ')
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPE', 'RS3_CELL_TYPE ')
                  EXIT RS3_CELL_LOOP
               ENDIF
            ENDDO RS3_CELL_LOOP
         ENDDO RS3_CYCLE_LOOP
         SL%NCW = 0
      ENDDO RS3_LOOP

   !! ---------------------------------------------------------------------------------------------
   !! set CELL_TYPEs for Falgout coarsening
   !! Care: ZPOINT is equal to FINE in this case !!!
   !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_FALGOUT)

      FALGOUT_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         !WRITE(LU_SCARC,*) 'A'
         !> copy column pointers for matrix
         DO IA = 1, SL%NA
            SL%S_COL(IA) = SL%A_COL(IA)
         ENDDO

         !WRITE(LU_SCARC,*) 'B'
         !> this part may be changed for other coarsening strategies (with respect to cfmarker != 1)
         ICOUNT = 1
         FALGOUT_CELL_LOOP: DO IC = 1, SL%NC
         !WRITE(LU_SCARC,*) 'PUH1:  IC=',IC
            IF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE) THEN
               IF (SL%MEASURE(IC) >= 1.0_EB .AND. SL%A_ROW(IC+1)-SL%A_ROW(IC)+1 > 0) THEN
                  SL%CELL_TYPE(IC)  = NSCARC_CELL_TYPE_NONE
         !WRITE(LU_SCARC,*) 'PUH2:  SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
                  SL%GRAPH(ICOUNT) = IC
         !WRITE(LU_SCARC,*) 'C:  SL%GRAPH(',ICOUNT,')=',SL%GRAPH(ICOUNT)
                  ICOUNT = ICOUNT + 1
                  IF (ICOUNT > SL%NC) THEN
                     WRITE(SCARC_MESSAGE,'(2A,I8)') TRIM(SCARC_ROUTINE),': Graph size must be increased'
                     CALL SHUTDOWN(SCARC_MESSAGE); RETURN
                  ENDIF
               ELSE
                  SL%CELL_TYPE(IC)  = NSCARC_CELL_TYPE_FINE
         !WRITE(LU_SCARC,*) 'D: SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
               ENDIF
            ELSEIF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE0) THEN
               SL%MEASURE(IC) = NSCARC_MEASURE_NONE
            ELSE
               SL%GRAPH(ICOUNT) = IC
         !WRITE(LU_SCARC,*) 'F:  SL%GRAPH(',ICOUNT,')=',SL%GRAPH(ICOUNT)
               ICOUNT = ICOUNT + 1
            ENDIF
         ENDDO FALGOUT_CELL_LOOP

!WRITE(LU_SCARC,*) '-----------------------------------------------'
!DO IC = 1, SL%NC
!   WRITE(LU_SCARC,*) 'E: SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
!ENDDO
!WRITE(LU_SCARC,*) '-----------------------------------------------'
!DO IC = 1, SL%NC
!   WRITE(LU_SCARC,*) 'E: SL%MEASURE(',IC,')=',SL%MEASURE(IC)
!ENDDO
!WRITE(LU_SCARC,*) '-----------------------------------------------'
!DO IC = 1, SL%NC
!   WRITE(LU_SCARC,*) 'E: SL%GRAPH(',IC,')=',SL%GRAPH(IC)
!ENDDO
!WRITE(LU_SCARC,*) '-----------------------------------------------'
!DO IC = 1, SL%NC
!   WRITE(LU_SCARC,*) 'E: SL%A_ROW(',IC,')=',SL%A_ROW(IC)
!ENDDO
!WRITE(LU_SCARC,*) '-----------------------------------------------'
!DO IC = 1, SL%NA
!   WRITE(LU_SCARC,*) 'E: SL%S_COL(',IC,')=',SL%S_COL(IC)
!ENDDO
!WRITE(LU_SCARC,*) '-----------------------------------------------'

         IGRAPH = ICOUNT
         ILOOP  = 0
         FALGOUT_COARSENING_LOOP: DO

!WRITE(LU_SCARC,*) 'iter =',ILOOP

            !> Care: Is graph_size always SL%NC ?
            IF (ILOOP > 0) THEN

!WRITE(LU_SCARC,*) 'HALLO2 ILOOP=',ILOOP

               DO IG = 1, IGRAPH

                  IC = SL%GRAPH(IG)
                  JC = IG

                  !> make IC a fine cell and look for all dependencies
                  IF ((SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_COARSE) .AND. (SL%MEASURE(IC) < 1)) THEN
                     SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
!WRITE(LU_SCARC,*) 'Setting %CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
                     DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                        IF (SL%S_COL(IROW) > -1) SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
!WRITE(LU_SCARC,*) ' jS: ',IROW, ' Setting SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
                     ENDDO
                  ENDIF

                  !> remove cell from graph
                  IF (SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_NONE) THEN
                     SL%MEASURE(JC) = NSCARC_MEASURE_NONE
                     IGRAPH = IGRAPH - 1
                     SL%GRAPH(JC) = SL%GRAPH(IGRAPH)
                     SL%GRAPH(IGRAPH) = IC
!WRITE(LU_SCARC,*) 'graph_size=',IGRAPH,' graph_array(',JC,')=',SL%GRAPH(JC),': graph_array(',IGRAPH,')=',SL%GRAPH(IGRAPH)
                     JC = JC - 1               !> correct ???
                  ENDIF

               ENDDO

            ENDIF
            ILOOP = ILOOP + 1

            FALGOUT_GRAPH_LOOP: DO IG = 1, IGRAPH

               IC = SL%GRAPH(IG)

!WRITE(LU_SCARC,*) 'SUSI ===== IC = ',IC,' ===================================================================='
!DO JC = 1, SL%NC
!   WRITE(LU_SCARC,*) 'E: SL%GRAPH(',JC,')=',SL%GRAPH(JC)
!ENDDO
!WRITE(LU_SCARC,*) '------------------------------------------------------------------------------'
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,10i3)') 'before i=',IC,&
!                                 ': CF_marker(',IC,')=',SL%CELL_TYPE(IC),&
!                                 ': graph_array(',ic,')=',SL%GRAPH(IC), &
!                                 ': S_diag_i(',IC,')=',SL%A_ROW(IC)
!WRITE(LU_SCARC,*) '------------------------------------------------------------------------------'

               !> C-points are not influenced by neighbors (no interpolation)
               IF (SL%CELL_TYPE(IC) > 0) THEN

!WRITE(LU_SCARC,*) '================== COARSE'
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE            !> define it as C-point

!WRITE(LU_SCARC,*) 'A0: SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)

                  DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                     JC = SL%S_COL(IROW)
!WRITE(LU_SCARC,*) 'A1: --------- IROW=',IROW,': JC=', JC
                     IF (JC > -1) THEN
                        SL%S_COL(IROW) = -SL%S_COL(IROW) - 1  !> remove edge
!WRITE(LU_SCARC,*) 'A2: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
!WRITE(LU_SCARC,*) 'A3: SL%CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
                        IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) THEN
!WRITE(LU_SCARC,*) 'A4: SL%MEASURE (',JC,')=',SL%MEASURE(JC)
                           SL%MEASURE(JC) = SL%MEASURE(JC) - 1.0_EB   !> decrement measure of unmarked neigh.
!WRITE(LU_SCARC,*) 'A5: SL%MEASURE (',JC,')=',SL%MEASURE(JC)
                        ENDIF
                     ENDIF
                  ENDDO

               ELSE

!WRITE(LU_SCARC,*) '================== FINE'
                  !> dependencies which have already been marked
                  DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                     JC = SL%S_COL(IROW)
!WRITE(LU_SCARC,*) 'C1: --------- IROW=',IROW,': JC=', JC
                     IF (JC < 0) JC = -JC-1
!WRITE(LU_SCARC,*) 'C2: --------- IROW=',IROW,': JC=', JC
!WRITE(LU_SCARC,*) 'C3: SL%CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
                     IF (SL%CELL_TYPE(JC) > 0) THEN
!WRITE(LU_SCARC,*) 'C4: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                        IF (SL%S_COL(IROW) > -1) THEN
                           SL%S_COL(IROW) = -SL%S_COL(IROW) -1   !> remove edge
!WRITE(LU_SCARC,*) 'C5: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                           SL%CELL_TYPE(JC)    = NSCARC_CELL_TYPE_COMMON   !> temporarily set CELL_TYPE to common
!WRITE(LU_SCARC,*) 'C6: SL%CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
                        ENDIF
                     ELSEIF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_FINE0) THEN
                        IF (SL%S_COL(IROW) > -1) THEN
                           SL%S_COL(IROW) = -SL%S_COL(IROW) -1   !> remove edge
!WRITE(LU_SCARC,*) 'C7: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                        ENDIF
                     ENDIF
                  ENDDO

                  !> dependencies which haven't been marked yet
                  DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                     IF (SL%S_COL(IROW) > -1) THEN
!WRITE(LU_SCARC,*) 'E: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                        BREMOVE = .TRUE.
                        JC = SL%S_COL(IROW)
                        DO JROW = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                           KC = SL%S_COL(JROW)
!WRITE(LU_SCARC,*) 'E1: SL%S_COL(',JROW,')=',SL%S_COL(JROW)
                           IF (KC < 0) KC = -KC-1                        !> check for all dependencies !!
                           IF (SL%CELL_TYPE(KC) == NSCARC_CELL_TYPE_COMMON) THEN
                              SL%S_COL(IROW) = -SL%S_COL(IROW)-1
!WRITE(LU_SCARC,*) 'E2: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                              SL%MEASURE(JC) = SL%MEASURE(JC) - 1.0_EB
!WRITE(LU_SCARC,*) 'E3: SL%MEASURE (',JC,')=',SL%MEASURE(JC)
                              BREMOVE = .FALSE.
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO

               ENDIF

               !> reset CELL_TYPES
               DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                  JC = SL%S_COL(IROW)
!WRITE(LU_SCARC,*) 'H1: jS = ', IROW,' j=',JC
                  IF (JC < 0) JC = -JC-1
!WRITE(LU_SCARC,*) 'H2: jS = ', IROW,' j=',JC
                  IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_COMMON) SL%CELL_TYPE(JC)=NSCARC_CELL_TYPE_COARSE
!WRITE(LU_SCARC,*) 'H3: SL%CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
               ENDDO

!WRITE(LU_SCARC,*) '------------------------------------------------------------------------------'
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,10i3)') 'after  i=',IC,&
!                                 ': CF_marker(',IC,')=',SL%CELL_TYPE(IC),&
!                                 ': graph_array(',ic,')=',SL%GRAPH(IC), &
!                                 ': S_diag_i(',IC,')=',SL%A_ROW(IC)
            ENDDO FALGOUT_GRAPH_LOOP

!WRITE(LU_SCARC,*) 'AFTER RESET'

!WRITE(LU_SCARC,*) '-----------------------------------------------'
!DO IC = 1, SL%NC
!   WRITE(LU_SCARC,*) 'Z: SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
!ENDDO
!WRITE(LU_SCARC,*) '-----------------------------------------------'
!DO IC = 1, SL%NC
!   WRITE(LU_SCARC,*) 'Z: SL%MEASURE(',IC,')=',SL%MEASURE(IC)
!ENDDO

         ENDDO FALGOUT_COARSENING_LOOP

      ENDDO FALGOUT_LOOP

   !! ---------------------------------------------------------------------------------------------
   !! set CELL_TYPEs for GMG-like interpolation
   !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_GMG)

      SELECT_GMG_DIMENSION: SELECT CASE (TYPE_DIMENSION)

         CASE (NSCARC_DIMENSION_TWO)

            GMG_LOOP: DO NM = 1, NMESHES

               IF (PROCESS(NM) /= MYID) CYCLE

               SLF => SCARC(NM)%LEVEL(NL)
               SLC => SCARC(NM)%LEVEL(NL+1)

               !> Analyze grid sizes
               INCRX = 2
               SLC%NX = SLF%NX/2
               IF (MOD(SLF%NX,2) == 0) THEN
                  BEVENX=.TRUE.
               ELSE
                  BEVENX=.FALSE.
               ENDIF

               INCRZ  = 2
               SLC%NZ = SLF%NZ/2
               IF (MOD(SLF%NZ,2) == 0) THEN
                  BEVENZ=.TRUE.
               ELSE
                  BEVENZ=.FALSE.
               ENDIF

               IF (BEVENX.AND.BEVENZ) THEN
                  ICASE = 0
               ELSE IF (BEVENX.AND..NOT.BEVENZ) THEN
                  ICASE = 1
               ELSE IF (.NOT.BEVENX.AND.BEVENZ) THEN
                  ICASE = 2
               ELSE
                  ICASE = 3
               ENDIF

               SLC%NY = 1
               SLC%NC = SLC%NX * SLC%NZ

               SELECT CASE(ICASE)
                  CASE (0)
                     DO IZ = 1, SLF%NZ ,2
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, 1,SLF%NX  ,1,1,IZ,2,1,2)
                     ENDDO
                  CASE (1)
                     DO IZ = 1, SLF%NZ-3,2
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, 1,SLF%NX  ,1,1,IZ,2,1,2)
                     ENDDO
                     CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, 1,SLF%NX  ,1,1,SLF%NZ-2,2,1,3)
                  CASE (2)
                     DO IZ = 1, SLF%NZ,2
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, 1       ,SLF%NX-3,1,1,IZ,2,1,2)
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, SLF%NX-2,SLF%NX-2,1,1,IZ,3,1,2)
                     ENDDO
                  CASE (3)
                     DO IZ = 1, SLF%NZ-3,2
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, 1       ,SLF%NX-3,1,1,IZ,2,1,2)
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, SLF%NX-2,SLF%NX-2,1,1,IZ,3,1,2)
                     ENDDO
                     CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, 1       ,SLF%NX-3,1,1,SLF%NZ-2,2,1,3)
                     CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, SLF%NX-2,SLF%NX-2,1,1,SLF%NZ-2,3,1,3)
               END SELECT

               IF (NMESHES > 1) THEN
               SLF%NCW = 0
               DO IC = 1, SLF%NC
                  IF (SLF%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     DO IOR0 = -3, 3
                        IW = SLF%WALL_INDEX(IC, IOR0)
                        IF (SLF%WALL(IW)%NOM /= 0) SLF%NCW = SLF%NCW + 1
                     ENDDO
                  ENDIF
               ENDDO
               ENDIF

            ENDDO GMG_LOOP

         CASE (NSCARC_DIMENSION_THREE)

            WRITE(*,*) 'GMG_COARSENING FOR 3D not yet implemented'

      END SELECT SELECT_GMG_DIMENSION

CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPE', 'GMG_COARSENING ')

   !! ---------------------------------------------------------------------------------------------
   !! set CELL_TYPEs for GMG-like interpolation
   !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_GMG3)

      SELECT_GMG3_DIMENSION: SELECT CASE (TYPE_DIMENSION)

         CASE (NSCARC_DIMENSION_TWO)

            GMG3_LOOP: DO NM = 1, NMESHES

               IF (PROCESS(NM) /= MYID) CYCLE

               SLF => SCARC(NM)%LEVEL(NL)
               SLC => SCARC(NM)%LEVEL(NL+1)

               !> Analyze grid sizes
               IF (MOD(SLF%NX,3) == 0) THEN
                  INCRX    = 3
                  SLC%NX   = SLF%NX/3
                  BTWO_X   = .FALSE.
                  BTHREE_X = .TRUE.
               ELSE IF (MOD(SLF%NX,2) == 0) THEN
                  INCRX    = 2
                  SLC%NX   = SLF%NX/2
                  BTWO_X   = .TRUE.
                  BTHREE_X = .FALSE.
               ELSE
                  INCRX    = 2
                  SLC%NX   = SLF%NX/2 + 1
                  BTWO_X   = .FALSE.
                  BTHREE_X = .FALSE.
               ENDIF

               IF (MOD(SLF%NZ,3) == 0) THEN
                  INCRZ    = 3
                  SLC%NZ   = SLF%NZ/3
                  BTWO_Z   = .FALSE.
                  BTHREE_Z = .TRUE.
               ELSE IF (MOD(SLF%NZ,2) == 0) THEN
                  INCRZ    = 2
                  SLC%NZ   = SLF%NZ/2
                  BTWO_Z   = .TRUE.
                  BTHREE_Z = .FALSE.
               ELSE
                  INCRZ    = 2
                  SLC%NZ   = SLF%NZ/2 + 1
                  BTWO_Z   = .FALSE.
                  BTHREE_Z = .FALSE.
               ENDIF

               IF (.NOT.(BTWO_X.OR.BTHREE_X).OR..NOT.(BTWO_Z.OR.BTHREE_Z)) THEN
                  WRITE(SCARC_MESSAGE,'(2A,I8)') TRIM(SCARC_ROUTINE),': Cell numbers not divisable by 2 or 3!'
                  CALL SHUTDOWN(SCARC_MESSAGE); RETURN
               ENDIF

               IF (BTHREE_X) THEN
                  ICASE = 0
               ELSE IF (BTWO_X) THEN
                  ICASE = 1
               ELSE
                  ICASE = 2
               ENDIF

               SLC%NY = 1
               SLC%NC = SLC%NX * SLC%NZ

               SELECT CASE(ICASE)
                  CASE (0)
                     DO IZ = 1, SLF%NZ, INCRZ
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, 1,SLF%NX  ,1,1,IZ,3,1,INCRZ)
                     ENDDO
                  CASE (1)
                     DO IZ = 1, SLF%NZ, INCRZ
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, 1,SLF%NX  ,1,1,IZ,2,1,INCRZ)
                     ENDDO
                  CASE (3)
                     DO IZ = 1, SLF%NZ, INCRZ
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, 1       ,SLF%NX-3,1,1,IZ,2,1,INCRZ)
                        CALL  SCARC_PATCH_CELL_TYPES(SLF%CELL_TYPE,SLF%NX,SLF%NZ, SLF%NX-2,SLF%NX-2,1,1,IZ,3,1,INCRZ)
                     ENDDO
               END SELECT

               IF (NMESHES > 1) THEN
               SLF%NCW = 0
               DO IC = 1, SLF%NC
                  IF (SLF%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     DO IOR0 = -3, 3
                        IW = SLF%WALL_INDEX(IC, IOR0)
                        IF (SLF%WALL(IW)%NOM /= 0) SLF%NCW = SLF%NCW + 1
                     ENDDO
                  ENDIF
               ENDDO
               ENDIF

            ENDDO GMG3_LOOP

         CASE (NSCARC_DIMENSION_THREE)

            WRITE(*,*) 'GMG3_COARSENING FOR 3D not yet implemented'

      END SELECT SELECT_GMG3_DIMENSION

CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPE', 'GMG3_COARSENING ')


   !! ---------------------------------------------------------------------------------------------
   !! first set CELL_TYPEs for cells in internal boundary layers (adjacent to neighboring meshes)
   !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_BDRY)

      BDRY_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SL => SCARC(NM)%LEVEL(NL)

         !> First define coarse cells on boundary layer(adjacent to neighboring meshes)
         BDRY_CYCLE_LOOP: DO

            !!!
            !> get maximum (remaining) measure for all cells
            !!!
            MEASURE_MAX = 0.0_EB
            BDRY_LOOP1: DO IW = 1, SL%NW
               IF (SL%WALL(IW)%NOM == 0) CYCLE BDRY_LOOP1
               IC = SL%WALL(IW)%ICW
               MEASURE_MAX = MAX(MEASURE_MAX, SL%MEASURE(IC))
            ENDDO BDRY_LOOP1
            IF (MEASURE_MAX <= EPS) EXIT BDRY_CYCLE_LOOP

            BDRY_LOOP2: DO IW = 1, SL%NW

               IF (SL%WALL(IW)%NOM == 0) CYCLE BDRY_LOOP2
               IC = SL%WALL(IW)%ICW

               !> Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, SL%MEASURE(IC))) THEN

                  SL%MEASURE(IC)  = 0.0_EB
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE

                  !> Determine set of fine cells
                  DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1

                   !> JC is set to be a fine cell which is no longer measured
                     JC = SL%A_COL(ICOL)

                     SL%MEASURE(JC)  = 0.0_EB
                     SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_SFINE

                     !>  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                     DO JCOL = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                        KC = SL%A_COL(JCOL)
                        IF (KC > 0 .AND. KC /= IC .AND. (SL%CELL_TYPE(KC)==NSCARC_CELL_TYPE_NONE)) THEN
                           SL%MEASURE(KC) = SL%MEASURE(KC) + 1.0_EB
                           MEASURE_MAX = MAX(MEASURE_MAX, SL%MEASURE(KC))
                        ENDIF
                     ENDDO

                  ENDDO

CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_MEASURES', 'BDRY_LOOP2 ')
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPE', 'BDRY_LOOP2 ')
                  EXIT BDRY_LOOP2
               ENDIF

            ENDDO BDRY_LOOP2
         ENDDO BDRY_CYCLE_LOOP
      ENDDO BDRY_LOOP



END SELECT


!> -----------------------------------------------------------------------------------------------------
!> Exchange CELL_TYPEs along internal boundaries
!> -----------------------------------------------------------------------------------------------------
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPES', 'CELL_TYPE FINAL0')
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPES', 'CELL_TYPE FINAL1')


CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPES', 'CELL_TYPE FINAL2')
!CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_MEASURES', 'MEASURE FINAL2')

END SUBROUTINE SCARC_SETUP_CELL_TYPES




!> ------------------------------------------------------------------------------------------------
!> Set cell types for A1 coarsening
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SET_TYPES_A1(MEASURE, CELL_TYPE, JC, OFFSET)
REAL(EB), DIMENSION(:), INTENT(OUT) :: MEASURE
INTEGER , DIMENSION(:), INTENT(OUT) :: CELL_TYPE
INTEGER , INTENT(IN) :: JC, OFFSET

MEASURE(JC)         = 0.0_EB
MEASURE(JC-OFFSET)  = 0.0_EB
MEASURE(JC+OFFSET)  = 0.0_EB

CELL_TYPE(JC)        = NSCARC_CELL_TYPE_SFINE
CELL_TYPE(JC-OFFSET) = NSCARC_CELL_TYPE_WFINE
CELL_TYPE(JC+OFFSET) = NSCARC_CELL_TYPE_WFINE

END SUBROUTINE SET_TYPES_A1


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

SELECT CASE (TYPE_DIMENSION)
   CASE (NSCARC_DIMENSION_TWO)
      IF (ABS(WALL(IW)%NOM) == 2) THEN
         IC = -1
      ELSE
         IC = (WALL(IW)%IZW-1)*NX + WALL(IW)%IXW
      ENDIF
   CASE (NSCARC_DIMENSION_THREE)
      IC = (WALL(IW)%IZW-1)*NX*NY + (WALL(IW)%IYW-1)*NX + WALL(IW)%IXW
END SELECT

GET_WALL_CELL = IC
RETURN

END FUNCTION GET_WALL_CELL


!> ------------------------------------------------------------------------------------------------
!> Perform interpolation to coarse level
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PROLONGATION(NL)
INTEGER , INTENT(IN) :: NL
INTEGER  :: NM, IP, IC, JC, KC, ICG, ICOL, JCOL, KCOL, ICA, SIGN0
INTEGER  :: IW, IW0, IW2, IDIAG, JDIAG, JCO, IC2, ICP, ICP2, ICOL0, ICOL_FIRST, ICOL_LAST
REAL(EB) :: SUM_COUPLED, SUM_CPOINTS, SCAL, SUM_COARSE, SUM_DIAG
REAL(EB) :: DATA_SUM, DATA_DIAG, DATA_INTERPOL
REAL(EB) :: VALUES(20), WEIGHTS(20)
INTEGER  :: NEIGHBOR(20), NWEIGHTS, NWEIGHTS2
INTEGER  :: COARSE_CELL(20), COARSE_INDEX(20)
LOGICAL  :: BFIRST
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL


!WRITE(LU_SCARC,*) 'TYPE_INTERPOL=',TYPE_INTERPOL

SELECT_INTERPOLATION: SELECT CASE (TYPE_INTERPOL)

   !! ------------------------------------------------------------------------------------------
   !! Standard interpolation
   !! ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_STANDARD)

      MESHES_LOOP_STANDARD: DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(*,*) 'Standard interpolatiion not yet implemented'
      ENDDO MESHES_LOOP_STANDARD
      stop

   !! ------------------------------------------------------------------------------------------
   !! GMG-like interpolation
   !! ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_GMG, NSCARC_INTERPOL_GMG3)
      IF (TYPE_DIMENSION == NSCARC_DIMENSION_TWO) THEN
         SCAL = 4.0_EB
      ELSE
         SCAL = 4.0_EB
      ENDIF
         SCAL = 2.0_EB

      MESHES_LOOP_GMG: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SL => SCARC(NM)%LEVEL(NL)

         ICP2 = 1
         DO IC = 1, SL%NCE
            BFIRST = .TRUE.
            DO ICP = 1, SL%NCCE

               ROW_LOOP: DO IC2 = SL%R_ROW(ICP),SL%R_ROW(ICP+1)-1
                  IF (SL%R_COL(IC2) == IC) THEN
                     SL%P(ICP2) = SCAL*SL%R(IC2)
                     IF (BFIRST) THEN
                        SL%P_ROW(IC) = ICP2
                     ENDIF
                     SL%P_COL(ICP2) = ICP
                     ICP2 = ICP2 + 1
                     BFIRST = .FALSE.
                     EXIT ROW_LOOP
                  ENDIF
               ENDDO ROW_LOOP

            ENDDO

         ENDDO
         SL%P_ROW(SL%NC+1)=ICP2

      ENDDO MESHES_LOOP_GMG

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)
    CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROLONG1')

   !! ------------------------------------------------------------------------------------------
   !! Classical interpolation
   !! ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_CLASSICAL)

      !WRITE(LU_SCARC,*) 'CLASSICAL'

      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NL, 'SETUP_PROLONGATION', 'Matrix structure ')
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION','CELL_TYPE INIT')

      MESHES_LOOP_CLASSICAL: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SL => SCARC(NM)%LEVEL(NL)

         IP = 1
         INTERNAL_CELL_LOOP_CLASSICAL: DO IC = 1, SL%NC

            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0

            !!!
            !> If IC is a coarse cell, its value is taken
            !!!
            IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               NEIGHBOR(1)= SL%CELL_TYPE(IC)
               WEIGHTS(1) = 1.0_EB
               NWEIGHTS   = 1

            !!!
            !> If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
            !!!
            ELSE

             !> Get main diagonal entry a_ii for that fine cell
               IDIAG = SL%A_ROW(IC)
               SUM_DIAG = SL%A(IDIAG)

               !> First search for all neighboring coarse grid cells (store them in NEIGHBOR)
               IW = 1
               DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                  JC = SL%A_COL(ICOL)
                  IF (SL%CELL_TYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     NEIGHBOR(IW) = JC
                     WEIGHTS(IW)  = -SL%A(ICOL)
                     IW = IW + 1
                  ENDIF
               ENDDO
               NWEIGHTS = IW - 1

               !> Then search for the strongly and weakly coupled fine grid cells
               DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1

                  IW2 = 1
                  SUM_COARSE = 0.0_EB

                  JC = SL%A_COL(ICOL)

                  SELECT CASE (SL%CELL_TYPE(JC))

                     CASE (NSCARC_CELL_TYPE_SFINE)

                        !> search for couplings KC of the strongly coupled JC which belong to the
                        !> coarse interpolatory set of IC
                        DO JCOL = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                           KC = SL%A_COL(JCOL)

                           NWEIGHTS_SFINE_LOOP: DO IW = 1, NWEIGHTS
                              IF (KC == NEIGHBOR(IW) ) THEN
                                 COARSE_CELL (IW2) = JCOL
                                 COARSE_INDEX(IW2) = IW
                                 SUM_COARSE = SUM_COARSE + SL%A(JCOL)
                                 IW2 = IW2 + 1
                                 EXIT NWEIGHTS_SFINE_LOOP
                              ENDIF
                           ENDDO NWEIGHTS_SFINE_LOOP
                        ENDDO
                        NWEIGHTS2 = IW2 - 1

                        DO IW2 = 1, NWEIGHTS2
                           JCOL = COARSE_CELL(IW2)
                           IW   = COARSE_INDEX(IW2)
                           WEIGHTS(IW) = WEIGHTS(IW) - SL%A(ICOL)*SL%A(JCOL)/REAL(SUM_COARSE,EB)
                        ENDDO

                     CASE (NSCARC_CELL_TYPE_WFINE)

                        SUM_DIAG = SUM_DIAG + SL%A(ICOL)

                  END SELECT
               ENDDO

               DO IW = 1, NWEIGHTS
                  NEIGHBOR(IW) = SL%CELL_TYPE(NEIGHBOR(IW))
                  WEIGHTS(IW)  = WEIGHTS(IW)/SUM_DIAG
               ENDDO

            ENDIF

            !>
            !> Define corresponding entry in interpolation matrix P by means of the upper weights
            !!!
            SL%P_ROW(IC) = IP
            DO IW = 1, NWEIGHTS
               IF  (NEIGHBOR(IW) /= -1) THEN
                  SL%P_COL(IP) = NEIGHBOR(IW)
                  SL%P(IP)     = WEIGHTS(IW)
                  IP = IP +1
               ENDIF
            ENDDO
            SL%P_ROW(SL%NCE+1) = IP

         ENDDO INTERNAL_CELL_LOOP_CLASSICAL

      ENDDO MESHES_LOOP_CLASSICAL

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)

   !! ------------------------------------------------------------------------------------------
   !! Classical interpolation2
   !! ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_CLASSICAL2)

      !WRITE(LU_SCARC,*) 'CLASSICAL2'
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NL, 'SETUP_PROLONGATION', 'TYPE2: Matrix structure ')
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION','TYPE2: CELL_TYPE INIT')

      MESHES_LOOP_CLASSICAL2: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE
         SL => SCARC(NM)%LEVEL(NL)

         !> Build Interpolation
         ICOL0 = 1
         INTERNAL_CELL_LOOP2_CLASSICAL2: DO IC = 1, SL%NC

            !> For a C-point IC, use identity and set mapping
            IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               SL%P_ROW(IC) = ICOL0
               SL%P_COL(ICOL0) = SL%CELL_TYPE(IC)
               SL%P(ICOL0) = 1.0_EB

               !WRITE(LU_SCARC,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)') &
               !  'C: IC=', IC,': P_ROW(',IC,')=', SL%P_ROW(IC),': P_COL(',ICOL0,')=',SL%P_COL(ICOL0), &
               !  ': P(',ICOL0,')==',SL%P(ICOL0)

               ICOL0 = ICOL0 + 1

            !> For a F-point IC, build interpolation
            ELSE

               SL%P_ROW(IC)= ICOL0               !> diagonal part of P
               ICOL_FIRST = ICOL0

               DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                  JC = SL%A_COL(ICOL)

                  !> If JC is a C-point, initialize interpolation to zero and set column number in P_COL
                  IF (SL%CELL_TYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     SL%CELL_MARKER(JC) = ICOL0
                     SL%P_COL(ICOL0) = SL%CELL_TYPE(JC)
                     SL%P(ICOL0) = 0.0_EB
                     !WRITE(LU_SCARC,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)') &
                     !  'D: IC=', IC,': CELL_MARKER(',JC,')=', SL%CELL_MARKER(JC),': P_COL(',ICOL0,')=',SL%P_COL(ICOL0), &
                     !  ': P(',ICOL0,')==',SL%P(ICOL0)
                     ICOL0 = ICOL0 + 1
                  !> If JC is a F-point, set it to be a strong F-point with relevant connections
                  ELSE IF (SL%CELL_TYPE(JC) /= NSCARC_CELL_TYPE_WFINE) THEN
                     SL%CELL_MARKER(JC) = NSCARC_CELL_TYPE_SFINE
                     !WRITE(LU_SCARC,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)') &
                     !  'E: IC=', IC,': CELL_MARKER(',JC,')=', SL%CELL_MARKER(JC)
                  ENDIF

               ENDDO
               ICOL_LAST = ICOL0 - 1

               !> consider the i-th row of A, start with diagonal part
               IDIAG = SL%A_ROW(IC)
               DATA_DIAG = SL%A(IDIAG)

               DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                  JC = SL%A_COL(ICOL)

                  !> If JC is a strongly coupled C-point to IC, the value a_(IC, JC) must be
                  !> incorporated to the interpolation weight
                  IF (SL%CELL_MARKER(JC) >= ICOL_FIRST) THEN
                     JCOL = SL%CELL_MARKER(JC)
                     SL%P(JCOL) = SL%P(JCOL) + SL%A(ICOL)

                  !WRITE(LU_SCARC,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,f12.6,a,i3,a,i3,a,i3)') &
                  !'G: IC=', IC,': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,': P(',JCOL,')=',SL%P(JCOL)

                  !> If JC is a strongly coupled F-point to IC, the value a_(IC, JC) must be distributed
                  !> to C-points which are strongly coupled to IC (no distribution to the diagonal part)
                  ELSE IF (SL%CELL_MARKER(JC) >= NSCARC_CELL_TYPE_SFINE) THEN
                     DATA_SUM  = 0.0_EB
                     SIGN0 = 1
                     IF (SL%A(IDIAG) > 0) SIGN0=-1

                     !> search JC in the row of A and get sum of the couplings to C-points strongly coupled to IC
                     DO JCOL = SL%A_ROW(JC), SL%A_ROW(JC+1)-1
                        KC = SL%A_COL(JCOL)
                        IF (SL%CELL_MARKER(KC) > ICOL_FIRST .AND. SIGN0*SL%A(JCOL) >0) THEN
                           DATA_SUM = DATA_SUM + SL%A(JCOL)
                           !WRITE(LU_SCARC,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,f12.6,a,i3,a,i3,a,i3)') &
                           !'H: IC=', IC,': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,': KC=',KC,': DATA_SUM=',DATA_SUM
                        ENDIF
                     ENDDO

                     IF (DATA_SUM .NE. 0) THEN

                        DATA_INTERPOL = SL%A(ICOL)/DATA_SUM

                        !> loop over row of A for JC and spread data
                        DO JCOL = SL%A_ROW(JC), SL%A_ROW(JC+1)-1
                           KC = SL%A_COL(JCOL)
                           IF (SL%CELL_MARKER(KC) >= ICOL_FIRST .AND. SIGN0*SL%A(JCOL) > 0) THEN
                              KCOL = SL%CELL_MARKER(KC)
                              SL%P(KCOL) = SL%P(KCOL) + DATA_INTERPOL * SL%A(JCOL)
                              !WRITE(LU_SCARC,'(a,i3, a, i3,a,i3,a,i3,a,i3,a,i3,a,f12.6,a,i3,a,i3)') &
                              !'I: IC=', IC,': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,': KC=',KC,': P(',KCOL,')=',SL%P(KCOL)
                           ENDIF
                        ENDDO


                     ENDIF

                  !> If JC is a weakly coupled F-point to IC, the value a_(IC, JC) must be
                  !> distributed to the diagonal
                  ELSE IF (SL%CELL_TYPE(JC) .NE. NSCARC_CELL_TYPE_WFINE) THEN
                     DATA_DIAG = DATA_DIAG + SL%A(ICOL)
                     !WRITE(LU_SCARC,'(a,i3, a, i3,a,i3,a,i3,a,i3,a,f12.6,a,i3,a,i3,a,i3)') &
                     !'J: IC=', IC,': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,': KC=',KC,': DATA_DIAG=',DATA_DIAG
                  ENDIF

               ENDDO

               IF (DATA_DIAG == 0.0_EB) THEN
                  WRITE(*,*) "SCARC_SETUP_PROLONGATION: WARNING - zero diagonal data !"
                  DO ICOL = ICOL_FIRST, ICOL_LAST
                     SL%P(ICOL) = 0.0_EB
                  ENDDO
               ELSE
                  DO ICOL = ICOL_FIRST, ICOL_LAST
                     SL%P(ICOL) = - SL%P(ICOL)/DATA_DIAG
                     !WRITE(LU_SCARC,'(a,i3,a,i3,a,f12.6,a,i3,a,i3)') &
                     !'K: IC=', IC,': P(',ICOL,')=',SL%P(ICOL)
                  ENDDO
               ENDIF
            ENDIF

         ENDDO INTERNAL_CELL_LOOP2_CLASSICAL2
         SL%P_ROW(IC) = ICOL0
         !WRITE(LU_SCARC,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)') &
         !'F: IC=', IC,': P_ROW(',IC,')=', SL%P_ROW(IC)



      ENDDO MESHES_LOOP_CLASSICAL2

      !WRITE(LU_SCARC,*) 'SIZE(P)=',SIZE(SL%P)
      !WRITE(LU_SCARC,*) 'SIZE(P_ROW)=',SIZE(SL%P_ROW)
      !WRITE(LU_SCARC,*) 'SIZE(P_COL)=',SIZE(SL%P_COL)
      !WRITE(LU_SCARC,*) 'SIZE(R)=',SIZE(SL%R)
      !WRITE(LU_SCARC,*) 'SIZE(R_ROW)=',SIZE(SL%R_ROW)
      !WRITE(LU_SCARC,*) 'SIZE(R_COL)=',SIZE(SL%R_COL)
      !WRITE(LU_SCARC,*) 'P:'
      !WRITE(LU_SCARC,'(f12.6)') (SL%P(IC), IC=1, SIZE(SL%P))
      !WRITE(LU_SCARC,*) 'P_ROW:'
      !WRITE(LU_SCARC,'(i4)') (SL%P_ROW(IC), IC=1, SIZE(SL%P_ROW))
      !WRITE(LU_SCARC,*) 'P_COL:'
      !WRITE(LU_SCARC,'(i4)') (SL%P_COL(IC), IC=1, SIZE(SL%P_COL))



      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'TYPE2: final prolongation ')

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_PROL, NL)
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'TYPE2: final prolongation ')

   !! ---------------------------------------------------------------------------------------------
   !! Direct interpolation
   !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_DIRECT)

      !WRITE(LU_SCARC,*) 'DIRECT'
      MESHES_LOOP_DIRECT: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SL => SCARC(NM)%LEVEL(NL)
         CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION','CELL_TYPE INIT')

         IP = 1
         INTERNAL_CELL_LOOP: DO IC = 1, SL%NC

            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0

            !!!
            !> If IC is a coarse cell, its value is taken
            !!!
            IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               NEIGHBOR(1)= SL%CELL_TYPE(IC)
               NEIGHBOR(1)= IC
               WEIGHTS(1) = 1.0_EB
               NWEIGHTS   = 1

            !!!
            !> If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
            !!!
            ELSE

               !> Get main diagonal entry a_ii for that fine cell
               IDIAG = SL%A_ROW(IC)

               !> Select type of fine cell (weakly/strongly coupled)
               SELECT_FPOINT_TYPE: SELECT CASE(SL%CELL_TYPE(IC))

                  !!!
                  !> Strongly coupled fine cell IC
                  !> approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                  !> Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                  !!!
                  CASE (NSCARC_CELL_TYPE_SFINE)

                     !> Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                     SUM_COUPLED = 0.0_EB
                     SUM_CPOINTS = 0.0_EB

                     DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                        JC = SL%A_COL(ICOL)
                        SUM_COUPLED = SUM_COUPLED + SL%A(ICOL)
                        IF (SL%CELL_TYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + SL%A(ICOL)
                     ENDDO

                     SCAL = - SUM_COUPLED/SUM_CPOINTS

                     !> for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                     IW = 1
                     DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                        JC = SL%A_COL(ICOL)
                        IF (SL%CELL_TYPE(JC) > 0 ) THEN
                           NEIGHBOR(IW) = SL%CELL_TYPE(JC)
                           NEIGHBOR(IW) = JC
                           WEIGHTS(IW)  = SCAL * SL%A(ICOL)/SL%A(IDIAG)
                           IW = IW +1
                        ENDIF
                     ENDDO
                     NWEIGHTS = IW - 1

                  !!!
                  !> Weakly coupled fine cell IC:
                  !> Determine strongly coupled fine cells JC surrounding IC and, in turn, replace
                  !> each of them by a mean value of their surrounding strongly coupled coarse cells
                  !!!
                  CASE (NSCARC_CELL_TYPE_WFINE)

                     IW = 1                                             !> weights counter
                     DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1         !> loop over couplings of weakly coupled fine cell

                        !> regard subdiagonal matrix entries a_ij of weakly coupled fine cell JC to IC
                        JC = SL%A_COL(ICOL)

                        !> Find all surrounding (coupled) points of JC and compute scaling factor
                        !> compute scaling factor SCAL = a_ij/a_ii * [sum_(k in N_j) a_jk]/[sum_(l in C_j) ajl]
                        SUM_COUPLED = 0.0_EB
                        SUM_CPOINTS = 0.0_EB

                        DO JCOL = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                           KC = SL%A_COL(JCOL)
                           SUM_COUPLED = SUM_COUPLED + SL%A(JCOL)
                           IF (SL%CELL_TYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + SL%A(JCOL)
                        ENDDO

                        IF (SUM_CPOINTS == 0.0_EB) THEN
                           WRITE(SCARC_MESSAGE,'(2A,I8,A,I8,A,I8,A,I8)') TRIM(SCARC_ROUTINE),&
                              ': Error in SCARC_INTERPOLATION on mesh ',NM,' for cell ',IC,' with neighbor ', JC, &
                              '  on level ', NL
                           CALL SHUTDOWN(SCARC_MESSAGE); RETURN
                        ENDIF
                        SCAL =  SL%A(ICOL)/SL%A(IDIAG) * SUM_COUPLED/SUM_CPOINTS

                        !> Get diagonal matrix a_jj for point JC
                        JDIAG = SL%A_ROW(JC)

                        !> Compute interpolation weights for all strong coarse cell couplings KC of JC
                        !> note that a coarse cell KC may be considered several times for different JC's
                        DO JCOL = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                           KC = SL%A_COL(JCOL)
                           IF (SL%CELL_TYPE(KC) > 0) THEN
                             NEIGHBOR(IW) = SL%CELL_TYPE(KC)
                             NEIGHBOR(IW) = KC
                             WEIGHTS(IW)  = SCAL * SL%A(JCOL)/SL%A(JDIAG)
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
            !!!
            SL%P_ROW(IC) = IP
            DO IW = 1, NWEIGHTS
               IF  (NEIGHBOR(IW) /= -1) THEN
                  SL%P_COL(IP) = NEIGHBOR(IW)
                  SL%P(IP)     = WEIGHTS(IW)
                  IP = IP +1
               ENDIF
            ENDDO

         ENDDO INTERNAL_CELL_LOOP

         SL%P_ROW(SL%NC+1) = IP
         !SL%NROW0 = IP                !> really neede - OUTDATED ?

         !> Determine new coarse number and store it in CELL_TYPE
         INTERNAL_CELL_LOOP2: DO IC = 1, SL%NC
            DO ICOL = SL%P_ROW(IC), SL%P_ROW(IC+1)-1
              JC = SL%P_COL(ICOL)
              IF (JC <= SL%NC) THEN
                 SL%P_COL(ICOL) = SL%CELL_TYPE(JC)
              ENDIF
            ENDDO
         ENDDO INTERNAL_CELL_LOOP2


         !> In the multi-mesh case determine the new coarse cell numbers for neighboring cells
         IF (NMESHES > 1) THEN
            GHOST_CELL_LOOP: DO IW = 1, SL%NW
               IF (SL%WALL(IW)%NOM==0) CYCLE GHOST_CELL_LOOP
               ICG = SL%WALL(IW)%ICG(1)
               ICA = SL%WALL(IW)%ICW
               DO ICOL = SL%P_ROW(ICA), SL%P_ROW(ICA+1)-1
                 JC = SL%P_COL(ICOL)
                 IF (JC > SL%NC) THEN
                    JCO = SL%CCVAR(SL%WALL(IW)%IXG,SL%WALL(IW)%IYG,SL%WALL(IW)%IZG,UNKH)  ! ACHTUNG: SCHAUEN!
                    SL%P_COL(ICOL) = - JCO
                 ENDIF
               ENDDO
            ENDDO GHOST_CELL_LOOP
         ENDIF

      ENDDO MESHES_LOOP_DIRECT

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)

 !! ---------------------------------------------------------------------------------------------
 !! Direct interpolation with special treatment of boundary layers
 !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_DIRECT_BDRY)

      MESHES_LOOP_DIRECT_BDRY: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SL => SCARC(NM)%LEVEL(NL)
         CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION','CELL_TYPE INIT')

         IP = 1
         INTERNAL_CELL_BDRY_LOOP: DO IC = 1, SL%NC

            DIRECT_BDRY_IF: IF (SL%INTERNAL_BDRY_CELL(IC)/=0) THEN

               VALUES   = 0.0_EB
               WEIGHTS  = 0.0_EB
               NEIGHBOR = 0

               !!!
               !> If IC is a coarse cell, its value is taken
               !!!
               IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

                  NEIGHBOR(1)= SL%CELL_TYPE(IC)
                  NEIGHBOR(1)= IC
                  WEIGHTS(1) = 1.0_EB
                  NWEIGHTS   = 1


               !!!
               !> If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
               !!!
               ELSE


                  !> Get main diagonal entry a_ii for that fine cell
                  IDIAG = SL%A_ROW(IC)

                  !> Strongly coupled fine cell IC
                  !> approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                  !> Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                  !> Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                  SUM_COUPLED = 0.0_EB
                  SUM_CPOINTS = 0.0_EB

                  DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                     JC = SL%A_COL(ICOL)
                     SUM_COUPLED = SUM_COUPLED + SL%A(ICOL)
                     IF (JC > SL%NC) CYCLE
                     IF (SL%INTERNAL_BDRY_CELL(JC)/=0.AND.SL%CELL_TYPE(JC)>0) SUM_CPOINTS = SUM_CPOINTS + SL%A(ICOL)
                  ENDDO

                  SCAL = - SUM_COUPLED/SUM_CPOINTS

                  !> for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                  IW = 1
                  DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                     JC = SL%A_COL(ICOL)
                     IF (JC > SL%NC) CYCLE
                     IF (SL%INTERNAL_BDRY_CELL(JC)/=0.AND.SL%CELL_TYPE(JC) > 0 ) THEN
                        NEIGHBOR(IW) = SL%CELL_TYPE(JC)
                        NEIGHBOR(IW) = JC
                        WEIGHTS(IW)  = SCAL * SL%A(ICOL)/SL%A(IDIAG)
                        IW = IW +1
                     ENDIF
                  ENDDO
                  NWEIGHTS = IW - 1

               ENDIF


             !> ----------------------------!!> real internal cell ------------------------
            ELSE

               VALUES   = 0.0_EB
               WEIGHTS  = 0.0_EB
               NEIGHBOR = 0

               !!!
               !> If IC is a coarse cell, its value is taken
               !!!
               IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

                  NEIGHBOR(1)= SL%CELL_TYPE(IC)
                  NEIGHBOR(1)= IC
                  WEIGHTS(1) = 1.0_EB
                  NWEIGHTS   = 1

               !!!
               !> If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
               !!!
               ELSE

                  !> Get main diagonal entry a_ii for that fine cell
                  IDIAG = SL%A_ROW(IC)

                  !> Select type of fine cell (weakly/strongly coupled)
                  SELECT_FPOINT_BDRY_TYPE: SELECT CASE(SL%CELL_TYPE(IC))

                     !!!
                     !> Strongly coupled fine cell IC
                     !> approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                     !> Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                     !!!
                     CASE (NSCARC_CELL_TYPE_SFINE)

                        !> Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                        SUM_COUPLED = 0.0_EB
                        SUM_CPOINTS = 0.0_EB

                        DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                           JC = SL%A_COL(ICOL)
                           SUM_COUPLED = SUM_COUPLED + SL%A(ICOL)
                           IF (SL%CELL_TYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + SL%A(ICOL)
                        ENDDO

                        SCAL = - SUM_COUPLED/SUM_CPOINTS

                        !> for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                        IW = 1
                        DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                           JC = SL%A_COL(ICOL)
                           IF (SL%CELL_TYPE(JC) > 0 ) THEN
                              NEIGHBOR(IW) = SL%CELL_TYPE(JC)
                              NEIGHBOR(IW) = JC
                              WEIGHTS(IW)  = SCAL * SL%A(ICOL)/SL%A(IDIAG)
                              IW = IW +1
                           ENDIF
                        ENDDO
                        NWEIGHTS = IW - 1

                     !!!
                     !> Weakly coupled fine cell IC:
                     !> Determine strongly coupled fine cells JC surrounding IC and, in turn, replace
                     !> each of them by a mean value of their surrounding strongly coupled coarse cells
                     !!!
                     CASE (NSCARC_CELL_TYPE_WFINE)

                        IW = 1                                             !> weights counter
                        DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1         !> loop over couplings of weakly coupled fine cell

                           !> regard subdiagonal matrix entries a_ij of weakly coupled fine cell JC to IC
                           JC = SL%A_COL(ICOL)

                           !> Find all surrounding (coupled) points of JC and compute scaling factor
                           !> compute scaling factor SCAL = a_ij/a_ii * [sum_(k in N_j) a_jk]/[sum_(l in C_j) ajl]
                           SUM_COUPLED = 0.0_EB
                           SUM_CPOINTS = 0.0_EB

                           DO JCOL = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                              KC = SL%A_COL(JCOL)
                              SUM_COUPLED = SUM_COUPLED + SL%A(JCOL)
                              IF (SL%CELL_TYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + SL%A(JCOL)
                           ENDDO

                           IF (SUM_CPOINTS == 0.0_EB) THEN
                              WRITE(SCARC_MESSAGE,'(2A,I8,A,I8,A,I8,A,I8)') TRIM(SCARC_ROUTINE),&
                                 ': Error in SCARC_INTERPOLATION on mesh ',NM,' for cell ',IC,' with neighbor ', JC, &
                                 '  on level ', NL
                              CALL SHUTDOWN(SCARC_MESSAGE); RETURN
                           ENDIF
                           SCAL =  SL%A(ICOL)/SL%A(IDIAG) * SUM_COUPLED/SUM_CPOINTS

                           !> Get diagonal matrix a_jj for point JC
                           JDIAG = SL%A_ROW(JC)

                           !> Compute interpolation weights for all strong coarse cell couplings KC of JC
                           !> note that a coarse cell KC may be considered several times for different JC's
                           DO JCOL = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                              KC = SL%A_COL(JCOL)
                              IF (SL%CELL_TYPE(KC) > 0) THEN
                                NEIGHBOR(IW) = SL%CELL_TYPE(KC)
                                NEIGHBOR(IW) = KC
                                WEIGHTS(IW)  = SCAL * SL%A(JCOL)/SL%A(JDIAG)
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
            !!!
            SL%P_ROW(IC) = IP
            DO IW = 1, NWEIGHTS
               IF  (NEIGHBOR(IW) /= -1) THEN
                  SL%P_COL(IP) = NEIGHBOR(IW)
                  SL%P(IP)     = WEIGHTS(IW)
                  IP = IP +1
               ENDIF
            ENDDO

         ENDDO INTERNAL_CELL_BDRY_LOOP

         SL%P_ROW(SL%NC+1) = IP
         !SL%NROW0 = IP                !> really neede - OUTDATED ?

         !> Determine new coarse number and store it in CELL_TYPE
         INTERNAL_CELL_BDRY_LOOP2: DO IC = 1, SL%NC
            DO ICOL = SL%P_ROW(IC), SL%P_ROW(IC+1)-1
              JC = SL%P_COL(ICOL)
              IF (JC <= SL%NC) THEN
                 SL%P_COL(ICOL) = SL%CELL_TYPE(JC)
              ENDIF
            ENDDO
         ENDDO INTERNAL_CELL_BDRY_LOOP2


      ENDDO MESHES_LOOP_DIRECT_BDRY

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)


   !! ---------------------------------------------------------------------------------------------
   !! Multipass interpolation
   !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_MULTIPASS)

      MESHES_LOOP_MULTIPASS: DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(*,*) 'Multipass interpolation not yet implemented'
      ENDDO MESHES_LOOP_MULTIPASS
      stop

END SELECT SELECT_INTERPOLATION


END SUBROUTINE SCARC_SETUP_PROLONGATION

!> ------------------------------------------------------------------------------------------------
!> Numbering of single patches
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PATCH_WEIGHTS(R, R_ROW, R_COL, CELL_TYPE, WEIGHTS, NX, NX1, NX2, NY1, NY2, IZ, INCRX, INCRY, INCRZ, IP, ICC)
INTEGER, INTENT(IN)    :: NX, NX1, NX2, NY1, NY2, IZ, INCRX, INCRY, INCRZ
INTEGER, INTENT(INOUT) :: IP, ICC
INTEGER,  DIMENSION(:)  , INTENT(OUT) :: R_ROW, R_COL, CELL_TYPE
REAL(EB), DIMENSION(:,:), INTENT(IN)  :: WEIGHTS
REAL(EB), DIMENSION(:)  , INTENT(OUT) :: R
INTEGER :: IX, INX, INZ, PX(3), PY(3), PZ(3), PC(3,3)


PY=1
DO INZ = 1, INCRZ
   PZ(INZ) = IZ + INZ - 1
ENDDO
WRITE(LU_SCARC,'(a,i4,a,7i4)') 'IZ=',IZ,': PZ=',PZ(1:INCRZ), NY1, NY2, INCRY

DO IX = NX1, NX2, INCRX

   DO INX = 1, INCRX
      PX(INX) = IX + INX - 1
   ENDDO

   R_ROW(ICC) = IP

   !WRITE(LU_SCARC,'(a,i4,a,4i4)') 'IX=',IX,': PX=',PX(1:INCRX)
   DO INZ = 1, INCRZ
      DO INX = 1, INCRX
         PC(INX,INZ) = (PZ(INZ)-1)*NX + PX(INX)
         IF (CELL_TYPE(PC(INX,INZ))== NSCARC_CELL_TYPE_COARSE) CELL_TYPE(PC(INX,INZ))=ICC
         R_COL(IP)   = PC(INX,INZ)
         !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4)') 'PC(',INX,',',INZ,')=',PC(INX,INZ)
         R(IP) = WEIGHTS(INX, INZ)
         IP = IP + 1
      ENDDO
   ENDDO
   ICC = ICC + 1

ENDDO

END SUBROUTINE SCARC_PATCH_WEIGHTS

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

!WRITE(LU_SCARC,*) 'STARTING MATRIX_TRANSPOSE'
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'NCE=',NCE
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'NCCE=',NCCE
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'M_ROW'
!WRITE(LU_SCARC,'(8i5)') (M_ROW(IC), IC=1, NCE+1)
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'M_COL'
!WRITE(LU_SCARC,'(8i5)') (M_COL(IC), IC=1, NCE)

!> identify the number of non-zero entries in every column of M (corresponds to a row in MT)
!> And store it in the MT_ROW-array (caution: starts from the second position)
MT_ROW(1) = 1
DO ICOL = 1, M_ROW(NCE+1)-1
   IC = M_COL(ICOL)
   MT_ROW(IC+1) = MT_ROW(IC+1) + 1
   !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,i3)') 'ICOL=',ICOL, ': MT_ROW(',IC+1,')=',MT_ROW(IC+1)
ENDDO

!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'MT_ROW 2:'
!WRITE(LU_SCARC,'(8i5)') (MT_ROW(IC), IC=1, NCE+1)

!> shift it to the first position while summing up the entries for M_ROW
DO IC = 2, NCCE+1
   MT_ROW(IC) = MT_ROW(IC) + MT_ROW(IC-1)
   !WRITE(LU_SCARC,'(a,i3,a,i3,i3,i3)') 'MT_ROW(',IC,')=',MT_ROW(IC), MT_ROW(IC-1)
ENDDO

!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'MT_ROW 3:'
!WRITE(LU_SCARC,'(8i5)') (MT_ROW(IC), IC=1, NCCE+1)

!> set correct entries for MT
DO IC = 1, NCE
   DO ICOL = M_ROW(IC), M_ROW(IC+1)-1
      JC   = M_COL(ICOL)
      JCOL = MT_ROW(JC)
      !WRITE(LU_SCARC,'(8i4,a,i3,a,i3)') (MT_ROW(I),I=1,8),'JCOL=',JCOL,': JC=',JC
      MT_COL(JCOL) = IC
      MT(JCOL)     = M(ICOL)
      MT_ROW(JC)   = MT_ROW(JC) + 1
      !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)') &
      !    'IC=',IC, ': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,&
      !    ': MT_COL(',JCOL,')=', MT_COL(JCOL), ': MT_ROW(',JCOL,')=',  MT_ROW(JC), ': MT(',JCOL,')=',MT(JCOL)
   ENDDO
ENDDO

!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'MT_ROW 4:'
!WRITE(LU_SCARC,'(8i5)') (MT_ROW(IC), IC=1, NCCE+1)

!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'IC=',NCCE+1, ': MT_ROW(',NCCE+1,')=',MT_ROW(NCCE+1)
DO IC = NCCE, 2, -1
   MT_ROW(IC) = MT_ROW(IC-1)
   !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'IC=',IC, ': MT_ROW(',IC,')=',MT_ROW(IC)
ENDDO
MT_ROW(1)=1
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'IC=',1, ': MT_ROW(',1,')=',MT_ROW(1)

END SUBROUTINE SCARC_MATRIX_TRANSPOSE


!> ------------------------------------------------------------------------------------------------
!> Build transpose R of prolongation matrix P
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM_AMG(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IERR = 0
INTEGER :: IROW, IROW_INIT=0, IROW_SAVE
INTEGER :: ICOL1, ICOL2, ICOL3, IC0, IC1, IC2, IC3
LOGICAL :: BSQUARE = .TRUE.
REAL(EB) :: R_VALUE, RA_VALUE, RAP_VALUE
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SLF, SLC

MESHES_LOOP_SYSTEM_AMG: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   SLF => SCARC(NM)%LEVEL(NL)
   SLC => SCARC(NM)%LEVEL(NL+1)

   !WRITE(LU_SCARC,*) 'SCARC_SETUP_SYSTEM_AMG: NC  =', SLF%NC
   !WRITE(LU_SCARC,*) 'SCARC_SETUP_SYSTEM_AMG: NCC =', SLF%NCC

   IROW = 1

!>! First look at all exterior entries
!>  IF (SLC%NW > 0) THEN
!>
!>     LOOP_NBR_CELLS: DO IC0 = 1, SL%NCC
!>
!>        IF (SLC%INTERNAL_BDRY(IC) < 0) CYCLE LOOP_NBR_CELLS
!>
!>        IROW_INIT = IROW
!>
!>      !> loop over all entries in row IC0 of R
!>        LOOP_NBR_R_ENTRIES: DO ICOL1 = SLF%R_ROW(IC0), SLF%R_ROW(IC0+1)-1
!>
!>           IC1 = SLF%R_COL(ICOL1)
!>           IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '--------- ICOL1 =', ICOL1,'   : IC1=',IC1
!>
!>         !> loop over all entries in row IC1 of A
!>           LOOP_NBR_A_ENTRIES: DO ICOL2 = SLF%A_ROW(IC1), SLF%A_ROW(IC1+1)-1
!>
!>              IC2 = SLF%A_COL(ICOL2)
!>              IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '   ------ ICOL2 =', ICOL2,'   : IC2=',IC2
!>
!>              IF (SLF%A_TAG(IC2) /= IC0) THEN
!>
!>               !> mark IC2 as already considered
!>                 SLF%A_TAG(IC2) = IC0
!>                 IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '          A_TAG(',IC2,') =', SLF%A_TAG(IC2)
!>
!>               !> loop over all entries in row IC2 of P
!>                 LOOP_NBR_P_ENTRIES: DO ICOL3 = SLF%P_ROW(IC2), SLF%P_ROW(IC2+1)-1
!>
!>                    IC3 = SLF%P_COL(ICOL3)
!>                    IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '    !> ---- ICOL3 =', ICOL3,'   : IC3=',IC3
!>                    !IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
!>                  !>   '          P_TAG(',IC3,') =', SLF%P_TAG(IC3),': IROW_INIT=',IROW_INIT,' : IROW=',IROW
!>
!>                  !> verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
!>                    IF (SLF%P_TAG(IC3) < IROW_INIT) THEN
!>                       SLF%P_TAG(IC3) = IROW
!>                       IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '          P_TAG(',IC3,') =', SLF%P_TAG(IC3)
!>                       IROW = IROW + 1
!>                       IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '          IROW =', IROW
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


   !! Determine size of matrix RAP
   !! loop over interior c-cells
   LOOP1_OWN_CELLS: DO IC0 = 1, SLF%NCCE

      !WRITE(LU_SCARC,*) '================ MARKMARK1:  IC0 =', IC0-1
      !WRITE(LU_SCARC,*) 'IROW     = ', IROW
      !WRITE(LU_SCARC,*) 'IROW_INIT= ', IROW_INIT
      !WRITE(LU_SCARC,*) '================ A_MARKER'
      !DO IC1 = 1, SLF%NCE
      !   WRITE(LU_SCARC,'(a,i4,a,i4)') 'A_MARKER(',IC1,')=', SLF%A_TAG(IC1)-1
      !ENDDO
      !WRITE(LU_SCARC,*) '================ P_MARKER'
      !DO IC1 = 1, SLF%NCCE
      !   WRITE(LU_SCARC,'(a,i4,a,i4)') 'P_MARKER(',IC1,')=', SLF%P_TAG(IC1)-1
      !ENDDO
      !WRITE(LU_SCARC,*) '================ A_ROW'
      !DO IC1 = 1, 12
      !   WRITE(LU_SCARC,'(a,i4,a,i4)') 'A_ROW(',IC1,')=', SLF%A_ROW(IC1)
      !ENDDO
      !WRITE(LU_SCARC,*) '================ A_COL'
      !DO IC1 = 1, 12
      !   WRITE(LU_SCARC,'(a,i4,a,i4)') 'A_COL(',IC1,')=', SLF%A_COL(IC1)
      !ENDDO
      !WRITE(LU_SCARC,*) '================ R_ROW'
      !DO IC1 = 1, 12
      !   WRITE(LU_SCARC,'(a,i4,a,i4)') 'R_ROW(',IC1,')=', SLF%R_ROW(IC1)
      !ENDDO
      !WRITE(LU_SCARC,*) '================ R_COL'
      !DO IC1 = 1, 12
      !   WRITE(LU_SCARC,'(a,i4,a,i4)') 'R_COL(',IC1,')=', SLF%R_COL(IC1)
      !ENDDO

      IROW_INIT = IROW

      IF (BSQUARE) THEN
         SLF%P_TAG(IC0) = IROW
         IROW = IROW + 1
         !WRITE(LU_SCARC,*) 'SQUARE P_TAG(',IC0,') =', SLF%P_TAG(IC0),'   : IROW=',IROW
      ENDIF

      !> loop over all entries in row IC0 of R
      LOOP1_R_ENTRIES: DO ICOL1 = SLF%R_ROW(IC0), SLF%R_ROW(IC0+1)-1

         IC1 = SLF%R_COL(ICOL1)
         !WRITE(LU_SCARC,*) '--------- ICOL1 =', ICOL1,'   : IC1=',IC1

         !> loop over all entries in row IC1 of A
         LOOP1_A_ENTRIES: DO ICOL2 = SLF%A_ROW(IC1), SLF%A_ROW(IC1+1)-1

            IC2 = SLF%A_COL(ICOL2)
            !WRITE(LU_SCARC,*) '   ------ ICOL2 =', ICOL2,'   : IC2=',IC2

            IF (SLF%A_TAG(IC2) /= IC0) THEN

               !> mark IC2 as already considered
               SLF%A_TAG(IC2) = IC0
               !WRITE(LU_SCARC,*) '          A_TAG(',IC2,') =', SLF%A_TAG(IC2)

               !> loop over all entries in row IC2 of P
               LOOP1_P_ENTRIES: DO ICOL3 = SLF%P_ROW(IC2), SLF%P_ROW(IC2+1)-1

                  IC3 = SLF%P_COL(ICOL3)
                  !WRITE(LU_SCARC,*) '    !> ---- ICOL3 =', ICOL3,'   : IC3=',IC3
                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
                  !>   '          P_TAG(',IC3,') =', SLF%P_TAG(IC3),': IROW_INIT=',IROW_INIT,' : IROW=',IROW

                  !> verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
                  IF (IC3 <= 0) WRITE (LU_SCARC,*) '------------!!!> CAUTION1: IC3 <=0 :' , IC3
                  IF (IC3 >0) THEN   !> ONLY TEMPORARILY
                  IF (SLF%P_TAG(IC3) < IROW_INIT) THEN
                     SLF%P_TAG(IC3) = IROW
                     !WRITE(LU_SCARC,*) '          P_TAG(',IC3,') =', SLF%P_TAG(IC3)
                     IROW = IROW + 1
                     !WRITE(LU_SCARC,*) '          IROW =', IROW
                  ENDIF
                  ENDIF   !> ONLY TEMPORARILY
               ENDDO LOOP1_P_ENTRIES

            ENDIF

         ENDDO LOOP1_A_ENTRIES
      ENDDO LOOP1_R_ENTRIES

      !> Store counters
      IROW_SAVE = IROW

   ENDDO LOOP1_OWN_CELLS

   !! Determine size of matrix RAP  == SLC%A
   SLC%NA = IROW
   SLC%NAE = SLC%NA + 100  !> ONLY TEMPORARILY

   !WRITE(LU_SCARC,*) 'AllOCATION RAP in length ', SLC%NA

   !! Allocate coarse matrix structures
   ALLOCATE (SLC%A(SLC%NAE+1), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLC%A', IERR)
   SLC%A = 0.0_EB

   ALLOCATE (SLC%A_COL(SLC%NAE+1), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLC%A_COL', IERR)
   SLC%A_COL = 0

   ALLOCATE (SLC%A_ROW(SLC%NCE+1), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SLC%A_ROW', IERR)
   SLC%A_ROW = 0
   SLC%A_ROW(SLF%NCC+1) = IROW

   SLF%A_TAG = 0
   SLF%P_TAG = 0

   IROW = 1

   !! loop over interior c-cells
   LOOP2_C_CELLS: DO IC0 = 1, SLF%NCC

      !WRITE(LU_SCARC,*) '================ MARKMARK2:  IC0 =', IC0-1
      !WRITE(LU_SCARC,*) 'IROW     = ', IROW
      !WRITE(LU_SCARC,*) 'IROW_INIT= ', IROW_INIT
      !WRITE(LU_SCARC,*) '================ A_MARKER'
      !DO IC1 = 1, SLF%NCE
      !   WRITE(LU_SCARC,'(a,i4,a,i4)') 'A_MARKER(',IC1,')=', SLF%A_TAG(IC1)-1
      !ENDDO
      !WRITE(LU_SCARC,*) '================ P_MARKER'
      !DO IC1 = 1, SLF%NCCE
      !   WRITE(LU_SCARC,'(a,i4,a,i4)') 'P_MARKER(',IC1,')=', SLF%P_TAG(IC1)-1
      !ENDDO

      IROW_INIT = IROW

      SLC%A_ROW(IC0) = IROW_INIT

      IF (BSQUARE) THEN
         SLF%P_TAG(IC0) = IROW
         SLC%A_COL(IROW) = IC0
         SLC%A(IROW) = 0.0_EB
         IROW = IROW + 1

         !WRITE(LU_SCARC,*) 'SQUARE P_TAG(',IC0,')         =', SLF%P_TAG(IC0)
         !WRITE(LU_SCARC,*) 'SQUARE RAP_diag_data(',IC0,') =', SLC%A(IROW)
         !WRITE(LU_SCARC,*) 'SQUARE RAP_diag_j(',IC0,')    =', SLC%A_COL(IROW)
         !WRITE(LU_SCARC,*) 'SQUARE IROW=',IROW
         !WRITE(LU_SCARC,*) 'NC =',SLF%NC
         !WRITE(LU_SCARC,*) 'NCC=',SLF%NCC

      ENDIF

      !WRITE(LU_SCARC,*) 'MYDEBUG: 0: ======== RAP 10: A_diag_i  ======================'
      !DO IC1 = 1, SLF%NC+1
      !   !WRITE(LU_SCARC,*) 'A_ROW(',IC1,')=',SLF%A_ROW(IC1)
      !   WRITE(LU_SCARC,'(a,i5)') 'MYDEBUG: 0: ', SLF%A_ROW(IC1)
      !ENDDO
      !WRITE(LU_SCARC,*) 'MYDEBUG: 0: ======== RAP 10: A_diag_i  ======================'
      !DO IC1 = 1, SLF%NA
      !   !WRITE(LU_SCARC,*) 'A_COL(',IC1,')=',SLF%A_COL(IC1)
      !   WRITE(LU_SCARC,'(a,i5)') 'MYDEBUG: 0: ', SLF%A_COL(IC1)
      !ENDDO
      !WRITE(LU_SCARC,*) 'MYDEBUG: 0: ======== RAP 10: A_diag_i  ======================'
      !DO IC1 = 1, SLF%NA
      !   !WRITE(LU_SCARC,*) 'A(',IC1,')=',SLF%A(IC1)
      !   WRITE(LU_SCARC,'(a,f12.6)') 'MYDEBUG: 0: ', SLF%A(IC1)
      !ENDDO
      !WRITE(LU_SCARC,*) '==============================================='
      !WRITE(LU_SCARC,*) 'MATRIX R_ROW:'
      !DO IC1 = 1, SLF%NCC+1
      !   WRITE(LU_SCARC,*) 'R_ROW(',IC1,')=',SLF%R_ROW(IC1)
      !ENDDO
      !WRITE(LU_SCARC,*) 'MATRIX R_COL:'
      !DO IC1 = 1, SLF%NR
      !   WRITE(LU_SCARC,*) 'R_COL(',IC1,')=',SLF%R_COL(IC1)
      !ENDDO
      !WRITE(LU_SCARC,*) 'MATRIX R:'
      !DO IC1 = 1, SLF%NR
      !   WRITE(LU_SCARC,*) 'R(',IC1,')=',SLF%R(IC1)
      !ENDDO

      !> loop over all entries in row IC0 of R
      LOOP2_R_ENTRIES: DO ICOL1 = SLF%R_ROW(IC0), SLF%R_ROW(IC0+1)-1

         IC1 = SLF%R_COL(ICOL1)
         R_VALUE = SLF%R(ICOL1)

         !WRITE(LU_SCARC,'(a,i3,a,i3,a,f12.6)') &
         !      '--------- ICOL1 =', ICOL1,'   : IC1 =',IC1,': R_VALUE=',R_VALUE

         !> loop over all entries in row IC1 of A
         LOOP2_A_ENTRIES: DO ICOL2 = SLF%A_ROW(IC1), SLF%A_ROW(IC1+1)-1

            IC2 = SLF%A_COL(ICOL2)
            RA_VALUE = R_VALUE * SLF%A(ICOL2)

            !WRITE(LU_SCARC,'(a,i3,a,i3,a,f12.6,a,i3,a,i3)') &
            !         '--------- ICOL2 =', ICOL2,'   : IC2 =',IC2,': RA_VALUE=',&
            !         RA_VALUE, ': A_TAG(',IC2,')=',SLF%A_TAG(IC2)

            !> Hasn't cell IC2 been considered before? (new values for RAP can only be done for unmarked cells)
            IF (SLF%A_TAG(IC2) /= IC0) THEN

               !> mark IC2 as already considered
               SLF%A_TAG(IC2) = IC0

               !> loop over all entries in row IC2 of P
               LOOP2_P_ENTRIES1: DO ICOL3 = SLF%P_ROW(IC2), SLF%P_ROW(IC2+1)-1

                  IC3 = SLF%P_COL(ICOL3)
                  RAP_VALUE = RA_VALUE * SLF%P(ICOL3)

                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,f12.6)') &
                  !   '---A----- ICOL3 =', ICOL3,'   : IC3 =',IC3,': RAP_VALUE=',RAP_VALUE

                    !> verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
                  IF (IC3 >0) THEN     !> ONLY TEMPORARILY

                  IF (SLF%P_TAG(IC3) < IROW_INIT) THEN
                     SLF%P_TAG(IC3)  = IROW
                     SLC%A_COL(IROW) = SLF%P_COL(ICOL3)
                     SLC%A(IROW)     = RAP_VALUE
                     !WRITE(LU_SCARC,*) '   A1     P_TAG(',IC3,')     =', SLF%P_TAG(IC3)
                     !WRITE(LU_SCARC,*) '   A1     SLC%A_COL(',IROW,') =', SLC%A_COL(IROW)
                     !WRITE(LU_SCARC,*) '   A1     SLC%A(',IROW,')     =', SLC%A(IROW)
                     IROW = IROW + 1
                     !WRITE(LU_SCARC,*) '   A1     IROW =', IROW
                  ELSE
                     SLC%A(SLF%P_TAG(IC3)) = SLC%A(SLF%P_TAG(IC3)) + RAP_VALUE
                     !WRITE(LU_SCARC,*) '   A2     SLC%A(',SLF%P_TAG(IC3),')     =', SLC%A(SLF%P_TAG(IC3))
                  ENDIF

                  ENDIF  !> ONLY TEMPORARILY
               ENDDO LOOP2_P_ENTRIES1

            !> or has it been already considered
            ELSE

               !> loop over all entries in row IC2 of P
               LOOP2_P_ENTRIES2: DO ICOL3 = SLF%P_ROW(IC2), SLF%P_ROW(IC2+1)-1

                  IC3 = SLF%P_COL(ICOL3)
                  RAP_VALUE = RA_VALUE * SLF%P(ICOL3)
                  SLC%A(SLF%P_TAG(IC3)) = SLC%A(SLF%P_TAG(IC3)) + RAP_VALUE

                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,f12.6,a,i3,a,f12.6)') &
                  !   '---B----- ICOL3 =', ICOL3,'   : IC3 =',IC3,': RAP_VALUE=',RAP_VALUE,&
                  !   ': A(',SLF%P_TAG(IC3),')=',SLC%A(SLF%P_TAG(IC3))

                  !WRITE(LU_SCARC,*) '    !> ---- ICOL3 =', ICOL3,'   : IC3=',IC3
                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
                  !   '          P_TAG(',IC3,') =', SLF%P_TAG(IC3),': IROW_INIT=',IROW_INIT,' : IROW=',IROW

               ENDDO LOOP2_P_ENTRIES2
            ENDIF

         ENDDO LOOP2_A_ENTRIES
      ENDDO LOOP2_R_ENTRIES

      !> Store counters
      IROW_SAVE = IROW

   ENDDO LOOP2_C_CELLS

   !WRITE(LU_SCARC,*) 'MYDEBUG: 0: ======== RAP 11: RAP_diag_i  ================= NCC:', SLF%NCC
   !DO IC1 = 1, SLF%NCC+1
   !   !WRITE(LU_SCARC,*) 'A_ROW(',IC1,')=',SLF%A_ROW(IC1)
   !   WRITE(LU_SCARC,'(a,i5)') 'MYDEBUG: 0: ', SLC%A_ROW(IC1)-1
   !ENDDO
   !WRITE(LU_SCARC,*) 'MYDEBUG: 0: ======== RAP 11: RAP_diag_i  ================= NA :', SLF%NA
   !DO IC1 = 1, SLC%NA
   !   !WRITE(LU_SCARC,*) 'A_COL(',IC1,')=',SLF%A_COL(IC1)
   !   WRITE(LU_SCARC,'(a,i5)') 'MYDEBUG: 0: ', SLC%A_COL(IC1)-1
   !ENDDO
   !WRITE(LU_SCARC,*) 'MYDEBUG: 0: ======== RAP 11: RAP_diag_i  ================= NA :=', SLF%NA
   !DO IC1 = 1, SLC%NA
   !   !WRITE(LU_SCARC,*) 'A(',IC1,')=',SLF%A(IC1)
   !   WRITE(LU_SCARC,'(a,f12.6)') 'MYDEBUG: 0: ', -SLC%A(IC1)
   !ENDDO


ENDDO MESHES_LOOP_SYSTEM_AMG

END SUBROUTINE SCARC_SETUP_SYSTEM_AMG

!> ------------------------------------------------------------------------------------------------
!> Define restriction matrix R (currently transpose of prolongation matrix P)
!>  - In spite of the additinal need for the storing of R, this is done to save computational time
!>  - during the later matrix transfer operations
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_RESTRICTION(NL)
INTEGER , INTENT(IN) :: NL
INTEGER  :: NM, NOM, ICP, ICP2, IC, IC2, IG, IG0, ICC, IP, INCRX, INCRZ, ICASE, IZ
LOGICAL  :: BFIRST, BEVENX, BEVENZ, BTWO_X, BTWO_Z, BTHREE_X, BTHREE_Z
REAL(EB) :: WEIGHTS2X2(2,2), WEIGHTS2X3(2,3), WEIGHTS3x2(3,2), WEIGHTS3X3(3,3)
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSL

SELECT CASE (TYPE_INTERPOL)

 !! ------------------------------------------------------------------------------------------
 !! GMG-like interpolation
 !! ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_GMG)

      WEIGHTS2X2 = RESHAPE ((/ 0.250_EB, 0.250_EB,    &
                               0.250_EB, 0.250_EB /), &
                   SHAPE(WEIGHTS2X2))

      WEIGHTS3X2 = RESHAPE ((/ 0.125_EB, 0.250_EB, 0.125_EB,    &
                               0.125_EB, 0.250_EB, 0.125_EB /), &
                   SHAPE(WEIGHTS3X2))

      WEIGHTS2X3 = RESHAPE ((/ 0.125_EB, 0.125_EB,    &
                               0.250_EB, 0.250_EB,    &
                               0.125_EB, 0.125_EB /), &
                   SHAPE(WEIGHTS2X3))

      WEIGHTS3X3 = RESHAPE ((/ 0.0625_EB, 0.1250_EB, 0.0625_EB,    &
                               0.1250_EB, 0.2500_EB, 0.1250_EB,    &
                               0.0625_EB, 0.1250_EB, 0.0625_EB /), &
                   SHAPE(WEIGHTS3X3))

      MESHES_LOOP_GMG: DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL => SCARC(NM)%LEVEL(NL)
         IP = 1
         ICC = 1

         !> Analyze grid sizes
         INCRX = 2
         IF (MOD(SL%NX,2) == 0) THEN
            BEVENX=.TRUE.
         ELSE
            BEVENX=.FALSE.
         ENDIF

         INCRZ  = 2
         IF (MOD(SL%NZ,2) == 0) THEN
            BEVENZ=.TRUE.
         ELSE
            BEVENZ=.FALSE.
         ENDIF

         IF (BEVENX.AND.BEVENZ) THEN
            ICASE = 0
         ELSE IF (BEVENX.AND..NOT.BEVENZ) THEN
            ICASE = 1
         ELSE IF (.NOT.BEVENX.AND.BEVENZ) THEN
            ICASE = 2
         ELSE
            ICASE = 3
         ENDIF

         SELECT CASE(ICASE)
            CASE (0)
               DO IZ = 1, SL%NZ ,2
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS2x2, SL%NX,&
                                           1      , SL%NX  , 1, 1, IZ , 2, 1, 2, IP, ICC)
               ENDDO
            CASE (1)
               DO IZ = 1, SL%NZ-3,2
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS2x2, SL%NX,&
                                           1      , SL%NX  , 1, 1, IZ , 2, 1, 2, IP, ICC)
               ENDDO
               CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS2x3, SL%NX,&
                                           1      , SL%NX  , 1, 1, SL%NZ-2, 2, 1, 3, IP, ICC)
            CASE (2)
               DO IZ = 1, SL%NZ,2
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS2x2, SL%NX,&
                                           1      , SL%NX-3, 1, 1, IZ, 2, 1, 2, IP, ICC)
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS3x2, SL%NX,&
                                           SL%NX-2, SL%NX-2, 1, 1, IZ, 3, 1, 2, IP, ICC)
               ENDDO
            CASE (3)
               DO IZ = 1, SL%NZ-3,2
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS2x2, SL%NX,&
                                           1      , SL%NX-3, 1, 1, IZ, 2, 1, 2, IP, ICC)
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS3x2, SL%NX,&
                                           SL%NX-2, SL%NX-2, 1, 1, IZ, 3, 1, 2, IP, ICC)
               ENDDO
               CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS2x3, SL%NX,&
                                        1      , SL%NX-3, 1, 1, SL%NZ-2, 2, 1, 3, IP, ICC)
               CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS3x3, SL%NX,&
                                        SL%NX-2, SL%NX-2, 1, 1, SL%NZ-2, 3, 1, 3, IP, ICC)
         END SELECT
         SL%R_ROW(SL%NCC+1) = IP

      ENDDO MESHES_LOOP_GMG

      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPES', 'CELL_TYPE FINAL3')


 !! ------------------------------------------------------------------------------------------
 !! GMG-like interpolation
 !! ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_GMG3)

      WEIGHTS2X2 = RESHAPE ((/ 0.250_EB, 0.250_EB,    &
                               0.250_EB, 0.250_EB /), &
                   SHAPE(WEIGHTS2X2))

      WEIGHTS3X2 = RESHAPE ((/ 0.125_EB, 0.250_EB, 0.125_EB,    &
                               0.125_EB, 0.250_EB, 0.125_EB /), &
                   SHAPE(WEIGHTS3X2))

      WEIGHTS2X3 = RESHAPE ((/ 0.125_EB, 0.125_EB,    &
                               0.250_EB, 0.250_EB,    &
                               0.125_EB, 0.125_EB /), &
                   SHAPE(WEIGHTS2X3))

      WEIGHTS3X3 = RESHAPE ((/ 0.0625_EB, 0.1250_EB, 0.0625_EB,    &
                               0.1250_EB, 0.2500_EB, 0.1250_EB,    &
                               0.0625_EB, 0.1250_EB, 0.0625_EB /), &
                   SHAPE(WEIGHTS3X3))

      MESHES_LOOP_GMG3: DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL => SCARC(NM)%LEVEL(NL)
         IP = 1
         ICC = 1

         !> Analyze grid sizes
         IF (MOD(SL%NX,3) == 0) THEN
            INCRX   = 3
            BTWO_X   = .FALSE.
            BTHREE_X = .TRUE.
         ELSE IF (MOD(SL%NX,2) == 0) THEN
            INCRX    = 2
            BTWO_X   = .TRUE.
            BTHREE_X = .FALSE.
         ELSE
            INCRX  = 2
            BTWO_X   = .FALSE.
            BTHREE_X = .FALSE.
         ENDIF

         IF (MOD(SL%NZ,3) == 0) THEN
            INCRZ   = 3
            BTWO_Z   = .FALSE.
            BTHREE_Z = .TRUE.
         ELSE IF (MOD(SL%NZ,2) == 0) THEN
            INCRZ    = 2
            BTWO_Z   = .TRUE.
            BTHREE_Z = .FALSE.
         ELSE
            INCRZ  = 2
            BTWO_Z   = .FALSE.
            BTHREE_Z = .FALSE.
         ENDIF

         IF (.NOT.(BTWO_X.OR.BTHREE_X).OR..NOT.(BTWO_Z.OR.BTHREE_Z)) THEN
            WRITE(SCARC_MESSAGE,'(2A)') TRIM(SCARC_ROUTINE),': Cell numbers not divisable by 2 or 3!'
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
         ENDIF

         IF (BTHREE_X) THEN
            ICASE = 0
         ELSE IF (BTWO_X) THEN
            ICASE = 1
         ELSE
            ICASE = 2
         ENDIF

         SELECT CASE(ICASE)
            CASE (0)
               DO IZ = 1, SL%NZ ,INCRZ
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS3x3, SL%NX,&
                                           1      , SL%NX  , 1, 1, IZ , 3, 1, INCRZ, IP, ICC)
               ENDDO
            CASE (1)
               DO IZ = 1, SL%NZ,INCRZ
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS2x3, SL%NX,&
                                           1      , SL%NX  , 1, 1, IZ , 2, 1, INCRZ, IP, ICC)
               ENDDO
            CASE (2)
               DO IZ = 1, SL%NZ,INCRZ
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS2x2, SL%NX,&
                                           1      , SL%NX-3, 1, 1, IZ, 3, 1, INCRZ, IP, ICC)
                  CALL SCARC_PATCH_WEIGHTS(SL%R, SL%R_ROW, SL%R_COL, SL%CELL_TYPE, WEIGHTS3x2, SL%NX,&
                                           SL%NX-2, SL%NX-2, 1, 1, IZ, 2, 1, INCRZ, IP, ICC)
               ENDDO
         END SELECT
         SL%R_ROW(SL%NCC+1) = IP

      ENDDO MESHES_LOOP_GMG3
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPES', 'CELL_TYPE FINAL3')


 !! ------------------------------------------------------------------------------------------
 !! DEFAULT
 !! ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_CLASSICAL2)

      MESHES_LOOP_CLASSICAL2: DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL => SCARC(NM)%LEVEL(NL)
         CALL SCARC_MATRIX_TRANSPOSE(SL%P, SL%P_ROW, SL%P_COL, SL%R, SL%R_ROW, SL%R_COL, &
                                     SL%NCE, SL%NCCE)
      ENDDO MESHES_LOOP_CLASSICAL2


 !! ------------------------------------------------------------------------------------------
 !! DEFAULT
 !! ------------------------------------------------------------------------------------------
   CASE DEFAULT

      MESHES_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE
         SL => SCARC(NM)%LEVEL(NL)

         IC2 = 1
         DO ICP = 1, SL%NCCE
            BFIRST = .TRUE.
            DO IC = 1, SL%NCE

               ROW_LOOP: DO ICP2 = SL%P_ROW(IC),SL%P_ROW(IC+1)-1
                  IF (SL%P_COL(ICP2) == ICP) THEN
                     SL%R(IC2) = SL%P(ICP2)
                     IF (BFIRST) THEN
                        SL%R_ROW(ICP) = IC2
                     ENDIF
                     SL%R_COL(IC2) = IC
                     IC2 = IC2 + 1
                     BFIRST = .FALSE.
                     EXIT ROW_LOOP
                  ENDIF
               ENDDO ROW_LOOP

            ENDDO

         ENDDO
         SL%R_ROW(SL%NCCE+1)=IC2

         OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

            IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
            OS  => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE OTHER_MESHES_LOOP
            OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

            IC2 = 1
            ICC = 0
            OTHER_GHOSTCELL_LOOP: DO IG0 = 1, OSL%NCG

               ICP = OSL%CELL_TYPE(IG0)
               IF (ICP < NSCARC_CELL_TYPE_COARSE) CYCLE OTHER_GHOSTCELL_LOOP
               ICC = ICC + 1

               BFIRST = .TRUE.
               DO IG = 1, OSL%NCG

                  OTHER_ROW_LOOP: DO ICP2 = OSL%P_ROW(IG),OSL%P_ROW(IG+1)-1
                     IF (OSL%P_COL(ICP2) == ICC) THEN
                        OSL%R(IC2) = OSL%P(ICP2)
                        IF (BFIRST) THEN
                           OSL%R_ROW(ICC) = IC2
                        ENDIF
                        OSL%R_COL(IC2) = IG
                        IC2 = IC2 + 1
                        BFIRST = .FALSE.
                        EXIT OTHER_ROW_LOOP
                     ENDIF
                  ENDDO OTHER_ROW_LOOP

               ENDDO

            ENDDO OTHER_GHOSTCELL_LOOP
            OSL%R_ROW(ICC+1)=IC2

         ENDDO OTHER_MESHES_LOOP

      ENDDO MESHES_LOOP

END SELECT


CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_RESTRICTION, NL, 'SETUP_RESTRICTION', 'FINAL RESTRICT')
END SUBROUTINE SCARC_SETUP_RESTRICTION

!> ------------------------------------------------------------------------------------------------
!> Define restriction matrix R (currently transpose of prolongation matrix P)
!>  - In spite of the additinal need for the storing of R, this is done to save computational time
!>  - during the later matrix transfer operations
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_RESTRICTION2(NL)
INTEGER , INTENT(IN) :: NL
INTEGER  :: NM, NOM, ICP, ICP2, IC2, IG, IG0, ICC, ICC0
LOGICAL  :: BFIRST
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSL

SELECT CASE (TYPE_INTERPOL)

 !! ------------------------------------------------------------------------------------------
 !! GMG-like interpolation
 !! ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_GMG, NSCARC_INTERPOL_GMG3)

      MESHES_LOOP: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         ICC0= 0
         OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

            IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
            OS  => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE OTHER_MESHES_LOOP
            OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

            IC2 = 1
            ICC = 0
            OTHER_GHOSTCELL_LOOP: DO IG0 = 1, OSL%NCG

               ICP = OSL%CELL_TYPE(IG0)
               IF (ICP < NSCARC_CELL_TYPE_COARSE) CYCLE OTHER_GHOSTCELL_LOOP
               ICC = ICC + 1

               BFIRST = .TRUE.
               DO IG = 1, OSL%NCG

!WRITE(LU_SCARC,*) '============= NM=',NM,': NOM=',NOM,': IG0=',IG0
                  OTHER_ROW_LOOP: DO ICP2 = OSL%P_ROW(IG),OSL%P_ROW(IG+1)-1
!WRITE(LU_SCARC,*) 'OSL%P_COL(',ICP2,')=',OSL%P_COL(ICP2)
                     IF (OSL%P_COL(ICP2) == ICC) THEN
                        OSL%R(IC2) = OSL%P(ICP2)
                        IF (BFIRST) THEN
                           ICC0 = ICC0 + 1
                           OSL%R_ROW(ICC0) = IC2
                        ENDIF
                        OSL%R_COL(IC2) = IG
                        IC2 = IC2 + 1
                        BFIRST = .FALSE.
                        EXIT OTHER_ROW_LOOP
                     ENDIF
                  ENDDO OTHER_ROW_LOOP

               ENDDO

            ENDDO OTHER_GHOSTCELL_LOOP
            OSL%R_ROW(ICC0+1)=IC2

         ENDDO OTHER_MESHES_LOOP

      ENDDO MESHES_LOOP

END SELECT


CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_RESTRICTION, NL, 'SETUP_RESTRICTION', 'FINAL RESTRICT')
END SUBROUTINE SCARC_SETUP_RESTRICTION2


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
INTEGER :: NM, IC, ICO2, IERR, ICP1, ICP2, IP, ICOL, IDIAG, IOR0, IW
REAL (EB), ALLOCATABLE, DIMENSION(:) :: VAL
!REAL(EB):: MAUX1(50,50)=0.0_EB!, MAUX2(30,30)=0.0_EB
REAL:: AUX1, AUX2, PW!, MATRIX(16,16)
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLF, SLC

IERR = 0
SELECT CASE (TYPE_COARSENING)

 !! ---------------------------------------------------------------------------------
 !! GMG-like coarsening
 !! ---------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_GMG, NSCARC_INTERPOL_GMG3)

      MESHES_GMG: DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE

         SLF => SCARC(NM)%LEVEL(NL)             !> pointer to fine level
         SLC => SCARC(NM)%LEVEL(NL+1)           !> pointer to coarse level

         IF (NMESHES==1) THEN
            ALLOCATE (VAL(SLF%NC), STAT=IERR)
         ELSE
            ALLOCATE (VAL(SLF%NCE), STAT=IERR)
         ENDIF
         CALL CHKMEMERR ('SCARC_TRANSFER_MATRIX', 'VAL', IERR)
         VAL = 0.0_EB

         IP = 1
         ICP1_GMG_LOOP: DO ICP1 = 1, SLF%NCC

!yyWRITE(*,*) 'ICP1=',ICP1, SLF%NCC

            ICP2_GMG_LOOP1: DO ICP2 = 1, SLF%NCC
!WRITE(*,*) 'ICP2=',ICP2

               AUX2 = 0.0_EB

               DO IC = 1, SLF%NC
!WRITE(*,*) 'IC=',IC, SLF%NC

                  IDIAG = SLF%A_ROW(IC)
                  !PW    = 0.5_EB* P_WEIGHT(IC, ICP2, NM, NL)
                  PW    = P_WEIGHT(IC, ICP2, NM, NL)
                  AUX1  = SLF%A(IDIAG) * PW

                  COUPLINGS_GMG_LOOP: DO ICOL = SLF%A_ROW(IC)+1, SLF%A_ROW(IC+1)-1
                     ICO2  = SLF%A_COL(ICOL)
                     !PW = 0.5_EB* P_WEIGHT(ICO2, ICP2, NM, NL)
                     PW =  P_WEIGHT(ICO2, ICP2, NM, NL)
                     AUX1 = AUX1 + SLF%A(ICOL) * PW
                  ENDDO COUPLINGS_GMG_LOOP
!MATRIX(IC, ICP2 ) = AUX1
!WRITE(LU_SCARC,*) 'MATRIX(',IC,',',ICP2,')=',MATRIX(IC,ICP2)

                  AUX2 = AUX2 + R_WEIGHT(IC, ICP1, NM, NL) * AUX1

               ENDDO
               VAL(ICP2) = AUX2

            ENDDO ICP2_GMG_LOOP1

!WRITE(LU_SCARC,'(4f12.6)') ((MATRIX(IC0,ICP0),ICP0=1,4),IC0=1,16)

            !> analyze new matrix line and store it corresponding to compact storage technique:
            !> (diagonal entry first)
            SLC%A(IP)       = VAL(ICP1)
            SLC%A_ROW(ICP1) = IP
            SLC%A_COL(IP)   = ICP1

            IP  = IP + 1
            ICP2_GMG_LOOP3: DO ICP2 = 1, SLF%NCC
               IF (ICP2 /= ICP1 .AND. ABS(VAL(ICP2)) >= 1.0E-12_EB) THEN
                  SLC%A(IP)     = VAL(ICP2)
                  SLC%A_COL(IP) = ICP2
                  IP  = IP + 1
               ENDIF
            ENDDO ICP2_GMG_LOOP3
            DO IOR0 = -3, 3
               IW = SLC%WALL_INDEX(ICP1, IOR0)
               IF (IW > 0) THEN
                  SLC%A(IP) = 0.0_EB
                  SLC%A_COL(IP) = -IW
                  IP  = IP + 1
               ENDIF
            ENDDO

            VAL = 0.0_EB
         ENDDO ICP1_GMG_LOOP

         SLC%A_ROW(SLC%NC+1) = IP
         SLC%NA = IP - 1

         DEALLOCATE(VAL)

      ENDDO MESHES_GMG

!WRITE(LU_SCARC,*) 'ENDING MATRIX'

 !! ---------------------------------------------------------------------------------
 !! Default case
 !! ---------------------------------------------------------------------------------
   CASE DEFAULT

      MESHES_LOOP: DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE

         SLF => SCARC(NM)%LEVEL(NL)             !> pointer to fine level
         SLC => SCARC(NM)%LEVEL(NL+1)           !> pointer to coarse level

         IF (NMESHES==1) THEN
            ALLOCATE (VAL(SLF%NC), STAT=IERR)
         ELSE
            ALLOCATE (VAL(SLF%NCE), STAT=IERR)
         ENDIF
         CALL CHKMEMERR ('SCARC_TRANSFER_MATRIX', 'VAL', IERR)
         VAL = 0.0_EB

         IP = 1
         ICP1_LOOP: DO ICP1 = 1, SLF%NCC


            !> ---------------------------------------------------------------------------------
            !> First: loop over the internal coarse cells
            !> ---------------------------------------------------------------------------------
            ICP2_LOOP1: DO ICP2 = 1, SLF%NCCE

               AUX2 = 0.0_EB

               DO IC = 1, SLF%NC

                  IDIAG = SLF%A_ROW(IC)
                  PW    = P_WEIGHT(IC, ICP2, NM, NL)
                  AUX1  = SLF%A(IDIAG) * PW

                  COUPLINGS_LOOP: DO ICOL = SLF%A_ROW(IC)+1, SLF%A_ROW(IC+1)-1
                     ICO2  = SLF%A_COL(ICOL)
                     PW = P_WEIGHT(ICO2, ICP2, NM, NL)
                     AUX1 = AUX1 + SLF%A(ICOL) * PW
                  ENDDO COUPLINGS_LOOP

                  AUX2 = AUX2 + P_WEIGHT(IC, ICP1, NM, NL) * AUX1

               ENDDO
               VAL(ICP2) = AUX2

            ENDDO ICP2_LOOP1

            !> analyze new matrix line and store it corresponding to compact storage technique:
            !> (diagonal entry first)
            SLC%A(IP)       = VAL(ICP1)
            SLC%A_ROW(ICP1) = IP
            SLC%A_COL(IP)   = ICP1

            IP  = IP + 1
            ICP2_LOOP2: DO ICP2 = 1, SLF%NCCE
            !DO ICP2 = 1, SLF%NCC
               IF (ICP2 /= ICP1 .AND. ABS(VAL(ICP2)) >= 1.0E-12_EB) THEN
                  SLC%A(IP)     = VAL(ICP2)
                  SLC%A_COL(IP) = ICP2
                  IP  = IP + 1
               ENDIF
            ENDDO ICP2_LOOP2

            VAL = 0.0_EB
         ENDDO ICP1_LOOP

         SLC%A_ROW(SLC%NC+1) = IP
         SLC%NA = IP - 1

         DEALLOCATE(VAL)

      ENDDO MESHES_LOOP

END SELECT

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NL+1, 'TRANSFER_MATRIX', 'Next coarser matrix1')
IF (NMESHES > 1) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_MATRIX_SUBDIAG, NL)

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NL+1, 'TRANSFER_MATRIX', 'Next coarser matrix2')


END SUBROUTINE SCARC_TRANSFER_MATRIX


!> ------------------------------------------------------------------------------------------------
!> Determine if corresponding coarse point contributes a non-zero interpolation weight or not
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION P_WEIGHT(IC, ICP, NM, NL)
INTEGER, INTENT(IN):: IC, ICP, NM, NL
INTEGER :: ICOL
REAL(EB) :: VAL
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

VAL = 0.0_EB
!IF (IC <= SL%NC) THEN
   P_WEIGHT_LOOP: DO ICOL = SL%P_ROW(IC), SL%P_ROW(IC+1)-1
      IF (ICOL > 0) THEN
        IF (SL%P_COL(ICOL) /= ICP) CYCLE P_WEIGHT_LOOP
      ENDIF
      VAL = SL%P(ICOL)
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
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

VAL = 0.0_EB
IF (ICP <= SL%NCC) THEN
   R_WEIGHT_LOOP: DO ICOL = SL%R_ROW(ICP), SL%R_ROW(ICP+1)-1
      IF (ICOL > 0) THEN
        IF (SL%R_COL(ICOL) /= IC) CYCLE R_WEIGHT_LOOP
      ENDIF
      VAL = SL%R(ICOL)
   ENDDO R_WEIGHT_LOOP
ELSE
  WRITE(*,*)
ENDIF

R_WEIGHT = VAL
RETURN
END FUNCTION R_WEIGHT



!> -----------------------------------------------------------------------------------------------
!> Interface for the CALL of ScaRC-solver with requested storage technique
!> -----------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SOLVER
REAL (EB) :: TNOW

TNOW = SECOND()

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !! Krylov method (CG/BICG)
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV: SELECT CASE (TYPE_KRYLOV)
         CASE (NSCARC_KRYLOV_CG)
            CALL SCARC_METHOD_CG  (NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)
         CASE (NSCARC_KRYLOV_BICG)
            CALL SCARC_METHOD_BICG(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)
      END SELECT SELECT_KRYLOV

   !! Multigrid method
   CASE (NSCARC_METHOD_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)

   !! MKL method
#ifdef WITH_MKL
   CASE (NSCARC_METHOD_MKL)

      SELECT_MKL: SELECT CASE (TYPE_MKL)
         CASE (NSCARC_MKL_GLOBAL) 
            CALL SCARC_METHOD_CLUSTER(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)
         CASE (NSCARC_MKL_LOCAL) 
            CALL SCARC_METHOD_PARDISO(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)
      END SELECT SELECT_MKL
#endif

END SELECT SELECT_METHOD

TUSED_SCARC(NSCARC_TIME_SOLVER,MYID+1)=TUSED_SCARC(NSCARC_TIME_SOLVER,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
END SUBROUTINE SCARC_SOLVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Set pointer to chosen vector for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!> Set pointer to chosen vector for compact storage technique
!> ------------------------------------------------------------------------------------------------
SUBROUTINE POINT_TO_VECTOR(NVECTOR, NM, NL, VEC)
REAL(EB), INTENT(OUT), POINTER, DIMENSION(:) :: VEC
INTEGER, INTENT(IN):: NVECTOR, NM, NL

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_F)
      VEC => SCARC(NM)%LEVEL(NL)%F
   CASE (NSCARC_VECTOR_X)
      VEC => SCARC(NM)%LEVEL(NL)%X
   CASE (NSCARC_VECTOR_D)
      VEC => SCARC(NM)%LEVEL(NL)%D
   CASE (NSCARC_VECTOR_Y)
      VEC => SCARC(NM)%LEVEL(NL)%Y
   CASE (NSCARC_VECTOR_G)
      VEC => SCARC(NM)%LEVEL(NL)%G
   CASE (NSCARC_VECTOR_W)
      VEC => SCARC(NM)%LEVEL(NL)%W
   CASE (NSCARC_VECTOR_Z)
      VEC => SCARC(NM)%LEVEL(NL)%Z
   CASE (NSCARC_VECTOR_X2)
      VEC => SCARC(NM)%LEVEL(NL)%X2
   CASE (NSCARC_VECTOR_D2)
      VEC => SCARC(NM)%LEVEL(NL)%D2
   CASE (NSCARC_VECTOR_Y2)
      VEC => SCARC(NM)%LEVEL(NL)%Y2
   CASE (NSCARC_VECTOR_G2)
      VEC => SCARC(NM)%LEVEL(NL)%G2
   CASE (NSCARC_VECTOR_W2)
      VEC => SCARC(NM)%LEVEL(NL)%W2
   CASE (NSCARC_VECTOR_Z2)
      VEC => SCARC(NM)%LEVEL(NL)%Z2
   CASE (NSCARC_VECTOR_MEASURE)
      VEC => SCARC(NM)%LEVEL(NL)%MEASURE
END SELECT

END SUBROUTINE POINT_TO_VECTOR


!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for compact storage technique for neighbor
!> ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_OVECTOR_REAL(NVECTOR, NM, NOM, NL)
REAL(EB), POINTER, DIMENSION(:) :: POINT_TO_OVECTOR_REAL
INTEGER, INTENT(IN):: NVECTOR, NM, NOM, NL

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_X)
      POINT_TO_OVECTOR_REAL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%X
   CASE (NSCARC_VECTOR_MEASURE)
      POINT_TO_OVECTOR_REAL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%MEASURE
END SELECT

RETURN
END FUNCTION POINT_TO_OVECTOR_REAL


!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for compact storage technique for neighbor
!> ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_OVECTOR_INT(NVECTOR, NM, NOM, NL)
INTEGER, POINTER, DIMENSION(:) :: POINT_TO_OVECTOR_INT
INTEGER, INTENT(IN):: NVECTOR, NM, NOM, NL

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_CELL_TYPE)
      POINT_TO_OVECTOR_INT => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%CELL_TYPE
END SELECT

RETURN
END FUNCTION POINT_TO_OVECTOR_INT


!> ------------------------------------------------------------------------------------------------
!> Set vector to a constant value - only for debugging reasons
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_VALUE(VAL, NVECTOR, NL)
INTEGER , INTENT(IN) :: NVECTOR, NL
REAL(EB), INTENT(IN) :: VAL
REAL(EB), POINTER, DIMENSION(:) :: V
INTEGER :: NM, IC
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE
   SL => SCARC(NM)%LEVEL(NL)
   CALL POINT_TO_VECTOR(NVECTOR, NM, NL, V)
   DO IC = 1, SL%NC
     !V(IC) = VAL*NM
     V(IC) = VAL
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_VALUE


!> ------------------------------------------------------------------------------------------------
!> Compute global matrix-vector product (including data exchange along internal boundaries)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATVEC_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), POINTER, DIMENSION(:)     :: V1, V2
REAL(EB), POINTER, DIMENSION(:)     :: AC
REAL(EB) :: VCO, TNOW
INTEGER , POINTER, DIMENSION(:)     :: AC_ROW, AC_COL
INTEGER , POINTER :: NC, NCE
INTEGER :: NM, IC, JC, ICOL

TNOW = SECOND()

CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_MATVEC_PRODUCT', 'V1 INIT', NL)
CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_MATVEC_PRODUCT', 'V2 INIT', NL)

!> ------------------------------------------------------------------------------------------------
!> Exchange internal boundary values of vector1 such that the ghost values contain the corresponding
!> overlapped values of adjacent neighbor
!> ------------------------------------------------------------------------------------------------
TYPE_VECTOR = NVECTOR1
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_MATVEC_PRODUCT', 'V1 AFTER', NL)
CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_MATVEC_PRODUCT', 'V2 AFTER', NL)

!> ------------------------------------------------------------------------------------------------
!> Perform global matrix-vector product:
!> Note: - matrix already contains subdiagonal values from neighbor along internal boundaries
!>       - if vector1 contains neighboring values, then correct values of global matvec are achieved
!> ------------------------------------------------------------------------------------------------
DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE

   CALL POINT_TO_VECTOR (NVECTOR1, NM, NL, V1)
   CALL POINT_TO_VECTOR (NVECTOR2, NM, NL, V2)

   NC => SCARC(NM)%LEVEL(NL)%NCS                                      !> number of cells
   NCE => SCARC(NM)%LEVEL(NL)%NCE                                      !> number of cells

   AC     => SCARC(NM)%LEVEL(NL)%A                                    !> system matrix
   AC_ROW => SCARC(NM)%LEVEL(NL)%A_ROW                                !> row pointer
   AC_COL => SCARC(NM)%LEVEL(NL)%A_COL                                !> column pointer

   DO IC = 1, NC

      !> diagonal entry
      ICOL = AC_ROW(IC)
      JC   = AC_COL(ICOL)

      V2 (IC) = AC(ICOL)* V1(JC)
      VCO = V2(IC)
!WRITE(LU_SCARC,'(a,i4,a,5e14.6,i4)') 'MATVEC 1: V2(',IC,')=',V2(IC), 0.0_EB, V1(JC), AC(ICOL), AC(ICOL)*V1(JC), JC

      !> subdiagonal entries
      DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1
         JC = AC_COL(ICOL)
         V2(IC) =  V2(IC) + AC(ICOL)* V1(JC)
!WRITE(LU_SCARC,'(a,i4,a,5e14.6,i4)') 'MATVEC 2: V2(',IC,')=',V2(IC), VCO, V1(JC), AC(ICOL), AC(ICOL)*V1(JC), JC
         VCO = V2(IC)
      ENDDO

   ENDDO

ENDDO

CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_MATVEC_PRODUCT', 'V1 FINAL', NL)
CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_MATVEC_PRODUCT', 'V2 FINAL', NL)

TUSED_SCARC(NSCARC_TIME_MATVEC,MYID+1)=TUSED_SCARC(NSCARC_TIME_MATVEC,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
END SUBROUTINE SCARC_MATVEC_PRODUCT


!> ------------------------------------------------------------------------------------------------
!> Compute global scalarproductt (including global data exchange)
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_SCALAR_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER, INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), DIMENSION(:)    , POINTER ::  V1, V2
INTEGER , POINTER :: NC
<<<<<<< HEAD
INTEGER  :: NM, IERR, NL0! , IC
=======
INTEGER  :: NM, IERR, NL0 , IC
>>>>>>> origin/scarc
#if defined(WITH_MKL)
REAL(EB) :: DDOT
EXTERNAL :: DDOT
#endif

SP_LOCAL = 0.0_EB
SP_GROUP = 0.0_EB

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE

   CALL POINT_TO_VECTOR (NVECTOR1, NM, NL, V1)
   CALL POINT_TO_VECTOR (NVECTOR2, NM, NL, V2)

   NC  => SCARC(NM)%LEVEL(NL)%NCS

#if defined(WITH_MKL)
   SP_LOCAL(NM) = DDOT(NC, V1, 1, V2, 1)
#else   
   SP_LOCAL(NM) = 0.0_EB
   DO IC = 1, NC
      SP_LOCAL(NM) = SP_LOCAL(NM) + V1(IC) * V2(IC)
   ENDDO
#endif

SP_GROUP(MYID+1) = SP_GROUP(MYID+1) + SP_LOCAL(NM)

ENDDO


!> ------------------------------------------------------------------------------------------------
!> Get global scalarproduct by a global summation of the local values
!> ------------------------------------------------------------------------------------------------
IERR = 0
NL0  = NL
SP_GLOBAL   = 0.0_EB

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLREDUCE (SP_GROUP(MYID+1), SP_GLOBAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   SP_GLOBAL = SP_GROUP(1)
ENDIF

SCARC_SCALAR_PRODUCT = SP_GLOBAL
RETURN

END FUNCTION SCARC_SCALAR_PRODUCT


!> ------------------------------------------------------------------------------------------------
!> Compute global L2-norm (including global data exchange)
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_L2NORM(NVECTOR1, NL)
INTEGER, INTENT(IN):: NVECTOR1, NL
REAL(EB), DIMENSION(:), POINTER ::  V
INTEGER , POINTER :: NC
<<<<<<< HEAD
INTEGER  :: NM, IERR!, IC
=======
INTEGER  :: NM, IERR, IC
>>>>>>> origin/scarc
#if defined(WITH_MKL)
REAL(EB) :: DDOT
EXTERNAL :: DDOT
#endif

SP_LOCAL = 0.0_EB
SP_GROUP = 0.0_EB

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE

   CALL POINT_TO_VECTOR (NVECTOR1, NM, NL, V)
   NC => SCARC(NM)%LEVEL(NL)%NCS

#if defined(WITH_MKL)
   SP_LOCAL(NM) = DDOT(NC, V, 1, V, 1)
#else
   SP_LOCAL(NM) = 0.0_EB
   DO IC = 1, NC
      SP_LOCAL(NM) = SP_LOCAL(NM) + V(IC) * V(IC)
   ENDDO
#endif

   SP_GROUP(MYID+1) = SP_GROUP(MYID+1) + SP_LOCAL(NM)

ENDDO


!> ------------------------------------------------------------------------------------------------
!> Get global scalarproduct by a global summation of the local values
!> ------------------------------------------------------------------------------------------------
IERR = 0
SP_GLOBAL = 0.0_EB
IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLREDUCE (SP_GROUP(MYID+1), SP_GLOBAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   SP_GLOBAL = SP_GROUP(1)
ENDIF

SP_GLOBAL = SQRT (SP_GLOBAL/REAL(NCS_GLOBAL(NL), EB))

SCARC_L2NORM = SP_GLOBAL
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
INTEGER, POINTER  :: NC
#if defined(WITH_MKL)
EXTERNAL :: DAXPBY
#endif

DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   CALL POINT_TO_VECTOR(NVECTOR1, NM, NL, V1)
   CALL POINT_TO_VECTOR(NVECTOR2, NM, NL, V2)

#if defined(WITH_MKL)
   NC  => SCARC(NM)%LEVEL(NL)%NCS
   CALL DAXPBY(NC, SCAL1, V1, 1, SCAL2, V2, 1)
#else
   V2 = SCAL1 * V1 + SCAL2 * V2
#endif

ENDDO


END SUBROUTINE SCARC_VECTOR_SUM


!> ------------------------------------------------------------------------------------------------
!> Copy one integer array to another (possible scaled)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_COPY_INTEGER(X, Y, SCAL, NLEN)
INTEGER,  INTENT(IN) , DIMENSION(:) :: X
INTEGER,  INTENT(OUT), DIMENSION(:) :: Y
INTEGER,  INTENT(IN) :: NLEN, SCAL
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

END SUBROUTINE SCARC_COPY_INTEGER


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
REAL(EB), DIMENSION(:)    , POINTER ::  V1, V2
INTEGER  :: NM

#if defined(WITH_MKL)
INTEGER , POINTER :: NC
EXTERNAL :: DCOPY, DSCAL
#endif

DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   CALL POINT_TO_VECTOR(NVECTOR1, NM, NL, V1)
   CALL POINT_TO_VECTOR(NVECTOR2, NM, NL, V2)

#if defined(WITH_MKL)
   NC  => SCARC(NM)%LEVEL(NL)%NCS
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
REAL(EB), DIMENSION(:)    , POINTER ::  VC
INTEGER  :: NM

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE
   CALL POINT_TO_VECTOR(NVECTOR, NM, NL, VC)
   VC =  0.0_EB
ENDDO

END SUBROUTINE SCARC_VECTOR_CLEAR

!> ------------------------------------------------------------------------------------------------
!> Define vector to a scalar value
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_DEFINE(NVECTOR, SCAL, NL)
INTEGER , INTENT(IN):: NVECTOR, NL
REAL(EB), INTENT(IN):: SCAL
REAL(EB), DIMENSION(:)    , POINTER ::  VC
INTEGER  :: NM

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE
   CALL POINT_TO_VECTOR(NVECTOR, NM, NL, VC)
   VC =  SCAL
ENDDO

END SUBROUTINE SCARC_VECTOR_DEFINE



!> ------------------------------------------------------------------------------------------------
!> Perform preconditioning
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRECONDITIONING (NVECTOR1, NVECTOR2, NL, NSCOPE)
USE POIS, ONLY: H2CZSS, H3CZSS
INTEGER, INTENT(IN):: NVECTOR1, NVECTOR2, NL, NSCOPE
REAL(EB), DIMENSION(:,:,:), POINTER ::  FFT
REAL(EB), DIMENSION(:)    , POINTER ::  V1, V2
REAL(EB), DIMENSION(:)    , POINTER ::  AC
INTEGER , DIMENSION(:)    , POINTER ::  AC_ROW, AC_COL
INTEGER , POINTER:: NC
INTEGER  :: NPRECON, NM, I, J, K, IC, ICOL, IROW, JCOL
REAL(EB) :: AUX, OMEGA=1.5_EB
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: SM
#endif

REAL (EB) :: TNOW

TNOW = SECOND()

SELECT CASE (NSCOPE)
   CASE (NSCARC_SCOPE_MAIN)
      NPRECON = TYPE_PRECON
   CASE (NSCARC_SCOPE_PRECON)
      NPRECON = TYPE_SMOOTH
   CASE (NSCARC_SCOPE_COARSE)
!>      NPRECON = NSCARC_PRECON_SSOR
     NPRECON = NSCARC_PRECON_JACOBI
END SELECT

SELECT CASE (NPRECON)

   !! ----------------------------------------------------------------------------------------
   !! Pardiso preconditioner
   !! ----------------------------------------------------------------------------------------
#ifdef WITH_MKL
   CASE (NSCARC_PRECON_CLUSTER)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
      
         SM => SCARC(NM)%MKL(NL)
         SL => SCARC(NM)%LEVEL(NL)
      
         !.. Back substitution and iterative refinement
         SM%PHASE = 33            ! only solving
      
         CALL POINT_TO_VECTOR (NVECTOR1, NM, NL, V1)
         CALL POINT_TO_VECTOR (NVECTOR2, NM, NL, V2)
         
         NC     => SCARC(NM)%LEVEL(NL)%NCS

         AC     => SCARC(NM)%MKL(NL)%ASYM
         AC_ROW => SCARC(NM)%MKL(NL)%ASYM_ROW
         AC_COL => SCARC(NM)%MKL(NL)%ASYM_COL

         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
            WRITE(LU_SCARC,*) 'PRECON_CLUSTER'
            !WRITE(LU_SCARC,*) 'PT=',SM%PT
            WRITE(LU_SCARC,*) 'NC    =',NC
            WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
            WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
            WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
            WRITE(LU_SCARC,*) 'PHASE=',SM%PHASE
            WRITE(LU_SCARC,*) 'PERM=',SM%PERM
            WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
            WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
            WRITE(LU_SCARC,*) 'ERROR2=',SM%ERROR
            DO IROW=1,NC
               DO JCOL=AC_ROW(IROW),AC_ROW(IROW)-1
                  WRITE(LU_SCARC,'(2I8,F18.12)') IROW,AC_COL(JCOL),AC(JCOL)
               ENDDO
               WRITE(LU_SCARC,*)
            ENDDO
            WRITE(LU_SCARC,*)
            DO IROW=1,NC
               DO JCOL=AC_ROW(IROW),AC_ROW(IROW)-1
                  IF (IROW == AC_COL(JCOL)) WRITE(LU_SCARC,'(2I8,F18.12)') IROW,AC_COL(JCOL),AC(JCOL)
               ENDDO
               WRITE(LU_SCARC,*)
            ENDDO
         ENDIF
         
         CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_PRECON_CLUSTER', 'X INIT', NL)
         CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_PRECON_CLUSTER', 'F INIT', NL)

         CALL CLUSTER_SPARSE_SOLVER(SM%CT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, NCS_GLOBAL(NL), &
                                    AC, AC_ROW, AC_COL, SM%PERM, SM%NRHS, SM%IPARM, &
                                    SM%MSGLVL, V1, V2, MPI_COMM_WORLD, SM%ERROR)

         IF (SM%ERROR /= 0) THEN
            WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR
         ENDIF
      
         CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_PRECON_CLUSTER', 'X FINAL', NL)
         CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_PRECON_CLUSTER', 'F FINAL', NL)

      ENDDO

   CASE (NSCARC_PRECON_PARDISO)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
      
         SM => SCARC(NM)%MKL(NL)
         SL => SCARC(NM)%LEVEL(NL)
      
         !.. Back substitution and iterative refinement
         SM%PHASE = 33            ! only solving
      
         CALL POINT_TO_VECTOR (NVECTOR1, NM, NL, V1)
         CALL POINT_TO_VECTOR (NVECTOR2, NM, NL, V2)
         
         NC     => SCARC(NM)%LEVEL(NL)%NCS

         AC     => SCARC(NM)%MKL(NL)%ASYM
         AC_ROW => SCARC(NM)%MKL(NL)%ASYM_ROW
         AC_COL => SCARC(NM)%MKL(NL)%ASYM_COL

         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
            WRITE(LU_SCARC,*) 'PRECON_PARDISO'
            !WRITE(LU_SCARC,*) 'PT=',SM%PT
            WRITE(LU_SCARC,*) 'NC    =',NC
            WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
            WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
            WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
            WRITE(LU_SCARC,*) 'PHASE=',SM%PHASE
            WRITE(LU_SCARC,*) 'SL%NC=',NC
            WRITE(LU_SCARC,*) 'PERM=',SM%PERM
            WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
            WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
            WRITE(LU_SCARC,*) 'ERROR2=',SM%ERROR
            DO IROW=1,NC
               DO JCOL=AC_ROW(IROW),AC_ROW(IROW+1)-1
                   WRITE(LU_SCARC,'(2I6,F18.12)') IROW,AC_COL(JCOL),AC(JCOL)
               ENDDO
               WRITE(LU_SCARC,*)
            ENDDO
            DO IROW=1,NC
               DO JCOL=AC_ROW(IROW),AC_ROW(IROW)-1
                  IF (IROW == AC_COL(JCOL)) WRITE(LU_SCARC,'(2I8,F18.12)') IROW,AC_COL(JCOL),AC(JCOL)
               ENDDO
               WRITE(LU_SCARC,*)
            ENDDO
         ENDIF
         
         CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_PRECON_PARDISO', 'X INIT', NL)
         CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_PRECON_PARDISO', 'F INIT', NL)
         
         CALL PARDISO_D(SM%PT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, NC, AC, AC_ROW, AC_COL, &
                        SM%PERM, SM%NRHS, SM%IPARM, SM%MSGLVL, V1, V2, SM%ERROR)
      
         IF (SM%ERROR /= 0) THEN
            WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR
         ENDIF
      
         CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_PRECON_PARDISO', 'X FINAL', NL)
         CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_PRECON_PARDISO', 'F FINAL', NL)

      ENDDO
#endif

   !! ----------------------------------------------------------------------------------------
   !! Pardiso preconditioner
   !! ----------------------------------------------------------------------------------------
#ifdef WITH_MKL
   CASE (NSCARC_PRECON_CLUSTER)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
      
         SM => SCARC(NM)%MKL(NL)
         SL => SCARC(NM)%LEVEL(NL)
      
         !.. Back substitution and iterative refinement
         SM%PHASE = 33            ! only solving
      
         CALL POINT_TO_VECTOR (NVECTOR1, NM, NL, V1)
         CALL POINT_TO_VECTOR (NVECTOR2, NM, NL, V2)
         
         NC     => SCARC(NM)%LEVEL(NL)%NCS

         AC     => SCARC(NM)%MKL(NL)%ASYM
         AC_ROW => SCARC(NM)%MKL(NL)%ASYM_ROW
         AC_COL => SCARC(NM)%MKL(NL)%ASYM_COL

         IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
            WRITE(LU_SCARC,*) 'PRECON_CLUSTER'
            !WRITE(LU_SCARC,*) 'PT=',SM%PT
            WRITE(LU_SCARC,*) 'NC    =',NC
            WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
            WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
            WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
            WRITE(LU_SCARC,*) 'PHASE=',SM%PHASE
            WRITE(LU_SCARC,*) 'PERM=',SM%PERM
            WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
            WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
            WRITE(LU_SCARC,*) 'ERROR2=',SM%ERROR
            DO IROW=1,NC
               DO JCOL=AC_ROW(IROW),AC_ROW(IROW)-1
                  WRITE(LU_SCARC,'(2I8,F18.12)') IROW,AC_COL(JCOL),AC(JCOL)
               ENDDO
               WRITE(LU_SCARC,*)
            ENDDO
            WRITE(LU_SCARC,*)
            DO IROW=1,NC
               DO JCOL=AC_ROW(IROW),AC_ROW(IROW)-1
                  IF (IROW == AC_COL(JCOL)) WRITE(LU_SCARC,'(2I8,F18.12)') IROW,AC_COL(JCOL),AC(JCOL)
               ENDDO
               WRITE(LU_SCARC,*)
            ENDDO
         ENDIF
         
         CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_PRECON_CLUSTER', 'X INIT', NL)
         CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_PRECON_CLUSTER', 'F INIT', NL)

         CALL CLUSTER_SPARSE_SOLVER(SM%CT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, NCS_GLOBAL(NL), &
                                    AC, AC_ROW, AC_COL, SM%PERM, SM%NRHS, SM%IPARM, &
                                    SM%MSGLVL, V1, V2, MPI_COMM_WORLD, SM%ERROR)

         IF (SM%ERROR /= 0) THEN
            WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR
         ENDIF
      
         CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_PRECON_CLUSTER', 'X FINAL', NL)
         CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_PRECON_CLUSTER', 'F FINAL', NL)

      ENDDO

   CASE (NSCARC_PRECON_PARDISO)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
      
         SM => SCARC(NM)%MKL(NL)
         SL => SCARC(NM)%LEVEL(NL)
      
         !.. Back substitution and iterative refinement
         SM%PHASE = 33            ! only solving
      
         CALL POINT_TO_VECTOR (NVECTOR1, NM, NL, V1)
         CALL POINT_TO_VECTOR (NVECTOR2, NM, NL, V2)
         
         NC     => SCARC(NM)%LEVEL(NL)%NCS

         AC     => SCARC(NM)%MKL(NL)%ASYM
         AC_ROW => SCARC(NM)%MKL(NL)%ASYM_ROW
         AC_COL => SCARC(NM)%MKL(NL)%ASYM_COL

         IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
            WRITE(LU_SCARC,*) 'PRECON_PARDISO'
            !WRITE(LU_SCARC,*) 'PT=',SM%PT
            WRITE(LU_SCARC,*) 'NC    =',NC
            WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
            WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
            WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
            WRITE(LU_SCARC,*) 'PHASE=',SM%PHASE
            WRITE(LU_SCARC,*) 'SL%NC=',NC
            WRITE(LU_SCARC,*) 'PERM=',SM%PERM
            WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
            WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
            WRITE(LU_SCARC,*) 'ERROR2=',SM%ERROR
            DO IROW=1,NC
               DO JCOL=AC_ROW(IROW),AC_ROW(IROW+1)-1
                   WRITE(LU_SCARC,'(2I6,F18.12)') IROW,AC_COL(JCOL),AC(JCOL)
               ENDDO
               WRITE(LU_SCARC,*)
            ENDDO
            DO IROW=1,NC
               DO JCOL=AC_ROW(IROW),AC_ROW(IROW)-1
                  IF (IROW == AC_COL(JCOL)) WRITE(LU_SCARC,'(2I8,F18.12)') IROW,AC_COL(JCOL),AC(JCOL)
               ENDDO
               WRITE(LU_SCARC,*)
            ENDDO
         ENDIF
         
         CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_PRECON_PARDISO', 'X INIT', NL)
         CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_PRECON_PARDISO', 'F INIT', NL)
         
         CALL PARDISO_D(SM%PT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, NC, AC, AC_ROW, AC_COL, &
                        SM%PERM, SM%NRHS, SM%IPARM, SM%MSGLVL, V1, V2, SM%ERROR)
      
         IF (SM%ERROR /= 0) THEN
            WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR
         ENDIF
      
         CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_PRECON_PARDISO', 'X FINAL', NL)
         CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_PRECON_PARDISO', 'F FINAL', NL)

      ENDDO
#endif

 !! ----------------------------------------------------------------------------------------
 !! Multigrid preconditioner
 !! ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID (NSCARC_SCOPE_PRECON, NVECTOR2)


 !! ----------------------------------------------------------------------------------------
 !! Jacobi preconditioner
 !! ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_JACOBI)

      DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE
         CALL POINT_TO_VECTOR(NVECTOR2, NM, NL, V2)

         NC     => SCARC(NM)%LEVEL(NL)%NCS
         AC     => SCARC(NM)%LEVEL(NL)%A
         AC_ROW => SCARC(NM)%LEVEL(NL)%A_ROW

         DO IC = 1, NC
            V2 (IC) = V2 (IC) / AC (AC_ROW(IC))
!WRITE(LU_SCARC,'(a,i4,a,2f16.6,i4)') 'JACOBI: V2(',IC,')=',V2(IC),AC(AC_ROW(IC)), AC_ROW(IC)
         ENDDO

      ENDDO

 !! ----------------------------------------------------------------------------------------
 !! SSOR preconditioner
 !! ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_SSOR)

      DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE
         CALL POINT_TO_VECTOR(NVECTOR2, NM, NL, V2)

         NC     => SCARC(NM)%LEVEL(NL)%NCS
         AC     => SCARC(NM)%LEVEL(NL)%A
         AC_ROW => SCARC(NM)%LEVEL(NL)%A_ROW
         AC_COL => SCARC(NM)%LEVEL(NL)%A_COL

         !> use only matrix superdiagonals
         FORWARD_CELL_LOOP: DO IC = 1, NC

            AUX = 0.0_EB
            LOWER_DIAG_LOOP: DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1
               IF (AC_COL(ICOL) >= IC) CYCLE LOWER_DIAG_LOOP
               IF (AC_COL(ICOL) <= NC) AUX = AUX + AC(ICOL) * V2(AC_COL(ICOL))
            ENDDO LOWER_DIAG_LOOP
            V2(IC) = (V2(IC) - AUX * OMEGA) / AC(AC_ROW(IC))
!WRITE(LU_SCARC,'(a,i4,a,2f16.6,i4)') 'SSOR1: V2(',IC,')=',V2(IC),AC(AC_ROW(IC)), AC_ROW(IC)

         ENDDO FORWARD_CELL_LOOP

         !> use only matrix subdiagonals
         BACKWARD_CELL_LOOP: DO IC = NC-1, 1, -1

            AUX = 0.0_EB
            UPPER_DIAG_LOOP: DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1
               IF (AC_COL(ICOL) <= IC) CYCLE UPPER_DIAG_LOOP
               IF (AC_COL(ICOL) <= NC) AUX = AUX + AC(ICOL) * V2(AC_COL(ICOL))
            ENDDO UPPER_DIAG_LOOP
            V2(IC) = V2(IC) - AUX * OMEGA / AC(AC_ROW(IC))
!WRITE(LU_SCARC,'(a,i4,a,2f16.6,i4)') 'SSOR2: V2(',IC,')=',V2(IC),AC(AC_ROW(IC)), AC_ROW(IC)

         ENDDO BACKWARD_CELL_LOOP

      ENDDO

 !! ----------------------------------------------------------------------------------------
 !! FFT preconditioner
 !! ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_FFT)

      DO NM = 1, NMESHES

         IF (PROCESS(NM) /= MYID) CYCLE
         CALL POINT_TO_VECTOR(NVECTOR1, NM, NL, V1)
         CALL POINT_TO_VECTOR(NVECTOR2, NM, NL, V2)

         M   => MESHES(NM)
         FFT => SCARC(NM)%PRECON(NLEVEL_MIN)%FFT

         DO K = 1, M%KBAR
            DO J = 1, M%JBAR
               DO I = 1, M%IBAR
                  IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                  FFT(I, J, K) = V1(IC)
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
                  V2(IC) = FFT(I, J, K)
               ENDDO
            ENDDO
         ENDDO

      ENDDO

END SELECT

TUSED_SCARC(NSCARC_TIME_PRECON,MYID+1)=TUSED_SCARC(NSCARC_TIME_PRECON,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW

END SUBROUTINE SCARC_PRECONDITIONING

!> ------------------------------------------------------------------------------------------------
!> Save settings of calling parent-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SAVE_PARENT(PARENT)
TYPE (SCARC_PARENT_TYPE), INTENT(OUT):: PARENT
PARENT%TYPE_METHOD   = TYPE_METHOD
PARENT%TYPE_SCOPE    = TYPE_SCOPE
PARENT%TYPE_PRECON   = TYPE_PRECON
PARENT%TYPE_SMOOTH   = TYPE_SMOOTH
PARENT%TYPE_CYCLE    = TYPE_CYCLE
PARENT%TYPE_ACCURACY = TYPE_ACCURACY
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) '------------------------------------------------'
   WRITE(LU_SCARC,*) 'SAVING PARENT:'
   WRITE(LU_SCARC,*) 'TYPE_METHOD=', PARENT%TYPE_METHOD
   WRITE(LU_SCARC,*) 'TYPE_SCOPE =', PARENT%TYPE_SCOPE
   WRITE(LU_SCARC,*) 'TYPE_PRECON=', PARENT%TYPE_PRECON
   WRITE(LU_SCARC,*) 'TYPE_SMOOTH=', PARENT%TYPE_SMOOTH
   WRITE(LU_SCARC,*) 'TYPE_CYCLE =', PARENT%TYPE_CYCLE
   WRITE(LU_SCARC,*) '------------------------------------------------'
ENDIF
END SUBROUTINE SCARC_SAVE_PARENT


!> ------------------------------------------------------------------------------------------------
!> Reset settings of calling parent-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESET_PARENT(PARENT)
TYPE (SCARC_PARENT_TYPE), INTENT(IN):: PARENT
IF ((TYPE_METHOD==NSCARC_METHOD_KRYLOV).AND.&
     TYPE_SCOPE ==NSCARC_SCOPE_MAIN) THEN
   TYPE_METHOD = TYPE_METHOD0
   TYPE_PRECON = TYPE_PRECON0
   TYPE_SCOPE  = TYPE_SCOPE0
ELSE
   TYPE_METHOD = PARENT%TYPE_METHOD
   TYPE_PRECON = PARENT%TYPE_PRECON
   TYPE_SCOPE  = PARENT%TYPE_SCOPE
ENDIF
TYPE_SMOOTH   = PARENT%TYPE_SMOOTH
TYPE_CYCLE    = PARENT%TYPE_CYCLE
TYPE_ACCURACY = PARENT%TYPE_ACCURACY
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,*) '------------------------------------------------'
   WRITE(LU_SCARC,*) 'RESETTING PARENT:'
   WRITE(LU_SCARC,*) 'TYPE_METHOD=', PARENT%TYPE_METHOD
   WRITE(LU_SCARC,*) 'TYPE_SCOPE =', PARENT%TYPE_SCOPE
   WRITE(LU_SCARC,*) 'TYPE_PRECON=', PARENT%TYPE_PRECON
   WRITE(LU_SCARC,*) 'TYPE_SMOOTH=', PARENT%TYPE_SMOOTH
   WRITE(LU_SCARC,*) 'TYPE_CYCLE =', PARENT%TYPE_CYCLE
   WRITE(LU_SCARC,*) '------------------------------------------------'
ENDIF
END SUBROUTINE SCARC_RESET_PARENT


!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
#ifdef WITH_MKL
SUBROUTINE SCARC_METHOD_CLUSTER(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER ::  NM, NL, IC, IP, IX, IY, IZ
REAL (EB) :: TNOW
REAL(EB), POINTER, DIMENSION(:) :: XS, FS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL
TYPE (SCARC_MKL_TYPE)   , POINTER :: SM
TYPE (SCARC_SCOPE_TYPE) :: MKL
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW = SECOND()

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'CALLING METHOD_CLUSTER'

CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(MKL, TYPE_MKL, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE

   SL => SCARC(NM)%LEVEL(NL)
   SM => SCARC(NM)%MKL(NL)


   !.. Back substitution and iterative refinement
   !CALL SCARC_SETUP_VALUE( 1.0_EB, MKL%F, NLEVEL_MIN)

   CALL POINT_TO_VECTOR (MKL%X, NM, NL, XS)
   CALL POINT_TO_VECTOR (MKL%F, NM, NL, FS)

   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      !WRITE(LU_SCARC,*) 'PT=',SM%PT
      WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
      WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
      WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
      WRITE(LU_SCARC,*) 'SL%NCS=',SL%NCS
      WRITE(LU_SCARC,*) 'PERM=',SM%PERM
      WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
      WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
      WRITE(LU_SCARC,*) 'ERROR3=',SM%ERROR
      WRITE(LU_SCARC,*) 'SL%NCS=', SL%NCS
      WRITE(LU_SCARC,*) 'LN4: SL%F='
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,CGSC)'
      DO IZ = SL%NZ+1,0, -1
      IF (SL%NX == 8) THEN
      WRITE(LU_SCARC,'(10I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=SL%NY+1,0, -1)
      ELSE
      WRITE(LU_SCARC,'(6I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=SL%NY+1,0, -1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,UNKH)'
      DO IZ = SL%NZ+1,0, -1
      IF (SL%NX == 8) THEN
      WRITE(LU_SCARC,'(10I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=SL%NY+1,0, -1)
      ELSE
      WRITE(LU_SCARC,'(6I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=SL%NY+1,0, -1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
      WRITE(LU_SCARC,'(4F18.6)') (FS(IC)/100.0_EB,IC=1,SL%NCS)
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') (XS(IC),IC=1,SL%NCS)
      DO IC=1,SL%NCS
         DO IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1
             WRITE(LU_SCARC,'(2I8,F18.12)') IC,SL%A_COL(IP),SL%A(IP)/100.0_EB
         ENDDO
         WRITE(LU_SCARC,*) 
      ENDDO
      WRITE(LU_SCARC,*) 'SIZE(SM%ASYM)=',SIZE(SM%ASYM)
      DO IC=1,SL%NCS
         DO IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1
             WRITE(LU_SCARC,'(2I8,F18.12)') IC,SM%ASYM_COL(IP),SM%ASYM(IP)/100.0_EB
         ENDDO
         WRITE(LU_SCARC,*) 
      ENDDO
   ENDIF


IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN

WRITE(81,*) 'CCVAR(.,.,.,CGSC)'
DO IZ = SL%NZ+1,0, -1
IF (SL%NX == 8) THEN
WRITE(81,'(10I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=1,SL%NY+1)
ELSE
WRITE(81,'(6I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=1,SL%NY+1)
ENDIF
WRITE(81,*) 
ENDDO

WRITE(82,*) 'CCVAR(.,.,.,UNKH)'
DO IZ = SL%NZ+1,0, -1
IF (SL%NX == 8) THEN
WRITE(82,'(10I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=1,SL%NY+1)
ELSE
WRITE(82,'(6I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=1,SL%NY+1)
ENDIF
WRITE(82,*) 
ENDDO
ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
WRITE(77,*) 'CALLING CLUSTER - inside'
WRITE(77,*) 'LN4: '
IF (SL%NX==8) THEN
WRITE(77,'(8E14.6)') (FS(IC)/1000.0_EB,IC=SL%NCS,1,-1)
ELSE
WRITE(77,'(4E14.6)') (FS(IC)/1000.0_EB,IC=SL%NCS,1,-1)
ENDIF
ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
WRITE(83,*) 'RHS'
DO IZ = 1, SL%NZ
DO IY = 1, SL%NY
DO IX = 1, SL%NX
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
   WRITE(LU_SCARC,*) 'CCVAR(',IX,',',IY,',',IZ,')=',SL%CCVAR(IX,IY,IZ,CGSC), SL%CCVAR(IX,IY,IZ,UNKH)
IF (SL%CCVAR(IX,IY,IZ,CGSC) /= IS_GASPHASE) CYCLE
IC = SL%CCVAR(IX,IY,IZ,UNKH)
WRITE(83,'(A,4I6,E14.6)') 'I,J,K,IC,F:', IX, IY, IZ, IC, SL%F(IC)/1000.0_EB
ENDDO
ENDDO
ENDDO
ENDIF
   
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
WRITE(79,*) 'D: IPARM: ', SM%MTYPE, SM%MAXFCT, SM%MNUM
WRITE(79,'(8i6)') SM%IPARM(1:42)
DO IC=1,SL%NCS
   DO IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1
       WRITE(80,'(2I8,F18.12)') IC,SM%ASYM_COL(IP),-SM%ASYM(IP)/1000.0_EB
   ENDDO
   WRITE(80,*) 
ENDDO
ENDIF
   
   CALL SCARC_DEBUG_LEVEL (MKL%X, 'SCARC_METHOD_CLUSTER', 'X INIT', NL)
   CALL SCARC_DEBUG_LEVEL (MKL%F, 'SCARC_METHOD_CLUSTER', 'F INIT', NL)
   
   SM%PHASE    = 33          ! only solving
   CALL CLUSTER_SPARSE_SOLVER(SM%CT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, NCS_GLOBAL(NL), &
                              SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, SM%PERM, SM%NRHS, SM%IPARM, &
                              SM%MSGLVL, FS, XS, MPI_COMM_WORLD, SM%ERROR)

   IF (SM%ERROR /= 0) THEN
      WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR
   ENDIF



IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
WRITE(79,*) 'E: IPARM: ', SM%MTYPE, SM%MAXFCT, SM%MNUM
WRITE(79,'(8i6)') SM%IPARM(1:42)
ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
WRITE(77,*) 'LN5'
IF (SL%NX==8) THEN
WRITE(77,'(8E14.6)') (-XS(IC),IC=SL%NCS,1,-1)
ELSE
WRITE(77,'(4E14.6)') (-XS(IC),IC=SL%NCS,1,-1)
ENDIF
ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
WRITE(84,*) 'SOL1'
DO IZ = 1, SL%NZ
DO IY = 1, SL%NY
DO IX = 1, SL%NX
IF (SL%CCVAR(IX,IY,IZ,CGSC) /= IS_GASPHASE) CYCLE
IC = SL%CCVAR(IX,IY,IZ,UNKH)
WRITE(84,'(A,4I6,E14.6)') 'I,J,K,IC,F:', IX, IY, IZ, IC, -SL%X(IC)
ENDDO
ENDDO
ENDDO
   
ENDIF


   CALL SCARC_DEBUG_LEVEL (MKL%X, 'SCARC_METHOD_CLUSTER', 'X FINAL', NL)
   CALL SCARC_DEBUG_LEVEL (MKL%F, 'SCARC_METHOD_CLUSTER', 'F FINAL', NL)

   IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
      CALL SCARC_FINISH_SOLVER  (NLEVEL_MIN)
      CALL SCARC_UPDATE_GHOSTCELLS (NLEVEL_MIN)
   ENDIF

   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) '------------------- FINISH ---------------------'
      WRITE(LU_SCARC,*) 'SL%F='
      WRITE(LU_SCARC,'(4F18.6)') FS
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') XS
   ENDIF


IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
WRITE(77,*) 'LN6'
IF (SL%NX==8) THEN
WRITE(77,'(8E14.6)') (-XS(IC),IC=SL%NCS,1,-1)
ELSE
WRITE(77,'(4E14.6)') (-XS(IC),IC=SL%NCS,1,-1)
ENDIF

WRITE(84,*) 'SOL2'
DO IZ = 1, SL%NZ
DO IY = 1, SL%NY
DO IX = 1, SL%NX
IF (SL%CCVAR(IX,IY,IZ,CGSC) /= IS_GASPHASE) CYCLE
IC = SL%CCVAR(IX,IY,IZ,UNKH)
WRITE(84,'(A,4I6,E14.6)') 'I,J,K,IC,F:', IX, IY, IZ, IC, -SL%X(IC)
ENDDO
ENDDO
ENDDO

ENDIF

ENDDO

CALL SCARC_RESET_PARENT(PARENT)

TUSED_SCARC(NSCARC_TIME_PARDISO,MYID+1)=TUSED_SCARC(NSCARC_TIME_PARDISO,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_CLUSTER
#endif


!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
#ifdef WITH_MKL
SUBROUTINE SCARC_METHOD_PARDISO(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER ::  NM, NL, IC, IP
REAL (EB) :: TNOW
REAL(EB), POINTER, DIMENSION(:) :: XS, FS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL
TYPE (SCARC_MKL_TYPE)   , POINTER :: SM
TYPE (SCARC_SCOPE_TYPE) :: MKL
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW = SECOND()

CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(MKL, TYPE_MKL, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   SL => SCARC(NM)%LEVEL(NL)
   SM => SCARC(NM)%MKL(NL)

   !.. Back substitution and iterative refinement
   SM%PHASE    = 33          ! only solving

   CALL POINT_TO_VECTOR (MKL%X, NM, NL, XS)
   CALL POINT_TO_VECTOR (MKL%F, NM, NL, FS)

   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      !WRITE(LU_SCARC,*) 'PT=',SM%PT
      WRITE(LU_SCARC,*) 'CALLING SCARC_METHOD_PARDISO'
      WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
      WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
      WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
      WRITE(LU_SCARC,*) 'PHASE=',SM%PHASE
      WRITE(LU_SCARC,*) 'SL%NCS=',SL%NCS
      WRITE(LU_SCARC,*) 'PERM=',SM%PERM
      WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
      WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
      WRITE(LU_SCARC,*) 'SOLVING FOR SYSTEM '
      WRITE(LU_SCARC,*) 'IPARM='
      WRITE(LU_SCARC,'(4I6)') SM%IPARM
      WRITE(LU_SCARC,*) 'SL%NCS=', SL%NCS
      WRITE(LU_SCARC,*) 'LN4: SL%F='
      WRITE(LU_SCARC,'(4F18.6)') (FS(IC),IC=1,SL%NCS)
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') (XS(IC),IC=1,SL%NCS)
      DO IC=1,SL%NCS
         DO IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1
             WRITE(LU_SCARC,'(2I8,F18.12)') IC,SM%ASYM_COL(IP),SM%ASYM(IP)
         ENDDO
         WRITE(LU_SCARC,*) 
      ENDDO
   ENDIF
   
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) 'LEN(XS) =',SIZE(XS)
      WRITE(LU_SCARC,*) 'LEN(FS) =',SIZE(FS)
   ENDIF
   
   CALL PARDISO_D(SM%PT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, SL%NCS, &
                  SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, SM%PERM, SM%NRHS, SM%IPARM, &
                  SM%MSGLVL, FS, XS, SM%ERROR)

   IF (SM%ERROR /= 0) THEN
      WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR
   WRITE(LU_SCARC,*) 'ERROR:3: ',SM%ERROR
   ENDIF

   IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
      WRITE(LU_SCARC,*) 'LN5: SL%F='
      WRITE(LU_SCARC,'(4F18.6)') (FS(IC),IC=1,SL%NCS)
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') (XS(IC),IC=1,SL%NCS)
   ENDIF

<<<<<<< HEAD
   CALL SCARC_DEBUG_LEVEL (MKL%X, 'SCARC_METHOD_PARDISO', 'X FINAL', NL)
   CALL SCARC_DEBUG_LEVEL (MKL%F, 'SCARC_METHOD_PARDISO', 'F FINAL', NL)
=======
TUSED_SCARC(NSCARC_TIME_PRECON,MYID+1)=TUSED_SCARC(NSCARC_TIME_PRECON,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
>>>>>>> origin/scarc

   IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
      CALL SCARC_FINISH_SOLVER  (NLEVEL_MIN)
      CALL SCARC_UPDATE_GHOSTCELLS (NLEVEL_MIN)
   ENDIF

<<<<<<< HEAD
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) '------------------- FINISH ---------------------'
      WRITE(LU_SCARC,*) 'SL%F='
      WRITE(LU_SCARC,'(4F18.6)') FS
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') XS
   ENDIF
=======
!> ------------------------------------------------------------------------------------------------
!> Save settings of calling parent-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SAVE_PARENT(PARENT)
TYPE (SCARC_PARENT_TYPE), INTENT(OUT):: PARENT
PARENT%TYPE_METHOD   = TYPE_METHOD
PARENT%TYPE_SCOPE    = TYPE_SCOPE
PARENT%TYPE_PRECON   = TYPE_PRECON
PARENT%TYPE_SMOOTH   = TYPE_SMOOTH
PARENT%TYPE_CYCLE    = TYPE_CYCLE
PARENT%TYPE_ACCURACY = TYPE_ACCURACY
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) '------------------------------------------------'
   WRITE(LU_SCARC,*) 'SAVING PARENT:'
   WRITE(LU_SCARC,*) 'TYPE_METHOD=', PARENT%TYPE_METHOD
   WRITE(LU_SCARC,*) 'TYPE_SCOPE =', PARENT%TYPE_SCOPE
   WRITE(LU_SCARC,*) 'TYPE_PRECON=', PARENT%TYPE_PRECON
   WRITE(LU_SCARC,*) 'TYPE_SMOOTH=', PARENT%TYPE_SMOOTH
   WRITE(LU_SCARC,*) 'TYPE_CYCLE =', PARENT%TYPE_CYCLE
   WRITE(LU_SCARC,*) '------------------------------------------------'
ENDIF
END SUBROUTINE SCARC_SAVE_PARENT
>>>>>>> origin/scarc

!WRITE(*,*) 'FINAL PARDISO: SCARC(1)%X:'
!WRITE(*,'(8f15.6)') (SCARC(1)%LEVEL(NLEVEL_MIN)%X(I), I=1, SCARC(1)%LEVEL(NL)%NC)

<<<<<<< HEAD
ENDDO

CALL SCARC_RESET_PARENT(PARENT)

TUSED_SCARC(NSCARC_TIME_PARDISO,MYID+1)=TUSED_SCARC(NSCARC_TIME_PARDISO,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_PARDISO
#endif
=======
!> ------------------------------------------------------------------------------------------------
!> Reset settings of calling parent-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESET_PARENT(PARENT)
TYPE (SCARC_PARENT_TYPE), INTENT(IN):: PARENT
IF ((TYPE_METHOD==NSCARC_METHOD_KRYLOV).AND.&
     TYPE_SCOPE ==NSCARC_SCOPE_MAIN) THEN
   TYPE_METHOD = TYPE_METHOD0
   TYPE_PRECON = TYPE_PRECON0
   TYPE_SCOPE  = TYPE_SCOPE0
ELSE
   TYPE_METHOD = PARENT%TYPE_METHOD
   TYPE_PRECON = PARENT%TYPE_PRECON
   TYPE_SCOPE  = PARENT%TYPE_SCOPE
ENDIF
TYPE_SMOOTH   = PARENT%TYPE_SMOOTH
TYPE_CYCLE    = PARENT%TYPE_CYCLE
TYPE_ACCURACY = PARENT%TYPE_ACCURACY
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) '------------------------------------------------'
   WRITE(LU_SCARC,*) 'RESETTING PARENT:'
   WRITE(LU_SCARC,*) 'TYPE_METHOD=', PARENT%TYPE_METHOD
   WRITE(LU_SCARC,*) 'TYPE_SCOPE =', PARENT%TYPE_SCOPE
   WRITE(LU_SCARC,*) 'TYPE_PRECON=', PARENT%TYPE_PRECON
   WRITE(LU_SCARC,*) 'TYPE_SMOOTH=', PARENT%TYPE_SMOOTH
   WRITE(LU_SCARC,*) 'TYPE_CYCLE =', PARENT%TYPE_CYCLE
   WRITE(LU_SCARC,*) '------------------------------------------------'
ENDIF
END SUBROUTINE SCARC_RESET_PARENT
>>>>>>> origin/scarc


!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
#ifdef WITH_MKL
SUBROUTINE SCARC_METHOD_CLUSTER(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER ::  NM, NL, IC, IP, IX, IY, IZ
REAL (EB) :: TNOW
REAL(EB), POINTER, DIMENSION(:) :: XS, FS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL
TYPE (SCARC_MKL_TYPE)   , POINTER :: SM
TYPE (SCARC_SCOPE_TYPE) :: MKL
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW = SECOND()

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'CALLING METHOD_CLUSTER'

CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(MKL, TYPE_MKL, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE

   SL => SCARC(NM)%LEVEL(NL)
   SM => SCARC(NM)%MKL(NL)


   !.. Back substitution and iterative refinement
   !CALL SCARC_SETUP_VALUE( 1.0_EB, MKL%F, NLEVEL_MIN)

   CALL POINT_TO_VECTOR (MKL%X, NM, NL, XS)
   CALL POINT_TO_VECTOR (MKL%F, NM, NL, FS)

   FS = 1.0_EB    
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      !WRITE(LU_SCARC,*) 'PT=',SM%PT
      WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
      WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
      WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
      WRITE(LU_SCARC,*) 'SL%NCS=',SL%NCS
      WRITE(LU_SCARC,*) 'PERM=',SM%PERM
      WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
      WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
      WRITE(LU_SCARC,*) 'ERROR3=',SM%ERROR
      WRITE(LU_SCARC,*) 'SL%NCS=', SL%NCS
      WRITE(LU_SCARC,*) 'LN4: SL%F='
      WRITE(LU_SCARC,'(4F18.6)') (FS(IC),IC=1,SL%NCS)
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') (XS(IC),IC=1,SL%NCS)
      WRITE(LU_SCARC,*) 'SIZE(SM%ASYM)=',SIZE(SM%ASYM)
      DO IC=1,SL%NCS
         DO IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1
             WRITE(LU_SCARC,'(2I8,F18.12)') IC,SM%ASYM_COL(IP),SM%ASYM(IP)/15.625_EB
         ENDDO
         WRITE(LU_SCARC,*) 
      ENDDO
   ENDIF

!IF (MYID == 0.AND.TYPE_DEBUG> NSCARC_DEBUG_NONE) THEN
IF (MYID == 1234) THEN
WRITE(77,*) '============================================================='
WRITE(77,*) '     SETUP_CCVAR_EXCHANGE ', NL
WRITE(77,*) '============================================================='
WRITE(77,*) 'CCVAR(.,.,.,CGSC)'
DO IZ = SL%NZ+1,1, -1
WRITE(77,'(10I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=SL%NY+1,0, -1)
WRITE(77,*) 
ENDDO
WRITE(77,*) 'CCVAR(.,.,.,UNKH)'
DO IZ = SL%NZ+1,1, -1
WRITE(77,'(10I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=SL%NY+1,0, -1)
WRITE(77,*) 
ENDDO
WRITE(77,*) 'CALLING CLUSTER - inside'
WRITE(77,*) 'LN4: '
WRITE(77,*) 'F ', SL%NCS
WRITE(77,'(8E14.6)') (FS(IC),IC=503,440,-1)
WRITE(77,*) 'F '
WRITE(77,'(8E14.6)') (FS(IC),IC=439,376,-1)
WRITE(77,*) 'F '
WRITE(77,'(8E14.6)') (FS(IC),IC=375,315,-1)
WRITE(77,*) 'F '
WRITE(77,'(8E14.6)') (FS(IC),IC=314,254,-1)
WRITE(77,*) 'F '
WRITE(77,'(8E14.6)') (FS(IC),IC=253,193,-1)
WRITE(77,*) 'F '
WRITE(77,'(8E14.6)') (FS(IC),IC=192,129,-1)
WRITE(77,*) 'F '
WRITE(77,'(8E14.6)') (FS(IC),IC=128,65,-1)
WRITE(77,*) 'F '
WRITE(77,'(8E14.6)') (FS(IC),IC=64,1,-1)
DO IC=1,SL%NCS
   DO IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1
       WRITE(77,'(2I8,F18.12)') IC,SM%ASYM_COL(IP),-SM%ASYM(IP)/15.625_EB
   ENDDO
   WRITE(77,*) 
ENDDO
ENDIF
   
   
   CALL SCARC_DEBUG_LEVEL (MKL%X, 'SCARC_METHOD_CLUSTER', 'X INIT', NL)
   CALL SCARC_DEBUG_LEVEL (MKL%F, 'SCARC_METHOD_CLUSTER', 'F INIT', NL)
   
   SM%PHASE    = 33          ! only solving
   CALL CLUSTER_SPARSE_SOLVER(SM%CT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, NCS_GLOBAL(NL), &
                              SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, SM%PERM, SM%NRHS, SM%IPARM, &
                              SM%MSGLVL, FS, XS, MPI_COMM_WORLD, SM%ERROR)

   IF (SM%ERROR /= 0) THEN
      WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR
   ENDIF

!IF (MYID == 0.AND.TYPE_DEBUG> NSCARC_DEBUG_NONE) THEN
IF (MYID == 1234) THEN
WRITE(77,*) 'LN5: SL%F='
WRITE(77,'(4E25.16)') (FS(IC),IC=1,SL%NCS)
WRITE(77,*) 'SL%X='
WRITE(77,'(4E25.16)') (XS(IC),IC=1,SL%NCS)
ENDIF

   CALL SCARC_DEBUG_LEVEL (MKL%X, 'SCARC_METHOD_CLUSTER', 'X FINAL', NL)
   CALL SCARC_DEBUG_LEVEL (MKL%F, 'SCARC_METHOD_CLUSTER', 'F FINAL', NL)

   IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
      CALL SCARC_FINISH_SOLVER  (NLEVEL_MIN)
      CALL SCARC_UPDATE_GHOSTCELLS (NLEVEL_MIN)
   ENDIF

   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) '------------------- FINISH ---------------------'
      WRITE(LU_SCARC,*) 'SL%F='
      WRITE(LU_SCARC,'(4F18.6)') FS
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') XS
   ENDIF

!IF (MYID == 0.AND.TYPE_DEBUG> NSCARC_DEBUG_NONE) THEN
IF (MYID == 1234) THEN
WRITE(77,*) 'LN6: SL%F='
WRITE(77,'(4E25.16)') FS
WRITE(77,*) 'SL%X='
WRITE(77,'(4E25.16)') XS
ENDIF

!WRITE(*,*) 'FINAL PARDISO: SCARC(1)%X:'
!WRITE(*,'(8f15.6)') (SCARC(1)%LEVEL(NLEVEL_MIN)%X(I), I=1, SCARC(1)%LEVEL(NL)%NC)

ENDDO

CALL SCARC_RESET_PARENT(PARENT)

TUSED_SCARC(NSCARC_TIME_PARDISO,MYID+1)=TUSED_SCARC(NSCARC_TIME_PARDISO,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_CLUSTER
#endif


!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
#ifdef WITH_MKL
SUBROUTINE SCARC_METHOD_PARDISO(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER ::  NM, NL, IC, IP
REAL (EB) :: TNOW
REAL(EB), POINTER, DIMENSION(:) :: XS, FS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL
TYPE (SCARC_MKL_TYPE)   , POINTER :: SM
TYPE (SCARC_SCOPE_TYPE) :: MKL
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW = SECOND()

CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(MKL, TYPE_MKL, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   SL => SCARC(NM)%LEVEL(NL)
   SM => SCARC(NM)%MKL(NL)

   !.. Back substitution and iterative refinement
   SM%PHASE    = 33          ! only solving

   CALL POINT_TO_VECTOR (MKL%X, NM, NL, XS)
   CALL POINT_TO_VECTOR (MKL%F, NM, NL, FS)

   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      !WRITE(LU_SCARC,*) 'PT=',SM%PT
      WRITE(LU_SCARC,*) 'CALLING SCARC_METHOD_PARDISO'
      WRITE(LU_SCARC,*) 'MAXFCT=',SM%MAXFCT
      WRITE(LU_SCARC,*) 'MNUM=',SM%MNUM
      WRITE(LU_SCARC,*) 'MTYPE=',SM%MTYPE
      WRITE(LU_SCARC,*) 'PHASE=',SM%PHASE
      WRITE(LU_SCARC,*) 'SL%NCS=',SL%NCS
      WRITE(LU_SCARC,*) 'PERM=',SM%PERM
      WRITE(LU_SCARC,*) 'NRHS=',SM%NRHS
      WRITE(LU_SCARC,*) 'MSGLVL=',SM%MSGLVL
      WRITE(LU_SCARC,*) 'SOLVING FOR SYSTEM '
      WRITE(LU_SCARC,*) 'IPARM='
      WRITE(LU_SCARC,'(4I6)') SM%IPARM
      WRITE(LU_SCARC,*) 'SL%NCS=', SL%NCS
      WRITE(LU_SCARC,*) 'LN4: SL%F='
      WRITE(LU_SCARC,'(4F18.6)') (FS(IC),IC=1,SL%NCS)
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') (XS(IC),IC=1,SL%NCS)
      DO IC=1,SL%NCS
         DO IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1
             WRITE(LU_SCARC,'(2I8,F18.12)') IC,SM%ASYM_COL(IP),SM%ASYM(IP)
         ENDDO
         WRITE(LU_SCARC,*) 
      ENDDO
   ENDIF
   
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) 'LEN(XS) =',SIZE(XS)
      WRITE(LU_SCARC,*) 'LEN(FS) =',SIZE(FS)
   ENDIF
   
   CALL PARDISO_D(SM%PT, SM%MAXFCT, SM%MNUM, SM%MTYPE, SM%PHASE, SL%NCS, &
                  SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, SM%PERM, SM%NRHS, SM%IPARM, &
                  SM%MSGLVL, FS, XS, SM%ERROR)

   IF (SM%ERROR /= 0) THEN
      WRITE(LU_SCARC,*) 'The following ERROR was detected: ', SM%ERROR
   WRITE(LU_SCARC,*) 'ERROR:3: ',SM%ERROR
   ENDIF

   IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) 'LN5: SL%F='
      WRITE(LU_SCARC,'(4F18.6)') (FS(IC),IC=1,SL%NCS)
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') (XS(IC),IC=1,SL%NCS)
   ENDIF

   CALL SCARC_DEBUG_LEVEL (MKL%X, 'SCARC_METHOD_PARDISO', 'X FINAL', NL)
   CALL SCARC_DEBUG_LEVEL (MKL%F, 'SCARC_METHOD_PARDISO', 'F FINAL', NL)

   IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
      CALL SCARC_FINISH_SOLVER  (NLEVEL_MIN)
      CALL SCARC_UPDATE_GHOSTCELLS (NLEVEL_MIN)
   ENDIF

   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) '------------------- FINISH ---------------------'
      WRITE(LU_SCARC,*) 'SL%F='
      WRITE(LU_SCARC,'(4F18.6)') FS
      WRITE(LU_SCARC,*) 'SL%X='
      WRITE(LU_SCARC,'(4F18.6)') XS
   ENDIF

!WRITE(*,*) 'FINAL PARDISO: SCARC(1)%X:'
!WRITE(*,'(8f15.6)') (SCARC(1)%LEVEL(NLEVEL_MIN)%X(I), I=1, SCARC(1)%LEVEL(NL)%NC)

ENDDO

CALL SCARC_RESET_PARENT(PARENT)

TUSED_SCARC(NSCARC_TIME_PARDISO,MYID+1)=TUSED_SCARC(NSCARC_TIME_PARDISO,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_PARDISO
#endif


!> ------------------------------------------------------------------------------------------------
!> Perform global CG-method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CG(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER   :: NL
INTEGER   :: ISTATE, ITE
REAL (EB) :: SIGMA0, SIGMA1, ALPHA0, GAMMA0
REAL (EB) :: TNOW
TYPE (SCARC_SCOPE_TYPE) :: CG
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW = SECOND()

!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Set environment variables and define working level
!>   - Initialize solution, right hand side vector and auxiliary vectors
!>   - Define iterations parameters
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(CG, TYPE_PRECON, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

!CALL SCARC_SETUP_VALUE( 1.0_EB, CG%X, NL)
!CALL SCARC_SETUP_VALUE(10.0_EB, CG%F, NL)

CALL SCARC_DEBUG_LEVEL (CG%X, 'SCARC_METHOD_CG', 'X INIT', NL)
CALL SCARC_DEBUG_LEVEL (CG%F, 'SCARC_METHOD_CG', 'F INIT', NL)

!CALL SCARC_DUMP_QUANTITY (CG%X, NL, ITE_MAIN, 0,'.sol')
!CALL SCARC_DUMP_QUANTITY (CG%F, NL, ITE_MAIN, 0,'.rhs')

!> ------------------------------------------------------------------------------------------------
!> Compute initial residual and perform initial preconditioning
!> ------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (CG%X, CG%W, NL)                                !>  W := A*X
CALL SCARC_VECTOR_SUM     (CG%F, CG%W, -1.0_EB, 1.0_EB, NL)               !>  W := W - F

CALL SCARC_DEBUG_LEVEL (CG%X, 'SCARC_METHOD_CG', 'X INIT2', NL)
CALL SCARC_DEBUG_LEVEL (CG%F, 'SCARC_METHOD_CG', 'F INIT2', NL)
CALL SCARC_DEBUG_LEVEL (CG%W, 'SCARC_METHOD_CG', 'W INIT2', NL)

CG%RES = SCARC_L2NORM (CG%W, NL)                                          !>  RESIN := ||W||
CG%RESIN = CG%RES

ISTATE = SCARC_CONVERGENCE_STATE (CG, 0, NL)                            !>  RES < TOL ??
CALL SCARC_VECTOR_COPY     (CG%W, CG%G, 1.0_EB, NL)                       !>  G := W
CALL SCARC_PRECONDITIONING (CG%W, CG%G, NL, NSCOPE)                       !>  G := PRECON(W)

CALL SCARC_DEBUG_LEVEL (CG%W, 'SCARC_METHOD_CG', 'W PRE_SCALPROD', NL)
CALL SCARC_DEBUG_LEVEL (CG%G, 'SCARC_METHOD_CG', 'G PRE_SCALPROD', NL)

SIGMA0 = SCARC_SCALAR_PRODUCT(CG%W, CG%G, NL)                             !>  SIGMA0 := (W,G)

CALL SCARC_VECTOR_COPY (CG%G, CG%D, -1.0_EB, NL)                          !>  D := -G

CALL SCARC_DEBUG_LEVEL (CG%W, 'SCARC_METHOD_CG', 'W PRECON INIT', NL)
CALL SCARC_DEBUG_LEVEL (CG%G, 'SCARC_METHOD_CG', 'G PRECON INIT', NL)
CALL SCARC_DEBUG_LEVEL (CG%D, 'SCARC_METHOD_CG', 'D PRECON INIT', NL)

!> ------------------------------------------------------------------------------------------------
!> Perform conjugate gradient looping
!> ------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, CG%NIT

   CALL SCARC_MATVEC_PRODUCT (CG%D, CG%Y, NL)                             !>  Y := A*D

CALL SCARC_DEBUG_LEVEL (CG%X, 'SCARC_METHOD_CG', 'X LOOP0', NL)
CALL SCARC_DEBUG_LEVEL (CG%D, 'SCARC_METHOD_CG', 'D LOOP0', NL)
CALL SCARC_DEBUG_LEVEL (CG%G, 'SCARC_METHOD_CG', 'G LOOP0', NL)
CALL SCARC_DEBUG_LEVEL (CG%Y, 'SCARC_METHOD_CG', 'Y LOOP0', NL)

   ALPHA0 = SCARC_SCALAR_PRODUCT (CG%D, CG%Y, NL)                         !>  ALPHA0 := (D,Y)
   ALPHA0 = SIGMA0/ALPHA0

!WRITE(LU_SCARC,*) 'ALPHA0=',ALPHA0,': SIGMA0=',SIGMA0, ': NL=',NL

   CALL SCARC_VECTOR_SUM (CG%D, CG%X, ALPHA0, 1.0_EB, NL)                 !>  X := ALPHA0*D + X
   CALL SCARC_VECTOR_SUM (CG%Y, CG%W, ALPHA0, 1.0_EB, NL)                 !>  W := ALPHA0*Y + W

CALL SCARC_DEBUG_LEVEL (CG%D, 'SCARC_METHOD_CG', 'D LOOP', NL)
CALL SCARC_DEBUG_LEVEL (CG%X, 'SCARC_METHOD_CG', 'X LOOP', NL)
CALL SCARC_DEBUG_LEVEL (CG%W, 'SCARC_METHOD_CG', 'W LOOP', NL)
!CALL SCARC_DUMP_QUANTITY (CG%X, NL, ITE_MAIN, ITE,'.sol')

   CG%RES = SCARC_L2NORM (CG%W, NL)                                       !>  RES := ||W||
   ISTATE = SCARC_CONVERGENCE_STATE (CG, ITE, NL)                         !>  RES < TOL ??
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP

!IF (TYPE_DEBUG>NSCARC_DEBUG_INFO1.AND.MYID==0) &
!   WRITE(*,'(a,i3,a,e14.5)') ' CG  -Iteration  =',ITE,': Residuum=',SCARC_RESIDUAL

   CALL SCARC_VECTOR_COPY     (CG%W, CG%G, 1.0_EB, NL)                    !>  G := W
   CALL SCARC_PRECONDITIONING (CG%W, CG%G, NL, NSCOPE)                    !>  G := PRECON(W)

CALL SCARC_DEBUG_LEVEL (CG%G, 'SCARC_METHOD_CG', 'G PRECON ', NL)
CALL SCARC_DEBUG_LEVEL (CG%W, 'SCARC_METHOD_CG', 'W PRECON ', NL)

   SIGMA1 = SCARC_SCALAR_PRODUCT (CG%W, CG%G, NL)                         !>  SIGMA1 := (W,G)
   GAMMA0 = SIGMA1/SIGMA0
   SIGMA0 = SIGMA1

   CALL SCARC_VECTOR_SUM (CG%G, CG%D, -1.0_EB, GAMMA0, NL)                !>  D := -G + GAMMA0*D

!   WRITE(LU_SCARC,*) 'SIGMA0=',SIGMA0,': SIGMA1=',SIGMA1,': GAMMA0=',GAMMA0
CALL SCARC_DEBUG_LEVEL (CG%W, 'SCARC_METHOD_CG', 'D VECTORSUM ', NL)

ENDDO CG_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of CG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(CG, ITE, ISTATE)
!IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS.AND.MYID==0) &
!>  WRITE(0,'(a,e14.5)') '                                       !> ---->  Konvergenzrate=',SCARC_CAPPA


IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_FINISH_SOLVER  (NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS (NLEVEL_MIN)
!ELSE IF (TYPE_SCOPE == NSCARC_SCOPE_COARSE) THEN
!>  CALL SCARC_VECTOR_COPY (CG%X, CG%F, 1.0_EB, NL)
ENDIF

CALL SCARC_DEBUG_LEVEL (CG%X, 'SCARC_METHOD_CG', 'X FINAL', NL)
CALL SCARC_DEBUG_LEVEL (CG%F, 'SCARC_METHOD_CG', 'F FINAL', NL)

!CALL SCARC_DUMP_QUANTITY (CG%X, NL, ITE_MAIN, ITE,'.sol')
!CALL SCARC_DUMP_QUANTITY (CG%F, NL, ITE_MAIN, ITE,'.rhs')

CALL SCARC_RESET_PARENT(PARENT)

!WRITE(*,*) 'FINAL CG: SCARC(1)%X:'
!WRITE(*,'(8f15.6)') (SCARC(1)%LEVEL(NLEVEL_MIN)%X(ITE), ITE=1, SCARC(1)%LEVEL(NL)%NC)

TUSED_SCARC(NSCARC_TIME_KRYLOV,MYID+1)=TUSED_SCARC(NSCARC_TIME_KRYLOV,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
END SUBROUTINE SCARC_METHOD_CG


!> ------------------------------------------------------------------------------------------------
!> Perform global BICGstab-method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_BICG(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER   :: NL = NSCARC_LEVEL_NONE
INTEGER   :: ISTATE, ITE
REAL (EB) :: ALPHA0, ALPHA1, ALPHA2, RHO0, RHO1, DTHETA, DBETA
REAL (EB) :: TNOW
TYPE (SCARC_SCOPE_TYPE)  :: BICG
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW = SECOND()

!> ------------------------------------------------------------------------------------------------
!> Initialization
!>   - Set environment variables and define working level
!>   - Initialize solution, right hand side vector and auxiliary vectors
!>   - Define iterations parameters
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(BICG, TYPE_PRECON, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

!> ------------------------------------------------------------------------------------------------
!> Compute initial defect and perform (double) initial preconditioning
!> ------------------------------------------------------------------------------------------------
ALPHA0 = 1.0_EB
RHO0   = 1.0_EB
RHO1   = 0.0_EB
DBETA  = 0.0_EB
DTHETA = 1.0_EB

!CALL SCARC_VECTOR_COPY    (BICG%F, BICG%W, 1.0_EB, NL)                         !>  W := F
!CALL SCARC_PRECONDITIONING (BICG%W, BICG%W, NL, NSCOPE)                        !>  W := PRECON(W)
CALL SCARC_MATVEC_PRODUCT  (BICG%X, BICG%W, NL)                                 !>  W := A*X
CALL SCARC_VECTOR_SUM      (BICG%F, BICG%W, 1.0_EB, -1.0_EB, NL)                !>  W := F - W
CALL SCARC_PRECONDITIONING (BICG%W, BICG%W, NL, NSCOPE)                         !>  W := PRECON(W)

BICG%RESIN = SCARC_L2NORM (BICG%W, NL)                                          !>  RESIN := ||W||
CALL SCARC_CONVERGENCE_INFO (BICG%RESIN, 0, NL, BICG%CROUTINE)

CALL SCARC_VECTOR_COPY (BICG%W, BICG%G, 1.0_EB, NL)                             !>  G := W

!> ------------------------------------------------------------------------------------------------
!> Perform bi-conjugate gradient looping:
!> ------------------------------------------------------------------------------------------------
BICG_LOOP: DO ITE = 1, BICG%NIT

   RHO1  = SCARC_SCALAR_PRODUCT (BICG%G, BICG%W, NL)                            !> RHO1 := (G,W)
   DBETA = (RHO1*DTHETA)/(RHO0*ALPHA0)
   RHO0  = RHO1

   CALL SCARC_VECTOR_SUM      (BICG%W, BICG%Z, 1.0_EB       , DBETA , NL)       !> Z := W + DBETA*Z
   CALL SCARC_VECTOR_SUM      (BICG%Y, BICG%Z, -DBETA*ALPHA0, 1.0_EB, NL)       !> Z := -DBETA*ALPHA0*Y + Z
   CALL SCARC_MATVEC_PRODUCT  (BICG%Z, BICG%Y, NL)                              !> Y := A*Z
   CALL SCARC_PRECONDITIONING (BICG%Y, BICG%Y, NL, NSCOPE)                      !> Z := PRECON(Z)

   DTHETA = SCARC_SCALAR_PRODUCT (BICG%G, BICG%Y, NL)                           !> DTHETA := (G,Y)
   DTHETA = RHO1/DTHETA

   CALL SCARC_VECTOR_SUM      (BICG%Y, BICG%W, -DTHETA, 1.0_EB, NL)             !> W := -DTHETA*Y + W
   CALL SCARC_MATVEC_PRODUCT  (BICG%W, BICG%D, NL)                              !> D := A*W
   CALL SCARC_PRECONDITIONING (BICG%D, BICG%D, NL, NSCOPE)                      !> D := PRECON(D)

   ALPHA1 = SCARC_SCALAR_PRODUCT (BICG%D, BICG%W, NL)                           !> ALPHA1 := (D,W)
   ALPHA2 = SCARC_SCALAR_PRODUCT (BICG%D, BICG%D, NL)                           !> ALPHA2 := (D,D)
   ALPHA0 = ALPHA1/ALPHA2

   CALL SCARC_VECTOR_SUM (BICG%Z, BICG%X,  DTHETA, 1.0_EB, NL)                  !> X :=  DTHETA*Z + X
   CALL SCARC_VECTOR_SUM (BICG%W, BICG%X,  ALPHA0, 1.0_EB, NL)                  !> X :=  ALPHA0*W + X
   CALL SCARC_VECTOR_SUM (BICG%D, BICG%W, -ALPHA0, 1.0_EB, NL)                  !> W := -ALPHA0*D + W

   BICG%RES = SCARC_L2NORM (BICG%W, NL)                                         !> RES := ||W||

   ISTATE = SCARC_CONVERGENCE_STATE(BICG, ITE, NL)                              !> RES < TOL ???
IF (TYPE_DEBUG>NSCARC_DEBUG_INFO1.AND.MYID==0) &
   WRITE(0,'(a,i3,a,e14.5,a,e14.5)') ' BICG-Iteration  =',ITE,': Residuum=',SCARC_RESIDUAL
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT BICG_LOOP

ENDDO BICG_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of BICG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(BICG, ITE, ISTATE)

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_FINISH_SOLVER(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
!ELSE IF (TYPE_SCOPE == NSCARC_SCOPE_COARSE) THEN
!>  CALL SCARC_VECTOR_COPY (BICG%X, BICG%F, 1.0_EB, NL)
ENDIF

CALL SCARC_DEBUG_LEVEL (BICG%X, 'SCARC_METHOD_BICG', 'X FINAL', NL)
CALL SCARC_DEBUG_LEVEL (BICG%F, 'SCARC_METHOD_BICG', 'F FINAL', NL)

CALL SCARC_RESET_PARENT(PARENT)

TUSED_SCARC(NSCARC_TIME_KRYLOV,MYID+1)=TUSED_SCARC(NSCARC_TIME_KRYLOV,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
END SUBROUTINE SCARC_METHOD_BICG


!> ------------------------------------------------------------------------------------------------
!> Perform geometric multigrid method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MULTIGRID(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER   :: NL = NSCARC_LEVEL_NONE
INTEGER   :: ISTATE, ICYCLE, ITE
REAL (EB) :: TNOW
TYPE (SCARC_SCOPE_TYPE)  :: MG
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW = SECOND()

<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) '==================================='
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'STARTING NEW MG METHOD'
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) '==================================='
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '==================================='
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'STARTING NEW MG METHOD'
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '==================================='
>>>>>>> origin/scarc
!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Set environment variables and define working level
!>   - Initialize solution, right hand side vector
!>   - Define iterations parameters (NL is set to finest level)
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(MG, TYPE_SMOOTH, NSCOPE, NVECTOR, NL)
CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MG', 'F init0 ', NL)
CALL SCARC_SETUP_WORKSPACE(NL)

CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MG', 'F init ', NL)
CALL SCARC_DEBUG_LEVEL (MG%X, 'SCARC_METHOD_MG', 'X init ', NL)


!> ------------------------------------------------------------------------------------------------
!> Compute initial defect:  RESIN := || F - A*X ||
!>   - Initialize cycle counts for MG-iteration
!>   - Perform initial matrix-vector product on finest level
!>   - calculate norm of initial residual on finest level
!> ------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (MG%X, MG%D, NL)                                       !>  D := A*X
CALL SCARC_VECTOR_SUM     (MG%F, MG%D, 1.0_EB, -1.0_EB, NL)                      !>  D := F - D

CALL SCARC_DEBUG_LEVEL (MG%D, 'SCARC_METHOD_MG', 'D init ', NL)

ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_SETUP, NL)
MG%RESIN = SCARC_L2NORM (MG%D, NL)                                               !>  RESIN := ||D||

CALL SCARC_CONVERGENCE_INFO(MG%RESIN, 0, NL, MG%CROUTINE)

MG%NIT=10    !!! ONLY TEMPORARILY

!> ------------------------------------------------------------------------------------------------

!> Perform multigrid-looping (start each iteration on finest level)
!> ------------------------------------------------------------------------------------------------
MULTIGRID_LOOP: DO ITE = 1, MG%NIT

   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLE_EXIT)

      !> presmoothing  (smoothing/restriction till coarsest level is reached)
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)
         CALL SCARC_SMOOTHING    (MG, NSCARC_CYCLE_PRESMOOTH, NL, NSCOPE)        !> D_fine   := smooth(defect)
         CALL SCARC_RESTRICTION  (MG%D, MG%F, NL, NL+1)                          !> F_coarse := rest(D_fine)
CALL SCARC_DEBUG_LEVEL (MG%D, 'SCARC_METHOD_MG', 'D after restrict', NL)
CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MG', 'F after restrict', NL+1)
         CALL SCARC_VECTOR_CLEAR (MG%X, NL+1)                                    !> X_coarse := 0.0
         NL = NL + 1                                                             !> set coarser level
      ENDDO PRESMOOTHING_LOOP

CALL SCARC_DEBUG_LEVEL (MG%X, 'SCARC_METHOD_MG', 'X before coarse', NL)
CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MG', 'F before coarse', NL)

      !> coarse grid solver
      SELECT CASE (TYPE_COARSE)
         CASE (NSCARC_COARSE_ITERATIVE)
IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) WRITE(LU_SCARC,*) 'CALLING ITERATIVE COARSE GRID SOLVER'
            CALL SCARC_METHOD_CG (NSCARC_SCOPE_COARSE, MG%F)                     !> X_coarse := exact_sol(.)
         CASE (NSCARC_COARSE_DIRECT)
IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) WRITE(LU_SCARC,*) 'CALLING DIRECT COARSE GRID SOLVER'
            !CALL SCARC_METHOD_DIRECT (MG%X, MG%F)
            IF (N_MPI_PROCESSES > 1) THEN
               CALL SCARC_METHOD_CLUSTER (NSCARC_SCOPE_COARSE, MG%F)             !> call CLUSTER_SPARSE_SOLVER
            ELSE
               CALL SCARC_METHOD_PARDISO (NSCARC_SCOPE_COARSE, MG%F)             !> call PARDISO
            ENDIF
      END SELECT

CALL SCARC_DEBUG_LEVEL (MG%X, 'SCARC_METHOD_MG', 'X after coarse', NL)
CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MG', 'F after coarse', NL)

      !> postsmoothing (smoothing/restriction till finest level is reached again)
      POSTSMOOTHING_LOOP: DO WHILE (NL > NLEVEL_MIN)
         NL=NL-1
         CALL SCARC_PROLONGATION (MG%X, MG%D, NL+1, NL)                          !> D_fine := prol(X_coarse)
CALL SCARC_DEBUG_LEVEL (MG%D, 'SCARC_METHOD_MG', 'D after prol', NL)
         CALL SCARC_VECTOR_SUM   (MG%D, MG%X, 1.0_EB, 1.0_EB, NL)                !> X_fine := D_fine + X_fine
         CALL SCARC_SMOOTHING    (MG, NSCARC_CYCLE_POSTSMOOTH, NL, NSCOPE)       !> D_fine := smooth(defect)
         ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_PROCEED, NL)                  !> perform requested cycle
         IF (ICYCLE /= NSCARC_CYCLE_POSTSMOOTH) CYCLE CYCLE_LOOP
      ENDDO POSTSMOOTHING_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) THEN
      WRITE(SCARC_MESSAGE,'(2A,I8,A)') TRIM(SCARC_ROUTINE),': Wrong level ', NL,' for multigrid method'
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   ENDIF

 !! ---------------------------------------------------------------------------------------------
 !! Compute norm of new residual on finest level and  leave loop correspondingly
 !! ---------------------------------------------------------------------------------------------
   CALL SCARC_MATVEC_PRODUCT (MG%X, MG%D, NL)                                    !> D := A*X
   CALL SCARC_VECTOR_SUM     (MG%F, MG%D, 1.0_EB, -1.0_EB, NL)                   !> D := F - D

   MG%RES = SCARC_L2NORM (MG%D, NL)                                              !> RES := ||D||
   ISTATE = SCARC_CONVERGENCE_STATE(MG, ITE, NL)                                 !> convergence ?
IF (TYPE_DEBUG>NSCARC_DEBUG_INFO1.AND.MYID==0) &
   WRITE(0,'(a,i3,a,e14.5,a,e14.5)') ' MG  -Iteration  =',ITE,': Residuum=',SCARC_RESIDUAL
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT MULTIGRID_LOOP

ENDDO MULTIGRID_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of MG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries (consistency!)
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(MG, ITE, ISTATE)

SELECT CASE (TYPE_SCOPE)
   CASE (NSCARC_SCOPE_MAIN)
      CALL SCARC_FINISH_SOLVER(NLEVEL_MIN)
      CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
   CASE (NSCARC_SCOPE_PRECON)
      CALL SCARC_VECTOR_COPY(MG%X, MG%F, 1.0_EB, NLEVEL_MIN)
END SELECT

CALL SCARC_DEBUG_LEVEL (MG%D, 'SCARC_METHOD_MULTIGRID', 'D FINAL', NLEVEL_MIN)
CALL SCARC_DEBUG_LEVEL (MG%X, 'SCARC_METHOD_MULTIGRID', 'X FINAL', NLEVEL_MIN)
CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MULTIGRID', 'F FINAL', NLEVEL_MIN)

CALL SCARC_RESET_PARENT(PARENT)

TUSED_SCARC(NSCARC_TIME_MULTIGRID,MYID+1)=TUSED_SCARC(NSCARC_TIME_MULTIGRID,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL    ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL    ,MYID+1)+SECOND()-TNOW
END SUBROUTINE SCARC_METHOD_MULTIGRID


!> ------------------------------------------------------------------------------------------------
!> Control multigrid cycling (F/V/W)
!> Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CYCLE_CONTROL(NSCOPE, NL)
INTEGER, INTENT(IN) :: NSCOPE, NL
INTEGER :: NM, NL0, ICYCLE

SELECT CASE (NSCOPE)

 !! ---------------------------------------------------------------------------------------------
 !! initialize cycle counts at beginning of multigrid method
 !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_SETUP)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
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

 !! ---------------------------------------------------------------------------------------------
 !! reset cycle counts at beginning of each new multigrid iteration
 !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_RESET)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         DO NL0 = NLEVEL_MIN, NLEVEL_MAX
            SCARC(NM)%CYCLE_COUNT(1,NL0)=SCARC(NM)%CYCLE_COUNT(2,NL0)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLE_NEXT

 !! ---------------------------------------------------------------------------------------------
 !! determine where to proceed with cycling
 !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_PROCEED)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE

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


!> ------------------------------------------------------------------------------------------------
!> Perform smoothing
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SMOOTHING(MG, NTYPE, NL, NSCOPE)
INTEGER , INTENT(IN) :: NTYPE, NL, NSCOPE
TYPE (SCARC_SCOPE_TYPE), INTENT(INOUT) :: MG
INTEGER :: ITE, ISTATE
REAL(EB):: TNOW
LOGICAL :: BMATVEC, BL2NORM
TYPE (SCARC_SCOPE_TYPE) :: SM

TNOW = SECOND()

!> ------------------------------------------------------------------------------------------------
!> Initialization
!> ------------------------------------------------------------------------------------------------
SM%X = MG%X
SM%F = MG%F
SM%D = MG%D

SELECT CASE (NTYPE)
   CASE (NSCARC_CYCLE_PRESMOOTH)
      SM%CROUTINE = 'SCARC_PRESMOOTHER'
   CASE (NSCARC_CYCLE_POSTSMOOTH)
      SM%CROUTINE = 'SCARC_POSTSMOOTHER'
END SELECT
SM%NIT   = SCARC_SMOOTH_ITERATIONS
SM%EPS   = SCARC_SMOOTH_ACCURACY
SM%OMEGA = SCARC_SMOOTH_OMEGA

BL2NORM  = .TRUE.
IF (NTYPE == NSCARC_CYCLE_PRESMOOTH.AND.NL==1) THEN
   BMATVEC = .FALSE.
ELSE
   BMATVEC = .TRUE.
ENDIF
BMATVEC = .TRUE.
BL2NORM = .TRUE.

CALL SCARC_DEBUG_LEVEL (SM%X, 'SCARC_SMOOTHER', 'X init ', NL)
CALL SCARC_DEBUG_LEVEL (SM%F, 'SCARC_SMOOTHER', 'F init ', NL)

!> ------------------------------------------------------------------------------------------------
!> Calculate initial defect on level NL (only if BMATVEC = .TRUE.)
!> Because initial vector is set to zero, this defect corresponds to F
!> ------------------------------------------------------------------------------------------------
IF (BMATVEC) THEN
   CALL SCARC_MATVEC_PRODUCT (SM%X, SM%D, NL)                                  !>  D := A*X
   CALL SCARC_VECTOR_SUM     (SM%F, SM%D, 1.0_EB, -1.0_EB, NL)                 !>  D := F - D
ENDIF

CALL SCARC_DEBUG_LEVEL (SM%D, 'SCARC_SMOOTHER', 'D init ', NL)

IF (BL2NORM.AND.BMATVEC) THEN
   SM%RESIN = SCARC_L2NORM (SM%D, NL)                                          !>  RESIN := ||D||
ELSE
   SM%RESIN = SCARC_RESIDUAL
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Smoothing loop
!> ------------------------------------------------------------------------------------------------
SMOOTH_LOOP: DO ITE=1, SM%NIT

   CALL SCARC_PRECONDITIONING (SM%D, SM%D, NL, NSCOPE)                         !>  D := PRECON (D)
CALL SCARC_DEBUG_LEVEL (SM%D, 'SCARC_SMOOTHER', 'D precon ', NL)
   CALL SCARC_VECTOR_SUM      (SM%D, SM%X, SM%OMEGA, 1.0_EB, NL)               !>  X := OMEGA*D + X
   CALL SCARC_MATVEC_PRODUCT  (SM%X, SM%D, NL)                                 !>  D := A*X
CALL SCARC_DEBUG_LEVEL (SM%D, 'SCARC_SMOOTHER', 'D matvec ', NL)
   CALL SCARC_VECTOR_SUM      (SM%F, SM%D, 1.0_EB, -1.0_EB, NL)                !>  D := F - D
CALL SCARC_DEBUG_LEVEL (SM%D, 'SCARC_SMOOTHER', 'D sum ', NL)

   IF (BL2NORM.OR.ITE==SM%NIT) THEN
      SM%RES = SCARC_L2NORM (SM%D, NL)                                         !>  RES := ||D||
!WRITE(LU_SCARC,*) 'SMOOTHER RESIDUUM =',SM%RES
      ISTATE = SCARC_CONVERGENCE_STATE(SM, ITE, NL)
      IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP                     !>  RES < TOL ?
   ENDIF
!WRITE(LU_SCARC,*) 'SMOOTHER LOOP ', ITE,': RESIDUAL =',SM%RES, BL2NORM

ENDDO SMOOTH_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> ------------------------------------------------------------------------------------------------
IF (BL2NORM) CALL SCARC_CONVERGENCE_RATE(SM, ITE, NL)

TUSED_SCARC(NSCARC_TIME_SMOOTH,MYID+1)=TUSED_SCARC(NSCARC_TIME_SMOOTH,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
END SUBROUTINE SCARC_SMOOTHING


!> ------------------------------------------------------------------------------------------------
!> Gaussian elimination for coarse grid solution
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_DIRECT(NVECTORX, NVECTORF)
INTEGER, INTENT(IN) :: NVECTORX, NVECTORF
INTEGER :: NM, IC, IOFFSET, IERR
REAL(EB), POINTER, DIMENSION(:)     :: VCX, VCF
REAL (EB) :: TNOW
TYPE (SCARC_TYPE), POINTER :: SM
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SCARC_ROUTINE = 'SCARC_METHOD_DIRECT'

TNOW = SECOND()

IERR=0
SM => SCARC(NMASTER)

!WRITE(LU_SCARC, *) 'STARTING DIRECT SOLVER'

!> ------------------------------------------------------------------------------------------------
!> Parallel version
!> ------------------------------------------------------------------------------------------------
IF (N_MPI_PROCESSES > 1) THEN

!WRITE(LU_SCARC, *) 'OLD VERSION'
   DO NM = 1, NMESHES
      SM%COUNTS1(NM-1) = NC_COARSE(NM)
      SM%DISPLS1(NM-1) = SM%OFFSET(NM)
   ENDDO

   DO NM = 1, NMESHES
      IF (PROCESS(NM) /= MYID) CYCLE

      SL => SCARC(NM)%LEVEL(NLEVEL_MAX)
      IOFFSET = SM%OFFSET(NM)

      DO IC = 1, SL%NZ
         SM%X_COARSE (IC) = VCF(IC)
      ENDDO

      WRITE(SCARC_MESSAGE,'(2A)') TRIM(SCARC_ROUTINE),': Method not implemented yet'
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN

   ENDDO

   WRITE(*,*) 'STILL MKL VERSION !!'
   !IF (MYID+1 == NMASTER) THEN
 !!   CALL DGETRS('N', NC_COARSE0, 1, SM%A_COARSE, NC_COARSE0, SM%PIVOT, &
 !!               SM%X_COARSE, NC_COARSE0, IERR)
   !ENDIF

   DO NM = 1, NMESHES
      IF (PROCESS(NM) /= MYID) CYCLE

      SL => SCARC(NM)%LEVEL(NLEVEL_MAX)
      IOFFSET = SM%OFFSET(NM)

      WRITE(*,*) "ACHTUNG; WAS IST HIER MIT BVECTOR???", NVECTORX, NVECTORF
      !VBX => POINT_TO_HVECTOR (NVECTORX, NM)

      DO IC = 1, SL%NC
         VCX(IC) = SM%X_COARSE (IC)
      ENDDO

   ENDDO

!> ------------------------------------------------------------------------------------------------
!> Serial version
!> ------------------------------------------------------------------------------------------------
ELSE

   DO NM = 1, NMESHES
      IF (PROCESS(NM) /= MYID) CYCLE

      SL => SCARC(NM)%LEVEL(NLEVEL_MAX)
      IOFFSET = SM%OFFSET(NM)

      DO IC = 1, SL%NZ
         SM%X_COARSE (IC) = VCF(IC)
      ENDDO

   ENDDO

   WRITE(*,*) 'STILL MKL VERSION !!'
   !IF (MYID+1 == NMASTER) THEN
 !!   CALL DGETRS('N', NC_COARSE0, 1, SM%A_COARSE, NC_COARSE0, SM%PIVOT, &
 !!               SM%X_COARSE, NC_COARSE0, IERR)
   !ENDIF

   DO NM = 1, NMESHES
      IF (PROCESS(NM) /= MYID) CYCLE

      SL => SCARC(NM)%LEVEL(NLEVEL_MAX)
      IOFFSET = SM%OFFSET(NM)

      WRITE(*,*) "ACHTUNG; WAS IST HIER MIT BVECTOR???"
      !VBX => POINT_TO_HVECTOR (NVECTORX, NM)

      DO IC = 1, SL%NC
         VCX(IC) = SM%X_COARSE (IC)
      ENDDO

   ENDDO

ENDIF

TUSED_SCARC(NSCARC_TIME_COARSE,MYID+1)=TUSED_SCARC(NSCARC_TIME_COARSE,MYID+1)+SECOND()-TNOW
TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)=TUSED_SCARC(NSCARC_TIME_TOTAL ,MYID+1)+SECOND()-TNOW
END SUBROUTINE SCARC_METHOD_DIRECT



!> ------------------------------------------------------------------------------------------------
!> Setup environement in every solver CALL (i.e. set pointers to used vectors)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SCOPE(SCOPE, NTYPE, NSCOPE, NRHS, NL)
INTEGER, INTENT(IN)  :: NSCOPE, NTYPE, NRHS
INTEGER, INTENT(OUT) :: NL
TYPE (SCARC_SCOPE_TYPE), INTENT(OUT):: SCOPE

IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
WRITE(LU_SCARC,*) 'SETUP_SCOPE:A: TYPE_METHOD=',TYPE_METHOD
WRITE(LU_SCARC,*) 'SETUP_SCOPE:A: NSCOPE=',NSCOPE
WRITE(LU_SCARC,*) 'SETUP_SCOPE:A: NTYPE=',NTYPE
WRITE(LU_SCARC,*) 'SETUP_SCOPE:A: NRHS=',NRHS
ENDIF

TYPE_SCOPE  = NSCOPE

SELECT CASE (NSCOPE)
   CASE (NSCARC_SCOPE_COARSE)
      SELECT CASE (TYPE_COARSE)
         CASE( NSCARC_COARSE_ITERATIVE)
            TYPE_METHOD = NSCARC_METHOD_KRYLOV
         CASE( NSCARC_COARSE_DIRECT)
            TYPE_METHOD = NSCARC_METHOD_MKL
      END SELECT
   CASE (NSCARC_SCOPE_PRECON)
      TYPE_METHOD = NSCARC_METHOD_MULTIGRID
   CASE DEFAULT
      TYPE_METHOD = TYPE_METHOD
END SELECT

IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
WRITE(LU_SCARC,*) 'SETUP_SCOPE:B: TYPE_METHOD=',TYPE_METHOD
WRITE(LU_SCARC,*) 'SETUP_SCOPE:B: NSCOPE=',NSCOPE
WRITE(LU_SCARC,*) 'SETUP_SCOPE:B: NTYPE=',NTYPE
WRITE(LU_SCARC,*) 'SETUP_SCOPE:B: NRHS=',NRHS
ENDIF

SCOPE%RESIN = 1.0_EB
SCOPE%RES   = 0.0_EB


SELECT CASE (TYPE_METHOD)

 !! ---------------------------------------------------------------------------------------------------
 !! Krylov method
 !! ---------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      TYPE_PRECON = NTYPE

      SELECT CASE (TYPE_KRYLOV)

         !> CG-method
         CASE (NSCARC_KRYLOV_CG)

            SCOPE%F = NRHS                                             !> set correct right hand side vector

            SELECT CASE (NSCOPE)
               CASE (NSCARC_SCOPE_MAIN)

                  NL = NLEVEL_MIN

                  SCOPE%CROUTINE = 'SCARC_GLOBAL_CG'
                  SCOPE%EPS      =  SCARC_KRYLOV_ACCURACY
                  SCOPE%NIT      =  SCARC_KRYLOV_ITERATIONS

                  SCOPE%X = NSCARC_VECTOR_X
                  SCOPE%D = NSCARC_VECTOR_D
                  SCOPE%G = NSCARC_VECTOR_G
                  SCOPE%Y = NSCARC_VECTOR_Y
                  SCOPE%W = NSCARC_VECTOR_W

                  !TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
               CASE (NSCARC_SCOPE_COARSE)

                  NL = NLEVEL_MAX

                  SCOPE%CROUTINE = 'SCARC_COARSE_CG'
                  SCOPE%EPS      =  SCARC_COARSE_ACCURACY
                  SCOPE%NIT      =  SCARC_COARSE_ITERATIONS

                  IF (TYPE_METHOD0 == NSCARC_METHOD_MULTIGRID) THEN              !>  pure MG-method
                     SCOPE%X = NSCARC_VECTOR_X                                   !>  take first-stage MG-vectors
                     SCOPE%D = NSCARC_VECTOR_D
                     SCOPE%G = NSCARC_VECTOR_G
                     SCOPE%Y = NSCARC_VECTOR_Y
                     SCOPE%W = NSCARC_VECTOR_W
                  ELSE                                                           !> CG-MG-method
                     SCOPE%X = NSCARC_VECTOR_X2                                  !> take second-stage MG-vectors
                     SCOPE%D = NSCARC_VECTOR_D2
                     SCOPE%G = NSCARC_VECTOR_G2
                     SCOPE%Y = NSCARC_VECTOR_Y2
                     SCOPE%W = NSCARC_VECTOR_W2
                  ENDIF

                  !TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
            END SELECT

         !> BICG-method
         CASE (NSCARC_KRYLOV_BICG)

            SCOPE%F = NRHS                                             !> set correct right hand side vector

            SELECT CASE (NSCOPE)
               CASE (NSCARC_SCOPE_MAIN)

                  NL = NLEVEL_MIN

                  SCOPE%CROUTINE = 'SCARC_GLOBAL_BICG'
                  SCOPE%EPS      =  SCARC_KRYLOV_ACCURACY
                  SCOPE%NIT      =  SCARC_KRYLOV_ITERATIONS

                  SCOPE%X = NSCARC_VECTOR_X
                  SCOPE%D = NSCARC_VECTOR_D
                  SCOPE%G = NSCARC_VECTOR_G
                  SCOPE%Y = NSCARC_VECTOR_Y
                  SCOPE%W = NSCARC_VECTOR_W
                  SCOPE%Z = NSCARC_VECTOR_Z
                  !TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE

               CASE (NSCARC_SCOPE_COARSE)

                  NL = NLEVEL_MAX

                  SCOPE%CROUTINE = 'SCARC_COARSE_BICG'
                  SCOPE%EPS      =  SCARC_COARSE_ACCURACY
                  SCOPE%NIT      =  SCARC_COARSE_ITERATIONS

                  IF (TYPE_METHOD0 == NSCARC_METHOD_MULTIGRID) THEN              !>  pure MG-method
                     SCOPE%X = NSCARC_VECTOR_X                                   !>  take first-stage MG-vectors
                     SCOPE%D = NSCARC_VECTOR_D
                     SCOPE%G = NSCARC_VECTOR_G
                     SCOPE%Y = NSCARC_VECTOR_Y
                     SCOPE%W = NSCARC_VECTOR_W
                     SCOPE%Z = NSCARC_VECTOR_Z
                  ELSE                                                           !> CG-MG-method
                     SCOPE%X = NSCARC_VECTOR_X2                                  !> take second-stage MG-vectors
                     SCOPE%D = NSCARC_VECTOR_D2
                     SCOPE%G = NSCARC_VECTOR_G2
                     SCOPE%Y = NSCARC_VECTOR_Y2
                     SCOPE%W = NSCARC_VECTOR_W2
                     SCOPE%Z = NSCARC_VECTOR_Z2
                  ENDIF

                  !TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE

            END SELECT

      END SELECT

   !! ---------------------------------------------------------------------------------------------------
   !! Multigrid method
   !! ---------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_MULTIGRID)

      TYPE_SMOOTH = NTYPE
      TYPE_PRECON = NTYPE

      SCOPE%EPS   = SCARC_MULTIGRID_ACCURACY
      SCOPE%NIT   = SCARC_MULTIGRID_ITERATIONS
      SCOPE%OMEGA = SCARC_SMOOTH_OMEGA

      NL = NLEVEL_MIN
      SCOPE%F = NRHS                                                   !> set correct right hand side vector

      !> select scope (multigrid as main solver or preconditioner)
      SELECT CASE (NSCOPE)
         CASE (NSCARC_SCOPE_MAIN)

            SCOPE%CROUTINE = 'SCARC_GLOBAL_MULTIGRID'

            SCOPE%X = NSCARC_VECTOR_X
            SCOPE%D = NSCARC_VECTOR_D

            !TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE

         CASE (NSCARC_SCOPE_PRECON)

            SCOPE%CROUTINE = 'SCARC_PRECON_MULTIGRID'

            SCOPE%X = NSCARC_VECTOR_X2
            SCOPE%D = NSCARC_VECTOR_D2

            !TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE

      END SELECT

   !! ---------------------------------------------------------------------------------------------------
   !! Pardiso method
   !! ---------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_MKL)

      SCOPE%CROUTINE = 'SCARC_COARSE_MKL'

      NL = NLEVEL_MAX

      SCOPE%F = NRHS                                                !> set correct right hand side vector
      SCOPE%X = NSCARC_VECTOR_X                                     !> set correct solution vector

END SELECT

END SUBROUTINE SCARC_SETUP_SCOPE


!> ----------------------------------------------------------------------------------------------------
!> Set initial solution corresponding to boundary data in BXS, BXF, ...
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WORKSPACE(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, I, J, K, IC
TYPE (MESH_TYPE), POINTER ::  M
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SCARC_ROUTINE = 'SCARC_SETUP_WORKSPACE'

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN

 !! Initialize solution and right hand side vector corresponding to boundary conditions
   SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)

      !> -------------- 2D --------------
      CASE (NSCARC_DIMENSION_TWO)

         DO NM = 1, NMESHES
            IF (PROCESS(NM) /= MYID) CYCLE

            M  => MESHES(NM)
            SL => SCARC(NM)%LEVEL(NL)

!WRITE(LU_SCARC,*) 'PRHS:'
!WRITE(LU_SCARC,'(4E14.6)') ((M%PRHS(I,1,K),I=1,4),K=4,1,-1)
!WRITE(LU_SCARC,*) 'H:'
!WRITE(LU_SCARC,'(4E14.6)') ((M%H(I,1,K),I=1,4),K=4,1,-1)
!WRITE(LU_SCARC,*) 'HS:'
!WRITE(LU_SCARC,'(4E14.6)') ((M%HS(I,1,K),I=1,4),K=4,1,-1)
!WRITE(LU_SCARC,*) 'CCVAR(CGSC):'
!WRITE(LU_SCARC,'(4i6)') ((SL%CCVAR(I,1,K,CGSC),I=1,4),K=4,1,-1)
!WRITE(LU_SCARC,*) 'CCVAR(UNKH):'
!WRITE(LU_SCARC,'(4i6)') ((SL%CCVAR(I,1,K,UNKH),I=1,4),K=4,1,-1)

            !> get right hand side and initial vector
            DO K = 1, M%KBAR
               DO I = 1, M%IBAR

<<<<<<< HEAD
                  IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.SL%CCVAR(I,1,K,CGSC) /= IS_GASPHASE) CYCLE
=======
                  IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC.AND.SL%CCVAR(I,1,K,CGSC) /= IS_GASPHASE) CYCLE
>>>>>>> origin/scarc
                  IC = SL%CCVAR(I,1,K,UNKH)

                  SL%F(IC) = M%PRHS (I, 1, K)
                  IF (PREDICTOR) THEN
                     SL%X(IC) = M%H (I, 1, K)
                  ELSE
                     SL%X(IC) = M%HS(I, 1, K)
                  ENDIF

               ENDDO
            ENDDO

<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'Starting SETUP_WORKSPACE:'
   WRITE(LU_SCARC,*) 'PRHS:'
   WRITE(LU_SCARC,'(4F18.6)') M%PRHS
   WRITE(LU_SCARC,*) 'BEFORE F:'
   WRITE(LU_SCARC,'(4F18.6)') SL%F
ENDIF
            !> set boundary conditions at exterior boundaries (corresponding to pois.f90)
<<<<<<< HEAD
            LEVEL_WALL_CELLS_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
=======
            LEVEL_WALLCELLS_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
>>>>>>> origin/scarc

               I    = SL%WALL(IW)%IXW
               J    = SL%WALL(IW)%IYW
               K    = SL%WALL(IW)%IZW

               IF (J /= 1) THEN
                  WRITE(SCARC_MESSAGE,'(2A,I8)') TRIM(SCARC_ROUTINE),': Wrong index for J=',J
                  CALL SHUTDOWN(SCARC_MESSAGE)
               ENDIF

               IOR0 = SL%WALL(IW)%IOR

<<<<<<< HEAD
               IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. SL%CCVAR(I,1,K,CGSC) /= IS_GASPHASE) CYCLE
               IC = SL%CCVAR(I,1,K,UNKH)

=======
               IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC .AND. SL%CCVAR(I,1,K,CGSC) /= IS_GASPHASE) CYCLE
               IC = SL%CCVAR(I,1,K,UNKH)

               !IC   = SL%WALL(IW)%ICW

               IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) &
                  WRITE(LU_SCARC,*) 'BDRY: IW=',IW,': I,J,K,IC,IOR0=',I,J,K,IC,IOR0

>>>>>>> origin/scarc
               SELECT CASE (IOR0)
                  CASE (1)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DXI2 * M%BXS(1,K)         !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,I4)') 'DIR: 1:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) + SL%DXI * M%BXS(1,K)                   !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,I4)') 'NEU: 1:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) + SL%DXI * M%BXS(1,K)                   !> Neumann
>>>>>>> origin/scarc
                     ENDIF
                  CASE (-1)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DXI2 *M%BXF(1,K)          !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'DIR:-1:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) - SL%DXI *M%BXF(1,K)                    !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'NEU:-1:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DXI *M%BXF(1,K)                    !> Neumann
>>>>>>> origin/scarc
                     ENDIF
                  CASE (3)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DZI2 * M%BZS(I,1)         !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'DIR: 3:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) + SL%DZI * M%BZS(I,1)                   !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'NEU: 3:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) + SL%DZI * M%BZS(I,1)                   !> Neumann
>>>>>>> origin/scarc
                     ENDIF
                  CASE (-3)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DZI2 * M%BZF(I,1)         !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'DIR: -3:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) - SL%DZI  * M%BZF(I,1)                  !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'NEU: -3:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DZI  * M%BZF(I,1)                  !> Neumann
>>>>>>> origin/scarc
                     ENDIF
               END SELECT

            ENDDO LEVEL_WALL_CELLS_LOOP2D



         ENDDO
      !> -------------- 3D --------------
      CASE (NSCARC_DIMENSION_THREE)

         DO NM = 1, NMESHES
            IF (PROCESS(NM) /= MYID) CYCLE

            M  => MESHES(NM)
            SL => SCARC(NM)%LEVEL(NL)

            !> get right hand side and initial vector
            DO K = 1, M%KBAR
               DO J = 1, M%JBAR
                  DO I = 1, M%IBAR
<<<<<<< HEAD
                     IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.SL%CCVAR(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
=======
                     IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC.AND.SL%CCVAR(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
>>>>>>> origin/scarc
                     IC = SL%CCVAR(I,J,K,UNKH)
                     SL%F(IC) = M%PRHS (I, J, K)
                     IF (PREDICTOR) THEN
                        SL%X(IC) = M%H (I, J, K)
                     ELSE
                        SL%X(IC) = M%HS(I, J, K)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

<<<<<<< HEAD
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
>>>>>>> origin/scarc
   WRITE(LU_SCARC,*) 'Starting SETUP_WORKSPACE:'
   WRITE(LU_SCARC,*) 'PRHS:'
   WRITE(LU_SCARC,'(4F18.6)') M%PRHS
   WRITE(LU_SCARC,*) 'BEFORE F:'
   WRITE(LU_SCARC,'(4F18.6)') SL%F
ENDIF
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
   WRITE(77,*) 'WALL(.)%BTYPE:'
   WRITE(77,'(4I4)') (SL%WALL(IW)%BTYPE, IW=1, M%N_EXTERNAL_WALL_CELLS)
ENDIF
            LEVEL_WALL_CELLS_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
=======
            LEVEL_WALLCELLS_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
>>>>>>> origin/scarc

               I    = SL%WALL(IW)%IXW
               J    = SL%WALL(IW)%IYW
               K    = SL%WALL(IW)%IZW

               IOR0 = SL%WALL(IW)%IOR

<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,*) 'Processing ',I,J,K,IOR0

               IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. SL%CCVAR(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
=======
               IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC .AND. SL%CCVAR(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
>>>>>>> origin/scarc
               IC = SL%CCVAR(I,J,K,UNKH)

               SELECT CASE (IOR0)
                  CASE (1)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DXI2 * M%BXS(J,K)         !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,I4)') 'DIR: 1:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) + SL%DXI * M%BXS(J,K)                   !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,I4)') 'NEU: 1:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,3f18.10)') 'A:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) + SL%DXI * M%BXS(J,K)                   !> Neumann
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,3f18.10)') 'B:F(',I,',',J,',',K,')=',SL%F(IC)
>>>>>>> origin/scarc
                     ENDIF
                  CASE (-1)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DXI2 *M%BXF(J,K)          !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'DIR:-1:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) - SL%DXI *M%BXF(J,K)                    !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'NEU:-1:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'C:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DXI *M%BXF(J,K)                    !> Neumann
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'D:F(',I,',',J,',',K,')=',SL%F(IC)
>>>>>>> origin/scarc
                     ENDIF
                  CASE (2)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DYI2 * M%BYS(I,K)         !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'DIR: 2:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) + SL%DYI * M%BYS(I,K)                   !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'NEU: 2:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'E:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) + SL%DYI * M%BYS(I,K)                   !> Neumann
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'F:F(',I,',',J,',',K,')=',SL%F(IC)
>>>>>>> origin/scarc
                     ENDIF
                  CASE (-2)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DYI2 *M%BYF(I,K)          !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'DIR:-2:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) - SL%DYI *M%BYF(I,K)                    !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'NEU:-2:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'G:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DYI *M%BYF(I,K)                    !> Neumann
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'H:F(',I,',',J,',',K,')=',SL%F(IC)
>>>>>>> origin/scarc
                     ENDIF
                  CASE (3)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DZI2 * M%BZS(I,J)         !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'DIR: 3:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) + SL%DZI * M%BZS(I,J)                   !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'NEU: 3:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'I:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) + SL%DZI * M%BZS(I,J)                   !> Neumann
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'J:F(',I,',',J,',',K,')=',SL%F(IC)
>>>>>>> origin/scarc
                     ENDIF
                  CASE (-3)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DZI2 * M%BZF(I,J)         !> Dirichlet
<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'DIR: -3:F(',I,',',J,',',K,')=',SL%F(IC),IC
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN.AND.SL%WALL(IW)%NOM==0) THEN
                        SL%F(IC) = SL%F(IC) - SL%DZI  * M%BZF(I,J)                  !> Neumann
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,E18.10,i4)') 'NEU: -3:F(',I,',',J,',',K,')=',SL%F(IC),IC
=======
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'K:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DZI  * M%BZF(I,J)                  !> Neumann
!WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,f18.10)') 'L:F(',I,',',J,',',K,')=',SL%F(IC)
>>>>>>> origin/scarc
                     ENDIF
               END SELECT

            ENDDO LEVEL_WALL_CELLS_LOOP3D

         ENDDO

   END SELECT SELECT_DIMENSION

<<<<<<< HEAD
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
   WRITE(LU_SCARC,*) 'AFTER F:'
   WRITE(LU_SCARC,'(4F18.6)') SL%F
ENDIF
IF (TYPE_DEBUG >= NSCARC_DEBUG_EXTREME) THEN
   WRITE(77,*) 'AFTER F:'
   WRITE(77,'(4F18.6)') SL%F
ENDIF
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'AFTER F:'
   WRITE(LU_SCARC,'(4F18.6)') SL%F
ENDIF
>>>>>>> origin/scarc

ELSE IF (TYPE_SCOPE == NSCARC_SCOPE_COARSE) THEN
   DO NM = 1, NMESHES
      IF (PROCESS(NM) /= MYID) CYCLE
      SCARC(NM)%LEVEL(NLEVEL_MAX)%X = 0.0_EB
   ENDDO
ENDIF

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !! In case of a Krylov method initialize auxiliary arrays
   CASE (NSCARC_METHOD_KRYLOV)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL => SCARC(NM)%LEVEL(NL)
         DO IC = SL%NCS+1, SL%NCE
            SL%D(IC) = 0.0_EB
            SL%G(IC) = 0.0_EB
            SL%Y(IC) = 0.0_EB
            SL%W(IC) = 0.0_EB
         ENDDO
      ENDDO

      IF (TYPE_KRYLOV == NSCARC_KRYLOV_BICG) THEN
         DO NM = 1, NMESHES
            IF (PROCESS(NM) /= MYID) CYCLE
            SL => SCARC(NM)%LEVEL(NL)
            DO IC = SL%NCS+1, SL%NCE
               SL%Z(IC) = 0.0_EB
            ENDDO
         ENDDO
      ENDIF

   !! In case of a multigrid method with coarse grid solution by CG, clear CG-vectors on max level
   CASE (NSCARC_METHOD_MULTIGRID)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL => SCARC(NM)%LEVEL(NL)
         DO IC = SL%NCS+1, SL%NCE
            SL%D = 0.0_EB
         ENDDO
         SL => SCARC(NM)%LEVEL(NLEVEL_MAX)
         DO IC = SL%NCS+1, SL%NCE
            SL%G = 0.0_EB
            SL%W = 0.0_EB
         ENDDO

      ENDDO

END SELECT SELECT_METHOD

END SUBROUTINE SCARC_SETUP_WORKSPACE


!> ------------------------------------------------------------------------------------------------
!> Print out residual information for loop ITE
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_INFO(RES, ITE, NL, CROUTINE)
INTEGER, INTENT(IN) :: ITE, NL
REAL(EB), INTENT(IN) :: RES
CHARACTER(*), INTENT(IN) :: CROUTINE
INTEGER:: NM

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE
   IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,1000) TRIM(CROUTINE), NL, ITE,  RES
   IF (TYPE_DEBUG>NSCARC_DEBUG_INFO1.AND.MYID==0) write(*,1000) TRIM(CROUTINE), NL, ITE,  RES
ENDDO

1000 FORMAT (5X,A30,': level=',i4,': #ite= ',i4,': res =',e25.16)
END SUBROUTINE SCARC_CONVERGENCE_INFO


!> ------------------------------------------------------------------------------------------------
!> Check if solver converges or diverges and print out residual information
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CONVERGENCE_STATE(SCOPE, ITE, NL)
INTEGER, INTENT(IN) :: ITE, NL
INTEGER :: ISTATE
TYPE (SCARC_SCOPE_TYPE), INTENT(IN) :: SCOPE

ISTATE = NSCARC_STATE_PROCEED
SCARC_RESIDUAL = SCOPE%RES

IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,1000) TRIM(SCOPE%CROUTINE), NL, ITE, SCOPE%RES
IF (TYPE_DEBUG>=NSCARC_DEBUG_INFO2.AND.MYID==0) WRITE(*,1000) TRIM(SCOPE%CROUTINE), NL, ITE, SCOPE%RES

SELECT CASE (TYPE_ACCURACY)
   CASE (NSCARC_ACCURACY_RELATIVE)
      IF (SCOPE%RES <= SCOPE%RESIN*SCOPE%EPS)  ISTATE = NSCARC_STATE_CONV
      IF (SCOPE%RES <= 1.0E-15)                ISTATE = NSCARC_STATE_CONV
   CASE (NSCARC_ACCURACY_ABSOLUTE)
      IF (SCOPE%RES <= SCOPE%EPS .AND. SCOPE%RES <= SCOPE%RESIN*SCARC_ACCURACY_RELATIVE) ISTATE = NSCARC_STATE_CONV
END SELECT
IF (SCOPE%RES > SCARC_ACCURACY_DIVERGENCE) ISTATE = NSCARC_STATE_DIVG

SCARC_CONVERGENCE_STATE = ISTATE
RETURN

1000 FORMAT (5X,A30,': level=',i4,': #ite= ',i4,': res =',e25.16)
END FUNCTION SCARC_CONVERGENCE_STATE


!> ------------------------------------------------------------------------------------------------
!> Compute convergence rate and print out resiual information for final loop
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_RATE(SCOPE, ITE, ISTATE)
INTEGER, INTENT(IN) :: ITE, ISTATE
INTEGER  :: ITERATIONS
REAL(EB) :: RESIDUAL, CAPPA
TYPE (SCARC_SCOPE_TYPE), INTENT(IN) :: SCOPE

RESIDUAL = SCOPE%RES
IF (ISTATE == NSCARC_STATE_DIVG) THEN
   ITERATIONS = - 1
   CAPPA      = 1.0_EB
ELSE
   IF (ISTATE == NSCARC_STATE_CONV) THEN
     ITERATIONS = ITE
   ELSE
     ITERATIONS = ITE-1
   ENDIF
   IF (SCOPE%RESIN >= 1.0E-70_EB) THEN
      CAPPA = (SCOPE%RES/SCOPE%RESIN) ** (1.0_EB/ITERATIONS)
   ELSE
      CAPPA = 0.0_EB
   ENDIF
ENDIF

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   SCARC_CAPPA      = CAPPA
   SCARC_RESIDUAL   = RESIDUAL
   SCARC_ITERATIONS = ITERATIONS
ENDIF

IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS.AND.TYPE_METHOD==TYPE_METHOD0.AND.TYPE_SCOPE==NSCARC_SCOPE_MAIN) THEN
   WRITE(LU_SCARC,2000) SCOPE%CROUTINE, ITERATIONS, CAPPA
   IF (MYID==0) WRITE(*,2000) SCOPE%CROUTINE, ITERATIONS, CAPPA
ELSE IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,2000) SCOPE%CROUTINE, ITERATIONS, CAPPA
   IF (MYID==0) WRITE(*,2000) SCOPE%CROUTINE, ITERATIONS, CAPPA
ENDIF

2000 FORMAT (5X,A30,': iterations: ',i6,':  convergence rate =',e14.6)
END SUBROUTINE SCARC_CONVERGENCE_RATE


!> ------------------------------------------------------------------------------------------------
!> Perform restriction from finer to coarser grid in multigrid method
!>    - 'FI' corresponds to finer   grid
!>    - 'CO' corresponds to coarser grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTRICTION (NVECTOR_FI, NVECTOR_CO, NL_FI, NL_CO)
INTEGER, INTENT(IN) :: NVECTOR_FI, NVECTOR_CO, NL_FI, NL_CO
INTEGER :: NM, ICOL, IC
REAL(EB), POINTER, DIMENSION(:)     :: FC_CO, DC_FI, R
INTEGER , POINTER, DIMENSION(:)     :: R_ROW, R_COL
INTEGER , POINTER :: NX_CO, NY_CO, NZ_CO, NC_CO
INTEGER  :: NX_FI, NY_FI, NZ_FI
INTEGER  :: IX_FI, IY_FI, IZ_FI, IC_FI(8)
INTEGER  :: IX_CO, IY_CO, IZ_CO, IC_CO
REAL(EB) :: AUX, AUX0
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLC, SLF

TYPE_VECTOR = NVECTOR_FI
IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC.AND.TYPE_COARSENING < NSCARC_COARSENING_GMG) &
   CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_VECTOR, NL_FI)

SELECT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

 !! ------------------------- Geometric multigrid -------------------------------------------
   CASE (NSCARC_MULTIGRID_GEOMETRIC)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE

         SLC => SCARC(NM)%LEVEL(NL_CO)
         SLF => SCARC(NM)%LEVEL(NL_FI)

         NX_CO => SLC%NX
         NY_CO => SLC%NY
         NZ_CO => SLC%NZ

         NX_FI = 2*NX_CO
         NY_FI = 2*NY_CO
         NZ_FI = 2*NZ_CO

         CALL POINT_TO_VECTOR(NVECTOR_FI, NM, NL_FI, DC_FI)
         CALL POINT_TO_VECTOR(NVECTOR_CO, NM, NL_CO, FC_CO)

<<<<<<< HEAD
         IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
=======
         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
>>>>>>> origin/scarc
            WRITE(LU_SCARC,*) 'NX_CO=',NX_CO
            WRITE(LU_SCARC,*) 'NY_CO=',NY_CO
            WRITE(LU_SCARC,*) 'NZ_CO=',NZ_CO
            WRITE(LU_SCARC,*) 'NX_FI=',NX_FI
            WRITE(LU_SCARC,*) 'NY_FI=',NY_FI
            WRITE(LU_SCARC,*) 'NZ_FI=',NZ_FI
            WRITE(LU_SCARC,*) 'SIZE(DC_FI)=',SIZE(DC_FI)
            WRITE(LU_SCARC,*) 'SIZE(FC_CO)=',SIZE(FC_CO)
         ENDIF

         SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)

            !> ----------------- 2D ---------------
            CASE (NSCARC_DIMENSION_TWO)

               DO IZ_CO = 1, NZ_CO
                  DO IX_CO = 1, NX_CO

<<<<<<< HEAD
                     IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. SLC%CCVAR(IX_CO,1,IZ_CO,CGSC)/=IS_GASPHASE) CYCLE
=======
                     IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC .AND. SLC%CCVAR(IX_CO,1,IZ_CO,CGSC)/=IS_GASPHASE) CYCLE
>>>>>>> origin/scarc

                     IX_FI = 2*IX_CO
                     IZ_FI = 2*IZ_CO

                     IC_CO = SLC%CCVAR(IX_CO, 1, IZ_CO, UNKH) 
<<<<<<< HEAD

         IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,*) 'PROCESSING: REST: ', IX_CO, IZ_CO, IX_FI, IZ_FI, IC_CO

=======

         IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) &
            WRITE(LU_SCARC,*) 'PROCESSING: REST: ', IX_CO, IZ_CO, IX_FI, IZ_FI, IC_CO

>>>>>>> origin/scarc
                     IC_FI(1) = SLF%CCVAR(IX_FI-1, 1, IZ_FI-1, UNKH) 
                     IC_FI(2) = SLF%CCVAR(IX_FI-1, 1, IZ_FI  , UNKH) 
                     IC_FI(3) = SLF%CCVAR(IX_FI  , 1, IZ_FI-1, UNKH) 
                     IC_FI(4) = SLF%CCVAR(IX_FI  , 1, IZ_FI  , UNKH) 

                     FC_CO(IC_CO) = 0.25_EB * (  DC_FI(IC_FI(1)) &
                                               + DC_FI(IC_FI(2)) &
                                               + DC_FI(IC_FI(3)) &
                                               + DC_FI(IC_FI(4)) )

                  ENDDO
               ENDDO

            !> ----------------- 3D ---------------
            CASE (NSCARC_DIMENSION_THREE)

               DO IZ_CO = 1, NZ_CO
                  DO IY_CO = 1, NY_CO
                     DO IX_CO = 1, NX_CO

<<<<<<< HEAD
                        IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. SLC%CCVAR(IX_CO,IY_CO,IZ_CO,CGSC)/=IS_GASPHASE) CYCLE
=======
                        IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC .AND. SLC%CCVAR(IX_CO,IY_CO,IZ_CO,CGSC)/=IS_GASPHASE) CYCLE
>>>>>>> origin/scarc

                        IX_FI = 2*IX_CO
                        IY_FI = 2*IY_CO
                        IZ_FI = 2*IZ_CO

                        IC_CO = SLC%CCVAR(IX_CO, IY_CO, IZ_CO, UNKH) 
                        IC_FI(1) = SLF%CCVAR(IX_FI-1, IY_FI-1, IZ_FI-1, UNKH) 
                        IC_FI(2) = SLF%CCVAR(IX_FI-1, IY_FI-1, IZ_FI  , UNKH) 
                        IC_FI(3) = SLF%CCVAR(IX_FI-1, IY_FI  , IZ_FI-1, UNKH) 
                        IC_FI(4) = SLF%CCVAR(IX_FI-1, IY_FI  , IZ_FI  , UNKH) 
                        IC_FI(5) = SLF%CCVAR(IX_FI  , IY_FI-1, IZ_FI-1, UNKH) 
                        IC_FI(6) = SLF%CCVAR(IX_FI  , IY_FI-1, IZ_FI  , UNKH) 
                        IC_FI(7) = SLF%CCVAR(IX_FI  , IY_FI  , IZ_FI-1, UNKH) 
                        IC_FI(8) = SLF%CCVAR(IX_FI  , IY_FI  , IZ_FI  , UNKH) 

                        !WRITE(LU_SCARC,'(A,I4,A,8I4)') 'REST:',IC_CO, ':', IC_FI

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

         END SELECT SELECT_DIMENSION

      ENDDO

 !! ------------------------- Algebraic multigrid -------------------------------------------
   CASE (NSCARC_MULTIGRID_ALGEBRAIC)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE

         SLC => SCARC(NM)%LEVEL(NL_CO)
         SLF => SCARC(NM)%LEVEL(NL_FI)

         CALL POINT_TO_VECTOR(NVECTOR_FI, NM, NL_FI, DC_FI)
         CALL POINT_TO_VECTOR(NVECTOR_CO, NM, NL_CO, FC_CO)

         NC_CO => SLC%NCS

         R     => SLF%R
         R_ROW => SLF%R_ROW
         R_COL => SLF%R_COL

         DO IC_CO = 1, NC_CO
            AUX = 0.0_EB
            DO ICOL = R_ROW(IC_CO), R_ROW(IC_CO+1)-1
               IC = R_COL(ICOL)
               AUX0 = AUX
               AUX = AUX + DC_FI(IC) * R(ICOL)
            ENDDO
            FC_CO(IC_CO) = AUX
            !WRITE(LU_SCARC,*) 'RESTRICT2: FC_CO(',IC_CO,')=',FC_CO(IC_CO)
         ENDDO

      ENDDO

END SELECT SELECT_MULTIGRID


END SUBROUTINE SCARC_RESTRICTION


!> ------------------------------------------------------------------------------------------------
!> Perform prolongation from coarser to finer grid in multigrid method
!>    - 'CO' corresponds to coarser grid
!>    - 'FI' corresponds to finer   grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PROLONGATION (NVECTOR_CO, NVECTOR_FI, NL_CO, NL_FI)
INTEGER, INTENT(IN) :: NVECTOR_CO, NVECTOR_FI, NL_CO, NL_FI
INTEGER :: NM, ICOL, IC, I
REAL(EB), POINTER, DIMENSION(:)     :: XC_CO, DC_FI, P
INTEGER , POINTER, DIMENSION(:)     :: P_ROW, P_COL
INTEGER , POINTER, DIMENSION(:,:)   :: P_PTR
INTEGER , POINTER :: NX_CO, NY_CO, NZ_CO, NC_FI, NCE_FI
INTEGER  :: NX_FI, NY_FI, NZ_FI
INTEGER  :: IX_FI, IY_FI, IZ_FI, IC_FI(8)
INTEGER  :: IX_CO, IY_CO, IZ_CO, IC_CO
REAL(EB) :: AUX, SCAL
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLC, SLF


SELECT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

 !! ------------------------- Geometric multigrid -------------------------------------------
   CASE (NSCARC_MULTIGRID_GEOMETRIC)

      DO NM = 1, NMESHES
      IF (PROCESS(NM) /= MYID) CYCLE

         SLC => SCARC(NM)%LEVEL(NL_CO)
         SLF => SCARC(NM)%LEVEL(NL_FI)

         NX_CO => SLC%NX
         NY_CO => SLC%NY
         NZ_CO => SLC%NZ

         NX_FI = 2*NX_CO
         NY_FI = 2*NY_CO
         NZ_FI = 2*NZ_CO

         CALL POINT_TO_VECTOR(NVECTOR_CO, NM, NL_CO, XC_CO)
         CALL POINT_TO_VECTOR(NVECTOR_FI, NM, NL_FI, DC_FI)

         SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)

            !> ----------------- 2D ---------------
            CASE (NSCARC_DIMENSION_TWO)

               DO IZ_CO = 1, NZ_CO
                  DO IX_CO = 1, NX_CO

<<<<<<< HEAD
                     IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.SLC%CCVAR(IX_CO,1,IZ_CO,CGSC)/=IS_GASPHASE) CYCLE
=======
                     IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC.AND.SLC%CCVAR(IX_CO,1,IZ_CO,CGSC)/=IS_GASPHASE) CYCLE
>>>>>>> origin/scarc

                     IX_FI = 2*IX_CO
                     IY_FI = 1
                     IZ_FI = 2*IZ_CO

                     IC_CO = SLC%CCVAR(IX_CO, 1, IZ_CO, UNKH) 

                     IC_FI(1) = SLF%CCVAR(IX_FI-1, 1, IZ_FI-1, UNKH) 
                     IC_FI(2) = SLF%CCVAR(IX_FI-1, 1, IZ_FI  , UNKH) 
                     IC_FI(3) = SLF%CCVAR(IX_FI  , 1, IZ_FI-1, UNKH) 
                     IC_FI(4) = SLF%CCVAR(IX_FI  , 1, IZ_FI  , UNKH) 

                     DO I = 1, 4
                        DC_FI(IC_FI(I)) = XC_CO(IC_CO)
                        !WRITE(LU_SCARC,*) 'PROL1:: DC_FI(',IC_FI(I),')=',DC_FI(IC_FI(I))
                     ENDDO


                  ENDDO
               ENDDO


            !> ----------------- 3D ---------------
            CASE (NSCARC_DIMENSION_THREE)

               DO IZ_CO = 1, NZ_CO
                  DO IY_CO = 1, NY_CO
                     DO IX_CO = 1, NX_CO

<<<<<<< HEAD
                        IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.SLC%CCVAR(IX_CO,IY_CO,IZ_CO,CGSC)/=IS_GASPHASE) CYCLE
=======
                        IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC.AND.SLC%CCVAR(IX_CO,IY_CO,IZ_CO,CGSC)/=IS_GASPHASE) CYCLE
>>>>>>> origin/scarc

                        IX_FI = 2*IX_CO
                        IY_FI = 2*IY_CO
                        IZ_FI = 2*IZ_CO

                        IC_CO = SLC%CCVAR(IX_CO, IY_CO, IZ_CO, UNKH) 

                        IC_FI(1) = SLF%CCVAR(IX_FI-1, IY_FI-1, IZ_FI-1, UNKH) 
                        IC_FI(2) = SLF%CCVAR(IX_FI-1, IY_FI-1, IZ_FI  , UNKH) 
                        IC_FI(3) = SLF%CCVAR(IX_FI-1, IY_FI  , IZ_FI-1, UNKH) 
                        IC_FI(4) = SLF%CCVAR(IX_FI-1, IY_FI  , IZ_FI  , UNKH) 
                        IC_FI(5) = SLF%CCVAR(IX_FI  , IY_FI-1, IZ_FI-1, UNKH) 
                        IC_FI(6) = SLF%CCVAR(IX_FI  , IY_FI-1, IZ_FI  , UNKH) 
                        IC_FI(7) = SLF%CCVAR(IX_FI  , IY_FI  , IZ_FI-1, UNKH) 
                        IC_FI(8) = SLF%CCVAR(IX_FI  , IY_FI  , IZ_FI  , UNKH) 
                        !WRITE(LU_SCARC,'(A,I4,A,8I4)') 'PROL:',IC_CO, ':', IC_FI

                        DO I = 1, 8
                           DC_FI(IC_FI(I)) = XC_CO(IC_CO)
                        ENDDO

                     ENDDO
                  ENDDO
               ENDDO

         END SELECT SELECT_DIMENSION

      ENDDO

 !! ------------------------- Algebraic multigrid -------------------------------------------
   CASE (NSCARC_MULTIGRID_ALGEBRAIC)

SCAL = 1.0_EB
IF (TYPE_COARSENING >= NSCARC_COARSENING_GMG) SCAL=2.0_EB

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE

         SLC => SCARC(NM)%LEVEL(NL_CO)
         SLF => SCARC(NM)%LEVEL(NL_FI)

         CALL POINT_TO_VECTOR(NVECTOR_CO, NM, NL_CO, XC_CO)
         CALL POINT_TO_VECTOR(NVECTOR_FI, NM, NL_FI, DC_FI)

         NC_FI  => SLF%NCS
         NCE_FI => SLF%NCE

         P      => SLF%P
         P_ROW  => SLF%P_ROW
         P_COL  => SLF%P_COL
         P_PTR  => SLF%P_PTR

         DO IC = 1, NC_FI
            AUX = 0.0_EB
            DO ICOL = P_ROW(IC), P_ROW(IC+1)-1
               IC_CO = P_COL(ICOL)
               AUX = XC_CO(IC_CO) * P(ICOL)
            ENDDO
            DC_FI(IC) = SCAL*AUX
            !WRITE(LU_SCARC,*) 'PROL2:: DC_FI(',IC,')=',DC_FI(IC)
         ENDDO
         DO IC = NC_FI+1, NCE_FI
            DC_FI(IC) = 0.0_EB
         ENDDO
      ENDDO

END SELECT SELECT_MULTIGRID


END SUBROUTINE SCARC_PROLONGATION



!> ------------------------------------------------------------------------------------------------
!> Finalize data
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_FINISH_SOLVER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, I, J, K, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE

   M  => MESHES(NM)
   SL => SCARC(NM)%LEVEL(NL)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

   DO K = 1, M%KBAR
      DO J = 1, M%JBAR
         DO I = 1, M%IBAR
<<<<<<< HEAD
            IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. SL%CCVAR(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
            IC = SL%CCVAR(I,J,K,UNKH)
            HP(I, J, K) = SL%X(IC)
IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) '1: HP(',I,',',J,',',K,')=',HP(I,J,K)
=======
            IF (.NOT.PRES_ON_WHOLE_DOMAIN_SCRC .AND. SL%CCVAR(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
            IC = SL%CCVAR(I,J,K,UNKH)
            HP(I, J, K) = SL%X(IC)
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '1: HP(',I,',',J,',',K,')=',HP(I,J,K)
>>>>>>> origin/scarc
         ENDDO
      ENDDO
   ENDDO

ENDDO

IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(LU_SCARC,'(6F18.10)') (((HP(I,J,K),I=0,5),J=1,1),K=5,0,-1)
ENDIF

END SUBROUTINE SCARC_FINISH_SOLVER


!> ------------------------------------------------------------------------------------------------
!> Set correct boundary values at external and internal boundaries
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_GHOSTCELLS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IXG, IYG, IZG, IXW, IYW, IZW
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL


DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE
   !! point to correct pressure vector on mesh 'NM'
   M  => MESHES(NM)
   SL => SCARC(NM)%LEVEL(NL)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

 !! compute ghost cell values
<<<<<<< HEAD
   WALL_CELLS_LOOP: DO IW = 1, SL%N_EXTERNAL_WALL_CELLS
=======
   WALLCELLS_LOOP: DO IW = 1, SL%N_EXTERNAL_WALL_CELLS
>>>>>>> origin/scarc

      IXG = SL%WALL(IW)%IXG
      IYG = SL%WALL(IW)%IYG
      IZG = SL%WALL(IW)%IZG

      IXW = SL%WALL(IW)%IXW
      IYW = SL%WALL(IW)%IYW
      IZW = SL%WALL(IW)%IZW

      IOR0 = SL%WALL(IW)%IOR

      SELECT CASE (IOR0)
         CASE ( 1)
            IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXG,IYG,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXS(IYW,IZW)
            ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXG,IYG,IZG) =  HP(IXW,IYW,IZW) - DXI *M%BXS(IYW,IZW)
            ENDIF
         CASE (-1)
            IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXG,IYG,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXF(IYW,IZW)
            ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXG,IYG,IZG) =  HP(IXW,IYW,IZW) + DXI *M%BXF(IYW,IZW)
            ENDIF
         CASE ( 2)
            IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXG,IYG,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYS(IXW,IZW)
            ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXG,IYG,IZG) =  HP(IXW,IYW,IZW) - DETA *M%BYS(IXW,IZW)
            ENDIF
         CASE (-2)
            IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXG,IYG,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYF(IXW,IZW)
            ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXG,IYG,IZG) =  HP(IXW,IYW,IZW) + DETA *M%BYF(IXW,IZW)
            ENDIF
         CASE ( 3)
            IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXG,IYG,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZS(IXW,IYW)
            ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXG,IYG,IZG) =  HP(IXW,IYW,IZW) - DZETA *M%BZS(IXW,IYW)
            ENDIF
         CASE (-3)
         IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXG,IYG,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZF(IXW,IYW)
            ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXG,IYG,IZG) =  HP(IXW,IYW,IZW) + DZETA *M%BZF(IXW,IYW)
            ENDIF
      END SELECT
   ENDDO WALL_CELLS_LOOP

ENDDO


!> -----------------------------------------------------------------------------------------------
!> Perform data exchange to achieve consistency of ghost values along internal boundaries
!> -----------------------------------------------------------------------------------------------
IF (N_MPI_PROCESSES > 1) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_PRESSURE, NL)

END SUBROUTINE SCARC_UPDATE_GHOSTCELLS


!> ------------------------------------------------------------------------------------------------
!>  Perform data exchange corresponding to requested exchange type (CALL receive and send-routines)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE (NTYPE, NL)
INTEGER, INTENT(IN):: NTYPE, NL

NREQ_SCARC = 0
TYPE_EXCHANGE = NTYPE

CALL SCARC_EXCHANGE_RECEIVE(NL)
CALL SCARC_EXCHANGE_SEND(NL)

END SUBROUTINE SCARC_EXCHANGE


!> ------------------------------------------------------------------------------------------------
!>  Receive data from neighbors (corresponds to POST_RECEIVES)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_RECEIVE (NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, IERR
TYPE (SCARC_TYPE)         , POINTER ::  S
TYPE (OSCARC_TYPE)        , POINTER ::  OS
TYPE (SCARC_LEVEL_TYPE), POINTER ::  OSL

IERR=0
RECEIVE_MESH_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE
   RECEIVE_OMESH_LOOP: DO NOM = 1, NMESHES

      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_OMESH_LOOP

      S   => SCARC(NM)
      OS  => SCARC(NM)%OSCARC(NOM)

      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE RECEIVE_OMESH_LOOP

      SELECT_EXCHANGE_TYPE: SELECT CASE (TYPE_EXCHANGE)


         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring step size along internal boundary
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER0(1),1,MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring wall data
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALLINFO)

            OSL => OS%LEVEL(NL)
            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),SIZE(OS%RECV_INTEGER),MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring wall data
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CCVAR)

            OSL => OS%LEVEL(NL)
            NREQ_SCARC = NREQ_SCARC+1
            OS%RECV_INTEGER = 0
            CALL MPI_IRECV(OS%RECV_INTEGER(1),SIZE(OS%RECV_INTEGER),MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring grid dimensions
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MESHINFO)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),8,MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring step size along internal boundary
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WIDTHINFO)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_REAL0(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange number of neighboring cells for AMG method (compact type only)
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_TRANSFER_SIZE)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),9,MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange number of neighboring cells for AMG method (compact type only)
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),1,MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange neighboring CELL_TYPEs
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_TYPE)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),SIZE(OS%RECV_INTEGER),MPI_INTEGER,&
                           SNODE,TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)


         !> ---------------------------------------------------------------------------------------
         !> Perform exchanges for
         !>    - internal values for matrix-vector multiplication
         !>    - internal boundariy values
         !>    - internal subdiagonal matrix values
         !>    - internal subdiagonal or ghost matrix values
         !>    - internal measure/CELL_TYPE values
         !> ---------------------------------------------------------------------------------------
         CASE DEFAULT

            NREQ_SCARC = NREQ_SCARC+1
            !WRITE(LU_SCARC,*) 'RECEIVING IN LENGTH ', SIZE(OS%RECV_REAL)
            CALL MPI_IRECV(OS%RECV_REAL(1),SIZE(OS%RECV_REAL),MPI_DOUBLE_PRECISION,&
                           SNODE,TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

      END SELECT SELECT_EXCHANGE_TYPE
   ENDDO RECEIVE_OMESH_LOOP
ENDDO RECEIVE_MESH_LOOP

END SUBROUTINE SCARC_EXCHANGE_RECEIVE


!> ------------------------------------------------------------------------------------------------
!> Send data to neighbors (corresponds to MESH_EXCHANGE)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_SEND (NL)
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, NOM, NCOL, NROW
INTEGER  :: IERR, IW, IWL, IWG, IWW, IW0, IROW, ICOL, IPTR, ICPL, ICELL_TYPE
<<<<<<< HEAD
INTEGER  :: IOR0, IC, JC, ICC, JCC, ICE, ICO, ICN, ICG, ICW
INTEGER  :: IX, IY, IZ
INTEGER  :: IXW, IYW, IZW
INTEGER  :: IXG, IYG, IZG
=======
INTEGER  :: IOR0, IC, JC, ICC, JCC, ICE, ICO, ICN, ICG, ICW, IX, IY, IZ
>>>>>>> origin/scarc
INTEGER  :: I, J, K, LL
REAL(EB) :: ZSUM, ZSUM1
INTEGER , POINTER, DIMENSION(:)     ::  RECV_INTEGER, RECV_INTEGER0
REAL(EB), POINTER, DIMENSION(:)     ::  RECV_REAL, RECV_REAL0
REAL(EB), POINTER, DIMENSION(:)     ::  VECTOR
REAL(EB), POINTER, DIMENSION(:,:,:) :: HVECTOR
TYPE (SCARC_TYPE)      , POINTER :: S
TYPE (OSCARC_TYPE)     , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, OSL

IERR = 0

!WRITE(LU_SCARC,*) 'STARTING EXCHANGE_SEND, TYPE_EXCHANGE=',TYPE_EXCHANGE

!> ------------------------------------------------------------------------------------------------
!> Collect data for sending corresponding to requested exchange type
!> ------------------------------------------------------------------------------------------------
MESH_PACK_LOOP: DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   IF (PROCESS(NM)/=MYID)  CYCLE MESH_PACK_LOOP

   S   => SCARC(NM)
   SL  => SCARC(NM)%LEVEL(NL)

   OMESH_PACK_LOOP: DO NOM = 1, NMESHES

      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      OS  => SCARC(NM)%OSCARC(NOM)

      IF (OS%NICMAX_S == 0 .AND. OS%NICMAX_R == 0) CYCLE OMESH_PACK_LOOP

      OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

      OMESH_PACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_BASIC: pack neighboring sizes for exchange
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC)

            !WRITE(LU_SCARC,'(a,i4,a,i4,a,3i4)') 'BASIC: NM=',NM,': NOM=',NOM,': SENDING ',OS%LEVEL(NL)%NWL
            OS%SEND_INTEGER0(1)=OS%LEVEL(NL)%NWL

            IF (RNODE /= SNODE) THEN
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER0(1),1,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_WIDTHINFO: pack neighboring grid resolution information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WIDTHINFO)

            SELECT CASE(OSL%IOR)
               CASE (1)
                  OSL%DH = SL%DXL(0)
               CASE (-1)
                  OSL%DH = SL%DXL(SL%NX)
                CASE (2)
                  OSL%DH = SL%DYL(0)
               CASE (-2)
                  OSL%DH = SL%DYL(SL%NY)
               CASE (3)
                  OSL%DH = SL%DZL(0)
               CASE (-3)
                  OSL%DH = SL%DZL(SL%NZ)
            END SELECT
            OS%SEND_REAL0(1) = OSL%DH

            IF (RNODE /= SNODE) THEN
               !WRITE(LU_SCARC,'(a,i4,a,i4,a,f12.6)') 'RESOLUTION: NM=',NM,': NOM=',NOM,': SENDING ',OS%SEND_REAL0(1)
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL0(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_WALLINFO: pack neighboring wallinfo information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALLINFO)

            !WRITE(LU_SCARC,*) 'NWL=',OSL%NWL
            IPTR=1
            DO IWL = 1, OSL%NWL
               IWG = OSL%IWL_TO_IWG(IWL)
               !WRITE(LU_SCARC,'(5(a,i4))') &
               !   'NL=',NL,': IWL=',IWL,': IWG=',IWG,': NCPLS=',OSL%NCPLS,': SIZE(SL%WALL)=', SIZE(SL%WALL)
               OS%SEND_INTEGER(IPTR   ) = SL%WALL(IWG)%IXG
               OS%SEND_INTEGER(IPTR+ 1) = SL%WALL(IWG)%IYG
               OS%SEND_INTEGER(IPTR+ 2) = SL%WALL(IWG)%IZG
               OS%SEND_INTEGER(IPTR+ 3) = SL%WALL(IWG)%IXW
               OS%SEND_INTEGER(IPTR+ 4) = SL%WALL(IWG)%IYW
               OS%SEND_INTEGER(IPTR+ 5) = SL%WALL(IWG)%IZW
               OS%SEND_INTEGER(IPTR+ 6) = SL%WALL(IWG)%IXN(1)
               OS%SEND_INTEGER(IPTR+ 7) = SL%WALL(IWG)%IXN(2)
               OS%SEND_INTEGER(IPTR+ 8) = SL%WALL(IWG)%IYN(1)
               OS%SEND_INTEGER(IPTR+ 9) = SL%WALL(IWG)%IYN(2)
               OS%SEND_INTEGER(IPTR+10) = SL%WALL(IWG)%IZN(1)
               OS%SEND_INTEGER(IPTR+11) = SL%WALL(IWG)%IZN(2)
               OS%SEND_INTEGER(IPTR+12) = SL%WALL(IWG)%NOM
               !WRITE(LU_SCARC,'(i4,i4)') (I, OS%SEND_INTEGER(I), I=IPTR,IPTR+11)
               IPTR = IPTR + 13
               DO ICPL=1,OSL%NCPLS
                  OS%SEND_INTEGER(IPTR)=SL%WALL(IWG)%ICE(ICPL)
                  !WRITE(LU_SCARC,'(i4,i4)') (I, OS%SEND_INTEGER(I), I=IPTR,IPTR)
                  IPTR = IPTR + 1
               ENDDO
            ENDDO

            IF (RNODE /= SNODE) THEN
               NREQ_SCARC = NREQ_SCARC+1
               !WRITE(LU_SCARC,*) 'SENDING IN SIZE ', SIZE(OS%SEND_INTEGER)
               CALL MPI_ISEND(OS%SEND_INTEGER(1),SIZE(OS%SEND_INTEGER),MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_CCVAR: pack neighboring CCVAR information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CCVAR)

            IPTR=1
            OS%SEND_INTEGER=0
            DO IWL = 1, OSL%NWL
               IWG = OSL%IWL_TO_IWG(IWL)
               OS%SEND_INTEGER(IPTR  ) = SL%CCVAR(SL%WALL(IWG)%IXW,SL%WALL(IWG)%IYW,SL%WALL(IWG)%IZW,CGSC)
               OS%SEND_INTEGER(IPTR+1) = SL%CCVAR(SL%WALL(IWG)%IXW,SL%WALL(IWG)%IYW,SL%WALL(IWG)%IZW,UNKH)
               !WRITE(LU_SCARC,*) 'SENDING CCVAR ', IWG, OS%SEND_INTEGER(iptr), OS%SEND_INTEGER(IPTR+1)
               IPTR = IPTR + 2
            ENDDO

            !WRITE(LU_SCARC,*) '================= >SENDING CCVAR '
            !WRITE(LU_SCARC,'(2i6)')  OS%SEND_INTEGER(1:2*OSL%NWL)

            IF (RNODE /= SNODE) THEN
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER(1),SIZE(OS%SEND_INTEGER),MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_CCVAR: pack neighboring CCVAR information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CCVAR)

            IPTR=1
            OS%SEND_INTEGER=0
            DO IWL = 1, OSL%NWL
               IWG = OSL%IWL_TO_IWG(IWL)
               OS%SEND_INTEGER(IPTR  ) = SL%CCVAR(SL%WALL(IWG)%IXW,SL%WALL(IWG)%IYW,SL%WALL(IWG)%IZW,CGSC)
               OS%SEND_INTEGER(IPTR+1) = SL%CCVAR(SL%WALL(IWG)%IXW,SL%WALL(IWG)%IYW,SL%WALL(IWG)%IZW,UNKH)
               !WRITE(LU_SCARC,*) 'SENDING CCVAR ', IWG, OS%SEND_INTEGER(iptr), OS%SEND_INTEGER(IPTR+1)
               IPTR = IPTR + 2
            ENDDO

            !WRITE(LU_SCARC,*) '================= >SENDING CCVAR '
            !WRITE(LU_SCARC,'(2i6)')  OS%SEND_INTEGER(1:2*OSL%NWL)

            IF (RNODE /= SNODE) THEN
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER(1),SIZE(OS%SEND_INTEGER),MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MESHINFO: pack neighboring mesh information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MESHINFO)

            SL => S%LEVEL(NL)

            OS%SEND_INTEGER(1)=SL%NX
            OS%SEND_INTEGER(2)=SL%NY
            OS%SEND_INTEGER(3)=SL%NZ
            OS%SEND_INTEGER(4)=SL%NC
            OS%SEND_INTEGER(5)=SL%NCS
            OS%SEND_INTEGER(6)=SL%NW
            OS%SEND_INTEGER(7)=SL%N_EXTERNAL_WALL_CELLS
            OS%SEND_INTEGER(8)=SL%N_INTERNAL_WALL_CELLS
            !WRITE(LU_SCARC,'(a,i4,a,i4,a,5i4)') 'MESHINFO: NM=',NM,': NOM=',NOM,': SEND ',OS%SEND_INTEGER(1:5)

            IF (RNODE /= SNODE) THEN
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER(1),8,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
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
            DO ICG=1, OSL%NCG
               IWG = OSL%ICG_TO_IWG(ICG)
               IX  = SL%WALL(IWG)%IXW
               IY  = SL%WALL(IWG)%IYW
               IZ  = SL%WALL(IWG)%IZW
               OS%SEND_REAL(LL) = HVECTOR(IX,IY,IZ)
               !WRITE(LU_SCARC,*) 'PACKING H(',IX,IY,IZ,')=',HVECTOR(IZ,IY,IZ)
               LL = LL+1
            ENDDO

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


         !! ---------------------------------------------------------------------------------------
         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_VECTOR: pack overlapping parts of a given vector
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_VECTOR)

            CALL POINT_TO_VECTOR(TYPE_VECTOR, NM, NL, VECTOR)

            LL = 1
            DO ICG= 1, OSL%NCG
               !ZSUM = 0.0_EB
               IWG = OSL%ICG_TO_IWG(ICG)
               IX  = SL%WALL(IWG)%IXW
               IY  = SL%WALL(IWG)%IYW
               IZ  = SL%WALL(IWG)%IZW

               !ZSUM = VECTOR(SL%CCVAR(IX, IY, IZ, UNKH))
               !OS%SEND_REAL(LL) = ZSUM/REAL(OSL%NCPLR,EB)

<<<<<<< HEAD
               IF (PRES_ON_WHOLE_DOMAIN.OR.SL%CCVAR(IX, IY, IZ, CGSC) == IS_GASPHASE) THEN
=======
               IF (PRES_ON_WHOLE_DOMAIN_SCRC.OR.SL%CCVAR(IX, IY, IZ, CGSC) == IS_GASPHASE) THEN
>>>>>>> origin/scarc
                  OS%SEND_REAL(LL) = VECTOR(SL%CCVAR(IX, IY, IZ, UNKH))
               ELSE
                  OS%SEND_REAL(LL) = -987654321.0_EB
               ENDIF
               !WRITE(LU_SCARC,'(a,i6,a,E20.10,a,i4,a,i4,a,4i4)') &
               !   'SEND_REAL(',LL,')=',OS%SEND_REAL(LL), &
               !  ' : IC=',SL%CCVAR(IX,IY,IZ,UNKH), ' : IWG=',IWG, ' : IX,IY,IZ,NOM',IX, IY, IZ

               LL = LL + 1
            ENDDO 

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               !WRITE(LU_SCARC, *) 'Sending in length ', OSL%NCG
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


         !! ---------------------------------------------------------------------------------------
         !! EXCHANGE_MEASURE_ADD: pack neighboring measure information (AMG only) with adding
         !! ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MEASURE_ADD)

            LL  = 1
            PACK_MEASURE_ADD: DO IWL=1, OSL%NWL
               ICO = OSL%IWL_TO_ICO(IWL)
               OS%SEND_REAL(LL) = SL%MEASURE(ICO)
               !WRITE(LU_SCARC,'(a,i4,a,f12.6)') 'PACKA VECTOR(',IWL,')=',SL%MEASURE(ICO)
               LL = LL+1
            ENDDO PACK_MEASURE_ADD

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> Send data along internal boundaries corresponding to requested exchange type
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_TYPE)

            LL  = 1
            PACK_CELL_TYPE: DO IWL=1, OSL%NWL
               ICW = OSL%IWL_TO_ICW(IWL)
               OS%SEND_INTEGER(LL) = SL%CELL_TYPE(ICW)
               !WRITE(LU_SCARC,*) 'EXCHANGE_CELL_TYPE: SENDING VECTOR(',ICW,')=',SL%CELL_TYPE(ICW)
               LL = LL+1
            ENDDO PACK_CELL_TYPE

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER(1), SIZE(OS%SEND_INTEGER), MPI_INTEGER, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_SIZE: Send sizes of overlapping matrix areas
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            OSL => OS%LEVEL(NL)
            OS%SEND_INTEGER(1) = OSL%NA0

            IF (RNODE /= SNODE) THEN
               !WRITE(LU_SCARC,*) 'SIZE_MATRIX: SENDING IN SIZE 1: ', OS%SEND_INTEGER(1)
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER(1),1,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_STENCIL: pack subdiagonal entries of system matrix on overlapping parts
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG)

            LL = 1
            IWW = 0

            WRITE(*,*) '!!!!! EXCHANGE_MATRIX_SUBDIAG, PACK: Achtung, hier nochmal checken  AMG !!!!'

            DO ICG=1,OSL%NCG
               IOR0 = ABS(OSL%WALL(ICG)%IOR)
               IWG = OSL%ICG_TO_IWG(ICG)
               IX  = SL%WALL(IWG)%IXG
               IY  = SL%WALL(IWG)%IYG
               IZ  = SL%WALL(IWG)%IZG
               DO ICPL = 1, OSL%NCPL
                  IC = SL%CCVAR(IX, IY, IZ, UNKH)
                  OS%SEND_REAL(LL) = SL%DI2(IOR0)
                  WRITE(LU_SCARC,*) 'PACK_MATRIX_SUBDIAG1: IC=',IC,': OS%SEND_REAL(',LL,')=', OS%SEND_REAL(LL)
                  LL = LL+1
               ENDDO
            ENDDO 


            !WRITE(LU_SCARC,*) 'NSCARC_EXCHANGE_MATRIX_SUBDIAG: SEND_REAL'
            !WRITE(LU_SCARC,'(8f12.6)') (OS%SEND_REAL(IC), IC=1, SIZE(OS%SEND_REAL))

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_STENCIL: pack matrix stencil information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_STENCIL)

             !> Pack first cell layer
            LL = 1
            OS%SEND_REAL = 0.0_EB

            DO ICG=1, OSL%NCG
  
            WRITE(*,*) '!!!!! EXCHANGE_MATRIX_STENCIL, PACK: Achtung, hier nochmal checken  AMG !!!!'
            !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4)') 'EXCHANGE_MATRIX_STENCIL: ICG=',ICG,&
            !    ': OSL%NCG=',OSL%NCG,': SL%NC=',SL%NC,': OSL%NCPLR=',OSL%NCPLR,&
            !    ': OSL%NCPLS=',OSL%NCPLS,': OSL%NCPL=',OSL%NCPL
 
               IWG = OSL%ICG_TO_IWG(ICG)
               IX  = SL%WALL(IWG)%IXG
               IY  = SL%WALL(IWG)%IYG
               IZ  = SL%WALL(IWG)%IZG
               DO ICPL = 1, OSL%NCPLR
                  IC = SL%CCVAR(IX, IY, IZ, UNKH)

                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') 'ICPL=',ICPL,': IC=',IC, ': ICG=',ICG

                  NCOL = 0
                  DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
                     JC = SL%A_COL(ICOL)
                     !WRITE(LU_SCARC,*) 'ICOL=',ICOL,': JC=',JC
                     IF (JC > SL%NC) THEN
                        IW0 = SL%ICE_TO_IWG(JC)
                        !WRITE(LU_SCARC,*) 'IW0=',IW0,': SL%NOM=',SL%WALL(IW0)%NOM,': NOM=',NOM,': NM=',NM
                        IF (SL%WALL(IW0)%NOM == NOM) NCOL = NCOL + 1
                     ELSE
                        NCOL = NCOL + 1
                     ENDIF
                     !WRITE(LU_SCARC,*) 'NCOL=',NCOL
                  ENDDO
                  OS%SEND_REAL(LL) = REAL(NCOL)
                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,f12.6)') &
                  !       'PACK_MATRIX_STENCIL: IC=',IC,': OS%SEND_REAL(',LL,')=', OS%SEND_REAL(LL)

                  LL = LL + 1
               ENDDO
            ENDDO 


            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_SYSTEM: pack overlapping parts of system matrix
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SYSTEM)

            LL = 1
            OS%SEND_REAL = 0.0_EB

            WRITE(*,*) '!!!!! EXCHANGE_MATRIX_SYSTEM, PACK: Achtung, hier nochmal checken  AMG !!!!'
            !> Pack first cell layer
            PACK_MATRIX_SYSTEM: DO ICG=1, OSL%NCG
               IF (OSL%WALL(ICG)%NOM/=NM) CYCLE PACK_MATRIX_SYSTEM
               IWG = OSL%ICG_TO_IWG(ICG)
               IX  = SL%WALL(IWG)%IXG
               IY  = SL%WALL(IWG)%IYG
               IZ  = SL%WALL(IWG)%IZG

            !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4)') 'EXCHANGE_MATRIX_STENCIL: ICG=',ICG,&
            !    ': OSL%NCG=',OSL%NCG,': SL%NC=',SL%NC,': OSL%NCPLR=',OSL%NCPLR,&
            !    ': OSL%NCPLS=',OSL%NCPLS,': OSL%NCPL=',OSL%NCPL

               DO ICPL = 1, OSL%NCPLR
                  IC = SL%CCVAR(IX, IY, IZ, UNKH)

                  MATRIX_COLUMN_LOOP: DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
                     JC = SL%A_COL(ICOL)
                     !WRITE(LU_SCARC,*) 'ICOL=',ICOL,': JC=',JC,': SL%NC=',SL%NC
                     IF (JC > SL%NC) THEN
                        IWL = SL%ICE_TO_IWL(JC)
                        !WRITE(LU_SCARC,*) 'A: IWL=',IWL,&
                        !  ': NOM=',OSL%WALL(IWL)%NOM,': NOM=',NOM,':NM=',NM
                        IF (OSL%WALL(IWL)%NOM /= NM) CYCLE MATRIX_COLUMN_LOOP
                        !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4)') &
                        !     'A: IWL=',IWL,':  ICGETO_IWG(',JC,')=',SL%ICE_TO_IWG(JC)
                        OS%SEND_REAL(LL) = - REAL(SL%ICE_TO_IWG(JC),EB)
                        OS%SEND_REAL(LL) = - REAL(JC)
                        !WRITE(LU_SCARC,'(a,i3,a,i3,a,f12.6)') &
                        !     'A: PACK_MATRIX_SYSTEM: IC=',IC,': OS%SEND_REAL(',LL,')=', OS%SEND_REAL(LL)
                     ELSE
                        OS%SEND_REAL(LL) =   REAL(JC,EB)
                        !WRITE(LU_SCARC,'(a,i3,a,i3,a,f12.6)') &
                        !     'B: PACK_MATRIX_SYSTEM: IC=',IC,': OS%SEND_REAL(',LL,')=', OS%SEND_REAL(LL)
                     ENDIF
                     OS%SEND_REAL(LL+1) = SL%A(ICOL)
                        !WRITE(LU_SCARC,'(a,i3,a,i3,a,f12.6)') &
                        !     'B: PACK_MATRIX_SYSTEM: IC=',IC,': OS%SEND_REAL(',LL+1,')=', OS%SEND_REAL(LL+1)
                     LL = LL + 2
                  ENDDO MATRIX_COLUMN_LOOP

               ENDDO

            ENDDO PACK_MATRIX_SYSTEM

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


         !! ---------------------------------------------------------------------------------------
         !! EXCHANGE_MATRIX_PROL: pack overlapping parts of prolongation matrix
         !! ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_PROL)

           ! WRITE(LU_SCARC,*) 'MATRIX_PROL:'

            LL = 1
            WRITE(*,*) '!!!!! EXCHANGE_MATRIX_PROL, PACK: Achtung, hier nochmal checken  AMG !!!!'

            PACK_MATRIX_PROL: DO ICG=1, OSL%NCG
               IF (OSL%WALL(ICG)%NOM/=NM) CYCLE PACK_MATRIX_PROL
               IWG = OSL%ICG_TO_IWG(ICG)
               IX  = SL%WALL(IWG)%IXG
               IY  = SL%WALL(IWG)%IYG
               IZ  = SL%WALL(IWG)%IZG

               !WRITE(LU_SCARC,'(a,2i4)') '========================== ICG=',ICG, OSL%NCPLR

               DO ICPL = 1, OSL%NCPLR
                  IC = SL%CCVAR(IX, IY, IZ, UNKH)

                  OS%SEND_REAL(LL)   = REAL(SL%CELL_TYPE(IC),EB)
                  OS%SEND_REAL(LL+1) = REAL(SL%P_ROW(IC+1)-SL%P_ROW(IC),EB)

                  !WRITE(LU_SCARC,'(a,i4,a,i4,a,2f12.6)') 'IC=',IC,': IC=',IC,': SEND ', OS%SEND_REAL(LL:LL+1)
                  LL = LL + 2

                  DO ICOL = SL%P_ROW(IC), SL%P_ROW(IC+1)-1
                     JC = SL%P_COL(ICOL)
                     IF (JC >= SL%NC) THEN
                        OS%SEND_REAL(LL) = - REAL(JC,EB)
                     ELSE
                        OS%SEND_REAL(LL) =   REAL(JC,EB)
                     ENDIF
                     OS%SEND_REAL(LL+1) = SL%P(ICOL)
                     !WRITE(LU_SCARC,'(a,i4,a,i4,a,2f12.6)') 'IC=',IC,': JC=',JC,': SEND ', OS%SEND_REAL(LL:LL+1)
                     LL = LL + 2
                  ENDDO

               ENDDO

            ENDDO PACK_MATRIX_PROL

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF

         !! ---------------------------------------------------------------------------------------
         !! EXCHANGE_MATRIX_REST: pack overlapping parts of restriction matrix
         !! ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_REST)

            WRITE(*,*) '!!!!! EXCHANGE_MATRIX_REST, PACK: Achtung, hier nochmal checken  AMG !!!!'
            LL = 1
            PACK_MATRIX_REST: DO ICG=1, OSL%NCG
               IF (OSL%WALL(IWG)%NOM/=NM) CYCLE PACK_MATRIX_REST
               
               WRITE(*,*) 'ACHTUNG: MATRIX_REST: Wegen CCVAR schauen!!'

               !WRITE(LU_SCARC,'(a,i4,a,i4,a,2f12.6)') '========================== ICG=',ICG,': NCPL=',NCPL
               ! ACHTUNG HIER WEGEN NCPL SCHAUEN!!
<<<<<<< HEAD

                  IX = SL%WALL(IWG)%IXG
                  IY = SL%WALL(IWG)%IYG
                  IZ = SL%WALL(IWG)%IZG
                  IC = SL%CCVAR(IX, IY, IZ, UNKH)

=======

                  IX = SL%WALL(IWG)%IXG
                  IY = SL%WALL(IWG)%IYG
                  IZ = SL%WALL(IWG)%IZG
                  IC = SL%CCVAR(IX, IY, IZ, UNKH)

>>>>>>> origin/scarc
                  ICC = SL%CELL_TYPE(IC)
                  IF (ICC > 0) THEN
                     OS%SEND_REAL(LL)   = REAL(ICC)
                     OS%SEND_REAL(LL+1) = REAL(SL%R_ROW(ICC+1)-SL%R_ROW(ICC),EB)

                     !WRITE(LU_SCARC,'(a,i4,a,i4,a,2f12.6)') 'IC=',IC,': ICC=',ICC,': SEND ', OS%SEND_REAL(LL:LL+1)
                     LL = LL + 2

                     DO ICOL = SL%R_ROW(ICC), SL%R_ROW(ICC+1)-1
                        JCC = SL%R_COL(ICOL)
                        IF (JCC >= SL%NCC) THEN
                           OS%SEND_REAL(LL) = - REAL(JCC,EB)
                        ELSE
                           OS%SEND_REAL(LL) =   REAL(JCC,EB)
                        ENDIF
                        OS%SEND_REAL(LL+1) = SL%R(ICOL)
                        !WRITE(LU_SCARC,'(a,i4,a,i4,a,2f12.6)') 'IC=',IC,': ICC=',ICC,': SEND ', OS%SEND_REAL(LL:LL+1)
                        LL = LL + 2
                     ENDDO

                  ENDIF

               !ENDDO REST_FIRST_LAYER_LOOP

            ENDDO PACK_MATRIX_REST

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_TRANSER_SIZE: pack sizes of transfer matrices (AMG only)
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_TRANSFER_SIZE)

            SL  =>  S%LEVEL(NL)
            OSL => OS%LEVEL(NL)
            OS%SEND_INTEGER(1) = SL%NC
            OS%SEND_INTEGER(2) = SL%NW
            OS%SEND_INTEGER(3) = SL%NCE
            OS%SEND_INTEGER(4) = SL%NW
            OS%SEND_INTEGER(5) = SL%NCCI
            OS%SEND_INTEGER(6) = OSL%NP0
            OS%SEND_INTEGER(7) = OSL%NR0
            OS%SEND_INTEGER(8) = OSL%NCC0
            OS%SEND_INTEGER(9) = OSL%NCF0
            !WRITE(LU_SCARC,*) 'SIZE_MATRIX: NM=',NM,': NOM=',NOM,': SENDING ',OS%SEND_INTEGER(1:9)

            IF (RNODE /= SNODE) THEN
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER(1),9,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF


      END SELECT OMESH_PACK_SELECT
   ENDDO OMESH_PACK_LOOP
ENDDO MESH_PACK_LOOP


!> ------------------------------------------------------------------------------------------------
!> Information from Mesh NM is received by Mesh NOM  (NOM receiver, NM sender)
!> ------------------------------------------------------------------------------------------------
IF (N_MPI_PROCESSES>1.AND.NREQ_SCARC/=0) CALL MPI_WAITALL(NREQ_SCARC,REQ_SCARC(1:NREQ_SCARC),MPI_STATUSES_IGNORE,IERR)

!> ------------------------------------------------------------------------------------------------
!> Extract communication data from corresponding RECEIVE-buffers
!> ------------------------------------------------------------------------------------------------
MESH_UNPACK_LOOP: DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE
   OMESH_UNPACK_LOOP: DO NOM=1,NMESHES

      SNODE  = PROCESS(NM)
      RNODE  = PROCESS(NOM)

      OS => SCARC(NM)%OSCARC(NOM)

      OMESH_UNPACK_IF: IF (OS%NICMAX_S/=0 .AND. OS%NICMAX_R/=0) THEN

         SL   => SCARC(NM)%LEVEL(NL)
         OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

         !>
         IF (RNODE/=SNODE) THEN
            RECV_INTEGER0  => SCARC(NM)%OSCARC(NOM)%RECV_INTEGER0
            RECV_INTEGER   => SCARC(NM)%OSCARC(NOM)%RECV_INTEGER
            RECV_REAL0     => SCARC(NM)%OSCARC(NOM)%RECV_REAL0
            RECV_REAL      => SCARC(NM)%OSCARC(NOM)%RECV_REAL
         ELSE
            RECV_INTEGER0  => SCARC(NOM)%OSCARC(NM)%SEND_INTEGER0
            RECV_INTEGER   => SCARC(NOM)%OSCARC(NM)%SEND_INTEGER
            RECV_REAL0     => SCARC(NOM)%OSCARC(NM)%SEND_REAL0
            RECV_REAL      => SCARC(NOM)%OSCARC(NM)%SEND_REAL
         ENDIF

         OMESH_UNPACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_BASIC: unpack basic neighboring sizes for exchange
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_BASIC)

               OSL%NCG  = RECV_INTEGER0(1)

               !WRITE(LU_SCARC,*) '======== NSCARC_EXCHANGE_BASIC:'
               !WRITE(LU_SCARC,*) 'NM=',NM,': NOM=',NOM,' : OSL%NCG =',OSL%NCG

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WIDTHINFO: unpack neighboring step size information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_WIDTHINFO)

               OSL%DH = RECV_REAL0(1)

               !WRITE(LU_SCARC,'(a,i4,a,2f12.6)') 'WIDTH1: FACE(', OSL%IOR,')%DH =',OSL%DH
               !WRITE(LU_SCARC,'(a,10f8.4)') 'WIDTH: DXL:',(SL%DXL(I), I=0,SL%NX)
               !WRITE(LU_SCARC,'(a,10f8.4)') 'WIDTH: DYL:',(SL%DYL(J), J=0,SL%NY)
               !WRITE(LU_SCARC,'(a,10f8.4)') 'WIDTH: DZL:',(SL%DZL(K), K=0,SL%NZ)

               SELECT CASE (OSL%IOR)
                  CASE ( 1)
                     SL%DXL(0)     = 0.5_EB*(OSL%DH + SL%DXL(0))
                  CASE (-1)
                     SL%DXL(SL%NX) = 0.5_EB*(OSL%DH + SL%DXL(SL%NX))
                  CASE ( 2)
                     SL%DYL(0)     = 0.5_EB*(OSL%DH + SL%DYL(0))
                  CASE (-2)
                     SL%DYL(SL%NY) = 0.5_EB*(OSL%DH + SL%DYL(SL%NY))
                  CASE ( 3)
                     SL%DZL(0)     = 0.5_EB*(OSL%DH + SL%DZL(0))
                  CASE (-3)
                     SL%DZL(SL%NZ) = 0.5_EB*(OSL%DH + SL%DZL(SL%NZ))
               END SELECT

               !WRITE(LU_SCARC,'(a,i4,a,2f12.6)') 'WIDTH2: FACE(', OSL%IOR,')%DH =',OSL%DH
               !WRITE(LU_SCARC,'(a,10f8.4)') 'WIDTH: DXL:',(SL%DXL(I), I=0,SL%NX)
               !WRITE(LU_SCARC,'(a,10f8.4)') 'WIDTH: DYL:',(SL%DYL(J), J=0,SL%NY)
               !WRITE(LU_SCARC,'(a,10f8.4)') 'WIDTH: DZL:',(SL%DZL(K), K=0,SL%NZ)


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WALLINFO: unpack neighboring wall information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_WALLINFO)
               IPTR=1
               DO ICG = 1, OSL%NCG
                  !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4)') 'RECEIVING ICG=',ICG,': NCPLR=',OSL%NCPLR
                  OSL%WALL(ICG)%IXG    = RECV_INTEGER(IPTR    )
                  OSL%WALL(ICG)%IYG    = RECV_INTEGER(IPTR + 1)
                  OSL%WALL(ICG)%IZG    = RECV_INTEGER(IPTR + 2)
                  OSL%WALL(ICG)%IXW    = RECV_INTEGER(IPTR + 3)
                  OSL%WALL(ICG)%IYW    = RECV_INTEGER(IPTR + 4)
                  OSL%WALL(ICG)%IZW    = RECV_INTEGER(IPTR + 5)
                  OSL%WALL(ICG)%IXN(1) = RECV_INTEGER(IPTR + 6)
                  OSL%WALL(ICG)%IXN(2) = RECV_INTEGER(IPTR + 7)
                  OSL%WALL(ICG)%IYN(1) = RECV_INTEGER(IPTR + 8)
                  OSL%WALL(ICG)%IYN(2) = RECV_INTEGER(IPTR + 9)
                  OSL%WALL(ICG)%IZN(1) = RECV_INTEGER(IPTR +10)
                  OSL%WALL(ICG)%IZN(2) = RECV_INTEGER(IPTR +11)
                  OSL%WALL(ICG)%NOM    = RECV_INTEGER(IPTR +12)
                  !WRITE(LU_SCARC,'(i4,i4)') (I, RECV_INTEGER(I), I=IPTR,IPTR+12)
                  IPTR = IPTR + 13
                  ALLOCATE (OSL%WALL(ICG)%ICE(OSL%NCPLR))
                  DO ICPL=1,OSL%NCPLR
                     OSL%WALL(ICG)%ICE(ICPL) = RECV_INTEGER(IPTR)
                     !WRITE(LU_SCARC,'(i4,i4)') (I, RECV_INTEGER(I), I=IPTR,IPTR)
                     IPTR = IPTR + 1
                  ENDDO
               ENDDO

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_CCVAR: unpack neighboring CCVAR information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_CCVAR)
               IPTR=1
               DO ICG = 1, OSL%NCG

                  OSL%WALL(ICG)%CGSC   = RECV_INTEGER(IPTR  )
                  OSL%WALL(ICG)%UNKH   = RECV_INTEGER(IPTR+1)
                  IPTR = IPTR + 2

                  IWG = OSL%ICG_TO_IWG(ICG)
<<<<<<< HEAD

                  IXG = SL%WALL(IWG)%IXG
                  IYG = SL%WALL(IWG)%IYG
                  IZG = SL%WALL(IWG)%IZG

                  IXW= SL%WALL(IWG)%IXW
                  IYW= SL%WALL(IWG)%IYW
                  IZW= SL%WALL(IWG)%IZW

                  !WRITE(LU_SCARC,'(A,i4,a,6i4)') 'OSL%WALL(',ICG,'):', &
                  !                   OSL%WALL(ICG)%CGSC, &
                  !                   OSL%WALL(ICG)%UNKH, IWG, IXG, IYG, IZG
                  
                  !IF (OSL%WALL(ICG)%CGSC == IS_GASPHASE.AND.SL%CCVAR(IXW,IYW,IZW,CGSC)/=IS_SOLID) THEN
                  IF (OSL%WALL(ICG)%CGSC == IS_GASPHASE) THEN
                     SL%CCVAR(IXG, IYG, IZG, CGSC) = OSL%WALL(ICG)%CGSC
                     SL%CCVAR(IXG, IYG, IZG, UNKH) = OSL%WALL(ICG)%UNKH
                  ENDIF

               ENDDO


=======
                  IX = SL%WALL(IWG)%IXG
                  IY = SL%WALL(IWG)%IYG
                  IZ = SL%WALL(IWG)%IZG

                  !WRITE(LU_SCARC,'(A,i4,a,6i4)') 'OSL%WALL(',ICG,'):', &
                  !                   OSL%WALL(ICG)%CGSC, &
                  !                   OSL%WALL(ICG)%UNKH, IWG, IX, IY, IZ
                  
                  SL%CCVAR(IX, IY, IZ, CGSC) = OSL%WALL(ICG)%CGSC
                  SL%CCVAR(IX, IY, IZ, UNKH) = OSL%WALL(ICG)%UNKH


               ENDDO


>>>>>>> origin/scarc
            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WALLINFO: unpack neighboring mesh information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MESHINFO)

               OSL%NX  = RECV_INTEGER(1)
               OSL%NY  = RECV_INTEGER(2)
               OSL%NZ  = RECV_INTEGER(3)
               OSL%NC  = RECV_INTEGER(4)
               OSL%NCS = RECV_INTEGER(5)
               OSL%NW  = RECV_INTEGER(6)
               OSL%N_EXTERNAL_WALL_CELLS  = RECV_INTEGER(7)
               OSL%N_INTERNAL_WALL_CELLS  = RECV_INTEGER(8)

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WALLINFO: unpack neighboring matrix size information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

               OSL%NA = RECV_INTEGER(1)

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WALLINFO: unpack neighboring transfer size information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_TRANSFER_SIZE)

               OSL%NC  =  RECV_INTEGER(1)
               OSL%NW  =  RECV_INTEGER(2)
               OSL%NCE =  RECV_INTEGER(3)
               OSL%NW  =  RECV_INTEGER(4)
               OSL%NCCI=  RECV_INTEGER(5)
               OSL%NP  =  RECV_INTEGER(6)
               OSL%NR  =  RECV_INTEGER(7)
               OSL%NCC =  RECV_INTEGER(8)
               OSL%NCF =  RECV_INTEGER(9)

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
               UNPACK_PRESSURE: DO IWL = 1, OSL%NWL

                  ZSUM=0.0_EB
                  IWG = OSL%IWL_TO_IWG(IWL)
                  DO ICPL = 1, OSL%NCPL
                     ZSUM=ZSUM+RECV_REAL(LL)
                     LL = LL+1
                  ENDDO

                  I=SL%WALL(IWG)%IXG
                  J=SL%WALL(IWG)%IYG
                  K=SL%WALL(IWG)%IZG

                  HVECTOR(I, J, K) = ZSUM/REAL(OSL%NCPL,EB)
                  !WRITE(LU_SCARC,'(a,3i4,a,f12.6)') 'UNPACKING H(', I, J, K,')=',HVECTOR(I, J, K)

               ENDDO UNPACK_PRESSURE

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_VECTOR: unpack overlapping parts of a given vector
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_VECTOR)

               CALL POINT_TO_VECTOR (TYPE_VECTOR, NM, NL, VECTOR)

               !WRITE(LU_SCARC,'(a,16E14.6)') 'RECV_REAL:',(RECV_REAL(ICG),ICG=1,4)

               LL = 1
               UNPACK_VECTOR: DO IWL = 1, OSL%NWL
                  ZSUM = 0.0_EB
                  !ICO = OSL%IWL_TO_ICO(IWL)

                  IWG = OSL%IWL_TO_IWG(IWL)
                  IX = SL%WALL(IWG)%IXG
                  IY = SL%WALL(IWG)%IYG
                  IZ = SL%WALL(IWG)%IZG

                  ICE = SL%WALL(IWG)%ICE(1)

                  VECTOR(ICE) = RECV_REAL(LL)
                  !WRITE(LU_SCARC,'(a,i4,a,i4,a,E20.10,a,3i4)') &
                  !  'IWL=',IWL,':  UNPACK VECTOR(',ICE,')=',VECTOR(ICE), ' : IWG=', IWG

                  !ICO = SL%WALL(IWG)%ICO
                  !DO ICPL = 1, OSL%NCPL
                  !   ZSUM = ZSUM + RECV_REAL(LL)
                  !   LL = LL + 1
                  !ENDDO
                  !VECTOR(ICO) = ZSUM/REAL(OSL%NCPL,EB)

                  LL = LL + 1
               ENDDO UNPACK_VECTOR


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MEASURE_ADD: unpack neighboring measure information (AMG only) with adding
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MEASURE_ADD)

               LL = 1
               UNPACK_MEASURE_ADD: DO IWL = 1, OSL%NWL
                  ICW = OSL%IWL_TO_ICW(IWL)
                  SL%MEASURE(ICW) = SL%MEASURE(ICW)  + RECV_REAL(LL)
                  !WRITE(LU_SCARC,'(a,i4,a,3f12.6)') &
                  !   'UNPACKA:  SL%MEASURE(',ICW,')=', SL%MEASURE(ICW),  RECV_REAL(LL)
                  LL = LL + 1
               ENDDO UNPACK_MEASURE_ADD

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_CELL_TYPE: unpack neighboring CELL_TYPE information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_CELL_TYPE)

!WRITE(LU_SCARC,*) 'RECEIVING RECV_REAL(1:8)', RECV_INTEGER(1:8), SIZE(RECV_INTEGER)
               LL = 1
               UNPACK_CELL_TYPE: DO IWL = 1, OSL%NWL

                  ZSUM1=0.0_EB
                  DO ICPL = 1, OSL%NCPL
                     ICG = OSL%IWL_TO_ICG(IWL, ICPL)
                     OSL%CELL_TYPE(ICG) = RECV_INTEGER(LL)
                  !WRITE(LU_SCARC,*) 'RECEIVING OSL%CELL_TYPE(',ICG,')=',OSL%CELL_TYPE(ICG)
                     LL = LL+1
                  ENDDO

               ENDDO UNPACK_CELL_TYPE


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_SUBDIAG: unpack subdiagonal information from !neighboring system matrix
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG)

               !WRITE(LU_SCARC,*) 'EXCHANGE_MATRIX_SUBDIAG, UNPACK: Achtung, hier nochmal checken f�r AMG!!'
               !WRITE(LU_SCARC,*) 'NSCARC_EXCHANGE_MATRIX_SUBDIAG: RECV_REAL'
               !WRITE(LU_SCARC,'(8f12.6)') (RECV_REAL(IC), IC=1, SIZE(RECV_REAL))

               LL = 1
               UNPACK_MATRIX_SUBDIAG: DO IWL = 1, OSL%NWL

                  ZSUM=0.0_EB
                  ! ACHTUNG WEGEN NCPL SCHAUEN 
                  IWG = OSL%IWL_TO_IWG(IWL)
                  IX = SL%WALL(IWG)%IXG
                  IY = SL%WALL(IWG)%IYG
                  IZ = SL%WALL(IWG)%IZG
                  IC = SL%CCVAR(IX, IY, IZ, UNKH)
                  ZSUM=ZSUM+RECV_REAL(LL)
                  LL = LL+1

                  IC=SL%WALL(IWG)%ICW

                  IROW = SL%A_ROW(IC)
                  DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                     IF (IW == -SL%A_COL(ICOL)) THEN
                        SL%A(ICOL)     = ZSUM/REAL(SL%WALL(IW)%NCPL,EB)
                        SL%A_COL(ICOL) = SL%WALL(IW)%ICE(1)
                        !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)') &
                        !'UNPACK_MATRIX_SUBDIAG1: NM=',NM,': NOM=',NOM,': A_COL(',ICOL,')=', &
                        !SL%A_COL(ICOL),' A(',ICOL,')=',SL%A(ICOL)
                     ENDIF
                  ENDDO

               ENDDO UNPACK_MATRIX_SUBDIAG

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_STENCIL: unpack stencil information from neighboring system matrix
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_STENCIL)

                    !WRITE(LU_SCARC,*) 'TRYING TO EXCHANGE MATRIX STENCIL'
               LL = 1
               UNPACK_MATRIX_STENCIL: DO IWL = 1, OSL%NWL
                  !IWG = OSL%IWL_TO_IWG(IWL)
                  !ICO = OSL%IWL_TO_ICO(IWL)
                  !ICO = SL%WALL(IWG)%ICO

                  OSL%A_SIZE(IWL) = NINT(RECV_REAL(LL))
                  !WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3)') &
                  !'UNPACK_MATRIX_STENCIL1: NM=',NM,': NOM=',NOM,': A_SIZE(',IWL,')=',OSL%A_SIZE(IWL)

                  LL = LL + 1

                ENDDO UNPACK_MATRIX_STENCIL

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_SYSTEM: extract neighboring system matrix information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SYSTEM)

               LL = 1
               UNPACK_MATRIX_SYSTEM: DO IWL = 1, OSL%NWL

                  IWG = OSL%IWL_TO_IWG(IWL)

                  ICE = SL%WALL(IWG)%ICE(1)
                  IX = SL%WALL(IWG)%IXG
                  IY = SL%WALL(IWG)%IYG
                  IZ = SL%WALL(IWG)%IZG
                  ICN = SL%CCVAR(IX, IY, IZ, UNKH)
                  ICG = SL%WALL(IWG)%ICG(1)
                  ICW = SL%WALL(IWG)%ICW

                  !WRITE(LU_SCARC,'(7(a,i3))') 'IWL=',IWL,': IWG=',IWG,': ICE=',ICE,': ICN=',ICN,': ICG=',ICG,': ICW=',ICW

                  DO ICOL = SL%A_ROW(ICE), SL%A_ROW(ICE+1)-1
                     SL%A_COL(ICOL)= NINT(RECV_REAL(LL))
                     SL%A(ICOL)    = RECV_REAL(LL+1)
                     !WRITE(LU_SCARC,'(6(a,i3),a,f12.6)') &
                     !'UNPACK_MATRIX_SYSTEM1: NM=',NM,': NOM=',NOM,': ICE=',ICE,&
                     !': A_COL(',ICOL,')=',SL%A_COL(ICOL),' A(',ICOL,')=',SL%A(ICOL)
                     LL = LL + 2
                  ENDDO

                ENDDO UNPACK_MATRIX_SYSTEM


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_PROL: Extract neighboring prolongation matrix information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_PROL)

               LL = 1
               NROW=1
               ICC = 1
               UNPACK_MATRIX_PROL: DO IWL = 1, OSL%NWL

                  IWG = OSL%IWL_TO_IWG(IWL)

                  ICG = SL%WALL(IWG)%ICG(1)
                  ICE = SL%WALL(IWG)%ICE(1)
                  IX = SL%WALL(IWG)%IXG
                  IY = SL%WALL(IWG)%IYG
                  IZ = SL%WALL(IWG)%IZG
                  ICN = SL%CCVAR(IX, IY, IZ, UNKH)
                  ICELL_TYPE = NINT(RECV_REAL(LL))
                  NCOL      = NINT(RECV_REAL(LL+1))

                  !WRITE(LU_SCARC,'(a,i4,a,f12.6,a,i4,a,f12.6)') &
                  !    'RECV_REAL(',LL,')=',RECV_REAL(LL), ': RECV_REAL(', LL+1,')=', RECV_REAL(LL+1)
                  !WRITE(LU_SCARC,'(a,i4,a,i4,a,i4,a,i4,a,i8)') &
                  !    '========================== IWG=',IWG,': NCOL=',NCOL,': ICG=',ICG,': ICE=',ICE,': ICN=',ICN
                  OSL%CELL_TYPE(ICG) = ICELL_TYPE
                  OSL%P_ROW(ICG) = NROW
                  !WRITE(LU_SCARC,*) 'UNPACK PROLONGATION: OSL%CELL_TYPE(',ICG,')=',OSL%CELL_TYPE(ICG)
                  !WRITE(LU_SCARC,*) 'UNPACK PROLONGATION: OSL%P_ROW(',ICG,')   =',OSL%P_ROW(ICG)

                  LL = LL + 2

                  DO ICOL = NROW, NROW+NCOL-1
                     OSL%P_COL(ICOL)= NINT(RECV_REAL(LL))
                     OSL%P(ICOL)    = RECV_REAL(LL+1)
                  !WRITE(LU_SCARC,*) 'UNPACK PROLONGATION: OSL%P_COL(',ICOL,')   =',OSL%P_COL(ICOL)
                  !WRITE(LU_SCARC,*) 'UNPACK PROLONGATION: OSL%P(',ICOL,')   =',OSL%P(ICOL)
                     LL = LL + 2
                  ENDDO
                  NROW = NROW + NCOL
                  OSL%P_ROW(ICG+1) = NROW

               ENDDO UNPACK_MATRIX_PROL


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_REST: exchange neighboring restriction matrix information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_REST)

               LL = 1
               NROW=1
               ICC = 1
               UNPACK_MATRIX_REST: DO

                  IW = NINT(RECV_REAL(LL))
                  LL = LL + 1
                  IF (IW==-888.OR.IW==-999) EXIT UNPACK_MATRIX_REST

                  WRITE(*,*) 'ACHTUNG, MUSS WIEDER AUSGEARBEITET WERDEN!'

                  !ICG = SL%WALL(IW)%ICG(1)
                  !ICE = SL%WALL(IW)%ICE(1)
                  !ICN = SL%WALL(IW)%ICN(1)
                  !ICELL_TYPE = NINT(RECV_REAL(LL))
                  !NCOL      = NINT(RECV_REAL(LL+1))
                 !>
                 !> OSL%CELL_TYPE(ICG) = ICELL_TYPE
                 !> OSL%P_ROW(ICG) = NROW
         !>
         !>        LL = LL + 2
         !>
         !>        DO ICOL = NROW, NROW+NCOL-1
         !>           OSL%P_COL(ICOL)= NINT(RECV_REAL(LL))
         !>           OSL%P(ICOL)    = RECV_REAL(LL+1)
         !>           LL = LL + 2
         !>        ENDDO
         !>        NROW = NROW + NCOL
         !>        OSL%P_ROW(ICG+1) = NROW

               ENDDO UNPACK_MATRIX_REST

          END SELECT OMESH_UNPACK_SELECT
      ENDIF OMESH_UNPACK_IF
   ENDDO OMESH_UNPACK_LOOP
ENDDO MESH_UNPACK_LOOP

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
!> Debug requested quantity
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_QUANTITY(NTYPE, NL, CROUTINE, CNAME)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM, IP, IC, IW, I, J, K, IOR0, INBR
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
CHARACTER (20) :: LINE
TYPE (MESH_TYPE), POINTER :: M
TYPE (OSCARC_TYPE) , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, OSL
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: SM
#endif

IF (TYPE_DEBUG < NSCARC_DEBUG_MEDIUM) RETURN


SELECT CASE (NTYPE)

   !! ------------------------------------------------------------------------------------------------
   !! Debug system matrix A (corresponding to system type)
   !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIX)


      IF (NL > NLEVEL_MIN) RETURN
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,*) '----------- SHOWING FULL MATRIX ENTRIES'
         WRITE(LU_SCARC,*) 'NA =',SL%NA
         WRITE(LU_SCARC,*) 'NC =',SL%NCS
         WRITE(LU_SCARC,*) 'SL%NCS =',SL%NCS
         WRITE(LU_SCARC,*) 'SIZE(A) =',SIZE(SL%A)
         WRITE(LU_SCARC,*) 'SIZE(A_COL) =',SIZE(SL%A_COL)
         WRITE(LU_SCARC,*) 'SIZE(A_ROW) =',SIZE(SL%A_ROW)
         WRITE(LU_SCARC,*) '---------------------- A_ROW:', SL%NCS
         WRITE(LU_SCARC,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NCS+1)
         WRITE(LU_SCARC,*) '---------------------- A_COL:'
         DO IC = 1, SL%NCS
            WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%A_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) '---------------------- A:', SL%A_ROW(IC+1)- SL%A_ROW(IC)
         DO IC = 1, SL%NCS
            WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SL%A(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
         ENDDO
#ifdef WITH_MKL
         IF ((TYPE_METHOD == NSCARC_METHOD_MKL) .AND. (TYPE_MKL == NSCARC_MKL_GLOBAL)) THEN
            SM => SCARC(NM)%MKL(NL)
            WRITE(LU_SCARC,*) '---------------------- AG_COL:'
            DO IC = 1, SL%NCS
               WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SM%AG_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
            ENDDO
         ENDIF
         DO IC=1,SL%NCS
            DO IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1
                WRITE(LU_SCARC,'(2I8,F18.12)') IC,SL%A_COL(IP),SL%A(IP)
            ENDDO
            WRITE(LU_SCARC,*)
         ENDDO
         DO IC=1,SL%NCS
            DO IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1
                IF (IC == SL%A_COL(IP)) WRITE(LU_SCARC,'(2I8,F18.12)') IC,SL%A_COL(IP),SL%A(IP)
            ENDDO
         ENDDO
     ENDDO
#endif

         !CALL SCARC_PRINT_MATRIX (SL%A, SL%A_ROW, SL%A_COL, SL%NCS, SL%NCS, NM, NL, 'A')
         !CALL SCARC_PRINT_MATRIX2(SL%A, SL%A_ROW, SL%A_COL, SL%NCS, SL%NCS, SL%NX, SL%NY, SL%NZ, NM, NL, 'A')
         !CALL SCARC_MATLAB_MATRIX(SL%A, SL%A_ROW, SL%A_COL, SL%NCS, SL%NCS, NM, NL, 'A')

#ifdef WITH_MKL
   !! ------------------------------------------------------------------------------------------------
   !! Debug symmetric system matrix AS
   !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIXS)

      IF (NL > NLEVEL_MIN) RETURN
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         SM => SCARC(NM)%MKL(NL)
         WRITE(LU_SCARC,*) '----------- SHOWING SYMMETRIC MATRIX ENTRIES'
         WRITE(LU_SCARC,*) 'NA =',SL%NA
         WRITE(LU_SCARC,*) 'NCS=',SL%NCS
         WRITE(LU_SCARC,*) 'SIZE(ASYM) =',SIZE(SM%ASYM)
         WRITE(LU_SCARC,*) 'SIZE(ASYM_COL) =',SIZE(SM%ASYM_COL)
         WRITE(LU_SCARC,*) 'SIZE(ASYM_ROW) =',SIZE(SM%ASYM_ROW)
         WRITE(LU_SCARC,*) '---------------------- ASYM_ROW:', SL%NCS
         WRITE(LU_SCARC,'(4i9)') (SM%ASYM_ROW(IC), IC=1,SL%NCS+1)
         WRITE(LU_SCARC,*) '---------------------- A_COL:'
         DO IC = 1, SL%NC
            WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SM%ASYM_COL(IP),IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1)
         ENDDO
         !WRITE(LU_SCARC,*) '---------------------- AS:', SM%ASYM_ROW(IC+1)- SM%ASYM_ROW(IC)
         DO IC = 1, SL%NC
            WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SM%ASYM(IP),IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1)
         ENDDO
         DO IC=1,SL%NC
            DO IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1
                WRITE(LU_SCARC,'(2I8,F18.12)') IC,SM%ASYM_COL(IP),SM%ASYM(IP)
            ENDDO
            WRITE(LU_SCARC,*)
         ENDDO
            WRITE(LU_SCARC,*)
            WRITE(LU_SCARC,*)
         DO IC=1,SL%NC
            DO IP=SM%ASYM_ROW(IC),SM%ASYM_ROW(IC+1)-1
                IF (IC == SM%ASYM_COL(IP)) WRITE(LU_SCARC,'(2I8,F18.12)') IC,SM%ASYM_COL(IP),SM%ASYM(IP)
            ENDDO
            WRITE(LU_SCARC,*)
         ENDDO

         !CALL SCARC_PRINT_MATRIX (SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, SL%NC, SL%NC, NM, NL, 'AS')
         !CALL SCARC_PRINT_MATRIX2(SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, SL%NC, SL%NC, SL%NX, SL%NY, SL%NZ, NM, NL, 'AS')
         !CALL SCARC_MATLAB_MATRIX(SM%ASYM, SM%ASYM_ROW, SM%ASYM_COL, SL%NC, SL%NC, NM, NL, 'AS')
      ENDDO
#endif

   !! ------------------------------------------------------------------------------------------------
   !! Debug extended system matrix A (corresponding to system type)
   !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIXE)

      IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,*) 'PRINTING EXTENDED MATRIX, NCE=',SL%NC
         WRITE(LU_SCARC,*) '---------------------- A_ROW:', SL%NC
         WRITE(LU_SCARC,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NC+1)
         WRITE(LU_SCARC,*) '---------------------- A_COL:'
         DO IC = 1, SL%NC
            WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(SL%A_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) '---------------------- A:'
         DO IC = 1, SL%NC
            WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(SL%A(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) 'SIZE(A) =',SIZE(SL%A)
         WRITE(LU_SCARC,*) 'SIZE(A_COL) =',SIZE(SL%A_COL)
         WRITE(LU_SCARC,*) 'SIZE(A_ROW) =',SIZE(SL%A_ROW)
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug prolongation matrix P (only for compact system)
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_PROLONGATION)

      IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         WRITE(LU_SCARC,*) '============= P_ROW:', NM
         WRITE(LU_SCARC,'(4i9)') (SL%P_ROW(IC), IC=1, SL%NCE+1)
         WRITE(LU_SCARC,*) '============= P_COL:', NM, SL%NCE
         DO IC=1,SL%NCE
            WRITE(LU_SCARC,'(i5,a,10i9)') IC,':',(SL%P_COL(IP), IP=SL%P_ROW(IC),SL%P_ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) '============= P:', NM, SL%NCE
         DO IC=1,SL%NCE
            WRITE(LU_SCARC,'(i5,a,10f9.2)') IC,':',(SL%P(IP), IP=SL%P_ROW(IC),SL%P_ROW(IC+1)-1)
         ENDDO
         DO NOM = 1, NMESHES
            IF (NOM == NM) CYCLE
            OS  => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE
            OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P_ROW:', OSL%NCG
            WRITE(LU_SCARC,'(4i9)') (OSL%P_ROW(IC), IC=1, OSL%NCG+1)
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P_COL:', OSL%NCG
            DO IC=1,OSL%NCG
               WRITE(LU_SCARC,'(i5,a,10i9)') IC,':',(OSL%P_COL(IP), IP=OSL%P_ROW(IC),OSL%P_ROW(IC+1)-1)
            ENDDO
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P:'
            DO IC=1,OSL%NCG
               WRITE(LU_SCARC,'(i5,a,10f9.2)') IC,':',(OSL%P(IP), IP=OSL%P_ROW(IC),OSL%P_ROW(IC+1)-1)
            ENDDO
         ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug restriction matrix R (only for compact system)
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_RESTRICTION)

      IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         WRITE(LU_SCARC,*) '============= R_ROW:', NM, SL%NCCE
         WRITE(LU_SCARC,'(4i9)') (SL%R_ROW(IC), IC=1, SL%NCCE)
         WRITE(LU_SCARC,*) '============= R_COL:', NM
         DO IC=1,SL%NCCE
            WRITE(LU_SCARC,'(i5,a,10i9)') IC,':',(SL%R_COL(IP), IP=SL%R_ROW(IC),SL%R_ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) '============= R:', NM
         DO IC=1,SL%NCCE
            WRITE(LU_SCARC,'(i5,a,10f9.2)') IC,':',(SL%R(IP), IP=SL%R_ROW(IC),SL%R_ROW(IC+1)-1)
         ENDDO
         !DO NOM = 1, NMESHES
         DO NOM = 1, NMESHES
            IF (NOM == NM) CYCLE
            OS  => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE
            OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%R_ROW:', OSL%NCC
            WRITE(LU_SCARC,'(4i9)') (OSL%R_ROW(IC), IC=1, OSL%NCC+1)
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%R_COL:', OSL%NCC
            DO IC=1,OSL%NCC
               WRITE(LU_SCARC,'(i5,a,10i9)') IC,':',(OSL%R_COL(IP), IP=OSL%R_ROW(IC),OSL%R_ROW(IC+1)-1)
            ENDDO
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P:', OSL%NCC
            DO IC=1,OSL%NCC
               WRITE(LU_SCARC,'(i5,a,10f9.2)') IC,':',(OSL%R(IP), IP=OSL%R_ROW(IC),OSL%R_ROW(IC+1)-1)
            ENDDO
         ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug FACEINFO
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_FACEINFO)

      IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         M  => MESHES(NM)
         DO IOR0 = -3, 3
            IF (IOR0 == 0) CYCLE
            WRITE(LU_SCARC,*) '========================================='
            WRITE(LU_SCARC,*) '============= DEBUGGING FACE(',IOR0,'): FOR LEVEL ', NL
            WRITE(LU_SCARC,*) '========================================='
            WRITE(LU_SCARC,*) 'FACE(.)%NFM:', SL%FACE(IOR0)%NFM
            WRITE(LU_SCARC,*) 'FACE(.)%NFX:', SL%FACE(IOR0)%NFX
            WRITE(LU_SCARC,*) 'FACE(.)%NFY:', SL%FACE(IOR0)%NFY
            WRITE(LU_SCARC,*) 'FACE(.)%NFZ:', SL%FACE(IOR0)%NFZ
            WRITE(LU_SCARC,*) 'FACE(.)%NWL:', SL%FACE(IOR0)%NFW
            WRITE(LU_SCARC,*) 'FACE(.)%IWG_MARKER:', SL%FACE(IOR0)%IWG_MARKER
            WRITE(LU_SCARC,*) '----------------------------------------------'
            WRITE(LU_SCARC,*) 'FACE(.)%BXY:', SL%FACE(IOR0)%BXY
            WRITE(LU_SCARC,*) 'FACE(.)%BXZ:', SL%FACE(IOR0)%BXZ
            WRITE(LU_SCARC,*) 'FACE(.)%BYZ:', SL%FACE(IOR0)%BYZ
            WRITE(LU_SCARC,*)
            !WRITE(LU_SCARC,*) 'FACE(.)%DH :', SL%FACE(IOR0)%DH
            WRITE(LU_SCARC,*) 'FACE(.)%IOFFSET_MATRIX:', SL%FACE(IOR0)%IOFF_MAT
            WRITE(LU_SCARC,*)
            IF (SL%FACE(IOR0)%NUM_NEIGHBORS /= 0) THEN
               DO INBR=1,SL%FACE(IOR0)%NUM_NEIGHBORS
                  NOM = SL%FACE(IOR0)%NEIGHBORS(INBR)
                  OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  WRITE(LU_SCARC,*) 'NUM_NEIGHBORS:', SL%FACE(IOR0)%NUM_NEIGHBORS
                  WRITE(LU_SCARC,*) 'NOM:', NOM
                  WRITE(LU_SCARC,*) 'SIZE(OSL%ICG_TO_IWG)=',SIZE(OSL%ICG_TO_IWG)
                  WRITE(LU_SCARC,*) 'SIZE(OSL%ICG_TO_ICE)=',SIZE(OSL%ICG_TO_ICE)
                  WRITE(LU_SCARC,*) 'SIZE(OSL%ICG_TO_ICO)=',SIZE(OSL%ICG_TO_ICO)
                  WRITE(LU_SCARC,*) 'SIZE(OSL%IWL_TO_ICW)=',SIZE(OSL%IWL_TO_ICW)
                  WRITE(LU_SCARC,*) 'SIZE(OSL%IWL_TO_ICO)=',SIZE(OSL%IWL_TO_ICO)
                  WRITE(LU_SCARC,*) 'SIZE(OSL%IWL_TO_ICG)=',SIZE(OSL%IWL_TO_ICG)
                  WRITE(LU_SCARC,'(a,i8,a,2f12.6)') '---OSL(',NOM,')%DH :',OSL%DH
                  WRITE(LU_SCARC,'(a,i8,a,2i8)') '---OSL(',NOM,')%IOR   :',OSL%IOR
                  WRITE(LU_SCARC,'(a,i8,a,2i8)') '---OSL(',NOM,')%NWL(.):',OSL%NWL
                  WRITE(LU_SCARC,'(a,i8,a,2i8)') '---OSL(',NOM,')%NCG(.):',OSL%NCG
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OSL(',NOM,')%ICG_TO_IWG:'
                  DO IW = 1, OSL%NCG
                     WRITE(LU_SCARC,'(32i8)') OSL%ICG_TO_IWG(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OSL(',NOM,')%ICG_TO_ICO:'
                  DO IW = 1, OSL%NCG
                     WRITE(LU_SCARC,'(32i8)') OSL%ICG_TO_ICO(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OSL(',NOM,')%ICG_TO_ICE:'
                  DO IW = 1, OSL%NCG
                     WRITE(LU_SCARC,'(32i8)') OSL%ICG_TO_ICE(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OSL(',NOM,')%IWL_TO_IWG:'
                  DO IW = 1, OSL%NWL
                     WRITE(LU_SCARC,'(32i8)') OSL%IWL_TO_IWG(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OSL(',NOM,')%IWL_TO_ICW:'
                  DO IW = 1, OSL%NWL
                     WRITE(LU_SCARC,'(32i8)') OSL%IWL_TO_ICW(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OSL(',NOM,')%IWL_TO_ICO:'
                  DO IW = 1, OSL%NWL
                     WRITE(LU_SCARC,'(32i8)') OSL%IWL_TO_ICO(IW)
                  ENDDO
                  !WRITE(LU_SCARC,*)
                  !WRITE(LU_SCARC,'(a,i4,a)') '------OSL(',NOM,')%IWL_TO_ICG:'
                  !DO IW = 1, OSL%NWL
                  !   WRITE(LU_SCARC,'(32i8)') OSL%IWL_TO_ICG(IW, 1:OSL%NCPL)
                  !ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug WALLINFO
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_WALLINFO)
      IF (TYPE_DEBUG < NSCARC_DEBUG_MUCH) RETURN

!WRITE(LU_SCARC,*) ': NL=',NL, ': CNAME=',CNAME,': CROUTINE=',CROUTINE
      !IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC) THEN
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         M  => MESHES(NM)
            WRITE(LU_SCARC,*) 'SIZE(SL%ICE_TO_IWG)=',SIZE(SL%ICE_TO_IWG)
            WRITE(LU_SCARC,*) 'SIZE(SL%ICE_TO_IWL)=',SIZE(SL%ICE_TO_IWL)
            WRITE(LU_SCARC,*) 'SIZE(SL%ICE_TO_ICG)=',SIZE(SL%ICE_TO_ICG)
            WRITE(LU_SCARC,*) 'NM  =',NM
            WRITE(LU_SCARC,*) 'NL  =',NL
            WRITE(LU_SCARC,*) 'NCE =',SL%NCE
            WRITE(LU_SCARC,*) 'NCS =',SL%NCS
            WRITE(LU_SCARC,*) 'NC  =',SL%NC
            WRITE(LU_SCARC,*) 'NW  =',SL%NW
            WRITE(LU_SCARC,*)
            WRITE(LU_SCARC,'(a,i8,a)') '------SL%ICE_TO_IWG:'
            DO IW = SL%NCS+1, SL%NCE
               WRITE(LU_SCARC,'(16i8)') SL%ICE_TO_IWG(IW)
            ENDDO
            WRITE(LU_SCARC,'(a,i8,a)') '------SL%ICE_TO_IWL:'
            DO IW = SL%NCS+1, SL%NCE
               WRITE(LU_SCARC,'(16i8)') SL%ICE_TO_IWL(IW)
            ENDDO
            WRITE(LU_SCARC,*)
            WRITE(LU_SCARC,'(a,i8,a)') '------SL%ICE_TO_ICG:'
            DO IW = SL%NCS+1, SL%NCE
               WRITE(LU_SCARC,'(16i8)') SL%ICE_TO_ICG(IW)
            ENDDO
            WRITE(LU_SCARC,*)
            WRITE(LU_SCARC,*)
            WRITE(LU_SCARC,*)
         IF (NL == 1) THEN
         WRITE(LU_SCARC,*) '============= WALL(.)%IXG:', NM
         WRITE(LU_SCARC,'(16i8)') (SL%WALL(IW)%IXG, IW=1,SL%NW)
         WRITE(LU_SCARC,*) '============= WALL(.)%IYG:', NM
         WRITE(LU_SCARC,'(16i8)') (SL%WALL(IW)%IYG, IW=1,SL%NW)
         WRITE(LU_SCARC,*) '============= WALL(.)%IZG:', NM
         WRITE(LU_SCARC,'(16i8)') (SL%WALL(IW)%IZG, IW=1,SL%NW)
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%IXW:', NM
         WRITE(LU_SCARC,'(16i8)') (SL%WALL(IW)%IXW, IW=1,SL%NW)
         WRITE(LU_SCARC,*) '============= WALL(.)%IYW:', NM
         WRITE(LU_SCARC,'(16i8)') (SL%WALL(IW)%IYW, IW=1,SL%NW)
         WRITE(LU_SCARC,*) '============= WALL(.)%IZW:', NM
         WRITE(LU_SCARC,'(16i8)') (SL%WALL(IW)%IZW, IW=1,SL%NW)
         ENDIF
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%IOR:', NM
         WRITE(LU_SCARC,'(16i8)') (SL%WALL(IW)%IOR, IW=1,SL%NW)
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%NOM:', NM
         WRITE(LU_SCARC,'(16i8)') (SL%WALL(IW)%NOM, IW=1,SL%NW)
         WRITE(LU_SCARC,*)
         IF (NL == 1) THEN
         WRITE(LU_SCARC,*) '============= WALL(.)%NCPL:', NM
         !WRITE(LU_SCARC,'(16i10)') (SL%WALL(IW)%NCPL, IW=1,SL%NW)
         WRITE(LU_SCARC,*) (SL%WALL(IW)%NCPL, IW=1,SL%NW)
         ENDIF
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%ICW:', NM
         DO IW=1,SL%NW
            WRITE(LU_SCARC,'(a,i8, a,i8)') 'IW=',IW,':',SL%WALL(IW)%ICW
         ENDDO
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%ICO:', NM
         DO IW=1,SL%NW
            IF (SL%WALL(IW)%NOM /=0) &
               WRITE(LU_SCARC,'(a,i8, a,i8)') 'IW=',IW,':',SL%WALL(IW)%ICO
         ENDDO
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%ICE:', NM
         DO IW=1,SL%NW
            IF (SL%WALL(IW)%NOM /=0) &
               WRITE(LU_SCARC,'(a,i8, a,i8)') 'IW=',IW,':',SL%WALL(IW)%ICE
         ENDDO
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%ICG:', NM
         DO IW=1,SL%NW
            IF (SL%WALL(IW)%NOM /=0) &
               WRITE(LU_SCARC,'(a,i8, a,i8)') 'IW=',IW,':',SL%WALL(IW)%ICG
         ENDDO
         WRITE(LU_SCARC,*) '====================================================='
         WRITE(LU_SCARC,*) ' Plotting out M%WALL-structure'
         WRITE(LU_SCARC,*) '====================================================='
         WRITE(LU_SCARC,*) ' M%WALL(.)%WALL_INDEX'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%WALL_INDEX, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%SURF_INDEX'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%SURF_INDEX, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%BACK_INDEX'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%BACK_INDEX, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%BOUNDARY_TYPE'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%BOUNDARY_TYPE, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%OBST_INDEX'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%OBST_INDEX, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%PRESSURE_BC_INDEX'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%PRESSURE_BC_INDEX, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%PRESSURE_ZONE'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%PRESSURE_ZONE, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%VENT_INDEX'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%VENT_INDEX, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%NOM'
         WRITE(LU_SCARC,'(16i8)') (M%EXTERNAL_WALL(IW)%NOM, IW=1,SL%N_EXTERNAL_WALL_CELLS)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%II'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%ONE_D%II, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%JJ'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%ONE_D%JJ, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%KK'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%ONE_D%KK, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%IIG'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%ONE_D%IIG, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%JJG'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%ONE_D%JJG, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%KKG'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%ONE_D%KKG, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%N_LAYER_CELLS'
         WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%ONE_D%IOR, IW=1,SL%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%IOR'
         !WRITE(LU_SCARC,'(16i8)') (M%WALL(IW)%ONE_D%N_LAYER_CELLS, IW=1,SL%NW)



      ENDDO
      !ENDIF


 !! ------------------------------------------------------------------------------------------------
 !! Debug PRESSURE_BC_INDEX
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_BCINDEX)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC, '(8i8)') (SL%WALL(J)%BTYPE, J=1,SL%NW)
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug WALL_CELL
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_ACELL)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC, '(8i8)') (SL%WALL(IW)%ICW, IW=1,SL%NW)
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug GHOST_CELL
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_GCELL)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,*) '-------- GHOST,  NGE=',SL%NCG
         DO IW=1,SL%NW
            !IF (SL%WALL(IW)%NOM/=0) &
            !>   WRITE(LU_SCARC, '(8i8)') (SL%WALL(IW)%ICG(I), I=1,SL%WALL(IW)%NCPL)
         ENDDO
         !DO IW=1,SL%NW
         !>   IF (SL%WALL(IW)%NOM/=0.AND.NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) &
         !>      WRITE(LU_SCARC, '(8i8)') (SL%WALL(IW)%ICG2(I), I=1,SL%WALL(IW)%NCPL)
         !ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug NOM_CELL
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_NCELL)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,*) '-------- NOM, NCE= ', SL%NCE
         WRITE(*,*) 'Achtung, hier wegen ICN schauen!'
         !DO IW=1,SL%NW
         !   !IF (SL%WALL(IW)%NOM/=0) &
         !   !>   WRITE(LU_SCARC, '(8i8)') (SL%WALL(IW)%ICN(I), I=1,SL%WALL(IW)%NCPL)
         !ENDDO
         !DO IW=1,SL%NW
         !>   IF (SL%WALL(IW)%NOM/=0.AND.NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) &
         !>      WRITE(LU_SCARC, '(8i8)') (SL%WALL(IW)%ICN2(I), I=1,SL%WALL(IW)%NCPL)
         !ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug SUBDIVISION
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_SUBDIVISION)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR= 1 '
         WRITE(LU_SCARC,'(3i8)') (SL%SUBDIVISION(I, 1), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR=-1 '
         WRITE(LU_SCARC,'(3i8)') (SL%SUBDIVISION(I,-1), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR= 2 '
         WRITE(LU_SCARC,'(3i8)') (SL%SUBDIVISION(I, 2), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR=-2 '
         WRITE(LU_SCARC,'(3i8)') (SL%SUBDIVISION(I,-2), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR= 3 '
         WRITE(LU_SCARC,'(3i8)') (SL%SUBDIVISION(I, 3), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR=-3 '
         WRITE(LU_SCARC,'(3i8)') (SL%SUBDIVISION(I,-3), I=1,3)
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug MEASURE
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MEASURE)

      WRITE(LU_SCARC,1000) CROUTINE, CNAME, 1, NL
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL  => SCARC(NM)%LEVEL(NL)
         IF (NL == 1) THEN
            DO K = SL%NZ, 1, -1
               WRITE(LU_SCARC, '(10f11.5)') (SL%MEASURE((K-1)*SL%NX*SL%NY+I), I=1, SL%NX)
            ENDDO
            WRITE(LU_SCARC,*) '--------------------'
            WRITE(LU_SCARC, '(4f11.5)') (SL%MEASURE(IC), IC=SL%NC+1, SL%NCE)
         ELSE
            WRITE(LU_SCARC, '(4f11.5)') (SL%MEASURE(IC), IC=1, SL%NC)
         ENDIF
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug CELL_TYPE
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_CELL_TYPE)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
        WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
        SL  => SCARC(NM)%LEVEL(NL)
        IF (NL == 1) THEN
           DO K = SL%NZ, 1, -1
              WRITE(LU_SCARC, '(8i6)') (SL%CELL_TYPE((K-1)*SL%NX*SL%NY+I), I=1, SL%NX)
           ENDDO
           WRITE(LU_SCARC,*) '--------------------'
           WRITE(LU_SCARC, '(4i6)') (SL%CELL_TYPE(IC), IC=SL%NC+1, SL%NCE)
        ELSE
           WRITE(LU_SCARC, '(4i6)') (SL%CELL_TYPE(IC), IC=1, SL%NC)
        ENDIF
      ENDDO
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL  => SCARC(NM)%LEVEL(NL)
         DO NOM = 1, NMESHES
            IF (NOM == NM) CYCLE
            OS  => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE
            OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%CELL_TYPE:', OSL%NCG
            DO IC=1,OSL%NCG
               WRITE(LU_SCARC,'(i8,a,i8)') IC,':',OSL%CELL_TYPE(IC)
            ENDDO
         ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug GRAPH
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_GRAPH)

      IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) THEN
       IF (NMESHES == 4 .OR. NMESHES==1) THEN                !> only temporarily
         DO NM = 1, NMESHES
            IF (PROCESS(NM) /= MYID) CYCLE
           WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
           SL  => SCARC(NM)%LEVEL(NL)
           IF (NL == 1) THEN
              DO K = SL%NZ, 1, -1
                 WRITE(LU_SCARC, '(8i6)') (SL%GRAPH((K-1)*SL%NX*SL%NY+I), I=1, SL%NX)
              ENDDO
              WRITE(LU_SCARC,*) '--------------------'
              WRITE(LU_SCARC, '(4i6)') (SL%GRAPH(IC), IC=SL%NC+1, SL%NCE)
           ELSE
              WRITE(LU_SCARC, '(4i6)') (SL%GRAPH(IC), IC=1, SL%NC)
           ENDIF
         ENDDO
       ELSEIF (NMESHES == 4.AND.NL==1) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(LU_SCARC,'(10i8)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(2)%LEVEL(NL)
         WRITE(LU_SCARC,'(10i8)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(3)%LEVEL(NL)
         WRITE(LU_SCARC,'(9i8)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(4)%LEVEL(NL)
         WRITE(LU_SCARC,'(9i8)') SL%GRAPH(1:SL%NC)
       ELSEIF (NMESHES == 4.AND.NL==2) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(LU_SCARC,'(5i8)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(2)%LEVEL(NL)
         WRITE(LU_SCARC,'(5i8)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(3)%LEVEL(NL)
         WRITE(LU_SCARC,'(4i8)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(4)%LEVEL(NL)
         WRITE(LU_SCARC,'(4i8)') SL%GRAPH(1:SL%NC)
       ELSEIF (NMESHES==1.AND.NL==1) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(LU_SCARC,'(9i8)') SL%GRAPH(1:SL%NC)
       ELSEIF (NMESHES==1.AND.NL==2) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(LU_SCARC,'(3i8)') SL%GRAPH(1:SL%NC)
       ENDIF
      ENDIF
      IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,1000) CROUTINE, CNAME, 1, NL
      IF (NMESHES == 1.OR.NL>100) THEN                !> only temporarily
      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         SL  => SCARC(NM)%LEVEL(NL)
         IF (NL == 1) THEN
            DO K = SL%NZ, 1, -1
               WRITE(LU_SCARC, '(10i8)') (SL%GRAPH((K-1)*SL%NX*SL%NY+I), I=1, SL%NX)
            ENDDO
         ELSE
            WRITE(LU_SCARC, '(4i8)') (SL%GRAPH(IC), IC=1, SL%NC)
         ENDIF
      ENDDO
      ENDIF
      ENDIF

 !! ------------------------------------------------------------------------------------------------
 !! Debug COARSE
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_COARSE)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE
         WRITE(LU_SCARC,1000) CROUTINE, CNAME, NM, NL
         SL  => SCARC(NM)%LEVEL(NL)
         DO K = SL%NZ,1,-1
            DO I = 1,SL%NX
               IC = (K-1)*SL%NX + I
               IF (SL%CELL_TYPE(IC) < 0) THEN
                  LINE(I:I) = 'O'
               ELSE IF (SL%CELL_TYPE(IC) > 0) THEN
                  LINE(I:I) = 'X'
               ENDIF
            ENDDO
            WRITE(LU_SCARC,'(12A4)') (LINE(I:I), I=1,SL%NX)
         ENDDO
      ENDDO

END SELECT

1000 FORMAT('======================================================================================',/, &
            '=== ',A20,' : ', A25,' for mesh ',i3,' on level ', i3, /, &
            '======================================================================================')

END SUBROUTINE SCARC_DEBUG_QUANTITY


!> ------------------------------------------------------------------------------------------------
!> Debug full vector information on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_VECTOR (NVECTOR, CROUTINE, CNAME)
INTEGER, INTENT(IN):: NVECTOR
REAL (EB), POINTER, DIMENSION(:)     :: VC
REAL (EB):: VALUES(10)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
TYPE (MESH_TYPE), POINTER :: M

IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN

DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE

   M  => MESHES(NM)
   CALL POINT_TO_VECTOR (NVECTOR, NM, NLEVEL_MIN, VC)

   IF (M%IBAR>8.OR.M%JBAR>8.OR.M%KBAR>8) CYCLE

   NX8=MIN(8,M%IBAR)
   NY8=MIN(8,M%JBAR)
   NZ8=MIN(8,M%KBAR)

   WRITE(LU_SCARC,*) '======================================================================'
   WRITE(LU_SCARC,2000) CROUTINE, CNAME, NM, NLEVEL_MIN
   WRITE(LU_SCARC,*) '======================================================================'
   SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)

      CASE (NSCARC_DIMENSION_TWO)
         DO KK = NZ8, 1, - 1
            DO II=1,NX8
               
<<<<<<< HEAD
               IF (PRES_ON_WHOLE_DOMAIN.OR.&
=======
               IF (PRES_ON_WHOLE_DOMAIN_SCRC.OR.&
>>>>>>> origin/scarc
                  SCARC(NM)%LEVEL(NLEVEL_MIN)%CCVAR(II,JJ,KK,CGSC) == IS_GASPHASE) THEN
                  IC=SCARC(NM)%LEVEL(NLEVEL_MIN)%CCVAR(II,JJ,KK,UNKH)
                  IF (ABS(VC(IC))<1.0E-14_EB) THEN
                     VALUES(II)=0.0_EB
                  ELSE
                     VALUES(II)=VC(IC)
                  ENDIF
               ELSE
                  VALUES(II)=0.0_EB
               ENDIF
            ENDDO
            WRITE(LU_SCARC, '(a,i8,a,8e13.5)') 'J= ',1,' : ',(VALUES(II), II=1, NX8)
         ENDDO
         !WRITE(LU_SCARC,*)  '------------------------------------------------',&
         !>                   '---------------------------------------------------'

      CASE (NSCARC_DIMENSION_THREE)
         DO KK = NZ8, 1, - 1
            WRITE(LU_SCARC,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                        '----------------------------------------'
            DO JJ = NY8, 1, - 1
               DO II=1,NX8
<<<<<<< HEAD
                  IF (PRES_ON_WHOLE_DOMAIN.OR.&
=======
                  IF (PRES_ON_WHOLE_DOMAIN_SCRC.OR.&
>>>>>>> origin/scarc
                     SCARC(NM)%LEVEL(NLEVEL_MIN)%CCVAR(II,JJ,KK,CGSC) == IS_GASPHASE) THEN
                     IF (ABS(VC(IC))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VC(IC)
                     ENDIF
                  ELSE
                     VALUES(II)=0.0_EB
                  ENDIF
               ENDDO
               WRITE(LU_SCARC, '(a,i8,a,8e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
            ENDDO
         ENDDO
         WRITE(LU_SCARC,*)  '------------------------------------------------',&
                            '---------------------------------------------------'
   END SELECT SELECT_DIMENSION

ENDDO


2000 FORMAT('=== ',A,' : ',A,' on mesh ',i8,' on level ',I4)
END SUBROUTINE SCARC_DEBUG_VECTOR


!> ------------------------------------------------------------------------------------------------
!> Only for debugging reasons: print out vector information on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_LEVEL (NVECTOR, CROUTINE, CNAME, NL)
INTEGER, INTENT(IN):: NVECTOR, NL
REAL (EB), POINTER, DIMENSION(:)     :: VC
REAL (EB):: VALUES(0:9)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

IF (TYPE_DEBUG < NSCARC_DEBUG_MUCH) RETURN


DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   SL => SCARC(NM)%LEVEL(NL)
   CALL POINT_TO_VECTOR (NVECTOR, NM, NL, VC)

   IF (SL%NX>8.OR.SL%NY>8.OR.SL%NZ>8) CYCLE

   NX8=MIN(8,SL%NX)
   NY8=MIN(8,SL%NY)
   NZ8=MIN(8,SL%NZ)

   WRITE(LU_SCARC,*) '=========================================================='
   WRITE(LU_SCARC,2001) CROUTINE, CNAME, NM, NL
   WRITE(LU_SCARC,2002) SL%NC, SL%NCE, NX8, NY8, NZ8, NVECTOR, SIZE(VC)
   WRITE(LU_SCARC,*) '=========================================================='
   !IF (NL == NLEVEL_MIN) THEN
         DO KK = NZ8, 1, - 1
            DO JJ = NY8, 1, - 1
               DO II=1,NX8
                  IF (SL%CCVAR(II,JJ,KK,CGSC) == IS_GASPHASE) THEN
                     IC=SL%CCVAR(II,JJ,KK,UNKH)
                     IF (ABS(VC(IC))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VC(IC)
                     ENDIF
                  ELSE
                     VALUES(II)=-999999.0_EB
                  ENDIF
               ENDDO
               WRITE(LU_SCARC, '(10e16.8)') (VALUES(II), II=1, NX8)
            ENDDO
         ENDDO
   !ENDIF

ENDDO

2000 FORMAT('=== ',A,' : ',A,' on mesh ',I8,' on level ',I8, ': NX, NY, NZ=',3I8,': NVECTOR=',I8)
2001 FORMAT('=== ',A,' : ',A,' on mesh ',I8,' on level ',I8)
2002 FORMAT('=== NC = ',I6,' NCE=',I6, ': NX, NY, NZ=',3I4,': NVECTOR=',I3,': Size=',I8)
END SUBROUTINE SCARC_DEBUG_LEVEL

!> ------------------------------------------------------------------------------------------------
!> Save dump of vector in dump-directory
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_QUANTITY (NVECTOR, NL, ITE_MAIN, ITE, CENDING)
INTEGER, INTENT(IN):: NVECTOR, NL, ITE_MAIN, ITE
CHARACTER(*), INTENT(IN):: CENDING
REAL (EB), POINTER, DIMENSION(:)     :: VC
CHARACTER(60) :: FN_DUMP
INTEGER :: LU_DUMP, IC, NM
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

DO NM = 1, NMESHES
   SL => SCARC(NM)%LEVEL(NL)
   CALL POINT_TO_VECTOR (NVECTOR, NM, NL, VC)
   WRITE (FN_DUMP, '(A,A,A,i3.3,A,i3.3,A,I3.3)') &
      'dump/',TRIM(CHID),TRIM(CENDING),NM,'_m',ITE_MAIN,'_cg',ITE
   LU_DUMP = GET_FILE_NUMBER()
   OPEN (LU_DUMP, FILE=FN_DUMP)
   DO IC = 1, SL%NC
      !WRITE(LU_DUMP,'(F25.16)') VC(IC)
      WRITE(LU_DUMP,*) VC(IC)
   ENDDO
   CLOSE(LU_DUMP)
ENDDO

END SUBROUTINE SCARC_DUMP_QUANTITY


!> ------------------------------------------------------------------------------------------------
!> Save dump of vector in dump-directory
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_LEVEL (NVECTOR, NL)
INTEGER, INTENT(IN):: NVECTOR, NL
REAL (EB), POINTER, DIMENSION(:)     :: VC
REAL (EB):: VALUES(0:10)
INTEGER :: NM, II, JJ, KK
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

IF (TYPE_DEBUG < NSCARC_DEBUG_INFO0) RETURN

DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE
   SL => SCARC(NM)%LEVEL(NL)
   CALL POINT_TO_VECTOR (NVECTOR, NM, NL, VC)

   IF (NL == NLEVEL_MIN) THEN
         DO KK = SL%NX, 1, - 1
            DO JJ = SL%NY, 1, - 1
               WRITE(LU_SCARC, '(16e25.16)') (VALUES((KK-1)*SL%NX*SL%NY+(JJ-1)*SL%NY+II), II=1,SL%NX)
            ENDDO
         ENDDO
   ENDIF

ENDDO

END SUBROUTINE SCARC_DUMP_LEVEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute current revision number
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_DUMP(T)
IMPLICIT NONE
REAL(EB), INTENT(IN) :: T
REAL(EB) :: DP, VORTX(200), VORTZ(200)
INTEGER :: II, JJ, KK, NM, III, KKK
INTEGER ::      MT, MH, MPRES, MUVEL, MVVEL, MWVEL, MVORTX, MVORTZ
CHARACTER(40):: CT, CH, CPRES, CUVEL, CVVEL, CWVEL, CVORTX, CVORTZ
TYPE (MESH_TYPE), POINTER :: M

DO NM = 1, NMESHES

  IF (NM .NE. MYID+1) CYCLE

M=>MESHES(NM)

MT     = GET_FILE_NUMBER()
MH     = GET_FILE_NUMBER()
MPRES  = GET_FILE_NUMBER()
MUVEL  = GET_FILE_NUMBER()
MVVEL  = GET_FILE_NUMBER()
MWVEL  = GET_FILE_NUMBER()
MVORTX = GET_FILE_NUMBER()
!MVORTY = GET_FILE_NUMBER()
MVORTZ = GET_FILE_NUMBER()

WRITE (CT    , '(a,a,a,i2.2)') 'save/time/' , TRIM(CHID), '.',NM
WRITE (CH    , '(a,a,a,i2.2)') 'save/h/'    , TRIM(CHID), '.',NM
WRITE (CPRES , '(a,a,a,i2.2)') 'save/pres/' , TRIM(CHID), '.',NM
WRITE (CUVEL , '(a,a,a,i2.2)') 'save/uvel/' , TRIM(CHID), '.',NM
WRITE (CVVEL , '(a,a,a,i2.2)') 'save/vvel/' , TRIM(CHID), '.',NM
WRITE (CWVEL , '(a,a,a,i2.2)') 'save/wvel/' , TRIM(CHID), '.',NM
WRITE (CVORTX, '(a,a,a,i2.2)') 'save/vortx/', TRIM(CHID), '.',NM
!WRITE (CVORTY, '(a,a,a,i2.2)') 'save/vorty/', TRIM(CHID), '.',NM
WRITE (CVORTZ, '(a,a,a,i2.2)') 'save/vortz/', TRIM(CHID), '.',NM


OPEN(MT    , file=CT    , FORM='FORMATTED', ACTION='READWRITE')
OPEN(MH    , file=CH    , FORM='FORMATTED', ACTION='READWRITE')
OPEN(MPRES , file=CPRES , FORM='FORMATTED', ACTION='READWRITE')
OPEN(MUVEL , file=CUVEL , FORM='FORMATTED', ACTION='READWRITE')
OPEN(MVVEL , file=CVVEL , FORM='FORMATTED', ACTION='READWRITE')
OPEN(MWVEL , file=CWVEL , FORM='FORMATTED', ACTION='READWRITE')
OPEN(MVORTX, file=CVORTX, FORM='FORMATTED', ACTION='READWRITE')
!OPEN(MVORTY, file=CVORTY, FORM='FORMATTED', ACTION='READWRITE')
OPEN(MVORTZ, file=CVORTZ, FORM='FORMATTED', ACTION='READWRITE')

WRITE(MT, 101) T

!WRITE(*,*) 'M%PRESSURE_ZONE:'
!WRITE(*,*) M%PRESSURE_ZONE(1,1,1), M%PRESSURE_ZONE(2,1,1)

DO KK=1,M%KBAR
   DO JJ=1,M%JBAR
      DO II=1,M%IBAR

         DP=PBAR(KK,M%PRESSURE_ZONE(II,JJ,KK)) + M%RHO(II,JJ,KK)*(M%H(II,JJ,KK)-M%KRES(II,JJ,KK)) - M%P_0(KK)
         WRITE(MH    ,102) M%H(II,JJ,KK)
         WRITE(MPRES ,102) DP

         IF (ABS(M%X(II)-M%Z(KK)).LT.1.0E-08) THEN

            ! VORTICITY X
            III = MAX(1,MIN(II,IBAR))
            VORTX(II) = (M%W(III,JJ+1,KK)-M%W(III,JJ,KK))*M%RDYN(JJ) - (M%V(III,JJ,KK+1)-M%V(III,JJ,KK))*M%RDZN(KK)
            WRITE(MVORTX,102) VORTX(II)

            ! VORTICITY Z
            !JJJ = MAX(1,MIN(JJ,JBAR))
            !VORTY(II) = (M%U(II,JJJ,KK+1)-M%U(II,JJJ,KK))*M%RDZN(KK) - (M%W(II+1,JJJ,KK)-M%W(II,JJJ,KK))*M%RDXN(II)
            !WRITE(MVORTY,102) VORTY(II)

            ! VORTICITY Z
            KKK = MAX(1,MIN(KK,KBAR))
            VORTZ(II) = (M%V(II+1,JJ,KKK)-M%V(II,JJ,KKK))*M%RDXN(II) - (M%U(II,JJ+1,KKK)-M%U(II,JJ,KKK))*M%RDYN(JJ)
            WRITE(MVORTZ,102) VORTZ(II)

         ENDIF

    ENDDO
 ENDDO
ENDDO

DO KK=1,M%KBAR
   DO JJ=1,M%JBAR
      DO II=1,M%IBAR+1
         WRITE(MUVEL ,102) M%U(II,JJ,KK)
      ENDDO
   ENDDO
ENDDO
DO KK=1,M%KBAR
   DO JJ=1,M%JBAR+1
      DO II=1,M%IBAR
         WRITE(MVVEL ,102) M%V(II,JJ,KK)
      ENDDO
   ENDDO
ENDDO
DO KK=1,M%KBAR+1
   DO JJ=1,M%JBAR
      DO II=1,M%IBAR
         WRITE(MWVEL ,102) M%W(II,JJ,KK)
      ENDDO
   ENDDO
ENDDO

CLOSE(MT)
CLOSE(MH)
CLOSE(MPRES)
CLOSE(MUVEL)
CLOSE(MVVEL)
CLOSE(MWVEL)
CLOSE(MVORTX)
!CLOSE(MVORTY)
CLOSE(MVORTZ)

ENDDO

101 FORMAT (f25.16)
102 FORMAT (4f25.16)

END SUBROUTINE SCARC_DUMP

!> ------------------------------------------------------------------------------------------------
!> Produce latex information about coarsening of grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_LATEX_GRID(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NL0, MLATEX, I, K, IC, IC1, IC2, IC3, IC4, IC5, IC6
CHARACTER (60) :: CLATEX

IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN

SELECT CASE(NTYPE)

 !! ---------------------------------------------------------------------------------------------
 !! Produce Latex information about grid coarsening with different node sizes for different levels
 !! ---------------------------------------------------------------------------------------------
   CASE(NSCARC_LATEX_STAGGERED)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE

         WRITE (CLATEX, '(A,A,A,i2.2,A,i2.2,A)') 'latex/',TRIM(CHID),'_mesh',NM,'_level',NL+1,'_s.tex'
         MLATEX=GET_FILE_NUMBER()
         OPEN(MLATEX,FILE=CLATEX)
         WRITE(MLATEX,1001)
         WRITE(MLATEX,2001)

         NL0 = NLEVEL_MIN
         DO K = SCARC(NM)%LEVEL(NL0)%NZ,1,-1
            DO I = 1,SCARC(NM)%LEVEL(NL0)%NX
               IC = (K-1)*SCARC(NM)%LEVEL(NL0)%NX + I
               WRITE(MLATEX,3002) 'L1','L1',I,K,I,K
               IC1=SCARC(NM)%LEVEL(NL0)%CELL_TYPE(IC)
               IF (IC1>0) THEN
                   WRITE(MLATEX,3002) 'L2','L2',I,K,I,K
!WRITE(LU_SCARC,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL0)%NCC
                   IF (SCARC(NM)%LEVEL(NL0)%NCC>=1) THEN
                   IF (NL==1) CYCLE
                   IC2=SCARC(NM)%LEVEL(NL0+1)%CELL_TYPE(IC1)
                   IF (IC2>0) THEN
                      WRITE(MLATEX,3002) 'L3','L3',I,K,I,K
!WRITE(LU_SCARC,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL0+1)%NCC
                      IF (SCARC(NM)%LEVEL(NL0+1)%NCC>=1) THEN
                      IF (NL==2) CYCLE
                      IC3=SCARC(NM)%LEVEL(NL0+2)%CELL_TYPE(IC2)
                      IF (IC3>0) THEN
                         WRITE(MLATEX,3002) 'L4','L4',I,K,I,K
!WRITE(LU_SCARC,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL0+2)%NCC
                         IF (SCARC(NM)%LEVEL(NL0+2)%NCC>=1) THEN
                         IF (NL==3) CYCLE
                         IC4=SCARC(NM)%LEVEL(NL0+3)%CELL_TYPE(IC3)
                         IF (IC4>0) THEN
                            WRITE(MLATEX,3002) 'L5','L5',I,K,I,K
!WRITE(LU_SCARC,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL0+3)%NCC
                            IF (SCARC(NM)%LEVEL(NL0+3)%NCC>=1) THEN
                            IF (NL==4) CYCLE
                            IC5=SCARC(NM)%LEVEL(NL0+4)%CELL_TYPE(IC4)
                            IF (IC5>0) THEN
                               WRITE(MLATEX,3002) 'L6','L6',I,K,I,K
!WRITE(LU_SCARC,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL+4)%NCC
                               IF (SCARC(NM)%LEVEL(NL0+4)%NCC>=1) THEN
                               IF (NL==5) CYCLE
                               IC6=SCARC(NM)%LEVEL(NL0+5)%CELL_TYPE(IC4)
                               IF (IC6>0) THEN
                                  WRITE(MLATEX,3002) 'L7','L7',I,K,I,K
                               ENDIF
                               ENDIF
                            ENDIF
                            ENDIF
                         ENDIF
                         ENDIF
                      ENDIF
                      ENDIF
                   ENDIF
                   ENDIF
               ENDIF
            ENDDO
         ENDDO

         WRITE(MLATEX,1002)
         CLOSE(MLATEX)

      ENDDO


 !! ---------------------------------------------------------------------------------------------
 !! Produce Latex information about grid coarsening with same node sizes for different levels
 !! ---------------------------------------------------------------------------------------------
   CASE(NSCARC_LATEX_EQUAL)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE

         WRITE (CLATEX, '(A,A,A,i2.2,A,i2.2,A)') 'latex/',TRIM(CHID),'_mesh',NM,'_level',NL+1,'_e.tex'
         MLATEX=GET_FILE_NUMBER()
         OPEN(MLATEX,FILE=CLATEX)
         WRITE(MLATEX,1001)
         WRITE(MLATEX,3001)

         NL0 = NLEVEL_MIN
         DO K = SCARC(NM)%LEVEL(NL0)%NZ,1,-1
            DO I = 1,SCARC(NM)%LEVEL(NL0)%NX
               IC = (K-1)*SCARC(NM)%LEVEL(NL0)%NX + I
               WRITE(MLATEX,3002) 'L1','L1',I,K,I,K
               IC1=SCARC(NM)%LEVEL(NL0)%CELL_TYPE(IC)
               IF (IC1>0) THEN
                   WRITE(MLATEX,4002) 'L2','L2',I,K,I,K, IC1
                   IF (SCARC(NM)%LEVEL(NL0)%NCC>=1) THEN
                   IF (NL==1) CYCLE
                   IC2=SCARC(NM)%LEVEL(NL0+1)%CELL_TYPE(IC1)
                   IF (IC2>0) THEN
                      WRITE(MLATEX,3002) 'L3','L3',I,K,I,K
                      IF (SCARC(NM)%LEVEL(NL0+1)%NCC>=1) THEN
                      IF (NL==2) CYCLE
                      IC3=SCARC(NM)%LEVEL(NL0+2)%CELL_TYPE(IC2)
                      IF (IC3>0) THEN
                         WRITE(MLATEX,3002) 'L4','L4',I,K,I,K
                         IF (SCARC(NM)%LEVEL(NL0+2)%NCC>=1) THEN
                         IF (NL==3) CYCLE
                         IC4=SCARC(NM)%LEVEL(NL0+3)%CELL_TYPE(IC3)
                         IF (IC4>0) THEN
                            WRITE(MLATEX,3002) 'L5','L5',I,K,I,K
                            IF (SCARC(NM)%LEVEL(NL0+3)%NCC>=1) THEN
                            IF (NL==4) CYCLE
                            IC5=SCARC(NM)%LEVEL(NL0+4)%CELL_TYPE(IC4)
                            IF (IC5>0) THEN
                               WRITE(MLATEX,3002) 'L6','L6',I,K,I,K
                               IF (SCARC(NM)%LEVEL(NL0+3)%NCC>=1) THEN
                               IF (NL==5) CYCLE
                               IC6=SCARC(NM)%LEVEL(NL0+4)%CELL_TYPE(IC4)
                               IF (IC6>0) THEN
                                  WRITE(MLATEX,3002) 'L7','L7',I,K,I,K
                               ENDIF
                               ENDIF
                            ENDIF
                            ENDIF
                         ENDIF
                         ENDIF
                      ENDIF
                      ENDIF
                   ENDIF
                   ENDIF
               ENDIF
            ENDDO
         ENDDO

         WRITE(MLATEX,1002)
         CLOSE(MLATEX)

      ENDDO

 !! ---------------------------------------------------------------------------------------------
 !! Produce Latex information about grid coarsening with different node sizes for different levels
 !! ---------------------------------------------------------------------------------------------
   CASE(NSCARC_LATEX_NUMBER)

      DO NM = 1, NMESHES
         IF (PROCESS(NM) /= MYID) CYCLE

         IF (NL == NLEVEL_MIN) THEN
            WRITE (CLATEX, '(A,A,A,i2.2,A,i2.2,A)') 'latex/',TRIM(CHID),'_mesh',NM,'_level',NL,'_n.tex'
            MLATEX=GET_FILE_NUMBER()
            OPEN(MLATEX,FILE=CLATEX)
            WRITE(MLATEX,1001)
            WRITE(MLATEX,4001)
            NL0 = NLEVEL_MIN
            DO K = SCARC(NM)%LEVEL(NL0)%NZ,1,-1
               DO I = 1,SCARC(NM)%LEVEL(NL0)%NX
                  IC = (K-1)*SCARC(NM)%LEVEL(NL0)%NX + I
!WRITE(LU_SCARC,*) 'IC=',IC,I,K
                  WRITE(MLATEX,4002) 'L1','L1',I,K,I,K, IC
               ENDDO
            ENDDO
            WRITE(MLATEX,1002)
            CLOSE(MLATEX)
         ENDIF

         WRITE (CLATEX, '(A,A,A,i2.2,A,i2.2,A)') 'latex/',TRIM(CHID),'_mesh',NM,'_level',NL+1,'_n.tex'
         MLATEX=GET_FILE_NUMBER()
         OPEN(MLATEX,FILE=CLATEX)
         WRITE(MLATEX,1001)
         WRITE(MLATEX,4001)

         NL0 = NLEVEL_MIN
         DO K = SCARC(NM)%LEVEL(NL0)%NZ,1,-1
            DO I = 1,SCARC(NM)%LEVEL(NL0)%NX
               IC = (K-1)*SCARC(NM)%LEVEL(NL0)%NX + I
               IC1=SCARC(NM)%LEVEL(NL0)%CELL_TYPE(IC)
               IF (IC1> 0) THEN
                  WRITE(MLATEX,4002) 'L2','L2',I,K,I,K, IC1
                  IF (SCARC(NM)%LEVEL(NL0)%NCC>=1) THEN
                  IF (NL==1) CYCLE
                  IC2=SCARC(NM)%LEVEL(NL0+1)%CELL_TYPE(IC1)
                  IF (IC2>0) THEN
                     WRITE(MLATEX,4002) 'L2','L2',I,K,I,K, IC2
                     IF (SCARC(NM)%LEVEL(NL0+1)%NCC>=1) THEN
                     IF (NL==2) CYCLE
                     IC3=SCARC(NM)%LEVEL(NL0+2)%CELL_TYPE(IC2)
                     IF (IC3>0) THEN
                        WRITE(MLATEX,4002) 'L3','L3',I,K,I,K, IC3
                        IF (SCARC(NM)%LEVEL(NL0+2)%NCC>=1) THEN
                        IF (NL==3) CYCLE
                        IC4=SCARC(NM)%LEVEL(NL0+3)%CELL_TYPE(IC3)
                        IF (IC4>0) THEN
                           WRITE(MLATEX,4002) 'L4','L4',I,K,I,K, IC4
                           IF (SCARC(NM)%LEVEL(NL0+3)%NCC>=1) THEN
                           IF (NL==4) CYCLE
                           IC5=SCARC(NM)%LEVEL(NL0+4)%CELL_TYPE(IC4)
                           IF (IC5>0) THEN
                              WRITE(MLATEX,4002) 'L5','L5',I,K,I,K, IC5
                              IF (SCARC(NM)%LEVEL(NL0+3)%NCC>=1) THEN
                              IF (NL==5) CYCLE
                              IC6=SCARC(NM)%LEVEL(NL0+4)%CELL_TYPE(IC4)
                              IF (IC6>0) THEN
                                 WRITE(MLATEX,4002) 'L6','L6',I,K,I,K, IC6
                                 WRITE(MLATEX,4002) 'L7','L7',I,K,I,K, IC6
                              ENDIF
                              ENDIF
                           ENDIF
                           ENDIF
                        ENDIF
                        ENDIF
                     ENDIF
                     ENDIF
                  ENDIF
                  ENDIF
               ELSE
                  WRITE(MLATEX,4002) 'L1','L1',I,K,I,K, IC1
               ENDIF
            ENDDO
         ENDDO

         WRITE(MLATEX,1002)
         CLOSE(MLATEX)

      ENDDO
END SELECT

1001 FORMAT('\documentclass[11pt]{minimal}' ,/,&
            '\usepackage{tikz}' ,/,&
            '\usetikzlibrary{positioning}' ,/,&
            '\begin{document}' ,/,&
            '\begin{tikzpicture}')
1002 FORMAT('\end{tikzpicture}',/,&
            '\end{document}')

2001 FORMAT ('[inner sep=1mm,',/,&
            ' scale=0.9,',/,&
            ' L7/.style={rectangle, draw=black!40, fill=yellow!100  , thin, minimum size=12mm},',/,&
            ' L6/.style={rectangle, draw=black!40, fill=black!80  , thin, minimum size=10mm},',/,&
            ' L5/.style={rectangle, draw=black!40, fill=blue!80 , thin, minimum size=8mm},',/,&
            ' L4/.style={rectangle, draw=black!40, fill=green!50 , thin, minimum size=6mm},',/,&
            ' L3/.style={rectangle, draw=black!40, fill=red!50   , thin, minimum size=4mm},',/,&
            ' L2/.style={rectangle, draw=black!40, fill=blue!40  , thin, minimum size=2mm},',/,&
            ' L1/.style={rectangle, draw=black!40, fill=black!02, thin, minimum size=10mm}],')
2002 FORMAT('\node [',A,']  (',A2,I3.3,I3.3,')  at (',I2,',',I2,')  {};')

3001 FORMAT ('[inner sep=0.1mm,',/,&
            ' scale=0.9,',/,&
            ' L6/.style={rectangle, draw=black!40, fill=black!80  , thin, minimum size=2mm},',/,&
            ' L5/.style={rectangle, draw=black!40, fill=blue!80 , thin, minimum size=3mm},',/,&
            ' L4/.style={rectangle, draw=black!40, fill=green!50 , thin, minimum size=5mm},',/,&
            ' L3/.style={rectangle, draw=black!40, fill=red!50   , thin, minimum size=7mm},',/,&
            ' L2/.style={rectangle, draw=black!40, fill=blue!40  , thin, minimum size=9mm},',/,&
            ' L1/.style={rectangle, draw=black!40, fill=black!02, thin, minimum size=10mm}],')
3002 FORMAT('\node [',A,']  (',A2,I3.3,I3.3,')  at (',I2,',',I2,')  {};')

4001 FORMAT ('[inner sep=1mm,',/,&
            ' scale=0.9,',/,&
            ' L7/.style={rectangle, draw=black!40, fill=yellow!100,thin, minimum size=10mm},',/,&
            ' L6/.style={rectangle, draw=black!40, fill=red!80  , thin, minimum size=10mm},',/,&
            ' L5/.style={rectangle, draw=black!40, fill=black!60 , thin, minimum size=10mm},',/,&
            ' L4/.style={rectangle, draw=black!40, fill=green!50 , thin, minimum size=10mm},',/,&
            ' L3/.style={rectangle, draw=black!40, fill=red!50   , thin, minimum size=10mm},',/,&
            ' L2/.style={rectangle, draw=black!40, fill=blue!40  , thin, minimum size=10mm},',/,&
            ' L1/.style={rectangle, draw=black!40, fill=black!02, thin, minimum size=10mm}],')
4002 FORMAT('\node [',A,']  (',A2,I3.3,I3.3,')  at (',I2,',',I2,')  {',i2,'};')

END SUBROUTINE SCARC_LATEX_GRID



!> ------------------------------------------------------------------------------------------------
!> Compute current revision number
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ANALYTICAL_SOLUTION(NM)
IMPLICIT NONE
!> Initialize flow variables with an analytical solution of the governing equations

INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K, FALL
REAL(EB) :: UU, WW,  RHO_R, RADIUS, VNU, THETA
REAL(EB) :: SHIFT_X, SHIFT_Y, SHIFT_Z, U_MAX, U_C, V_C, W_C, RHO_C, SCALE_RADIUS, CENTER_X, CENTER_Y, CENTER_Z

CALL POINT_TO_MESH(NM)

SELECT CASE(TYPE_CASE)

!> ---------------------------------------------------------------------------------
!> CD_NSA_2D
!> ---------------------------------------------------------------------------------
   CASE(NSCARC_CASE_CD_NSA_2D)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               U(I,J,K) = 1._EB - 2._EB*COS(X(I))*SIN(ZC(K))
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               W(I,J,K) = 1._EB + 2._EB*SIN(XC(I))*COS(Z(K))
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               UU = 1._EB - 2._EB*COS(XC(I))*SIN(ZC(K))
               WW = 1._EB + 2._EB*SIN(XC(I))*COS(ZC(K))
               !H(I,J,K) = -( COS(2._EB*XC(I)) + COS(2._EB*ZC(K)) ) + 0.5_EB*(UU**2+WW**2)
            ENDDO
         ENDDO
      ENDDO

!> ---------------------------------------------------------------------------------
!> VD_NSA_2D
!> Alles noch mal nachpruefen 
!> ---------------------------------------------------------------------------------
   CASE(NSCARC_CASE_VD_NSA_2D)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               U(I,J,K) = 1._EB - 2._EB*COS(X(I))*SIN(ZC(K))
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               W(I,J,K) = 1._EB + 2._EB*SIN(XC(I))*COS(Z(K))
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               UU = 1._EB - 2._EB*COS(XC(I))*SIN(ZC(K))
               WW = 1._EB + 2._EB*SIN(XC(I))*COS(ZC(K))
               !H(I,J,K) = -( COS(2._EB*XC(I)) + COS(2._EB*ZC(K)) ) + 0.5_EB*(UU**2+WW**2)
               RHO(I,J,K) = 1._EB + 10._EB * COS(XC(I))**2 * COS(ZC(K))**2
            ENDDO
         ENDDO
      ENDDO



!> ---------------------------------------------------------------------------------
!> CD_VA_2D
!> Bitte alles noch mal nachpruefen --------
!> ---------------------------------------------------------------------------------
   CASE(NSCARC_CASE_CD_VA_2D)
      WRITE(*,*) '------------------- Case 2 Testfall --------------------------'
      FALL = 1
      IF (FALL==0) THEN
 !! CD_VA_2D Urspruenglicher Testfall
              WRITE(*,*) '------------------- alter Testfall --------------------------'
         U_MAX = 1.0_EB
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=0,IBAR
             SHIFT_X= X(I) - 0.5_EB
                  SHIFT_Y= 0.0_EB
                  SHIFT_Z= ZC(K) - 0.5_EB
                  RADIUS = SQRT( (SHIFT_X)**2 + (SHIFT_Z)**2 )

                  IF (0.0_EB.LE.RADIUS.AND.RADIUS.LT.0.2_EB) THEN
                     VNU = 5.0_EB * RADIUS
                  ELSE IF (RADIUS.LT.0.4_EB) THEN
                     VNU = 2.0_EB - 5.0_EB * RADIUS
                  ELSE IF (RADIUS.GE.0.4_EB) THEN
                     VNU = 0.0_EB
                  ENDIF
                  U(I,J,K) = U_MAX - VNU * SHIFT_Z / RADIUS
               ENDDO
            ENDDO
         ENDDO
         DO K=0,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
             SHIFT_X= XC(I) - 0.5_EB
                  SHIFT_Y= 0.0_EB
                  SHIFT_Z= Z(K) - 0.5_EB
                  RADIUS = SQRT( (SHIFT_X)**2 + (SHIFT_Z)**2 )

                  IF (0.0_EB.LE.RADIUS.AND.RADIUS.LT.0.2_EB) THEN
                     VNU=5.0_EB*RADIUS
                  ELSE IF (RADIUS.LT.0.4_EB) THEN
                     VNU=2.0_EB-5.0_EB*RADIUS
                  ELSE IF (RADIUS.GE.0.4_EB) THEN
                     VNU=0.0_EB
                  ENDIF
                  W(I,J,K) = VNU * SHIFT_X / RADIUS
               ENDDO
            ENDDO
         ENDDO
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  RHO(I,J,K) = 1.0_EB
               ENDDO
            ENDDO
         ENDDO
   ENDIF
   IF (FALL==1) THEN
     !> Traveling vortex mit glatter Verteilung, Gebiet [0,1]x[0,1], Wirbel im Mittelpunkt 0.5,0.5
      WRITE(*,*) '------------------- neuer Testfall --------------------------'
     CENTER_X    = 0.5_EB
          CENTER_Y    = 0.0_EB
       CENTER_Z    = 0.5_EB
     SCALE_RADIUS= 0.4_EB
          RHO_C       = 0.5_EB
          U_C         = 0.0_EB
          V_C         = 0.0_EB

     !> Dichte, zellzentriert
     DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
         SHIFT_X= XC(I) - CENTER_X
              SHIFT_Y=         CENTER_Y
              SHIFT_Z= ZC(K) - CENTER_Z
              RADIUS = SQRT( (SHIFT_X)**2 + (SHIFT_Z)**2 ) / SCALE_RADIUS

                   IF (RADIUS .LT. 1.0_EB) THEN
             RHO(I,J,K) = RHO_C + 0.5_EB*(1-RADIUS*RADIUS)**6
                   ELSE
                 RHO(I,J,K) = RHO_C
         ENDIF
         RHO(I,J,K) = 1.0_EB
             ENDDO
          ENDDO
      ENDDO

     !> u-Geschwindigkeit, flaechenzentriert
     DO K=1,KBAR
        DO J=1,JBAR
           DO I=0,IBAR
         SHIFT_X= X(I)  - CENTER_X
              SHIFT_Y=         CENTER_Y
              SHIFT_Z= ZC(K) - CENTER_Z
              RADIUS = SQRT( (SHIFT_X)**2 + (SHIFT_Z)**2 ) / SCALE_RADIUS
                   THETA  = ATAN( SHIFT_Z / SHIFT_X )

              IF (RADIUS .LT. 1.0_EB) THEN
                 U(I,J,K) = -1024.0_EB * SIN( THETA ) * (1.0_EB - RADIUS)**6 * RADIUS**6 + U_C
                   ELSE
                      U(I,J,K) = U_C
                   ENDIF
           ENDDO
        ENDDO
     ENDDO

     !> v-Geschwindigkeit, flaechenzentriert
     DO K=1,KBAR
        DO J=0,JBAR
           DO I=1,IBAR
                   V(I,J,K) = 0.0_EB
           ENDDO
        ENDDO
     ENDDO

     !> w-Geschwindigkeit, flaechenzentriert
     DO K=0,KBAR
        DO J=1,JBAR
           DO I=1,IBAR
         SHIFT_X= XC(I) - CENTER_X
              SHIFT_Y=         CENTER_Y
              SHIFT_Z= Z(K)  - CENTER_Z
              RADIUS = SQRT( (SHIFT_X)**2 + (SHIFT_Z)**2 ) / SCALE_RADIUS
                   THETA  = ATAN( SHIFT_Z / SHIFT_X )

              IF (RADIUS .LT. 1.0_EB) THEN
                 W(I,J,K) = 1024.0_EB * SIN( THETA ) * (1.0_EB - RADIUS)**6 * RADIUS**6 + W_C
                   ELSE
                      W(I,J,K) = W_C
                   ENDIF
           ENDDO
        ENDDO
     ENDDO
   ENDIF


!> ---------------------------------------------------------------------------------
!> ZM_GRAV_ADVECTED_2D
!> Alles noch mal nachpruefen !!
!> ---------------------------------------------------------------------------------
   CASE(NSCARC_CASE_ZM_GRAV_ADV_2D)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               U(I,J,K) = 1.00_EB
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               W(I,J,K) = 1.25_EB
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               RADIUS=SQRT((XC(I)-0.75_EB)**2+(ZC(K)-0.75_EB)**2)
               IF (0.0_EB.LE.RADIUS.AND.RADIUS.LT.0.5_EB) THEN
                  RHO_R=1.0_EB+0.05_EB*(SIN(RADIUS/0.5_EB*PI+0.5_EB*PI)+1.0_EB)
               ELSE IF (0.5_EB.LE.RADIUS) THEN
                  RHO_R=1.0_EB
               ENDIF
               RHO(I,J,K) = RHO_R
            ENDDO
         ENDDO
      ENDDO
   END SELECT


END SUBROUTINE SCARC_ANALYTICAL_SOLUTION


!> ------------------------------------------------------------------------------------------------
!> Print out timings for ScaRC - not updated at the moment
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_TIMINGS
INTEGER:: IERR, I
INTEGER , ALLOCATABLE, DIMENSION(:)   :: COUNTS_SCARC_TIMERS, DISPLS_SCARC_TIMERS
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: BUFFER
CHARACTER(40) :: NAME_SCARC(0:NSCARC_TIMER)

!> only temporarily - use routine only in debug mode
!IF (TYPE_DEBUG == NSCARC_DEBUG_NONE) RETURN

IERR=0

ALLOCATE(COUNTS_SCARC_TIMERS(0:N_MPI_PROCESSES-1))
ALLOCATE(DISPLS_SCARC_TIMERS(0:N_MPI_PROCESSES-1))
ALLOCATE(BUFFER(NSCARC_TIMER, NMESHES))

COUNTS_SCARC_TIMERS = COUNTS_SCARC*NSCARC_TIMER
DISPLS_SCARC_TIMERS = DISPLS_SCARC*NSCARC_TIMER

BUFFER = TUSED_SCARC
IF (N_MPI_PROCESSES>1) CALL MPI_GATHERV(TUSED_SCARC(1,DISPLS_SCARC(MYID)+1),COUNTS_SCARC_TIMERS(MYID),&
                                        MPI_DOUBLE_PRECISION,TUSED_SCARC, COUNTS_SCARC_TIMERS,DISPLS_SCARC_TIMERS,&
                                        MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

!> Printout subroutine timings
!> outdated version, must be revised (not used at the moment)

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
   NAME_SCARC(NSCARC_TIME_PARDISO_SETUP)= 'ScaRC Pardiso setup'
   NAME_SCARC(NSCARC_TIME_PARDISO)      = 'ScaRC Pardiso method'

   !IF (N_MPI_PROCESSES > 1) THEN
   !  CALL MPI_ALLREDUCE(NCS_BLOCK(MYID+1,NL),TUSED_SCARC_MAX(NL),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
   !  WRITE(LU_OUTPUT,443) NM
   !ELSE
   !ENDIF
   WRITE(LU_OUTPUT,*) 
   DO I = 1, NSCARC_TIMER
      !WRITE(LU_OUTPUT,444) NAME_SCARC(I),TUSED_SCARC(I,MYID+1)
      WRITE(LU_OUTPUT,*) NAME_SCARC(I),TUSED_SCARC(I,MYID+1)
   ENDDO
ENDIF

   IF (TYPE_DEBUG /= NSCARC_DEBUG_NONE ) CLOSE(LU_SCARC)

443 FORMAT(//' ScaRC: CPU Time Usage, Mesh ',I3// &
         47X,' CPU (s)        %  '/ &
         7X,' -----------------------------------------------------------------')
440 FORMAT(7X,A40,1F15.6)
444 FORMAT(7X,A40,1F16.8)

END SUBROUTINE SCARC_TIMINGS


SUBROUTINE SCARC_HANDLE_MPI_ERROR(ERROR_MESSAGE,ERROR_CODE)

INTEGER :: ERROR_CODE, IERR = 0
CHARACTER(*) :: ERROR_MESSAGE

WRITE(LU_ERR,'(A,A,I2)') TRIM(ERROR_MESSAGE),', ERROR_CODE=',ERROR_CODE
CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)

END SUBROUTINE SCARC_HANDLE_MPI_ERROR


!> -----------------------------------------------------------------------------
!> MV: SCARC_SETUP_CCVAR 
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CCVAR(NL)

INTEGER, INTENT(IN) :: NL
INTEGER :: NM
INTEGER :: I, J, K, IX, IY, IZ
INTEGER, PARAMETER :: IMPADD = 1
INTEGER, PARAMETER :: SHFTM(1:3,1:6) = RESHAPE((/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL


MESHES_LOOP1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   SL => SCARC(NM)%LEVEL(NL)

   IF (.NOT. ALLOCATED(SL%CCVAR)) ALLOCATE(SL%CCVAR(0:SL%NX+1,0:SL%NY+1,0:SL%NZ+1,NCVARS))
   SL%CCVAR = 0

   SL%CCVAR(:,:,:,CGSC) = IS_UNDEFINED
   SL%CCVAR(1:SL%NX, 1:SL%NY, 1:SL%NZ, CGSC) = IS_GASPHASE 

ENDDO MESHES_LOOP1


! Now add OBST SOLID cells to CCVAR CGSC array:
MESHES_LOOP2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M  => MESHES(NM)
   SL => SCARC(NM)%LEVEL(NL)

   DO K = 1, SL%NZ
      DO J = 1, SL%NY
         DO I = 1, SL%NX
            IF (M%SOLID(M%CELL_INDEX(I, J, K))) SL%CCVAR(I, J, K, CGSC) = IS_SOLID
         ENDDO
      ENDDO
   ENDDO


<<<<<<< HEAD
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
=======
   IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) '     SET_CCVAR - 1 - ', SL%NX, SL%NY, SL%NZ
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) 'TYPE_DISCRET=',TYPE_DISCRET
      WRITE(LU_SCARC,*) 'SL%NCS=',SL%NCS
<<<<<<< HEAD
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,CGSC)'
      DO IZ = SL%NZ+1,0,-1
      IF (SL%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=SL%NY+1,0,-1)
      ELSE IF (SL%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=SL%NY+1,0,-1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=SL%NY+1,0,-1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,UNKH)'
      DO IZ = SL%NZ+1,0,-1
      IF (SL%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=SL%NY+1,0,-1)
      ELSE IF (SL%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=SL%NY+1,0,-1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=SL%NY+1,0,-1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
=======
      IF (SL%NX == 16) THEN
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(18I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(18I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ELSE IF (SL%NX == 8) THEN
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(10I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(10I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ELSE
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(6I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(6I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ENDIF
>>>>>>> origin/scarc
   ENDIF

ENDDO MESHES_LOOP2

! Define local number of cut-cell:
! Cell numbers for Poisson equation:
MESHES_LOOP3 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   SL => SCARC(NM)%LEVEL(NL)

   IF (ALLOCATED(SL%NUNKH_LOC)) DEALLOCATE(SL%NUNKH_LOC)
   ALLOCATE(SL%NUNKH_LOC(1:NMESHES)) 
   SL%NUNKH_LOC = 0

   SL%CCVAR(:,:,:,UNKH) = IS_UNDEFINED

<<<<<<< HEAD
   IF (PRES_ON_WHOLE_DOMAIN) THEN ! Classic IBM.
=======
   IF (PRES_ON_WHOLE_DOMAIN_SCRC) THEN ! Classic IBM.
>>>>>>> origin/scarc

      ! Loop on Cartesian cells, define cut cells and solid cells CGSC:
      DO K=1,SL%NZ
         DO J=1,SL%NY
            DO I=1,SL%NX
               SL%NUNKH_LOC(NM) = SL%NUNKH_LOC(NM) + 1
               SL%CCVAR(I,J,K,UNKH) = SL%NUNKH_LOC(NM)
            ENDDO
         ENDDO
      ENDDO

<<<<<<< HEAD
=======
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR2(.,.,.,UNKH)'
   WRITE(LU_SCARC,'(6I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
ENDIF

>>>>>>> origin/scarc
   ELSE

      ! Loop on Cartesian cells, define cut cells and solid cells CGSC:
      DO K=1,SL%NZ
         DO J=1,SL%NY
            DO I=1,SL%NX
               IF (SL%CCVAR(I,J,K,CGSC) == IS_GASPHASE ) THEN
                  SL%NUNKH_LOC(NM)     = SL%NUNKH_LOC(NM) + 1
                  SL%CCVAR(I,J,K,UNKH) = SL%NUNKH_LOC(NM)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

   ENDIF 

   SL%NCS = SL%NUNKH_LOC(NM)

<<<<<<< HEAD
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
=======
   IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
>>>>>>> origin/scarc
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) '     SET_CCVAR - 2 - ', SL%NX, SL%NY, SL%NZ
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) 'TYPE_DISCRET=',TYPE_DISCRET
      WRITE(LU_SCARC,*) 'SL%NCS=',SL%NCS
<<<<<<< HEAD
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,CGSC)'
      DO IZ = SL%NZ+1,0,-1
      IF (SL%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ELSE IF (SL%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,UNKH)'
      DO IZ = SL%NZ+1,0,-1
      IF (SL%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ELSE IF (SL%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
=======
      IF (SL%NX == 16) THEN
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(18I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(18I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ELSE IF (SL%NX == 8) THEN
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(10I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(10I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ELSE
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(6I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(6I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ENDIF
>>>>>>> origin/scarc
   ENDIF

ENDDO MESHES_LOOP3

END SUBROUTINE SCARC_SETUP_CCVAR


!> -----------------------------------------------------------------------------
!> MV: Get information about global numbers of unknowns
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CCVAR_GLOBAL(NL)

INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2, IERR
INTEGER :: IX, IY, IZ
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

MESHES_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   SL => SCARC(NM)%LEVEL(NL)

   ! Define total number of unknowns and global unknow index start per MESH:
   IF (ALLOCATED(SL%NUNKH_TOT)) DEALLOCATE(SL%NUNKH_TOT)
   ALLOCATE(SL%NUNKH_TOT(1:NMESHES))
   SL%NUNKH_TOT = 0

   IF (N_MPI_PROCESSES > 1) THEN
      CALL MPI_ALLREDUCE(SL%NUNKH_LOC, SL%NUNKH_TOT, NMESHES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
   ELSE
      SL%NUNKH_TOT = SL%NUNKH_LOC
   ENDIF

   ! Define global start indexes for each mesh:
   IF (ALLOCATED(SL%UNKH_IND)) DEALLOCATE(SL%UNKH_IND)
   ALLOCATE(SL%UNKH_IND(1:NMESHES)) 
   SL%UNKH_IND = 0

   DO NM2=2,NMESHES
      SL%UNKH_IND(NM2) = SL%UNKH_IND(NM2-1) + SL%NUNKH_TOT(NM2-1)
   ENDDO

<<<<<<< HEAD
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) '     SETUP_CCVAR_GLOBAL ', NL
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,CGSC)'
      DO IZ = SL%NZ+1,0,-1
      IF (SL%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ELSE IF (SL%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,UNKH)'
      DO IZ = SL%NZ+1,0,-1
      IF (SL%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ELSE IF (SL%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
=======
   IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) '     SETUP_CCVAR_GLOBAL ', NL
      WRITE(LU_SCARC,*) '============================================================='
      IF (SL%NX == 16) THEN
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(18I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(18I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ELSE IF (SL%NX == 8) THEN
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(10I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(10I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ELSE IF (SL%NX == 4) THEN
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(6I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(6I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ELSE IF (SL%NX == 2) THEN
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,CGSC)'
         WRITE(LU_SCARC,'(4I6)') (((SL%CCVAR(IX, IY, IZ, CGSC), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
         WRITE(LU_SCARC,*) 'NL=',NL,': CCVAR(.,.,.,UNKH)'
         WRITE(LU_SCARC,'(4I6)') (((SL%CCVAR(IX, IY, IZ, UNKH), IX=0, SL%NX+1), IY=0, SL%NY+1), IZ=0, SL%NZ+1)
      ENDIF
>>>>>>> origin/scarc
      WRITE(LU_SCARC,*) 'NUNKH_LOC=', SL%NUNKH_LOC
      WRITE(LU_SCARC,*) 'NUNKH_TOT=', SL%NUNKH_TOT
      WRITE(LU_SCARC,*) 'NCS=', SL%NCS
      !WRITE(77,*) SL%UNKH_IND
   ENDIF

ENDDO MESHES_LOOP

! Achtung: Hier nochmal schauen!
!SL%NUNKH_LOCAL = sum(SL%NUNKH_LOC(1:NMESHES)) ! only nonzeros are for meshes that belong to this process.

END SUBROUTINE SCARC_SETUP_CCVAR_GLOBAL


!> -----------------------------------------------------------------------------
!> MV: Setup CCVAR on coarser levels
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CCVAR_LEVEL(NL)

INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IXF, IYF, IZF, IX, IY, IZ, NSTEP
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLF, SLC

SCARC_ROUTINE = 'SCARC_SETUP_CCVAR_LEVEL'

MAIN_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   SLF => SCARC(NM)%LEVEL(NLEVEL_MIN)
   SLC => SCARC(NM)%LEVEL(NL)

   IF (.NOT. ALLOCATED(SLC%CCVAR)) &
   ALLOCATE(SLC%CCVAR(0:SLC%NX+1,0:SLC%NY+1,0:SLC%NZ+1,NCVARS))
   SLC%CCVAR(:,:,:,CGSC) = IS_UNDEFINED
   SLC%CCVAR(:,:,:,UNKH) = IS_UNDEFINED

   IF (ALLOCATED(SLC%NUNKH_LOC)) DEALLOCATE(SLC%NUNKH_LOC)
   ALLOCATE(SLC%NUNKH_LOC(1:NMESHES)); SLC%NUNKH_LOC = 0

   NSTEP = 2**(NL - NLEVEL_MIN)
<<<<<<< HEAD
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'NSTEP=',NSTEP


   IF (PRES_ON_WHOLE_DOMAIN) THEN
=======
   IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) 'NSTEP=',NSTEP


   IF (PRES_ON_WHOLE_DOMAIN_SCRC) THEN
>>>>>>> origin/scarc
      DO IZ = 1, SLC%NZ
         IZF = (IZ-1)*NSTEP + 1
         DO IY = 1, SLC%NY
            IYF = (IY-1)*NSTEP + 1
            DO IX = 1, SLC%NX
               IXF = (IX-1)*NSTEP + 1
               SLC%CCVAR(IX,IY,IZ,CGSC) = SLF%CCVAR(IXF, IYF, IZF, CGSC)
               SLC%NUNKH_LOC(NM) = SLC%NUNKH_LOC(NM) + 1
               SLC%CCVAR(IX,IY,IZ,UNKH) = SLC%NUNKH_LOC(NM)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO IZ = 1, SLC%NZ
         IZF = (IZ-1)*NSTEP + 1
         DO IY = 1, SLC%NY
            IYF = (IY-1)*NSTEP + 1
            DO IX = 1, SLC%NX
               IXF = (IX-1)*NSTEP + 1
               SLC%CCVAR(IX,IY,IZ,CGSC) = SLF%CCVAR(IXF, IYF, IZF, CGSC)
               IF (SLF%CCVAR(IXF, IYF, IZF, CGSC) == IS_GASPHASE) THEN
                  SLC%NUNKH_LOC(NM) = SLC%NUNKH_LOC(NM) + 1
                  SLC%CCVAR(IX,IY,IZ,UNKH) = SLC%NUNKH_LOC(NM)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   SLC%NCS = SLC%NUNKH_LOC(NM)

<<<<<<< HEAD
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) '     SETUP_CCVAR_LEVEL ', NL, SLC%NCS
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,CGSC)'
      DO IZ = SLC%NZ+1,0,-1
      IF (SLC%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((SLC%CCVAR(IX, IY, IZ, CGSC), IX=0, SLC%NX+1), IY=0, SLC%NY+1)
      ELSE IF (SLC%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((SLC%CCVAR(IX, IY, IZ, CGSC), IX=0, SLC%NX+1), IY=0, SLC%NY+1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((SLC%CCVAR(IX, IY, IZ, CGSC), IX=0, SLC%NX+1), IY=0, SLC%NY+1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,UNKH)'
      DO IZ = SLC%NZ+1,0,-1
      IF (SLC%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((SLC%CCVAR(IX, IY, IZ, UNKH), IX=0, SLC%NX+1), IY=0, SLC%NY+1)
      ELSE IF (SLC%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((SLC%CCVAR(IX, IY, IZ, UNKH), IX=0, SLC%NX+1), IY=0, SLC%NY+1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((SLC%CCVAR(IX, IY, IZ, UNKH), IX=0, SLC%NX+1), IY=0, SLC%NY+1)
      ENDIF
      WRITE(LU_SCARC,*) 
      ENDDO
=======
   IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) '     SETUP_CCVAR_LEVEL ', NL
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,CGSC)'
      WRITE(LU_SCARC,'(6I6)') (((SLC%CCVAR(IX, IY, IZ, CGSC), IX=0, SLC%NX+1), IY=0, SLC%NY+1), IZ=0, SLC%NZ+1)
      WRITE(LU_SCARC,*) 'CCVAR(.,.,.,UNKH)'
      WRITE(LU_SCARC,'(6I6)') (((SLC%CCVAR(IX, IY, IZ, UNKH), IX=0, SLC%NX+1), IY=0, SLC%NY+1), IZ=0, SLC%NZ+1)
>>>>>>> origin/scarc
      WRITE(LU_SCARC,*) 'NUNKH_LOC=', SLC%NUNKH_LOC
      WRITE(LU_SCARC,*) 'NCS=', SLC%NCS
      WRITE(LU_SCARC,*) 'NX =', SLC%NX
      WRITE(LU_SCARC,*) 'NY =', SLC%NY
      WRITE(LU_SCARC,*) 'NZ =', SLC%NZ
   ENDIF


ENDDO MAIN_MESH_LOOP

END SUBROUTINE SCARC_SETUP_CCVAR_LEVEL

END MODULE SCRC
!FORMAT( 3F<2*N+M>.1 )
