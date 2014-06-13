MODULE SCRC  
            
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: SECOND, GET_FILE_NUMBER
USE MPI
USE TYPES, ONLY: MULTIPLIER_TYPE
 
IMPLICIT NONE

CHARACTER(255), PARAMETER :: scrchid='$Id$'
CHARACTER(255), PARAMETER :: scrcrev='$Revision$'
CHARACTER(255), PARAMETER :: scrcdate='$Date$'

PRIVATE

! Public structures (needed in main, read, divg, dump)

!> ------------------------------------------------------------------------------------------------
!> Public subroutines (initialization, solver, time measurement and revisioning)
!> ------------------------------------------------------------------------------------------------
PUBLIC SCARC_SETUP                             !< Setup routine for ScaRC, needed in main.f90
PUBLIC SCARC_SOLVER                            !< Call of basic ScaRC solver
PUBLIC SCARC_TIMINGS                           !< Call of time measurements for ScaRC
PUBLIC GET_REV_SCRC                            !< Get revision number of ScaRC

!> ------------------------------------------------------------------------------------------------
!> Public variables   (explanations in declaration part below)
!> Note: For input parameters in character format corresponding INTEGER type-parameters will
!> Be introduced later to simplify inquiries
!> ------------------------------------------------------------------------------------------------
PUBLIC SCARC_METHOD                            !< Selected ScaRC method
PUBLIC SCARC_INITIAL                           !< Initial solution for ScaRC

PUBLIC SCARC_RESIDUAL                          !< Residual of iterative solver
PUBLIC SCARC_ITERATIONS                        !< Number of iterations
PUBLIC SCARC_CAPPA                             !< Convergence rate
PUBLIC SCARC_ACCURACY                          !< Chosen accuracy type (relative/absolute)
PUBLIC SCARC_ACCURACY_DIVERGENCE               !< Divergence accuracy
PUBLIC SCARC_ACCURACY_RELATIVE                 !< Relative accuracy

PUBLIC SCARC_KRYLOV                            !< Type of Krylov method
PUBLIC SCARC_KRYLOV_ITERATIONS                 !< Maximum number of iterations for Krylov method
PUBLIC SCARC_KRYLOV_ACCURACY                   !< Requested accuracy for Krylov method

PUBLIC SCARC_MULTIGRID                         !< Type of multigrid method
PUBLIC SCARC_MULTIGRID_LEVEL                   !< Multigrid level
PUBLIC SCARC_MULTIGRID_CYCLE                   !< Type of multigrid cycle (V/W/F)
PUBLIC SCARC_MULTIGRID_COARSENING              !< Coarsening method for multigrid (AMG only)
PUBLIC SCARC_MULTIGRID_INTERPOL                !< Interpolation method for multigrid (AMG only)
PUBLIC SCARC_MULTIGRID_ITERATIONS              !< Maximum number of iterations for multigrid method
PUBLIC SCARC_MULTIGRID_ACCURACY                !< Requested accuracy for multigrid method

PUBLIC SCARC_SMOOTH                            !< Smoother for multigrid method
PUBLIC SCARC_SMOOTH_ITERATIONS                 !< Maximum number of iterations for smoothing method
PUBLIC SCARC_SMOOTH_ACCURACY                   !< Requested accuracy for smoothing method
PUBLIC SCARC_SMOOTH_OMEGA                      !< Damping parameter for smoothing method

PUBLIC SCARC_PRECON                            !< Preconditioner of defect correction method
PUBLIC SCARC_PRECON_ITERATIONS                 !< Maximum number of iterations for perconditioning method
PUBLIC SCARC_PRECON_ACCURACY                   !< Requested accuracy for preconditioning method
PUBLIC SCARC_PRECON_OMEGA                      !< Damping parameter for perconditioning method

PUBLIC SCARC_COARSE                            !< Coarse grid solver for multigrid method
PUBLIC SCARC_COARSE_ITERATIONS                 !< Maximum number of iterations for coarse grid solver
PUBLIC SCARC_COARSE_ACCURACY                   !< Requested accuracy for coarse grid solver
PUBLIC SCARC_COARSE_OMEGA                      !< Damping parameter coarse grid solver
PUBLIC SCARC_COARSE_PRECON                     !< Preconditioner for coarse grid solver
 
PUBLIC SCARC_CASE                              !< Predefined ScaRC case
PUBLIC SCARC_DEBUG                             !< Debugging parameter
PUBLIC SCARC_LAYER                             !< Number of layers for data exchange
PUBLIC SCARC_LU                                !< Debug file unit
!PUBLIC SCARC_MKL                              !< MKL code used ?

!> ------------------------------------------------------------------------------------------------
!> Corresponding declarations (with default settings)
!> ------------------------------------------------------------------------------------------------

!! General definitions
CHARACTER(20) :: SCARC_METHOD    = 'null'                   !< Requested solver method (KRYLOV/MULTIGRID)
CHARACTER(20) :: SCARC_INITIAL   = 'null'                   !< Initial solution (currently only default is used)

!! General iteration parameters
REAL (EB)     :: SCARC_RESIDUAL             = -1.0_EB       !< Residual of global selected solver
INTEGER       :: SCARC_ITERATIONS           =  0            !< Number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_CAPPA                =  1.0_EB       !< Convergence rate of selected ScarC solver
REAL (EB)     :: SCARC_ACCURACY_DIVERGENCE  = 1.E+6_EB      !< Divergence epsilon for all solvers
REAL (EB)     :: SCARC_ACCURACY_RELATIVE    = 1.E-2_EB      !< Minimum relative accuracy for all solvers
CHARACTER(20) :: SCARC_ACCURACY             = 'ABSOLUTE'    !< Accuracy type (ABSOLUTE/RELATIVE)

!! Parameters for multigrid-type methods
CHARACTER(20) :: SCARC_MULTIGRID               = 'GEOMETRIC'!! Type of MG-method (GEOMETRIC/ALGEBRAIC)
INTEGER       :: SCARC_MULTIGRID_LEVEL         = -1         !< User defined number of MG-levels (optionally)
CHARACTER(1)  :: SCARC_MULTIGRID_CYCLE         = 'V'        !< Cycling type  (F/V/W)
CHARACTER(20) :: SCARC_MULTIGRID_COARSENING    = 'FALGOUT'  !< Coarsening strategy  (RS3/A1/A2/PMIS/FDS...)
CHARACTER(20) :: SCARC_MULTIGRID_INTERPOL      = 'DIRECT'   !< Interpolation strategy (DIRECT/RS/STANDARD)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS    = 1000       !< Max number of iterations
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY      = 1.E-8_EB  !< Requested accuracy for convergence

!! Parameters for Krylov-type methods
CHARACTER(20) :: SCARC_KRYLOV            = 'CG'             !< Type of Krylov-method (CG/BICG)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 1000             !< Max number of iterations 
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-12_EB        !< Requested accuracy for convergence

!! Parameters for smoothing method (used in multigrids-methods)
CHARACTER(20) :: SCARC_SMOOTH            = 'SSOR'           !< Smoother for MG (JACOBI/SSOR)
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 10               !< Max number of iterations 
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-12_EB        !< Requested accuracy for convergence
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.90E+0_EB       !< Relaxation parameter 

!! Parameters for preconditioning method (used in Krylov-methods)
CHARACTER(20) :: SCARC_PRECON            = 'SSOR'           !< Preconditioner for CG/BICG (JACOBI/SSOR/MG)
INTEGER       :: SCARC_PRECON_ITERATIONS = 1000             !< Max number of iterations 
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-12_EB        !< Requested accuracy for convergence
REAL (EB)     :: SCARC_PRECON_OMEGA      = 0.9E+0_EB        !< Relaxation parameter 

!! Parameters for coarse grid method
CHARACTER(20) :: SCARC_COARSE            = 'ITERATIVE'      !< Coarse grid solver (iterative/direct)
INTEGER       :: SCARC_COARSE_ITERATIONS = 100              !< Max number of iterations for iterative solver
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-12_EB        !< Requested accuracy for iterative solver
REAL (EB)     :: SCARC_COARSE_OMEGA      = 1.5E+0_EB        !< Relaxation parameter for iterative solver
CHARACTER(20) :: SCARC_COARSE_PRECON     = 'SSOR'           !< Preconditioner for iterative solver
 
!! Debugging parameters
CHARACTER(20) :: SCARC_DEBUG = 'NONE'                       !< Debugging level (NONE/LESS/MEDIUM/MUCH)
CHARACTER(40) :: SCARC_FN, SCARC_FN_DUMP                    !< File name for ScaRC debug messages
INTEGER       :: SCARC_LU, SCARC_LU_DUMP                    !< Unit number for ScaRC debug file

CHARACTER(20) :: SCARC_CASE = 'NONE'                        !< Debugging level (NONE/LESS/MEDIUM/MUCH)
INTEGER       :: SCARC_LAYER = 1                            !< Unit number for ScaRC debug file

REAL (EB)     :: SCARC_WALL_TIME = -1.0_EB                  !< Time measurement parameter
REAL (EB)     :: SCARC_WALL_TIME_START = -1.0_EB            !< Initial time for time measurment

!MKL library
!LOGICAL       :: SCARC_MKL = .FALSE.

!! order for the treatment of the single mesh faces
INTEGER :: FACE_ORDER_STD(6) = (/-3,-1,-1,1,2,3/)           !< Standard order of mesh faces
INTEGER :: FACE_ORDER_XYZ(6) = (/1,-1,2,-2,3,-3/)           !< Coordinate direction related order of mesh faces


! Public declarations

!> ------------------------------------------------------------------------------------------------
!> Global constants
!> ------------------------------------------------------------------------------------------------
INTEGER, PARAMETER :: NSCARC_DIMENSION_NONE         = -1, &
                      NSCARC_DIMENSION_TWO          =  2, &    !< two-dimensional problem
                      NSCARC_DIMENSION_THREE        =  3       !< three-dimensional problem

INTEGER, PARAMETER :: NSCARC_SCOPE_NONE             = -1, &         
                      NSCARC_SCOPE_MAIN             =  1, &    !< method used as main solver
                      NSCARC_SCOPE_SMOOTH           =  2, &    !< method used as smoother
                      NSCARC_SCOPE_PRECON           =  3, &    !< method used as preconditiner
                      NSCARC_SCOPE_COARSE           =  4       !< method used as coarse grid solver

INTEGER, PARAMETER :: NSCARC_METHOD_NONE            = -1, &
                      NSCARC_METHOD_KRYLOV          =  1, &    !< Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID       =  2       !< multigrid-method as global solver

INTEGER, PARAMETER :: NSCARC_KRYLOV_NONE            = -1, &
                      NSCARC_KRYLOV_CG              =  1, &    !< CG   as Krylov solver
                      NSCARC_KRYLOV_BICG            =  2       !< BICG as Krylov solver

INTEGER, PARAMETER :: NSCARC_MULTIGRID_NONE         = -1, &
                      NSCARC_MULTIGRID_GEOMETRIC    =  1, &    !< geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC    =  2       !< algebraic multigrid

INTEGER, PARAMETER :: NSCARC_EXCHANGE_NONE          = -1, &
                      NSCARC_EXCHANGE_BASIC         =  1, &    !< initialize wall information
                      NSCARC_EXCHANGE_WALLINFO      =  2, &    !< initialize wall information
                      NSCARC_EXCHANGE_VECTOR        =  3, &    !< matrix-vector communication 
                      NSCARC_EXCHANGE_PRESSURE      =  4, &    !< vector values along internal boundaries
                      NSCARC_EXCHANGE_MEASURE       =  5, &    !< measure values along internal boundaries
                      NSCARC_EXCHANGE_MEASURE_ADD   =  6, &    !< measure values along internal boundaries
                      NSCARC_EXCHANGE_CELL_TYPE     =  7, &    !< cell types along internal boundaries
                      NSCARC_EXCHANGE_CELL_TYPE2    =  8, &    !< cell types II along internal boundaries
                      NSCARC_EXCHANGE_PROLONGATION  =  9, &    !< internal transfer weights
                      NSCARC_EXCHANGE_RESTRICTION   = 10, &    !< internal transfer weights
                      NSCARC_EXCHANGE_STENCIL       = 11, &    !< internal subdiagonal matrix values 
                      NSCARC_EXCHANGE_MATRIX_SIZE   = 12, &    !< neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_SUBDIAG= 13, &    !< neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_STENCIL= 14, &    !< neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_SYSTEM = 15, &    !< neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_PROL   = 16, &    !< neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_REST   = 17, &    !< neighboring matrix size
                      NSCARC_EXCHANGE_TRANSFER_SIZE = 18, &    !< neighboring transfer matrix size
                      NSCARC_EXCHANGE_WIDTHINFO     = 19, &    !< neighboring grid resolution
                      NSCARC_EXCHANGE_MESHINFO      = 20, &    !< neighboring mesh information
                      NSCARC_EXCHANGE_GRAPH         = 21       !< measure values along internal boundaries

INTEGER, PARAMETER :: NSCARC_SMOOTH_NONE            = -1, &
                      NSCARC_SMOOTH_JACOBI          =  1, &    !< smoothing by JACOBI-method
                      NSCARC_SMOOTH_SSOR            =  2, &    !< smoothing by SSOR-method
                      NSCARC_SMOOTH_FFT             =  4       !< smoothing by FFT-method

INTEGER, PARAMETER :: NSCARC_PRECON_NONE            = -1, &
                      NSCARC_PRECON_JACOBI          =  1, &    !< preconditioning by JACOBI-method
                      NSCARC_PRECON_SSOR            =  2, &    !< preconditioning by SSOR-method
                      NSCARC_PRECON_FFT             =  4, &    !< preconditioning by FFT-method
                      NSCARC_PRECON_MULTIGRID       =  5       !< preconditioning by MG-method

INTEGER, PARAMETER :: NSCARC_CYCLE_NONE             = -1, &
                      NSCARC_CYCLE_F                =  0, &    !< F-cycle for mg-method
                      NSCARC_CYCLE_V                =  1, &    !< V-cycle for mg-method
                      NSCARC_CYCLE_W                =  2, &    !< W-cycle for mg-method
                      NSCARC_CYCLE_SETUP            =  3, &    !< initialize cycle counts
                      NSCARC_CYCLE_RESET            =  4, &    !< reset cycle counts
                      NSCARC_CYCLE_PROCEED          =  5, &    !< proceed cycle counts
                      NSCARC_CYCLE_PRESMOOTH        =  6, &    !< presmoothing cycle
                      NSCARC_CYCLE_POSTSMOOTH       =  7, &    !< postsmoothing cycle
                      NSCARC_CYCLE_NEXT             =  8, &    !< perform next cycling loop
                      NSCARC_CYCLE_EXIT             =  9       !< exit cycling loop

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED           =  0, &    !< proceed loop
                      NSCARC_STATE_CONV              =  1, &    !< convergence
                      NSCARC_STATE_DIVG              =  2       !< divergence

INTEGER, PARAMETER :: NSCARC_DEBUG_NONE             = -1, &    !< no debugging requested
                      NSCARC_DEBUG_INFO0            =  0, &    !< info1  level of debugging requested
                      NSCARC_DEBUG_INFO1            =  1, &    !< info1  level of debugging requested
                      NSCARC_DEBUG_INFO2            =  2, &    !< info2  level of debugging requested
                      NSCARC_DEBUG_LESS             =  3, &    !< low    level of debugging requested
                      NSCARC_DEBUG_MEDIUM           =  4, &    !< medium level of debugging requested
                      NSCARC_DEBUG_MUCH             =  5, &    !< strong level of debugging requested
                      NSCARC_DEBUG_MATRIX           =  6, &    !< show matrix
                      NSCARC_DEBUG_MATRIXE          =  7, &    !< show matrix
                      NSCARC_DEBUG_IJKW             =  8, &    !< show IJKW
                      NSCARC_DEBUG_WALLINFO         =  9, &    !< show WALLINFO
                      NSCARC_DEBUG_FACEINFO         = 10, &    !< show FACEINFO
                      NSCARC_DEBUG_BCINDEX          = 11, &    !< show PRESSURE_BC_INDEX
                      NSCARC_DEBUG_ACELL            = 12, &    !< show WALL_CELL
                      NSCARC_DEBUG_GCELL            = 13, &    !< show GHOST_CELL
                      NSCARC_DEBUG_NCELL            = 14, &    !< show NOM_CELL
                      NSCARC_DEBUG_SUBDIVISION      = 15, &    !< show SUBDIVISION
                      NSCARC_DEBUG_MEASURE          = 16, &    !< show MEASURE
                      NSCARC_DEBUG_CELL_TYPE        = 17, &    !< show CELL_TYPE
                      NSCARC_DEBUG_GRAPH            = 18, &    !< show CELL_TYPE
                      NSCARC_DEBUG_COARSE           = 19, &    !< show coarse grid
                      NSCARC_DEBUG_PROLONGATION     = 29, &    !< show prolongation matrix
                      NSCARC_DEBUG_RESTRICTION      = 21       !< show restriction matrix 

INTEGER, PARAMETER :: NSCARC_COARSENING_NONE        = -1, &
                      NSCARC_COARSENING_BASIC       =  1, &    !< basic coarsening
                      NSCARC_COARSENING_FALGOUT     =  2, &    !< parallel Falgout
                      NSCARC_COARSENING_RS3         =  3, &    !< parallel RS3
                      NSCARC_COARSENING_A1          =  4, &    !< aggressive 1 (path=1, length=2)
                      NSCARC_COARSENING_A2          =  5, &    !< aggressive 2 (path=2, length=2)
                      NSCARC_COARSENING_PMIS        =  6, &    !< PMIS 
                      NSCARC_COARSENING_PMISG       =  7, &    !< PMIS 
                      NSCARC_COARSENING_FDSRS3      =  8, &    !< FDSRS3 : FDS variant similar to RS3
                      NSCARC_COARSENING_FDSA1       =  9, &    !< FDSA1  : FDS variant similar to A1
                      NSCARC_COARSENING_FDSA2       = 10, &    !< FDSA2  : FDS variant similar to A2
                      NSCARC_COARSENING_BDRY        = 11, &    !< FDSA2  : FDS variant similar to A2
                      NSCARC_COARSENING_GMG         = 12, &    !< GMG    : GMG as AMG-variang
                      NSCARC_COARSENING_GMG3        = 13       !< GMG3   : for cell numbers divisable by 3

INTEGER, PARAMETER :: NSCARC_COARSE_NONE            = -1, &
                      NSCARC_COARSE_ITERATIVE       =  1, &    !< iterative solution of coarse grid problem
                      NSCARC_COARSE_DIRECT          =  2       !< direct solution of coarse grid problem

INTEGER, PARAMETER :: NSCARC_DIRECT_NONE            = -1, &
                      NSCARC_DIRECT_GE              =  1, &    !< direct solution by Gaussian elmination
                      NSCARC_DIRECT_LU              =  2       !< direct solution by LU-decomposition 

INTEGER, PARAMETER :: NSCARC_SIZE_NONE              = -1, &
                      NSCARC_SIZE_MATRIX            =  2, &    !< size of system matrix for compact system
                      NSCARC_SIZE_TRANSFER          =  3       !< size of transfer matrices for compact system

INTEGER, PARAMETER :: NSCARC_VECTOR_NONE            = -1, &
                      NSCARC_VECTOR_X               =  1, &    !< selection parameter for vector X
                      NSCARC_VECTOR_F               =  2, &    !< selection parameter for vector F
                      NSCARC_VECTOR_Y               =  3, &    !< selection parameter for vector Y
                      NSCARC_VECTOR_G               =  4, &    !< selection parameter for vector G
                      NSCARC_VECTOR_W               =  5, &    !< selection parameter for vector R
                      NSCARC_VECTOR_D               =  6, &    !< selection parameter for vector D
                      NSCARC_VECTOR_Z               =  7, &    !< selection parameter for vector Z
                      NSCARC_VECTOR_X2              =  8, &    !< selection parameter for vector X2
                      NSCARC_VECTOR_G2              =  9, &    !< selection parameter for vector D2
                      NSCARC_VECTOR_D2              = 10, &    !< selection parameter for vector D2
                      NSCARC_VECTOR_W2              = 11, &    !< selection parameter for vector R2
                      NSCARC_VECTOR_Y2              = 12, &    !< selection parameter for vector Y2
                      NSCARC_VECTOR_Z2              = 13, &    !< selection parameter for vector Y2
                      NSCARC_VECTOR_H               = 14, &    !< selection parameter for vector Y2
                      NSCARC_VECTOR_HS              = 15, &    !< selection parameter for vector Y2
                      NSCARC_VECTOR_MEASURE         = 16, &    !< selection parameter for vector MEASURE
                      NSCARC_VECTOR_CELL_TYPE       = 17, &    !< selection parameter for vector CELL_TYPE
                      NSCARC_VECTOR_GRAPH           = 18, &    !< selection parameter for vector GRAPH
                      NSCARC_VECTOR_PTR             = 19       !< selection parameter for vector CELL_TYPE

INTEGER, PARAMETER :: NSCARC_MATRIX_NONE            = -9999999, &
                      NSCARC_MATRIX_SUBDIAG         =  1, &    !< exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SUBDIAG_LOWER   =  2, &    !< exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SUBDIAG_UPPER   =  3, &    !< exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SYSTEM          =  4, &    !< exchange subdiagonal matrix entries
                      NSCARC_MATRIX_TRANSFER        =  5, &    !< exchange prolongation matrix 
                      NSCARC_MATRIX_STENCIL         =  6, &    !< exchange prolongation matrix 
                      NSCARC_MATRIX_PROLONGATION    =  7, &    !< exchange prolongation matrix
                      NSCARC_MATRIX_RESTRICTION     =  8, &    !< exchange restriction matrix 
                      NSCARC_MATRIX_GMG             =  9       !< exchange prolongation matrix

INTEGER, PARAMETER :: NSCARC_ACCURACY_NONE          = -1, &
                      NSCARC_ACCURACY_ABSOLUTE      =  1, &    !< absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE      =  2       !< relative accuracy must be reached

REAL(EB), PARAMETER:: NSCARC_MEASURE_NONE           =  0.0_EB, &
                      NSCARC_MEASURE_ONE            =  1.0_EB, &  !< coarse-grid cell
                      NSCARC_MEASURE_COARSE         =  6.0_EB, &  !< coarse-grid cell
                      NSCARC_MEASURE_FINE           =  5.0_EB, &  !< fine-grid cell
                      NSCARC_MEASURE_SFINE          =  5.0_EB, &  !< strongly coupled fine-grid cell
                      NSCARC_MEASURE_WFINE          =  4.0_EB, &  !< weakly   coupled fine-grid cell
                      NSCARC_MEASURE_BDRY           =  1.0_EB     !< boundry weight

REAL(EB), PARAMETER:: NSCARC_GRAPH_NONE             =  0

INTEGER, PARAMETER :: NSCARC_CELL_TYPE_NONE          =  0, &
                      NSCARC_CELL_TYPE_COARSE        =  1, &   !< coarse-grid cell
                      NSCARC_CELL_TYPE_COARSE0       =  3, &   !< special coarse-grid cell
                      NSCARC_CELL_TYPE_COMMON        =  3, &   !< common cell
                      NSCARC_CELL_TYPE_FINE          = -1, &   !< fine-grid cell
                      NSCARC_CELL_TYPE_FINE0         = -3, &   !< special fine-grid cell
                      NSCARC_CELL_TYPE_SFINE         = -1, &   !< strongly coupled fine-grid cell
                      NSCARC_CELL_TYPE_WFINE         = -2, &   !< weakly   coupled fine-grid cell
                      NSCARC_CELL_TYPE_FPNT          = -1, &   !< special f-point
                      NSCARC_CELL_TYPE_ZPNT          = -2, &   !< special z-point
                      NSCARC_CELL_TYPE_SFPNT         = -3, &   !< special sf-point
                      NSCARC_CELL_TYPE_CPNT          =  2      !< special c-point

INTEGER, PARAMETER :: NSCARC_INTERPOL_NONE          =  0, &
                      NSCARC_INTERPOL_STANDARD      =  1, &    !< standard interpolation
                      NSCARC_INTERPOL_CLASSICAL     =  2, &    !< classical interpolation
                      NSCARC_INTERPOL_CLASSICAL2    =  3, &    !< classical interpolation
                      NSCARC_INTERPOL_DIRECT        =  4, &    !< direct interpolation
                      NSCARC_INTERPOL_DIRECT_BDRY   =  5, &    !< direct interpolation with special boundary
                      NSCARC_INTERPOL_MULTIPASS     =  6, &    !< multipass interpolation
                      NSCARC_INTERPOL_GMG           =  7, &    !< GMG-like interpolation
                      NSCARC_INTERPOL_GMG3          =  8       !< GMG3-like interpolation

INTEGER, PARAMETER :: NSCARC_LATEX_NONE             = -1, &    !< no latex information requested
                      NSCARC_LATEX_STAGGERED        =  1, &    !< show staggered latex information
                      NSCARC_LATEX_EQUAL            =  2, &    !< show equal latex information
                      NSCARC_LATEX_NUMBER           =  3, &    !< show number latex information
                      NSCARC_LATEX_ALL              =  4, &    !< show all latex information
                      NSCARC_LATEX_TABLE            =  5       !< show latex table

INTEGER, PARAMETER :: NSCARC_TIME_NONE              = -1, &
                      NSCARC_TIME_TOTAL             =  1, &    !< time for complete ScaRC part of FDS
                      NSCARC_TIME_SETUP             =  2, &    !< time for setup phase
                      NSCARC_TIME_SOLVER            =  3, &    !< time for ScaRC solver 
                      NSCARC_TIME_KRYLOV            =  4, &    !< time for Krylov solver
                      NSCARC_TIME_MULTIGRID         =  5, &    !< time for multigrid solver
                      NSCARC_TIME_PRECON            =  6, &    !< time for preconditioner
                      NSCARC_TIME_SMOOTH            =  7, &    !< time for smoother
                      NSCARC_TIME_COARSE            =  8, &    !< time for coarse grid solver
                      NSCARC_TIME_MATVEC            =  9, &    !< time for matrix-vector product
                      NSCARC_TIME_SCALPROD          = 10, &    !< time for scalar product
                      NSCARC_TIME_L2NORM            = 11, &    !< time for l2norm
                      NSCARC_TIME_EXCH_INIT         = 12, &    !< time for exchange initialization
                      NSCARC_TIME_EXCH_VECTOR       = 13, &    !< time for exchange of internal boundary
                      NSCARC_TIME_EXCH_MATRIX       = 14, &    !< time for exchange of internal matrix
                      NSCARC_TIME_EXCH_MEASURE      = 15       !< time for exchange of internal measure

INTEGER, PARAMETER :: NSCARC_CASE_NONE              = -1, &    !< no predefined initial solution used
                      NSCARC_CASE_CD_NSA_2D         =  1, &    !< CD_NSA_2D-case
                      NSCARC_CASE_VD_NSA_2D         =  2, &    !< VD_NSA_2D-case
                      NSCARC_CASE_CD_VA_2D          =  3, &    !< CD_VA_2D-case
                      NSCARC_CASE_ZM_GRAV_ADV_2D    =  4       !< ZM_GRAV_ADVECTED_2D

INTEGER, PARAMETER :: NSCARC_IOR_TOP_Z              = -3, &    !< top z-face
                      NSCARC_IOR_BACK_Y             = -2, &    !< back y-face
                      NSCARC_IOR_RIGHT_X            =  1, &    !< right x-face
                      NSCARC_IOR_DUMMY              =  0, &    !< dummy value
                      NSCARC_IOR_LEFT_X             =  1, &    !< left x-face
                      NSCARC_IOR_FRONT_Y            =  2, &    !< front y-face
                      NSCARC_IOR_BOTTOM_Z           =  3       !< bottom z-face

INTEGER, PARAMETER :: NSCARC_LEVEL_NONE             = -1, &    !< no predefined initial solution used
                      NSCARC_LEVEL_MIN              =  0, &    !< minimum multigrid level 
                      NSCARC_LEVEL_MAX              =  15      !< maximum multigrid level 

INTEGER, PARAMETER :: NSCARC_LAYER_NONE             = -1, &    !< no neighbor
                      NSCARC_LAYER_ONE              =  1, &    !< communication of one abutting layer
                      NSCARC_LAYER_TWO              =  2       !< communication of two abutting layers

INTEGER, PARAMETER :: NSCARC_DUMP_NONE              = -1, &    !< no dumping
                      NSCARC_DUMP_RHS               =  1, &    !< dump rhs
                      NSCARC_DUMP_PRES              =  2       !< dump pressure

INTEGER, PARAMETER :: NSCARC_STENCIL_NONE           = -1, &    !< no matrix stencil available
                      NSCARC_STENCIL_CENTRAL        =  1, &    !< standard 5- or 7-point stencil
                      NSCARC_STENCIL_AMG            =  2       !< arbitrary AMG-stencil

INTEGER, PARAMETER :: NSCARC_DIAG_MAIN              =  1, &    !< standard 5- or 7-point stencil
                      NSCARC_DIAG_LOWER             =  2, &    !< standard 5- or 7-point stencil
                      NSCARC_DIAG_UPPER             =  3       !< standard 5- or 7-point stencil

INTEGER, PARAMETER :: NSCARC_INITIAL_NONE           = -1       !< another initial function ?

INTEGER, PARAMETER :: NSCARC_COUPLING_MAX           = 10       !< maximum of possible couplings in stencil
 
INTEGER, PARAMETER :: NSCARC_DUMMY                  = -1       !< dummy variable (needed at several places)

INTEGER, PARAMETER :: NSCARC_DIRECTION_X            =  1, &    !< x-direction indicator
                      NSCARC_DIRECTION_Y            =  2, &    !< y-direction indicator
                      NSCARC_DIRECTION_Z            =  3       !< z-direction indicator

INTEGER, PARAMETER :: NSCARC_NUM_FACES              =  6       !< number of faces per mesh
INTEGER, PARAMETER :: NSCARC_MAX_FACE_NEIGHBORS     = 20       !< max number neighbors per mesh face
INTEGER, PARAMETER :: NSCARC_MAX_MESH_NEIGHBORS     =  6*NSCARC_MAX_FACE_NEIGHBORS

!> --------------------------------------------------------------------------------------------
!> Global variables 
!> --------------------------------------------------------------------------------------------
!> use integer types for the user defined input data (based on SCARC_TYPE_... variables)
INTEGER :: TYPE_DIMENSION  = NSCARC_DIMENSION_NONE        !< Dimension type (2D/3D)
INTEGER :: TYPE_SCOPE      = NSCARC_SCOPE_NONE            !< Type of surrounding solver scopy
INTEGER :: TYPE_SCOPE0     = NSCARC_SCOPE_MAIN            !< Type of surrounding solver scope II
INTEGER :: TYPE_METHOD     = NSCARC_METHOD_NONE           !< Type of global ScaRC method
INTEGER :: TYPE_METHOD0    = NSCARC_METHOD_NONE           !< Type of local ScaRC method
INTEGER :: TYPE_KRYLOV     = NSCARC_KRYLOV_NONE           !< Type of Krylov method
INTEGER :: TYPE_MULTIGRID  = NSCARC_MULTIGRID_NONE        !< Type of multigrid method
INTEGER :: TYPE_MATRIX     = NSCARC_MATRIX_NONE           !< Type of matrix for data exchange
INTEGER :: TYPE_ACCURACY   = NSCARC_ACCURACY_NONE         !< Type of requested accuracy
INTEGER :: TYPE_SMOOTH     = NSCARC_SMOOTH_SSOR           !< Type of smoother for multigrid method
INTEGER :: TYPE_PRECON     = NSCARC_PRECON_SSOR           !< Type of preconditioner for global iterative solver
INTEGER :: TYPE_PRECON0    = NSCARC_PRECON_SSOR           !< Type of preconditioner for local iterative solver
INTEGER :: TYPE_CYCLE      = NSCARC_CYCLE_V               !< Type of cycling for multigrid method
INTEGER :: TYPE_COARSENING = NSCARC_COARSENING_NONE       !< Type of coarsening algorithm for AMG
INTEGER :: TYPE_INTERPOL   = NSCARC_INTERPOL_NONE         !< Type of interpolation for AMG
INTEGER :: TYPE_COARSE     = NSCARC_COARSE_NONE           !< Type of coarse grid solver for multigrid method
INTEGER :: TYPE_CASE       = NSCARC_CASE_NONE             !< Type of predefined test case
INTEGER :: TYPE_DEBUG      = NSCARC_DEBUG_NONE            !< Type of debugging level
INTEGER :: TYPE_INITIAL    = NSCARC_INITIAL_NONE          !< Type of initial solution
INTEGER :: TYPE_EXCHANGE   = NSCARC_EXCHANGE_NONE         !< Type of data exchange
INTEGER :: TYPE_VECTOR     = NSCARC_VECTOR_NONE           !< Type of vector 
INTEGER :: TYPE_DIRECT     = NSCARC_DIRECT_NONE           !< Type of direct solver for multigrid method
INTEGER :: TYPE_LATEX      = NSCARC_LATEX_ALL             !< Type of Latex output level
INTEGER :: TYPE_DUMP       = NSCARC_DUMP_RHS              !< Type of Dump level
INTEGER :: TYPE_LAYER      = NSCARC_LAYER_ONE             !< Type of layers (for overlap)

INTEGER, PARAMETER :: N_TIMERS_SCARC=18                   !< Number of timers for ScaRC
!>
INTEGER :: NMESHES_MIN, NMESHES_MAX                       !< Range of meshes which must be processed for MYID
INTEGER :: NLEVEL, NLEVEL_MAX, NLEVEL_MIN                 !< Total, minimum and maximum number of multigrid levels
INTEGER :: NMASTER, NC_COARSE0                            !< Settings for global coarse solver
INTEGER :: NREQ_SCARC, N_EXCHANGES, TAG_SCARC             !< Variables for data exchange
INTEGER :: SNODE, RNODE                                   !< Process identifier for data exchange
INTEGER :: STATUS2_SCARC(MPI_STATUS_SIZE)                 !< Status array for data excange
!>
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: NC_GLOBAL        !< Global number of cells (over all meshes)
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: NC_LOCAL         !< Local number of cells (related to single meshes)
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: NC_COARSE        !< Global number of cells for coarse solver
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: REQ_SCARC        !< Request array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: COUNTS_SCARC     !< Counter array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: DISPLS_SCARC     !< Displacement array for data exchange
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: SP_LOCAL         !< Local scalar procucts
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TUSED_SCARC      !< Time measurement array

 
!! Private type declarations

!> --------------------------------------------------------------------------------------------
!> Main administration type for ScaRC-method
!> --------------------------------------------------------------------------------------------
TYPE SCARC_TYPE
! 
INTEGER :: NUM_NEIGHBORS = 0                              !< Number of adjacent neighbors of whole mesh
INTEGER, POINTER, DIMENSION(:) :: NEIGHBORS               !< List of adjacent neighbors of whole mesh
! 
INTEGER :: CYCLE_COUNT(2, NSCARC_LEVEL_MAX) = 0           !< Counter for multigrid cycling
REAL (EB), POINTER, DIMENSION (:,:) :: A_COARSE           !< System matrix for coarse grid solver
REAL (EB), POINTER, DIMENSION (:)   :: X_COARSE, X_BUF    !< Solution and RHS vecors for coarse grid solver
INTEGER  , POINTER, DIMENSION (:,:) :: PIVOT              !< Pivot vector for coarse grid solver
! 
INTEGER, POINTER, DIMENSION (:)   :: OFFSET               !< Offset vector for global data exchange
INTEGER, POINTER, DIMENSION (:)   :: COUNTS1, COUNTS2     !< Counters for data exchange
INTEGER, POINTER, DIMENSION (:)   :: DISPLS1, DISPLS2     !< Displacement vectors for data exchange

INTEGER  :: NX, NY, NZ                                    !< number of cells in x-, y- and z-direction on finest level
REAL(EB) :: XS, XF, YS, YF, ZS, ZF                        !< x-, y- and z-bounds of mesh (same for all levels)
! 
TYPE (SCARC_PRECON_TYPE), POINTER, DIMENSION(:) :: PRECON !< Preconditioning type
TYPE (SCARC_LEVEL_TYPE) , POINTER, DIMENSION(:) :: LEVEL  !< Grid level type
TYPE (OSCARC_TYPE)      , POINTER, DIMENSION(:) :: OSCARC !< ScaRC type on other mesh
! 
END TYPE SCARC_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration other mesh data needed for the coupling of adjacent neighbors
!> --------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE
REAL (EB), POINTER, DIMENSION (:) :: SEND_REAL, RECV_REAL       !< main real send and receive buffers 
INTEGER  , POINTER, DIMENSION (:) :: SEND_INTEGER, RECV_INTEGER !< main integer send and receive buffers 
REAL(EB) :: SEND_REAL0(050), RECV_REAL0(50)                     !< initial real send and receive buffers 
INTEGER  :: SEND_INTEGER0(50)  , RECV_INTEGER0(50)              !< initial integer send and receive buffers 
INTEGER  :: NICMAX_R=0, NICMAX_S=0
INTEGER  :: I_MIN_R=-10,I_MAX_R=-10,J_MIN_R=-10,J_MAX_R=-10,K_MIN_R=-10,K_MAX_R=-10,NIC_R=0, &
            I_MIN_S=-10,I_MAX_S=-10,J_MIN_S=-10,J_MAX_S=-10,K_MIN_S=-10,K_MAX_S=-10,NIC_S=0
TYPE (SCARC_LEVEL_TYPE), POINTER, DIMENSION(:) :: LEVEL
END TYPE OSCARC_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of face information related to wall cells and neighbors
!> --------------------------------------------------------------------------------------------
TYPE SCARC_FACE_TYPE
!>
LOGICAL :: BXY, BXZ, BYZ                       !< type of face
!>
INTEGER :: NFM                                 !< own mesh number 
INTEGER :: NFC, NFW                            !< number of cells and wall cells along face
INTEGER :: NFX, NFY, NFZ                       !< local face dimensions 
INTEGER :: NCPL = 1                            !< number of adjacent couplings
INTEGER :: NUM_NEIGHBORS = 0                   !< number of adjacent neighbors
INTEGER :: IWG_MARKER                          !< first IW number (corresponding to global numbering)
INTEGER :: IOFFSET_WALL   = 0                  !< counter for wall cells over all faces
INTEGER :: IOFFSET_MATRIX = 0                  !< counter for wall cells over all faces
!>
INTEGER , POINTER, DIMENSION(:) :: NEIGHBORS   !< adjacent neighbors
REAL(EB), POINTER, DIMENSION(:) :: DH          !< adjacent grid sizes
!>
END TYPE SCARC_FACE_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of wall information related to neighbors and BC's
!> --------------------------------------------------------------------------------------------
TYPE SCARC_WALL_TYPE
!>
!> properties of the single wall cells
!>
INTEGER :: BTYPE                        !< boundary type of wall cell (1:Dirichlet/2:Neumann/3:Internal)
INTEGER :: IOR                          !< orientation of wall cell
INTEGER :: NOM                          !< neighbor at wall cell
INTEGER :: NCPL                         !< neighbor at wall cell
INTEGER :: IWL                          !< corresponding local wall cell number for neighbor NOM
!>
!> Cell and coordinate numberings
!>
INTEGER:: ICW, ICW2                     !< internal wall cell for IW
INTEGER:: ICO, ICO2                     !< overlapping cell for IW
INTEGER:: IXG, IYG, IZG                 !< x-, y- and z-indices of ghost cells
INTEGER:: IXW, IYW, IZW                 !< x-, y- and z-indices of (internal) wall cells
INTEGER:: IXN(2), IYN(2), IZN(2)        !< x-, y- and z-indices of neighboring cells 
!>
!> pointer vecors to ghost, extended and neighboring cells
!>
INTEGER, POINTER, DIMENSION(:) :: ICE, ICE2         !< extended cell for IW
INTEGER, POINTER, DIMENSION(:) :: ICG, ICG2         !< ghost cell for IW
INTEGER, POINTER, DIMENSION(:) :: ICN, ICN2         !< neighboring cell for IW
!>
END TYPE SCARC_WALL_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of data related to single grid level
!> --------------------------------------------------------------------------------------------
TYPE SCARC_LEVEL_TYPE
!>
INTEGER :: SUBDIVISION(3,-3:3)=0
INTEGER :: NLEN_MATRIX_SYSTEM = 0
!>
!> Different numbers, lengths and settings
!>
INTEGER :: NX, NY, NZ           !< number of grid cells in x-, y- and z-direction
INTEGER :: NW, NWL              !< number of global and local wall cells
INTEGER :: NC                   !< number of cells
INTEGER :: NCG=0, NCG0=0        !< number of ghost cells 
INTEGER :: NCE=0, NCE0=0        !< number of extended cells plus auxiliary variable
INTEGER :: NCO=0, NCO0=0        !< number of overlapping cells plus auxiliary variable
INTEGER :: NA, NAE, NA0         !< number of elements in internal and extended system matrix
INTEGER :: NP, NPE, NP0         !< number of elements in internal and extended prolongation matrix
INTEGER :: NR, NRE, NR0         !< number of elements in internal and extended restriction matrix
INTEGER :: NLAYER               !< number of layers 
INTEGER :: NCPL, NCPL_MAX=-10   !< number of couplings
INTEGER :: NCPLS, NCPLR         !< number of couplings to send and read
INTEGER :: NCW                  !< number of coarse wall cells
INTEGER :: NCF, NCFE, NCF0      !< number of fine cells in internal and extended grid (AMG)
INTEGER :: NCC, NCCE, NCC0      !< number of coarse cells in internal and extended grid (AMG)
INTEGER :: NCCI                 !< number of internal coarse cells
INTEGER :: NOM                  !< number of adjacent neighbor
INTEGER :: NPOINTS              !< number of points in matrix stencil
!>
!> Different pointers
!>
INTEGER :: IWL = 0              !< local wall cell numbers
INTEGER :: ICG = 0              !< ghost cell counter
INTEGER :: ICG0 = 0             !< auxiliary ghost cell counter
!>
INTEGER :: ICG_MARKER = 0       !< ghost cell marker
INTEGER :: ICO_MARKER = 0       !< overlapping cell marker
INTEGER :: ICE_MARKER = 0       !< extended cell marker
INTEGER :: IWL_MARKER = 0       !< local wall cell marker
INTEGER :: IWG_MARKER = 0       !< global wall cell marker
!>
!> Coordinate data
!>
REAL(EB) :: DI2(3)              !< sample vector of step sizes in x-, y- and z-direction
REAL(EB) :: DX , DY , DZ        !< step sizes in x-, y- and z-direction
REAL(EB) :: DXI, DYI, DZI       !< inversed of step sizes in x-, y- and z-direction
REAL(EB) :: DXI2, DYI2, DZI2    !< squared and inversed step sizes in x-, y- and z-direction
!>
REAL (EB), POINTER, DIMENSION (:)  :: DXL, DYL, DZL        !< step size vectors in x-, y- and z-direction
REAL (EB), POINTER, DIMENSION (:)  :: XCOR, YCOR, ZCOR     !< coordinate vectors in x-, y- and z-direction
REAL (EB), POINTER, DIMENSION (:)  :: XMID, YMID, ZMID     !< midpoint vectors in x-, y- and z-direction
!>
!> matrices with corresponding storage vectors
!>
REAL (EB), POINTER, DIMENSION (:) :: A , P , R, S     !< system, prolongation, restriction and strength matrix
INTEGER,   POINTER, DIMENSION (:) :: A_ROW , A_COL    !< row and column pointers for system matrix A
INTEGER,   POINTER, DIMENSION (:) :: P_ROW , P_COL    !< row and column pointers for prolongation matrix P
INTEGER,   POINTER, DIMENSION (:) :: R_ROW , R_COL    !< row and column pointers for restriction matrix A
INTEGER,   POINTER, DIMENSION (:) :: S_ROW , S_COL    !< row and column pointers for strength matrix S
INTEGER,   POINTER, DIMENSION (:) :: ST_ROW, ST_COL   !< row and column pointers for transpose of strength matrix S
INTEGER,   POINTER, DIMENSION (:) :: A_TAG , P_TAG    !< auxiliary arrays for mark positions in A and P
INTEGER,   POINTER, DIMENSION (:) :: A_SIZE, AUX      !< some more auxiliary arrays
!>
!> Different vector definitions for global and local solvers
!>
REAL (EB), POINTER, DIMENSION (:) :: X , X2           !< solution vectors
REAL (EB), POINTER, DIMENSION (:) :: F , F2           !< right hand side vectors
REAL (EB), POINTER, DIMENSION (:) :: D , D2           !< defect vectors
REAL (EB), POINTER, DIMENSION (:) :: Y , Y2           !< auxiliary vectors
REAL (EB), POINTER, DIMENSION (:) :: G , G2           !< auxiliary vectors
REAL (EB), POINTER, DIMENSION (:) :: W , W2           !< auxiliary vectors
REAL (EB), POINTER, DIMENSION (:) :: Z , Z2           !< auxiliary vectors
!>
!> Different pointer vectors for wall and communication related data
!>
INTEGER, POINTER, DIMENSION(:,:):: WALL_INDEX         !< wall index 
INTEGER, POINTER, DIMENSION(:)  :: INTERNAL_BDRY_CELL !! index of internal boundary cells
INTEGER, POINTER, DIMENSION(:)  :: EXT_PTR            !< index of external cells
!>
!> neighbor related settings for IOR, NCPL and NWL
!>
REAL(EB):: DH = 0.0_EB                                !< local step sizes
INTEGER :: IOR = 0                                    !< local orientations

!>
!> Different mappings of quantities
!>
INTEGER, POINTER, DIMENSION (:)  :: ICE_TO_IWG        !< mapping from ICO to IWG
INTEGER, POINTER, DIMENSION (:)  :: ICE_TO_IWL        !< mapping from ICO to IWL
INTEGER, POINTER, DIMENSION (:)  :: ICE_TO_ICW        !< mapping from ICO to ICW
INTEGER, POINTER, DIMENSION (:)  :: ICE_TO_ICG        !< mapping from ICO to ICG
INTEGER, POINTER, DIMENSION (:)  :: ICE_TO_ICN        !< mapping from ICO to ICN
!>
INTEGER, POINTER, DIMENSION (:)  :: ICG_TO_IWG        !< mapping from ICG to IWG
INTEGER, POINTER, DIMENSION (:)  :: ICG_TO_ICW        !< mapping from ICG to ICW
INTEGER, POINTER, DIMENSION (:)  :: ICG_TO_ICE        !< mapping from ICG to ICE
INTEGER, POINTER, DIMENSION (:)  :: ICG_TO_ICN        !< mapping from ICG to ICN
!>
INTEGER, POINTER, DIMENSION (:)  :: ICN_TO_IWG        !< mapping from ICN to IWG
INTEGER, POINTER, DIMENSION (:)  :: ICN_TO_ICW        !< mapping from ICN to ICW
INTEGER, POINTER, DIMENSION (:)  :: ICN_TO_ICE        !< mapping from ICN to ICE
INTEGER, POINTER, DIMENSION (:)  :: ICN_TO_ICG        !< mapping from ICN to ICG
!>
INTEGER, POINTER, DIMENSION (:)  :: IWL_TO_IWG        !< mapping from IWL to IWG
INTEGER, POINTER, DIMENSION (:)  :: IWL_TO_ICW        !< mapping from IWL to ICW
INTEGER, POINTER, DIMENSION (:)  :: IWL_TO_ICO        !< mapping from IWL to ICO
INTEGER, POINTER, DIMENSION (:,:):: IWL_TO_ICG        !< mapping from IWL to ICG
INTEGER, POINTER, DIMENSION (:,:):: IWL_TO_ICN        !< mapping from IWL to ICN
!>
!> special administration vectors/arrays for Algebraic multigrid (AMG)
!>
REAL(EB):: MAX_ROWSUM=0.9_EB
REAL(EB), POINTER, DIMENSION (:)  :: MEASURE          !< mapping for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:,:):: CELL_MAP         !< cell mapping for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:)  :: CELL_TYPE        !< cell types for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:)  :: CELL_MARKER      !< cell markers for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:)  :: GRAPH            !< graph vector for AMG-coarsening
INTEGER,  POINTER, DIMENSION (:,:):: P_PTR            !< pointer vector for AMG-prolongation
!>
!> local type definitions
!>
TYPE (SCARC_WALL_TYPE), POINTER, DIMENSION(:) :: WALL !! description of single wall cells
TYPE (SCARC_FACE_TYPE), POINTER, DIMENSION(:) :: FACE !! description of single faces
!>
END TYPE SCARC_LEVEL_TYPE


!> --------------------------------------------------------------------------------------------
!> Administration of scopes for solution methods
!> --------------------------------------------------------------------------------------------
TYPE SCARC_SCOPE_TYPE
INTEGER :: X, F, Y, G, W, D, Z              !< local references of solver vectors
INTEGER :: ITE, NIT                         !< local references to iteration counters
REAL(EB) :: EPS, RES, RESIN, OMEGA          !< local references to solver parameters
CHARACTER(30) :: CROUTINE = 'null'           
END TYPE SCARC_SCOPE_TYPE


!> --------------------------------------------------------------------------------------------
!> Save parent settings when starting new solution method
!> --------------------------------------------------------------------------------------------
TYPE SCARC_PARENT_TYPE
INTEGER :: TYPE_METHOD           !< type of solution method
INTEGER :: TYPE_SCOPE            !< tyoe of related scope
INTEGER :: TYPE_PRECON           !< type of preconditioning method
INTEGER :: TYPE_SMOOTH           !< type of smoothing
INTEGER :: TYPE_ACCURACY         !< type of accuracy requirements (relative/absolute)
INTEGER :: TYPE_CYCLE            !< type of multigrid cycle
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
REAL(EB):: TNOW_SETUP

!> ------------------------------------------------------------------------------------------------
!> Initialize time measurement 
!> ------------------------------------------------------------------------------------------------
IERR = 0
TNOW_SETUP = SECOND()

ALLOCATE(TUSED_SCARC(0:N_TIMERS_SCARC,NMESHES),STAT=IERR)
CALL ChkMemErr('SCARC_SETUP','TUSED_SCARC',IERR)

TUSED_SCARC = 0._EB
TUSED_SCARC(NSCARC_TIME_TOTAL,:) = SECOND()


!> ------------------------------------------------------------------------------------------------
!> Parse input parameters
!> ------------------------------------------------------------------------------------------------
CALL SCARC_INPUT_PARSER                  !< parse input parameters (corresponding to FDS-file)


!> ------------------------------------------------------------------------------------------------
!> Setup different components of ScaRC 
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_DIMENSION               !< define dimension of underlying problem
CALL SCARC_SETUP_DEBUGGING               !< open debug file if requested
CALL SCARC_SETUP_PROCESSES               !< determine set of meshes which must be processed on MYID
CALL SCARC_SETUP_LEVELS                  !< define number of necessary grid levels 
CALL SCARC_SETUP_STRUCTURES              !< allocate requested ScaRC-types for all necessary grid levels
CALL SCARC_SETUP_MESHES                  !< set mesh information
CALL SCARC_SETUP_TYPES                   !< allocate wall and face types for all requested grid levels
CALL SCARC_SETUP_WALLINFO                !< set wall cell information
CALL SCARC_SETUP_EXCHANGE                !< set information for data exchange
CALL SCARC_SETUP_SYSTEM                  !< assemble system matrix with boundary conditions and solver vectors
CALL SCARC_SETUP_COARSENING              !< perform coarsening on different grid levels if requested (AMG only)
CALL SCARC_SETUP_VECTORS                 !< allocate solution and auxiliary vectors on all needed grid levels
CALL SCARC_SETUP_GLOBAL_DIMENSIONS       !< define some global variables
!CALL SCARC_SETUP_NEIGHBORS              !< compute information about abutting neighbors on coarser levels


!> ------------------------------------------------------------------------------------------------
!> Measure time for setup routine
!> ------------------------------------------------------------------------------------------------
TUSED_SCARC(NSCARC_TIME_SETUP,:)=TUSED_SCARC(NSCARC_TIME_SETUP,:)+SECOND()-TNOW_SETUP
TUSED_SCARC(NSCARC_TIME_TOTAL,:)=TUSED_SCARC(NSCARC_TIME_TOTAL,:)+SECOND()-TNOW_SETUP

END SUBROUTINE SCARC_SETUP


!> ----------------------------------------------------------------------------------------------------
!> Define dimension of problem
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DIMENSION
IF (TWO_D) THEN
   TYPE_DIMENSION = NSCARC_DIMENSION_TWO
ELSE
   TYPE_DIMENSION = NSCARC_DIMENSION_THREE
ENDIF
END SUBROUTINE SCARC_SETUP_DIMENSION


!> ----------------------------------------------------------------------------------------------------
!> Determine types of input parameters
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INPUT_PARSER
CHARACTER(300):: CMESSAGE

!> ------------------------------------------------------------------------------------------------
!> set type of global solver
!> ------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_METHOD))

   CASE ('KRYLOV')           

      TYPE_METHOD  = NSCARC_METHOD_KRYLOV
      TYPE_METHOD0 = NSCARC_METHOD_KRYLOV

      !< set type of Krylov-method (CG/BICG)
      SELECT CASE (TRIM(SCARC_KRYLOV))
         CASE ('CG')
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
         CASE ('BICG')
            TYPE_KRYLOV = NSCARC_KRYLOV_BICG
         CASE DEFAULT
            WRITE(CMESSAGE,1002) 'Krylov-solver',TRIM(SCARC_KRYLOV),'ScaRC','CG','BICG'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT 

      !< set type of preconditioner (JACOBI/SSOR/MG)
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
         CASE DEFAULT
            WRITE(CMESSAGE,1005) 'preconditioner',TRIM(SCARC_PRECON),&
                                 'Krylov-method','JACOBI','SSOR','MG','FFT'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT
      TYPE_PRECON0 = TYPE_PRECON

   CASE ('MULTIGRID')

      TYPE_METHOD  = NSCARC_METHOD_MULTIGRID
      TYPE_METHOD0 = NSCARC_METHOD_MULTIGRID

      !< set type of multigrid method (GEOMETRIC/ALGEBRAIC)
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

      !< set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_COARSE))
         CASE ('ITERATIVE')
            TYPE_COARSE = NSCARC_COARSE_ITERATIVE
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
         CASE ('DIRECT')
            TYPE_COARSE = NSCARC_COARSE_DIRECT
            TYPE_DIRECT = NSCARC_DIRECT_GE
      END SELECT 

      !< set type of smoother (JACOBI/SSOR)
      SELECT CASE (TRIM(SCARC_SMOOTH))                        !< use same parameters as for preconditioner
         CASE ('JACOBI')
            TYPE_PRECON = NSCARC_PRECON_JACOBI
            TYPE_SMOOTH = NSCARC_SMOOTH_JACOBI
         CASE ('SSOR')
            TYPE_PRECON = NSCARC_PRECON_SSOR
            TYPE_SMOOTH = NSCARC_SMOOTH_SSOR
         CASE DEFAULT
            WRITE(CMESSAGE,1003) 'smoother',TRIM(SCARC_SMOOTH),'multigrid','JACOBI','SSOR'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT

   CASE DEFAULT
      WRITE(CMESSAGE,1002) 'solver',TRIM(SCARC_METHOD),'ScaRC','KRYLOV','MULTIGRID'
      CALL SCARC_SHUTDOWN(CMESSAGE)

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

         !< set type of coarsening strategy (STANDARD/AGGRESSIVE)
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
               WRITE(CMESSAGE,1005) 'coarsening',TRIM(SCARC_MULTIGRID_COARSENING),&
                                    'algebraic multigrid','RS3'   ,'A1'   ,'A2'   ,'PMIS',&
                                                          'PMISG','TEST','FDSA2','FDSRS3'
               CALL SCARC_SHUTDOWN(CMESSAGE)
         END SELECT

      CASE DEFAULT
         WRITE(CMESSAGE,1002) 'multigrid',TRIM(SCARC_MULTIGRID),&
                              'Krylov-preconditioner','GEOMETRIC','ALGEBRAIC'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

   !! set type of coarse grid solver (CG/GE)
   SELECT CASE (TRIM(SCARC_COARSE))
      CASE ('ITERATIVE')
         TYPE_COARSE = NSCARC_COARSE_ITERATIVE
         TYPE_KRYLOV = NSCARC_KRYLOV_CG
      CASE ('DIRECT')
         TYPE_COARSE = NSCARC_COARSE_DIRECT
         TYPE_DIRECT = NSCARC_DIRECT_GE
      CASE DEFAULT
         WRITE(CMESSAGE,1002) 'coarse grid solver',TRIM(SCARC_COARSE),'multigrid','ITERATIVE','DIRECT'
         CALL SCARC_SHUTDOWN(CMESSAGE)
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
         WRITE(CMESSAGE,1003) 'cycling ',TRIM(SCARC_MULTIGRID_CYCLE),'multigrid','F','V','W'
         CALL SCARC_SHUTDOWN(CMESSAGE)
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
         WRITE(CMESSAGE,1004) 'cycling ',TRIM(SCARC_MULTIGRID_INTERPOL), &
                              'multigrid','STANDARD','DIRECT','DIRECT_BDRY','MULTIPASS'
         CALL SCARC_SHUTDOWN(CMESSAGE)
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
      WRITE(CMESSAGE,1002) 'accuracy',TRIM(SCARC_ACCURACY),'ScaRC','ABSOLUTE','RELATIVE'
      CALL SCARC_SHUTDOWN(CMESSAGE)
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
   CASE DEFAULT
      WRITE(CMESSAGE,1006) 'debugging',TRIM(SCARC_DEBUG),'ScaRC','NONE','INFO0','INFO1','INFO2','LESS','MEDIUM','MUCH'
      CALL SCARC_SHUTDOWN(CMESSAGE)
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
1006 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
1007 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
!1008 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
!>            'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
!>            'ScaRC: Aborting program ...')
2002 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, ' in case of ',A,'=',A,/, &
             'ScaRC: Redefining ',A,' to ',A,' !',/)
END SUBROUTINE SCARC_INPUT_PARSER


!> ------------------------------------------------------------------------------------------------
!> Setup debug file if requested
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DEBUGGING
INTEGER:: NM

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
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


!> ------------------------------------------------------------------------------------------------
!> Serial or parallel version ? Decide which meshes must be processed
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PROCESSES
IF (USE_MPI) THEN
   NMESHES_MIN = MYID+1
   NMESHES_MAX = MYID+1
ELSE
   NMESHES_MIN = 1
   NMESHES_MAX = NMESHES
ENDIF
END SUBROUTINE SCARC_SETUP_PROCESSES

 
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
      IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
         NLEVEL_MAX  = NLEVEL + 1
      ELSE IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
         NLEVEL_MAX  = NLEVEL_MIN + SCARC_MULTIGRID_LEVEL - 1
      ELSE
         NLEVEL_MAX  = NLEVEL
      ENDIF
      NLEVEL = NLEVEL_MAX

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
   STOP
ENDIF

!> Divide by 2 as often as possible or till user defined max-level is reached
NC0=NC
DO NL=1,NSCARC_LEVEL_MAX
   NC0=NC0/2
   IF (MOD(NC0,2)/=0) EXIT                !< NC no longer divisable by two
   IF (NL==SCARC_MULTIGRID_LEVEL) EXIT    !< max number of levels defined by user
   IF (NC0==1) EXIT                       !< NC is power of two, minimum has been reached
ENDDO

SCARC_GET_MAXLEVEL=NL
RETURN
1000 FORMAT(A,'=',I3,' must be divisable by two for ScaRC-Multigrid!')
END FUNCTION SCARC_GET_MAXLEVEL


!> ------------------------------------------------------------------------------------------------
!> Allocate ScaRC-structures for all needed levels
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STRUCTURES
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
MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   S => SCARC(NM)

WRITE(*,*) NM, ': ALLOCATING OSCARC(',NMESHES,')'
   ALLOCATE (S%OSCARC(NMESHES), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_STRUCTURES', 'OSCARC', IERR)
      
WRITE(*,*) NM, ': ALLOCATING LEVEL(',NLEVEL_MIN,':',NLEVEL_MAX,')'
   ALLOCATE (SCARC(NM)%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_STRUCTURES', 'LEVEL', IERR)

   IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
      ALLOCATE (SCARC(NM)%PRECON(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_STRUCTURES', 'PRECON', IERR)
   ENDIF

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_STRUCTURES


!> ----------------------------------------------------------------------------------------------------
!> Setup geometry information for mesh NM
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MESHES
INTEGER  :: IERR, NL, NM, IX, IY, IZ
INTEGER  :: NX0, NY0, NZ0
TYPE (MESH_TYPE)       , POINTER :: M
TYPE (SCARC_TYPE)      , POINTER :: S
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

IERR=0

LEVEL_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   M => MESHES(NM)
   S => SCARC(NM)
   
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      DO IX=1,M%N_EXTERNAL_WALL_CELLS
         WRITE(SCARC_LU,'(a,7i4)') '0:  IX=',IX, M%WALL(IX)%NOM_IB(1:6)
      ENDDO
   ENDIF

   !<  define hierarchy of meshes depending on the chosen method
   NX0=MESHES(NM)%IBAR
   NY0=MESHES(NM)%JBAR
   NZ0=MESHES(NM)%KBAR
   
   !< store bounds of mesh in SCARC-structure
   S%XS = M%XS                          !< x-bound left
   S%XF = M%XF                          !< x-bound right
   S%YS = M%YS                          !< y-bound front
   S%YF = M%YF                          !< y-bound back
   S%ZS = M%ZS                          !< z-bound bottom
   S%ZF = M%ZF                          !< z-bound top

   S%NX = M%IBAR
   S%NY = M%JBAR
   S%NZ = M%KBAR

   LEVEL_LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX
   
      IF (NL > NLEVEL_MIN .AND. TYPE_MULTIGRID /= NSCARC_MULTIGRID_GEOMETRIC) EXIT LEVEL_LEVEL_LOOP

      !< let SM point to SCARC(NM)%MESHES(NL)
      SL => S%LEVEL(NL)
      
      !< numbers of cells in x-, y- and z-direction for level 'NL'
      SL%NX = NX0
      SL%NY = NY0
      SL%NZ = NZ0
      
      SL%NC = SL%NX * SL%NY * SL%NZ 
      SL%NW = 2*SL%NX*SL%NY + 2*SL%NX*SL%NZ + 2*SL%NY*SL%NZ

      !< step widths in x-, y- and z-direction for level 'NL'
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

      NX0=NX0/2
      IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) NY0=NY0/2
      NZ0=NZ0/2

      !< Allocate vectors for coordinate information
      IF (NL == NLEVEL_MIN) THEN

         SL%XCOR => M%X
         SL%YCOR => M%Y
         SL%ZCOR => M%Z

         SL%XMID => M%XC
         SL%YMID => M%YC
         SL%ZMID => M%ZC

      ELSE

         !< compute coordinates in x-, y- and z-direction
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

         !< compute midpoints in x-, y- and z-direction
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

      !< Allocate vectors for step sizes in different directions 
      ALLOCATE(SL%DXL(0:SL%NX), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','DXL',IERR)
      SL%DXL = 0.0_EB
      ALLOCATE(SL%DYL(0:SL%NY), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','DYL',IERR)
      SL%DYL = 0.0_EB
      ALLOCATE(SL%DZL(0:SL%NZ), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','DZL',IERR)
      SL%DZL = 0.0_EB

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

      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,*) '============ SETUP_MESHES ==============='
         WRITE(SCARC_LU,*) 'NX     =',SL%NX
         WRITE(SCARC_LU,*) 'NY     =',SL%NY
         WRITE(SCARC_LU,*) 'NZ     =',SL%NZ
         WRITE(SCARC_LU,*) 'NC     =',SL%NC
         WRITE(SCARC_LU,*) 'NW     =',SL%NW
         WRITE(SCARC_LU,*) 'XS     =',S%XS
         WRITE(SCARC_LU,*) 'XF     =',S%XF
         WRITE(SCARC_LU,*) 'YS     =',S%YS
         WRITE(SCARC_LU,*) 'YF     =',S%YF
         WRITE(SCARC_LU,*) 'ZS     =',S%ZS
         WRITE(SCARC_LU,*) 'ZF     =',S%ZF
         WRITE(SCARC_LU,*) 'DX     =',SL%DX
         WRITE(SCARC_LU,*) 'DY     =',SL%DY
         WRITE(SCARC_LU,*) 'DZ     =',SL%DZ
         WRITE(SCARC_LU,*) 'DXI    =',SL%DXI
         WRITE(SCARC_LU,*) 'DYI    =',SL%DYI
         WRITE(SCARC_LU,*) 'DZI    =',SL%DZI
         WRITE(SCARC_LU,*) 'DXI2   =',SL%DXI2
         WRITE(SCARC_LU,*) 'DYI2   =',SL%DYI2
         WRITE(SCARC_LU,*) 'DZI2   =',SL%DZI2
         WRITE(SCARC_LU,*) 'DI2(1) =',SL%DI2(1)
         WRITE(SCARC_LU,*) 'DI2(2) =',SL%DI2(2)
         WRITE(SCARC_LU,*) 'DI2(3) =',SL%DI2(3)
         WRITE(SCARC_LU,'(a,20f8.3)') 'XCOR=',SL%XCOR(0:SL%NX)
         WRITE(SCARC_LU,'(a,20f8.3)') 'YCOR=',SL%YCOR(0:SL%NY)
         WRITE(SCARC_LU,'(a,20f8.3)') 'ZCOR=',SL%ZCOR(0:SL%NZ)
         WRITE(SCARC_LU,'(a,20f8.3)') 'XMID=',SL%XMID(0:SL%NX+1)
         WRITE(SCARC_LU,'(a,20f8.3)') 'YMID=',SL%YMID(0:SL%NY+1)
         WRITE(SCARC_LU,'(a,20f8.3)') 'ZMID=',SL%ZMID(0:SL%NZ+1)
         WRITE(SCARC_LU,'(a,20f8.3)') 'DXL =',SL%DXL(0:SL%NX)
         WRITE(SCARC_LU,'(a,20f8.3)') 'DYL =',SL%DYL(0:SL%NY)
         WRITE(SCARC_LU,'(a,20f8.3)') 'DZL =',SL%DZL(0:SL%NZ)
      ENDIF
   ENDDO LEVEL_LEVEL_LOOP

ENDDO LEVEL_MESHES_LOOP


END SUBROUTINE SCARC_SETUP_MESHES
 

!> ------------------------------------------------------------------------------------------------
!> Setup communication structure for data exchange 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TYPES 
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
MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   
   M => MESHES(NM)
   S => SCARC(NM)
   
   SLF => S%LEVEL(NLEVEL_MIN)

   !! Allocate basic wall and face type for finest level
   ALLOCATE(SLF%WALL(SLF%NW), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_TYPES','WALL',IERR)

   ALLOCATE(SLF%FACE(-3:3), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_TYPES','FACE',IERR)


   !> Allocate basic wall and face type for coarser levels if requested
   IF (NLEVEL_MAX > NLEVEL_MIN .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN                   

      DO NL=NLEVEL_MIN+1,NLEVEL_MAX
  
         SLC => S%LEVEL(NL)           

         ALLOCATE(SLC%WALL(SLC%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_TYPES','WALL',IERR)

         ALLOCATE(SLC%FACE(-3:3), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_TYPES','WALL',IERR)
      ENDDO

   ENDIF

   !> 
   !> Allocate basic WALL and FACE types on neighboring meshes NOM for all requested grid levels
   !> 
   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES
      
      IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
   
      OM => M%OMESH(NOM)
      OS => S%OSCARC(NOM)
   
      OS%I_MIN_S  = OM%I_MIN_S
      OS%I_MAX_S  = OM%I_MAX_S
      OS%J_MIN_S  = OM%J_MIN_S
      OS%J_MAX_S  = OM%J_MAX_S
      OS%K_MIN_S  = OM%K_MIN_S
      OS%K_MAX_S  = OM%K_MAX_S
       
      OS%I_MIN_R  = OM%I_MIN_R
      OS%I_MAX_R  = OM%I_MAX_R
      OS%J_MIN_R  = OM%J_MIN_R
      OS%J_MAX_R  = OM%J_MAX_R
      OS%K_MIN_R  = OM%K_MIN_R
      OS%K_MAX_R  = OM%K_MAX_R

      OS%NIC_S    = OM%NIC_S
      OS%NICMAX_S = OM%NIC_S

      OS%NIC_R    = OM%NIC_R
      OS%NICMAX_R = OM%NIC_R
   
      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0)  CYCLE OTHER_MESHES_LOOP

      N_EXCHANGES  = N_EXCHANGES+1
   
   ENDDO OTHER_MESHES_LOOP
ENDDO MESHES_LOOP


!> ------------------------------------------------------------------------------------------------
!> Initialize level structures on neighboring meshes
!> ------------------------------------------------------------------------------------------------
LEVEL_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   
   M => MESHES(NM)
   S => SCARC(NM)
   
   LEVEL_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES
      
      IF (NOM == NM) CYCLE LEVEL_OTHER_MESHES_LOOP

      OS => S%OSCARC(NOM)
      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0)  CYCLE LEVEL_OTHER_MESHES_LOOP

      !> Allocate OSCARC level structure for mesh NM
WRITE(*,*) NM,': ALLOCATING SCARC(',NM,')%OSCARC(',NOM,')%LEVEL(',NLEVEL_MIN,':',NLEVEL_MAX,')'

      ALLOCATE (S%OSCARC(NOM)%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OS%LEVEL', IERR)
   
      !> point to grid structure of OSCARC(NM) on finest level
      OSLF => S%OSCARC(NOM)%LEVEL(NLEVEL_MIN)

      OSLF%NX = MESHES(NOM)%IBAR 
      OSLF%NY = MESHES(NOM)%JBAR 
      OSLF%NZ = MESHES(NOM)%KBAR 

      OSLF%NW  = 2*OSLF%NX*OSLF%NY + 2*OSLF%NX*OSLF%NZ + 2*OSLF%NY*OSLF%NZ  
      OSLF%NC  = OSLF%NX*OSLF%NY*OSLF%NZ 
      OSLF%NA  = 0
      OSLF%NCG = 0

      IF (S%OSCARC(NOM)%NICMAX_S == 0 .AND. S%OSCARC(NOM)%NICMAX_R == 0) CYCLE LEVEL_OTHER_MESHES_LOOP

      !< In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN                   

         DO NL=NLEVEL_MIN+1,NLEVEL_MAX
   
            OSLC => S%OSCARC(NOM)%LEVEL(NL)                          !< pointer to coarser level
            OSLF => S%OSCARC(NOM)%LEVEL(NL-1)                        !< pointer to finer level
   
            !< get number of internal cells and external wall cells on neighbor NOM for level NL
            OSLC%NX = OSLF%NX/2
            SELECT CASE (TYPE_DIMENSION)
               CASE (NSCARC_DIMENSION_TWO)
                  OSLC%NY = 1
               CASE (NSCARC_DIMENSION_THREE)
                  OSLC%NY = OSLF%NY/2
            END SELECT
            OSLC%NZ = OSLF%NZ/2

            OSLC%NC  = OSLC%NX * OSLC%NY * OSLC%NZ
            OSLC%NW  = 2*OSLC%NX*OSLC%NY + 2*OSLC%NX*OSLC%NZ + 2*OSLC%NY*OSLC%NZ  
            OSLC%NA  = 0
            OSLC%NCG = 0
   
         ENDDO
      ENDIF
         
   ENDDO LEVEL_OTHER_MESHES_LOOP
ENDDO LEVEL_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_TYPES
   

!> ------------------------------------------------------------------------------------------------
!> Initialize arrays for data exchange 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE 
INTEGER :: NL, NM, NOM, NP, NLEN_SEND, NLEN_RECV
INTEGER :: NLMIN, NLMAX
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
!> Allocate counter and displacement vector for global data exchanges
!>
ALLOCATE(COUNTS_SCARC(0:NUMPROCS-1), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'COUNTS_SCARC', IERR)

ALLOCATE(DISPLS_SCARC(0:NUMPROCS-1), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'COUNTS_SCARC', IERR)

COUNTS_SCARC = 0
DO NP = 0, NUMPROCS-1
   DO NM = 1, NMESHES
      IF (PROCESS(NM)==NP) COUNTS_SCARC(NP) = COUNTS_SCARC(NP) + 1
   ENDDO
ENDDO

DISPLS_SCARC(0) = 0
DO NP = 1, NUMPROCS-1
   DISPLS_SCARC(NP) = COUNTS_SCARC(NP-1) + DISPLS_SCARC(NP-1)
ENDDO


!>
!> Exchange basic information about wall sizes (needed for the dimensioning of !the exchange buffers)
!>
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_BASIC , NLEVEL_MIN)


!>
!> Allocate send and receive buffers (real and integer) in correct lengths
!>
DO NM = NMESHES_MIN, NMESHES_MAX
   S => SCARC(NM)
   DO INBR = 1, S%NUM_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)

      SL  => SCARC(NM)%LEVEL(NLEVEL_MIN)
      OS  => SCARC(NM)%OSCARC(NOM)               
      OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)

      !< -------------------------- real buffers
      !NLEN_SEND = MAX(OS%NIC_R, OS%NIC_S)*TYPE_LAYER*SL%NCPL_MAX*2+10
      NLEN_SEND = MAX(OSL%NWL, OSL%NCG)+10
      NLEN_RECV = NLEN_SEND

      ALLOCATE (OS%SEND_REAL(NLEN_SEND), STAT = IERR)
      CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'SEND_REAL', IERR)
      OS%SEND_REAL = 0.0_EB

      ALLOCATE (OS%RECV_REAL(NLEN_RECV), STAT = IERR)
      CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'RECV_REAL', IERR)
      OS%RECV_REAL = 0.0_EB

      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,*) 'NOM=',NOM,': ALLOCATING SEND-VECTORS in length ',NLEN_SEND
         WRITE(SCARC_LU,*) 'NOM=',NOM,': ALLOCATING RECV-VECTORS in length ',NLEN_RECV
         WRITE(SCARC_LU,*) 'NOM=',NOM,': OSL%NWL=',OSL%NWL
         WRITE(SCARC_LU,*) 'NOM=',NOM,': OSL%NCG=',OSL%NCG
         WRITE(SCARC_LU,*) 'NOM=',NOM,': SL%NCPL_MAX=',SL%NCPL_MAX
      ENDIF

      !< -------------------------- integer buffers
      OSL%NCPLS = OSL%NCPL
      IF (OSL%NCG == OSL%NWL) THEN
         NLEN_SEND = 14*OSL%NWL 
         NLEN_RECV = 14*OSL%NCG 
         OSL%NCPLR = OSL%NCPL
      ELSE IF (OSL%NCG == 2*OSL%NWL) THEN
         NLEN_SEND = 12*OSL%NWL + 2*OSL%NCG 
         NLEN_RECV = 14*OSL%NCG 
         OSL%NCPLR =  1
      ELSE IF (OSL%NWL == 2*OSL%NCG) THEN
         NLEN_SEND = 14*OSL%NWL 
         NLEN_RECV = 12*OSL%NCG + 2*OSL%NWL
         OSL%NCPLR =  2
      ENDIF

      ALLOCATE (OS%SEND_INTEGER(NLEN_SEND), STAT = IERR)
      CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'SEND_INTEGER', IERR)
      OS%SEND_INTEGER = 0
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'NOM=',NOM,': ALLOCATING SEND-INT in length ',NLEN_SEND

      ALLOCATE (OS%RECV_INTEGER(NLEN_RECV), STAT = IERR)
      CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'RECV_INTEGER', IERR)
      OS%RECV_INTEGER = 0

      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'NOM=',NOM,': ALLOCATING RECV-INT in length ',NLEN_RECV
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'NOM=',NOM,': ALLOCATING WALL in length ',OSL%NCG

      !< -------------------- neighboring wall (and face?) structures for common wall cells
      ALLOCATE (OSL%WALL(OSL%NCG), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSL%WALL', IERR)

      !ALLOCATE (OSL%FACE(-3:3), STAT=IERR)
      !CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSL%FACE', IERR)

      !< In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN                   
         DO NL=NLEVEL_MIN+1,NLEVEL_MAX
   
            OSL => S%OSCARC(NOM)%LEVEL(NL)                          !< pointer to coarser level

            ALLOCATE (OSL%WALL(OSL%NCG), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSL%WALL', IERR)

            !ALLOCATE (OSL%FACE(-3:3), STAT=IERR)
            !CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSL%FACE', IERR)

         ENDDO
      ENDIF

   ENDDO
ENDDO

!> 
!> Initialize communication structures on finest level (if there is more than 1 mesh) 
!> 
IF (NMESHES>1) THEN

   NLMIN = NLEVEL_MIN
   NLMAX = NLEVEL_MIN
   IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) NLMAX = NLEVEL_MAX

   DO NL = NLMIN, NLMAX
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WIDTHINFO, NL)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WALLINFO, NL)
   ENDDO

ENDIF


END SUBROUTINE SCARC_SETUP_EXCHANGE


!> ----------------------------------------------------------------------------------------------------
!> Find the cell index corresponding to cell (II,JJ,KK) 
!> ----------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CELL_INDEX(II, JJ, KK, NM, NL)
INTEGER, INTENT(IN) :: II, JJ, KK, NM, NL
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

SCARC_CELL_INDEX = (KK-1)*SL%NX*SL%NY + (JJ-1)*SL%NX + II

END FUNCTION SCARC_CELL_INDEX


!> ----------------------------------------------------------------------------------------------------
!> Setup neighborship structures and boundary conditions 
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLINFO
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES
INTEGER :: NM, NL
INTEGER :: IWC, IREFINE, IERR, IFACE
INTEGER :: INBR_FACE, INBR_MESH
INTEGER :: IWG, IWL
INTEGER :: NOM_LAST , NOM
INTEGER :: NCPL_LAST, NCPL
INTEGER :: IOR0_LAST, IOR0
INTEGER :: MLATEX
INTEGER :: FACE_NEIGHBORS(NSCARC_MAX_FACE_NEIGHBORS, -3:3) 
INTEGER :: MESH_NEIGHBORS(NSCARC_MAX_MESH_NEIGHBORS)       
INTEGER :: NUM_FACE_NEIGHBORS(-3:3) 
INTEGER :: NUM_MESH_NEIGHBORS       
CHARACTER (40) :: CLATEX
TYPE (MESH_TYPE)        , POINTER :: M
TYPE (SCARC_TYPE)       , POINTER :: S
TYPE (OSCARC_TYPE)      , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SLF , SLC
TYPE (SCARC_LEVEL_TYPE) , POINTER :: OSLF, OSLC

IERR=0


MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   M => MESHES(NM)
   S => SCARC(NM)
      
   FACE_NEIGHBORS = -1
   MESH_NEIGHBORS = -1

   NUM_FACE_NEIGHBORS = 0
   NUM_MESH_NEIGHBORS = 0

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      DO IWG=1,M%N_EXTERNAL_WALL_CELLS
         WRITE(SCARC_LU,'(a,7i4)') 'A:  IWG=',IWG, M%WALL(IWG)%NOM_IB(1:6)
      ENDDO
   ENDIF

   !! For all solvers: 
   !! Determine array WALL and PRESSURE_BC_INDEX on finest level SLF
   SLF => S%LEVEL(NLEVEL_MIN)
   
   SLF%NCE = SLF%NC
   SLF%NCO = SLF%NC
   SLF%NLAYER = 1

   !> -----------------------------------------------------------------------
   !> Determine basic data for single faces (orientation, dimensions)
   !> -----------------------------------------------------------------------
   FACES_OF_MESH_LOOP1: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE FACES_OF_MESH_LOOP1
      SLF%FACE(IOR0)%NFM = NM

      !< information about face orientation and local dimensions
      SELECT CASE (ABS(IOR0))
         CASE (1)
            SLF%FACE(IOR0)%BXY  =  .FALSE.                 !< face in x-y-direction ?
            SLF%FACE(IOR0)%BXZ  =  .FALSE.                 !< face in x-z-direction ?
            SLF%FACE(IOR0)%BYZ  =  .TRUE.                  !< face in y-z-direction ?
            SLF%FACE(IOR0)%NFC  =  SLF%NX                  !< number of cells between opposite mesh faces
            SLF%FACE(IOR0)%NFX  =  1                       !< number of cells in x-direction
            SLF%FACE(IOR0)%NFY  =  SLF%NY                  !< number of cells in y-direction
            SLF%FACE(IOR0)%NFZ  =  SLF%NZ                  !< number of cells in z-direction
            SLF%FACE(IOR0)%NFW  =  SLF%NY*SLF%NZ           !< number of wall cells at that face
            SLF%FACE(IOR0)%DH   => SLF%DXL                 !< step size vector between opposite mesh faces
            SLF%FACE(IOR0)%IOFFSET_MATRIX = 1              !< related offset for matrix diagonals
         CASE (2)
            SLF%FACE(IOR0)%BXY  =  .FALSE.                 !< same as above ...
            SLF%FACE(IOR0)%BXZ  =  .TRUE.           
            SLF%FACE(IOR0)%BYZ  =  .FALSE.
            SLF%FACE(IOR0)%NFC  =  SLF%NY
            SLF%FACE(IOR0)%NFX  =  SLF%NX
            SLF%FACE(IOR0)%NFY  =  1
            SLF%FACE(IOR0)%NFZ  =  SLF%NZ
            SLF%FACE(IOR0)%NFW  =  SLF%NX*SLF%NZ
            SLF%FACE(IOR0)%DH   => SLF%DYL
            SLF%FACE(IOR0)%IOFFSET_MATRIX = SLF%NX
         CASE (3)
            SLF%FACE(IOR0)%BXY  =  .TRUE.                  !< same as above ...
            SLF%FACE(IOR0)%BXZ  =  .FALSE.
            SLF%FACE(IOR0)%BYZ  =  .FALSE.
            SLF%FACE(IOR0)%NFC  =  SLF%NZ
            SLF%FACE(IOR0)%NFX  =  SLF%NX
            SLF%FACE(IOR0)%NFY  =  SLF%NY
            SLF%FACE(IOR0)%NFZ  =  1
            SLF%FACE(IOR0)%NFW  =  SLF%NX*SLF%NY
            SLF%FACE(IOR0)%DH   => SLF%DZL
            SLF%FACE(IOR0)%IOFFSET_MATRIX = SLF%NX * SLF%NY
      END SELECT
      IF (IOR0 > 0) SLF%FACE(IOR0)%IOFFSET_MATRIX = -SLF%FACE(IOR0)%IOFFSET_MATRIX      !< take negative offset for IOR>0
   
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
   !> store basic data and determine numberof adajacent neighbors to each
   !> face with corresponding number of IW's
   !> -----------------------------------------------------------------------
   IOR0_LAST =  0
   NOM_LAST  = -1
   NCPL_LAST = -1
   IWL = 0
   GLOBAL_WALLCELLS_LOOP1: DO IWG = 1, SLF%NW

      !< Determine and store neighbors, orientation and number of couplings for a single wall cell
      NOM  =  M%WALL(IWG)%NOM
      IOR0 =  M%WALL(IWG)%ONE_D%IOR
      NCPL = (M%WALL(IWG)%NOM_IB(4) - M%WALL(IWG)%NOM_IB(1) + 1) * &
             (M%WALL(IWG)%NOM_IB(5) - M%WALL(IWG)%NOM_IB(2) + 1) * &
             (M%WALL(IWG)%NOM_IB(6) - M%WALL(IWG)%NOM_IB(3) + 1) 

      SLF%WALL(IWG)%NOM  = NOM                              !< store data in SCARC-structure
      SLF%WALL(IWG)%IOR  = IOR0
      SLF%WALL(IWG)%NCPL = NCPL

      IF (IOR0_LAST /= IOR0) IWL = 0
      IWL = IWL + 1                                         !< count local wall cells for that face

      !< for wall cells with a neighbor build extended and overlapping structures
      IF (NOM /= 0) THEN

         SLF%NCE = SLF%NCE + NCPL                           !< increase number of extended grid cells
         SLF%NCO = SLF%NCO + 1                              !< increase number of overlapping grid cells

         OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)       
         OSLF%NCPL = NCPL                                   !< initialize own counter for local wall cells
         OSLF%IOR  = -IOR0                                  !< initialize own orientation variable

         IF (NOM_LAST  /= NOM ) THEN                        !< new neighbor ?????
            OSLF%NWL = 1                                    !< initialize own counter for local wall cells
            OSLF%NCG = NCPL                                 !< initialize counter for local ghost cells
            SLF%NCPL_MAX  = MAX(SLF%NCPL_MAX, NCPL)         !< get max NCPL ever used on this mesh

            IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
               WRITE(SCARC_LU,'(10(a,i4))') &
                     'A: IWG=',IWG,': NOM=',NOM,': OSLF%NWL=', OSLF%NWL, ': SLF%NCO=', SLF%NCO, &
                     ': SLF%NCE=',SLF%NCE,': IWL=',IWL, ': IOR0=',IOR0,': NCPL=',OSLF%NCPL
            ENDIF
         ELSE                                               !< already known neighbor ???
            OSLF%NWL = OSLF%NWL + 1                         !< increase own counter for local wall cells
            OSLF%NCG = OSLF%NCG + NCPL                      !< increase counter for local ghost cells

            IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
               WRITE(SCARC_LU,'(10(a,i4))') &
                     'B: IWG=',IWG,': NOM=',NOM,': OSLF%NWL=', OSLF%NWL, ': SLF%NCO=', SLF%NCO, &
                     ': SLF%NCE=',SLF%NCE,': IWL=',IWL, ': IOR0=',IOR0,': NCPL=',OSLF%NCPL
            ENDIF
         ENDIF
      ENDIF

      IF (NOM /= 0 .AND. NOM_LAST /= NOM) THEN                            !< new neighbor at that face?

         NUM_FACE_NEIGHBORS(IOR0) = NUM_FACE_NEIGHBORS(IOR0) + 1          !< increase neighbor counter
         NUM_MESH_NEIGHBORS       = NUM_MESH_NEIGHBORS       + 1          !< increase neighbor counter

         FACE_NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0), IOR0) = NOM             !< store number of neighbor
         MESH_NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0))       = NOM             !< store number of neighbor

WRITE(*,'(5(a,i5))') 'NM=',NM,': NOM=',NOM,':IOR0=',IOR0,': #NBR=',NUM_FACE_NEIGHBORS(IOR0),': NBR=',NOM
      ENDIF

      IOR0_LAST = IOR0                                                    !< save former values
      NCPL_LAST = NCPL
      NOM_LAST  = NOM

   ENDDO GLOBAL_WALLCELLS_LOOP1


   !! allocate array which stores numbers of all neighboring meshes
   IF (NUM_MESH_NEIGHBORS /= 0) THEN
      ALLOCATE(S%NEIGHBORS(NUM_MESH_NEIGHBORS), STAT=IERR)   
      CALL ChkMemErr('SCARC_SETUP_WALLINFO','MESH_NEIGHBORS',IERR)
   ENDIF
   S%NUM_NEIGHBORS = NUM_MESH_NEIGHBORS


   !> -----------------------------------------------------------------------
   !> Store information about adjacent neighbors on different faces
   !> Allocate corresponding indexing arrays in OSCARC-structures
   !> -----------------------------------------------------------------------
   !! allocate administrative mapping arrays for own mesh 
   CALL SCARC_SETUP_MAPPINGS(NM, NLEVEL_MIN)

   INBR_MESH = 1
   NEIGHBORS_OF_FACE_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE NEIGHBORS_OF_FACE_LOOP

      !< if there are neighbors at face IOR0 store information about them
      IF (NUM_FACE_NEIGHBORS(IOR0) /= 0) THEN

         SLF%FACE(IOR0)%NUM_NEIGHBORS = NUM_FACE_NEIGHBORS(IOR0)        !< store number of neighbors on face

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(a,i4,a,i4,a,i4)') 'IOR0=',IOR0,': NUM_FACE_NEIGHBORS=',NUM_FACE_NEIGHBORS(IOR0)

         !< store numbers of neighbors
         ALLOCATE(SLF%FACE(IOR0)%NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0)), STAT=IERR)   
         CALL ChkMemErr('SCARC_SETUP_WALLINFO','NOMS',IERR)

         !< store every neighbor and allocate corresponding administration arrays 
         DO INBR_FACE = 1, NUM_FACE_NEIGHBORS(IOR0)

            NOM = FACE_NEIGHBORS(INBR_FACE,IOR0)

            SLF%FACE(IOR0)%NEIGHBORS(INBR_FACE) = NOM                   !< store NOM as a neighbor of that face
            S%NEIGHBORS(INBR_MESH)              = NOM                   !< store NOM as a neighbor of that mesh                
            INBR_MESH = INBR_MESH + 1

WRITE(*,'(4(a,i4))') 'NM=',NM,': IOR0=',IOR0,': NUM_FACE_NBRS=',NUM_FACE_NEIGHBORS(IOR0), &
           ': NBR=',SLF%FACE(IOR0)%NEIGHBORS(INBR_FACE)

            !! allocate administrative arrays for neighboring meshes
            CALL SCARC_SETUP_OMAPPINGS(NM, NOM, NLEVEL_MIN)

         ENDDO

      ENDIF

   ENDDO NEIGHBORS_OF_FACE_LOOP

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      DO IWG=1,SLF%NW
         WRITE(SCARC_LU,'(a,7i4)') 'B: IWG=',IWG, M%WALL(IWG)%NOM_IB(1:6)
      ENDDO
   ENDIF

   !> -----------------------------------------------------------------------
   !> Second loop over global IW's:
   !> store detailed coordinate and cell data
   !> and determine kind of boundary condition
   !> -----------------------------------------------------------------------
   SLF%ICE_MARKER = SLF%NC
   SLF%ICO_MARKER = SLF%NC

   IOR0_LAST = 0
   WALLCELLS_LOOP2: DO IWG = 1, SLF%NW

      !< Determine neighbors and orientation again
      NOM  = SLF%WALL(IWG)%NOM
      IOR0 = SLF%WALL(IWG)%IOR
      NCPL = SLF%WALL(IWG)%NCPL

      !< Determine boundary type for IW
      IF (M%WALL(IWG)%NOM /= 0) THEN 
         SLF%WALL(IWG)%BTYPE = INTERNAL
         IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
            WRITE(SCARC_LU,*) NM,': INTERNAL : WALL(',IWG,')=',SLF%WALL(IWG)%BTYPE
      ELSE IF (M%WALL(IWG)%PRESSURE_BC_INDEX == DIRICHLET) THEN
         SLF%WALL(IWG)%BTYPE = DIRICHLET
         IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
            WRITE(SCARC_LU,*) NM,': DIRICHLET: WALL(',IWG,')=',SLF%WALL(IWG)%BTYPE
      ELSE
         SLF%WALL(IWG)%BTYPE = NEUMANN
         IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
            WRITE(SCARC_LU,*) NM,': NEUMANN  : WALL(',IWG,')=',SLF%WALL(IWG)%BTYPE
      ENDIF

      !< setup coordinate information for that wall cell (this is independent of a possible neighbor)
      SLF%WALL(IWG)%ICW = (M%WALL(IWG)%ONE_D%KKG-1)*M%IBAR*M%JBAR + &           !< adjacent (internal) cell number
                          (M%WALL(IWG)%ONE_D%JJG-1)*M%IBAR        + &
                           M%WALL(IWG)%ONE_D%IIG

      SLF%WALL(IWG)%IXG =  M%WALL(IWG)%ONE_D%II                                 !< ghost cell indices
      SLF%WALL(IWG)%IYG =  M%WALL(IWG)%ONE_D%JJ
      SLF%WALL(IWG)%IZG =  M%WALL(IWG)%ONE_D%KK

      SLF%WALL(IWG)%IXW =  M%WALL(IWG)%ONE_D%IIG                                !< (internal) wall cell indices
      SLF%WALL(IWG)%IYW =  M%WALL(IWG)%ONE_D%JJG
      SLF%WALL(IWG)%IZW =  M%WALL(IWG)%ONE_D%KKG

      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
         WRITE(SCARC_LU,*) 'CALLING WALLCELL_NEIGHBOR FOR IWG=',IWG
         WRITE(SCARC_LU,*) ' NOM=',NOM
         WRITE(SCARC_LU,*) ' IOR=',IOR0
         WRITE(SCARC_LU,*) ' NOM_IB(1)=',M%WALL(IWG)%NOM_IB(1)
         WRITE(SCARC_LU,*) ' NOM_IB(2)=',M%WALL(IWG)%NOM_IB(2)
         WRITE(SCARC_LU,*) ' NOM_IB(3)=',M%WALL(IWG)%NOM_IB(3)
         WRITE(SCARC_LU,*) ' NOM_IB(4)=',M%WALL(IWG)%NOM_IB(4)
         WRITE(SCARC_LU,*) ' NOM_IB(5)=',M%WALL(IWG)%NOM_IB(5)
         WRITE(SCARC_LU,*) ' NOM_IB(6)=',M%WALL(IWG)%NOM_IB(6)
      ENDIF

      !< If there exists a neighbor for that wall cell, setup corresponding neighborship information
      IF (NOM /= 0) THEN
         CALL SCARC_SETUP_WALLCELL_NEIGHBOR(M%WALL(IWG)%NOM_IB(1), &
                                            M%WALL(IWG)%NOM_IB(4), &
                                            M%WALL(IWG)%NOM_IB(2), &
                                            M%WALL(IWG)%NOM_IB(5), &
                                            M%WALL(IWG)%NOM_IB(3), &
                                            M%WALL(IWG)%NOM_IB(6), &
                                            IWG, IOR0, NM, NOM, NLEVEL_MIN)
      ENDIF

      NOM_LAST  = NOM
      IOR0_LAST = IOR0

   ENDDO WALLCELLS_LOOP2

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '---> AFTER FINEST LEVEL SETTINGS '

   !! Only in case of MG-method:
   !! Determine WALL, FACE and OSCARC types for coarser levels
   LEVEL_GMG_IF: IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
   
      IREFINE=1
      LEVEL_GMG_LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX
      
         !< point to SCARC grid structures on coarser and finer level
         SLC => S%LEVEL(NL)
         SLF => S%LEVEL(NL-1)
   
         IREFINE=IREFINE*2

         CALL SCARC_CHECK_DIVISIBILITY(SLF%NCE-SLF%NC, 'SLC%NCE', 'SCARC_SETUP_WALLINFO')
         CALL SCARC_CHECK_DIVISIBILITY(SLF%NCO-SLF%NC, 'SLC%NCO', 'SCARC_SETUP_WALLINFO')

         SLC%NCE = SLC%NC + (SLF%NCE-SLF%NC)/2
         SLC%NCO = SLC%NC + (SLF%NCO-SLF%NC)/2

         SLC%ICE_MARKER = SLC%NC
         SLC%ICO_MARKER = SLC%NC

         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
            WRITE(SCARC_LU,*) 'SLC%NCE=',SLC%NCE
            WRITE(SCARC_LU,*) 'SLC%NCO=',SLC%NCO
         ENDIF

         !! allocate administrative mapping arrays for own level
         CALL SCARC_SETUP_MAPPINGS(NM, NL)

         !< compute FACE and WALL information for all faces of coarser level
         IWC = 1
         IWG = 1
         DO IFACE = 1, 6

            IOR0 = FACE_ORDER_XYZ(IFACE)

            !< compute mesh dimensions of coarser mesh level 
            CALL SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)
   
            !< for every neighbor do:
            IF (SLF%FACE(IOR0)%NUM_NEIGHBORS /= 0) THEN
               DO INBR_FACE = 1, SLF%FACE(IOR0)%NUM_NEIGHBORS

                  NOM = SLF%FACE(IOR0)%NEIGHBORS(INBR_FACE)
 
                  OSLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL-1)

                  OSLC%NCPL = OSLF%NCPL

                  CALL SCARC_CHECK_DIVISIBILITY(OSLF%NWL, 'OSLF%NWL', 'SCARC_SETUP_WALLINFO')
                  OSLC%NWL = OSLF%NWL/2

                  CALL SCARC_CHECK_DIVISIBILITY(OSLF%NCG, 'OSLF%NCG', 'SCARC_SETUP_WALLINFO')
                  OSLC%NCG = OSLF%NCG/2

                  IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
                     WRITE(SCARC_LU,'(4(a,i4))') 'NOM=',NOM,': NCPL=',OSLC%NCPL,': NWL=',OSLC%NWL,': NL=',NL
                  ENDIF

                  !< compute exchange dimensions of coarser mesh level (analogues to NIC, I_MIN ... )
                  CALL SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)

                  !< allocate administrative mapping arrays 
                  CALL SCARC_SETUP_OMAPPINGS(NM, NOM, NL)

               ENDDO
            ENDIF

            !< setup complete face information for coarser mesh
            CALL SCARC_SETUP_FACE(IOR0, IWC, IREFINE, NM, NL)
         ENDDO
   
   
      ENDDO LEVEL_GMG_LEVEL_LOOP
   ENDIF LEVEL_GMG_IF

ENDDO MESHES_LOOP


IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'BEFORE EXCHANGE'
DO NL=NLEVEL_MIN, NLEVEL_MAX
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL, 'SETUP_WALL', 'WALL')
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACEINFO , NL, 'SETUP_WALL', 'FACE')
ENDDO


IF (TYPE_DEBUG  > NSCARC_DEBUG_NONE      .AND. &
    TYPE_LATEX  > NSCARC_LATEX_ALL)   THEN
!>   TYPE_LATEX  == NSCARC_LATEX_ALL)   THEN
   MLATEX = 98
   DO NM = NMESHES_MIN, NMESHES_MAX
      CLATEX = "tables/wall  .tex"
      WRITE(CLATEX(12:13),'(i2.2)') NM
      OPEN (MLATEX, FILE=CLATEX)
      SLF => SCARC(NM)%LEVEL(NLEVEL_MIN)
      CALL SCARC_LATEX_WALL(MLATEX, NM, NLEVEL_MIN)
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

END SUBROUTINE SCARC_SETUP_WALLINFO

!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLCELL_NEIGHBOR(NX1, NX2, NY1, NY2, NZ1, NZ2, IWG, IOR0, NM, NOM, NL)
INTEGER, INTENT(IN) :: NX1, NX2, NY1, NY2, NZ1, NZ2
INTEGER, INTENT(IN) :: IWG, IOR0, NM, NOM, NL
INTEGER :: NOMX, NOMY, NOMZ
INTEGER :: ICG, ICO, ICE, ICN, IWL, ICPL, IERR=0, IX, IY, IZ, JL
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, OSL

SL  => SCARC(NM)%LEVEL(NL)
OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

!< store basic information about neighbor and orientation
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

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) '========== CALLING WALLCELL_NEIGHBOR FOR ', IOR0
   WRITE(SCARC_LU,*) 'XXX: OSL%NCPL=',OSL%NCPL
   WRITE(SCARC_LU,*) 'XXX: OSL%IOR =',OSL%IOR 
   WRITE(SCARC_LU,*) 'XXX: OSL%NCG =',OSL%NCG 
   WRITE(SCARC_LU,*) 'XXX: OSL%NWL =',OSL%NWL 
   WRITE(SCARC_LU,*) 'XXX: OSL%NCPLR =',OSL%NCPLR 
   WRITE(SCARC_LU,*) 'XXX: OSL%NCPLS =',OSL%NCPLS 
   WRITE(SCARC_LU,*) 'XXX: ICE     =', ICE
   WRITE(SCARC_LU,*) 'XXX: ICO     =', ICO
   WRITE(SCARC_LU,*) 'XXX: ICG     =', ICG
   WRITE(SCARC_LU,*) 'XXX: IWL     =', IWL
   WRITE(SCARC_LU,*) 'XXX: IWG     =', IWG
   WRITE(SCARC_LU,*) 'XXX: NX1, NX2     =', NX1, NX2
   WRITE(SCARC_LU,*) 'XXX: NY1, NY2     =', NY1, NY2
   WRITE(SCARC_LU,*) 'XXX: NZ1, NZ2     =', NZ1, NZ2
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'1: ICE_MARKER=',SL%ICE_MARKER
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'1: ICO_MARKER=',SL%ICO_MARKER
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'1: ICG_MARKER=',OSL%ICG_MARKER
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'1: IWL_MARKER=',OSL%IWL_MARKER
ENDIF

!< set neighboring coordinates
SL%WALL(IWG)%IXN(1) = NX1
SL%WALL(IWG)%IXN(2) = NX2
SL%WALL(IWG)%IYN(1) = NY1
SL%WALL(IWG)%IYN(2) = NY2
SL%WALL(IWG)%IZN(1) = NZ1
SL%WALL(IWG)%IZN(2) = NZ2

!< allocate pointer arrays for extended, ghost and neighboring cells
ALLOCATE(SL%WALL(IWG)%ICE(OSL%NCPL), STAT=IERR)   
CALL ChkMemErr('SCARC_SETUP_WALLINFO','ICE',IERR)
SL%WALL(IWG)%ICE = NSCARC_DUMMY

ALLOCATE(SL%WALL(IWG)%ICG(OSL%NCPL), STAT=IERR)   
CALL ChkMemErr('SCARC_SETUP_WALLINFO','ICG',IERR)
SL%WALL(IWG)%ICG = NSCARC_DUMMY

ALLOCATE(SL%WALL(IWG)%ICN(OSL%NCPL), STAT=IERR)   
CALL ChkMemErr('SCARC_SETUP_WALLINFO','ICN',IERR)
SL%WALL(IWG)%ICN = NSCARC_DUMMY

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

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'MESHES(',NM,')%IBAR=',MESHES(NM)%IBAR
   WRITE(SCARC_LU,*) 'MESHES(',NM,')%JBAR=',MESHES(NM)%JBAR
   WRITE(SCARC_LU,*) 'MESHES(',NM,')%KBAR=',MESHES(NM)%KBAR
   WRITE(SCARC_LU,*) 'MESHES(',NOM,')%IBAR=',MESHES(NOM)%IBAR
   WRITE(SCARC_LU,*) 'MESHES(',NOM,')%JBAR=',MESHES(NOM)%JBAR
   WRITE(SCARC_LU,*) 'MESHES(',NOM,')%KBAR=',MESHES(NOM)%KBAR
   WRITE(SCARC_LU,*) 'NOMX=',NOMX
   WRITE(SCARC_LU,*) 'NOMY=',NOMY
   WRITE(SCARC_LU,*) 'NOMZ=',NOMZ
ENDIF

!< store information about overlapped cells and set mapping arrays
ICPL = 0
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ICPL = ICPL + 1
         ICG  = ICG  + 1
         ICE  = ICE  + 1

         ICN  = (IZ-1)*NOMX*NOMY + (IY-1)*NOMX + IX

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,'(8(a,i4))') 'IX=',IX, ' :IY=',IY, ' :IZ=',IZ, ' :IWL=',IWL,' :ICE=', ICE, ' :ICN=',ICN, ' :ICG=',ICG
ENDIF

         SL%WALL(IWG)%ICN(ICPL) = ICN                    !< number of neighboring grid cell
         SL%WALL(IWG)%ICE(ICPL) = ICE                    !< number of extended grid cell 
         SL%WALL(IWG)%ICG(ICPL) = ICG                    !< number of ghost grid cell

         SL%ICE_TO_IWG(ICE)     = IWG                    !< map extended cell to global wall cell
         SL%ICE_TO_IWL(ICE)     = IWL                    !< map extended cell to local wall cell
         SL%ICE_TO_ICG(ICE)     = ICG                    !< map extended cell to ghost cell
         SL%ICE_TO_ICN(ICE)     = ICN                    !< map extended cell to neighboring cell
         SL%ICE_TO_ICW(ICE)     = SL%WALL(IWG)%ICW       !< map extended cell to ghost cell

         OSL%ICG_TO_IWG(ICG)    = IWG                    !< map ghost cell to global wall cell
         OSL%ICG_TO_ICN(ICG)    = ICN                    !< map ghost cell to neighboring cell
         OSL%ICG_TO_ICE(ICG)    = SL%WALL(IWG)%ICE(ICPL) !< map ghost cell to extended grid cell
         OSL%ICG_TO_ICW(ICG)    = SL%WALL(IWG)%ICW       !< map ghost cell to internal grid cell

         OSL%IWL_TO_ICG(IWL, ICPL)  = ICG                !< map extended cell to ghost cell
         OSL%IWL_TO_ICN(IWL, ICPL)  = ICN                !< map extended cell to neighboring cell
      ENDDO
   ENDDO
ENDDO

SL%WALL(IWG)%ICO = ICO                                   !< number of overlapping cell
SL%WALL(IWG)%IWL = IWL                                   !< number of local wall cell

OSL%IWL_TO_IWG(IWL) = IWG                                !< map local wall cell to global wall cell
OSL%IWL_TO_ICO(IWL) = ICO                                !< map local wall cell to internal grid cell
OSL%IWL_TO_ICW(IWL) = SL%WALL(IWG)%ICW                   !< map local wall cell to internal grid cell

OSL%IWL_MARKER = IWL                                     !< store local wall cell marker
OSL%ICG_MARKER = ICG                                     !< store ghost cell counter

SL%ICO_MARKER  = ICO                                     !< store overlapping cell counter
SL%ICE_MARKER  = ICE                                     !< store extended cell counter

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'NM=',NM,': NOM=',NOM,': IWL          =',IWL
   WRITE(SCARC_LU,*) 'NM=',NM,': NOM=',NOM,': ICE          =',ICE
   WRITE(SCARC_LU,*) 'NM=',NM,': NOM=',NOM,': ICO          =',ICO
   WRITE(SCARC_LU,*) 'NM=',NM,': NOM=',NOM,': OSL%NWL     =',OSL%NWL
   WRITE(SCARC_LU,*) 'NM=',NM,': NOM=',NOM,': OSL%IOR     =',OSL%IOR
   WRITE(SCARC_LU,*) 'NM=',NM,': NOM=',NOM,': OSL%NCG     =',OSL%NCG 
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,' ICO=',SL%WALL(IWG)%ICO
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,' ICE=',SL%WALL(IWG)%ICE
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,' ICE=',SL%WALL(IWG)%ICG
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,' ICN=',SL%WALL(IWG)%ICN
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'2: ICE_MARKER=',SL%ICE_MARKER
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'2: ICO_MARKER=',SL%ICO_MARKER
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'2: ICG_MARKER=',OSL%ICG_MARKER
   WRITE(SCARC_LU,'(a,i4,a,i4,a,4i4)') 'NM=',NM,': NOM=',NOM,'2: IWL_MARKER=',OSL%IWL_MARKER
ENDIF

END SUBROUTINE SCARC_SETUP_WALLCELL_NEIGHBOR


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MAPPINGS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: NC, NCO, NCE, NW, IERR
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL  => SCARC(NM)%LEVEL(NL)
NC  = SL%NC
NW  = SL%NW
NCE = SL%NCE
NCO = SL%NCO

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,'(8(a,i4))') 'SETUP MAPPINGS FOR NC=',NC, ': NCE=',NCE,': NCO=',NCO,': NW=',NW,': NL=',NL
   WRITE(SCARC_LU,'(a,2i4)') 'ALLOCATING ICE_TO_IWG in length ', NC+1, NCE
ENDIF
      
!> Allocate mappings for ICE
ALLOCATE(SL%ICE_TO_IWG(NC+1:NCE), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','ICE_TO_IWG',IERR)
SL%ICE_TO_IWG=0

ALLOCATE(SL%ICE_TO_IWL(NC+1:NCE), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','ICE_TO_IWL',IERR)
SL%ICE_TO_IWL=0

ALLOCATE(SL%ICE_TO_ICW(NC+1:NCE), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','ICE_TO_ICW',IERR)
SL%ICE_TO_ICW=0

ALLOCATE(SL%ICE_TO_ICG(NC+1:NCE), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','ICE_TO_ICG',IERR)
SL%ICE_TO_ICG=0

ALLOCATE(SL%ICE_TO_ICN(NC+1:NCE), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','SL%ICE_TO_ICN',IERR)
SL%ICE_TO_ICN = 0
      

!> Allocate remaining mappings (outdated)
ALLOCATE(SL%INTERNAL_BDRY_CELL(NC), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','INTERNAL_BDRY_CELL',IERR)
SL%INTERNAL_BDRY_CELL=0

ALLOCATE(SL%WALL_INDEX(NC, -3:3), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','WALL_INDEX',IERR)
SL%WALL_INDEX = 0


END SUBROUTINE SCARC_SETUP_MAPPINGS


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OMAPPINGS(NM, NOM, NL)
INTEGER, INTENT(IN) :: NM, NOM, NL
INTEGER:: NWL, NCPL, IERR
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSL

!> Get corresponding dimensions
OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

NWL  = OSL%NCG
NCPL = OSL%NCPL

!> initialize some counters
OSL%IWL = 0       
OSL%ICG = 0

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,'(8(a,i4))') 'SETUP OMAPPINGS FOR NWL=',NWL, ': NCPL=',NCPL,': NOM=',NOM,': NL=',NL
ENDIF

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

ALLOCATE(OSL%IWL_TO_ICN(NWL, NCPL), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','SL%IWL_TO_ICN',IERR)
OSL%IWL_TO_ICN = 0
      
ALLOCATE(OSL%IWL_TO_ICG(NWL, NCPL), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_MAPPINGS','SL%IWL_TO_ICG',IERR)
OSL%IWL_TO_ICG = 0
      

!> Allocate mapping arrays for ghost cells 
ALLOCATE(OSL%ICG_TO_IWG(NWL), STAT=IERR)              
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_IWG',IERR)
OSL%ICG_TO_IWG = 0

ALLOCATE(OSL%ICG_TO_ICE(NWL), STAT=IERR)              
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_ICE',IERR)
OSL%ICG_TO_ICE = 0

ALLOCATE(OSL%ICG_TO_ICW(NWL), STAT=IERR)              
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_ICW',IERR)
OSL%ICG_TO_ICW = 0

ALLOCATE(OSL%ICG_TO_ICN(NWL), STAT=IERR)              
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_ICN',IERR)
OSL%ICG_TO_ICN = 0



!ALLOCATE(OSL%ICN_TO_ICE(1:NCO), STAT=IERR)   !< kann verbessert werden??? riesig!
!CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICN_TO_ICE',IERR)
!OSL%IWL_TO_IWG = 0

      
END SUBROUTINE SCARC_SETUP_OMAPPINGS


!> ----------------------------------------------------------------------------------------------------
!> Check divisibility by 2 of a given number of elements (in one grid direction) 
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_DIVISIBILITY(NN, CDIR, CROUTINE)
INTEGER, INTENT(IN) :: NN
CHARACTER(*) , INTENT(IN) :: CDIR, CROUTINE
IF (MOD(NN,2) /= 0) THEN
   WRITE(SCARC_LU,1000) CROUTINE, CDIR
   STOP
ENDIF
1000 FORMAT('ERROR in ', a20,' number of wall cells in ',a1,' not divisible by 2 !!')
END SUBROUTINE SCARC_CHECK_DIVISIBILITY


!> ----------------------------------------------------------------------------------------------------
!> Set wall cell information on coarse level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)
INTEGER, INTENT(IN)    :: IOR0, NM, NL
INTEGER, INTENT(INOUT) :: IWG
INTEGER:: INBR_FACE, IERR
TYPE (SCARC_FACE_TYPE) , POINTER               :: FC, FF              !< coarse and fine FACE types
TYPE (SCARC_LEVEL_TYPE), POINTER               :: SLC, SLF            !< coarse and fine LEVEL types 

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
      FC%IOFFSET_MATRIX = 1

      FC%NFX = 1                                                         !< no extension in x-direction
      FC%NFC = SLC%NX                                                    !< number of cells between opposite faces

      IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) THEN                 !< only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y', 'SCARC_SETUP_FACE')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF

      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z', 'SCARC_SETUP_FACE')     !< number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2

   CASE (2)

      FC%DH => SLC%DYL
      FC%IOFFSET_MATRIX = SLC%NX

      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X', 'SCARC_SETUP_FACE')     !< number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2

      FC%NFY = 1                                                         !< no extension in y-direction   
      FC%NFC = SLC%NY                                                    !< number of cells between opposite faces

      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z', 'SCARC_SETUP_FACE')     !< number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2

   CASE (3)

      FC%DH => SLC%DZL
      FC%IOFFSET_MATRIX = SLC%NX * SLC%NY

      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X', 'SCARC_SETUP_FACE')     !< number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2

      IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) THEN                 !< only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y', 'SCARC_SETUP_FACE')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF

      FC%NFZ = 1                                                         !< no extension in y-direction   
      FC%NFC = SLC%NZ                                                    !< number of cells between opposite faces
END SELECT
IF (IOR0 > 0) FC%IOFFSET_MATRIX = -FC%IOFFSET_MATRIX

FC%NFW = FC%NFX * FC%NFY * FC%NFZ                                        !< get number of wall cells for that face

IWG = IWG + FC%NFW                                                       !< increase global wall cell counter

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) '========================== SETUP_FACE_DIMENSIONS:  IOR0=',IOR0
   WRITE(SCARC_LU,*) 'FC%NFW=',FC%NFW
   WRITE(SCARC_LU,*) 'FC%NFX=',FC%NFX
   WRITE(SCARC_LU,*) 'FC%NFY=',FC%NFY
   WRITE(SCARC_LU,*) 'FC%NFZ=',FC%NFZ
   WRITE(SCARC_LU,*) 'FF%IWG_MARKER=',FF%IWG_MARKER
   WRITE(SCARC_LU,*) 'IWG=',IWG
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
TYPE (SCARC_WALL_TYPE) , POINTER, DIMENSION(:) :: WC, WF              !< coarse and fine WALL types
TYPE (SCARC_FACE_TYPE) , POINTER               :: FF, FC              !< coarse and fine FACE types
TYPE (SCARC_LEVEL_TYPE), POINTER               :: SLC, SLF, OSLC      !< coarse and fine LEVEL types 

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
      IF (IOR0 > 0) THEN                                           !< set dimensions for wall cell counting
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
      NX1 = 1                                                      !< set dimensions for wall cell counting
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
      NX1 = 1                                                      !< set dimensions for wall cell counting
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

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) '========================== SETUP_FACE:  IOR0=',IOR0
   WRITE(SCARC_LU,*) 'NX1=',NX1
   WRITE(SCARC_LU,*) 'NX2=',NX2
   WRITE(SCARC_LU,*) 'NY1=',NY1
   WRITE(SCARC_LU,*) 'NY2=',NY2
   WRITE(SCARC_LU,*) 'NZ1=',NZ1
   WRITE(SCARC_LU,*) 'NZ2=',NZ2
ENDIF

!>
!> Loop over all wall cells of face IOR0
!>
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2
         
         !< Set orientation of neiboring face, indices of ghost and adjacent cell for coarse IW
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
         
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
            WRITE(SCARC_LU,*) 'WC(',IWC,')%IOR0=',WC(IWC)%IOR
            WRITE(SCARC_LU,*) 'WC(',IWC,')%ICW =',WC(IWC)%ICW
            WRITE(SCARC_LU,*) 'WC(',IWC,')%IXG =',WC(IWC)%IXG
            WRITE(SCARC_LU,*) 'WC(',IWC,')%IYG =',WC(IWC)%IYG
            WRITE(SCARC_LU,*) 'WC(',IWC,')%IZG =',WC(IWC)%IZG
            WRITE(SCARC_LU,*) 'WC(',IWC,')%IXW =',WC(IWC)%IXW
            WRITE(SCARC_LU,*) 'WC(',IWC,')%IYW =',WC(IWC)%IYW
            WRITE(SCARC_LU,*) 'WC(',IWC,')%IZW =',WC(IWC)%IZW
         ENDIF

         SELECT CASE (TYPE_DIMENSION)
         
            !< ------------------------------------------------------------------------------------
            !< 2D-version
            !< ------------------------------------------------------------------------------------
            CASE (NSCARC_DIMENSION_TWO)
         
               !< determine fine IW's, which must be merged to one coarse IW
               SELECT CASE (ABS(IOR0))
                  CASE ( 1)
                     IWF(1) = FF%IWG_MARKER + 2*(IZ-1)
                  CASE ( 2)
                     IWF(1) = FF%IWG_MARKER + 2*(IZ-1)*SLF%NX + 2*(IX - 1)
                  CASE ( 3)
                     IWF(1) = FF%IWG_MARKER + 2*(IX-1)
               END SELECT
               IWF(2) = IWF(1)+1
            
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
                  WRITE(SCARC_LU,*) 'IWF(1)=',IWF(1)
                  WRITE(SCARC_LU,*) 'IWF(2)=',IWF(2)
               ENDIF

               !< set fine cell neighbors (they must be the same for all fine IW's)
               NOMF(1) = WF(IWF(1))%NOM
               NOMF(2) = WF(IWF(2))%NOM
               IF (NOMF(1) /= NOMF(2)) THEN
                  WRITE(*,*) 'SCARC_SETUP_FACE: Inconsistent neighbors on IOR=', IOR0,' not allowed!'
                  STOP
               ENDIF
            
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
                  WRITE(SCARC_LU,*) 'NOMF(1)=',NOMF(1)
                  WRITE(SCARC_LU,*) 'NOMF(2)=',NOMF(2)
               ENDIF

               WC(IWC)%NOM = NOMF(1) 
            
               !< set corresponding pressure_bc_index on coarser level
               IBCF(1) = WF(IWF(1))%BTYPE
               IBCF(2) = WF(IWF(2))%BTYPE
               IF (IBCF(1) == INTERNAL .OR. IBCF(2) == INTERNAL) THEN
                  WC(IWC)%BTYPE = INTERNAL
               ELSE IF (IBCF(1) == DIRICHLET .OR. IBCF(2) == DIRICHLET) THEN
                  WC(IWC)%BTYPE = DIRICHLET
               ELSE
                  WC(IWC)%BTYPE = NEUMANN
               ENDIF
            
               !< in case of an internal boundary set neighboring WALL cells
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
                           WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                           STOP
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
                           WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                           STOP
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
                  !< Allocate and specify ICN and ICE arrays for OSLC
                  !!!
                  NCPL = (IZ2-IZ1+1)*(IX2-IX1+1)
                  OSLC%NCPL = NCPL
         
                  CALL SCARC_SETUP_WALLCELL_NEIGHBOR(IX1, IX2, 1, 1, IZ1, IZ2, IWC, IOR0, NM, NOMF(1), NL)

               ENDIF
                  
         
            !< ------------------------------------------------------------------------------------
            !< 3D-version
            !< ------------------------------------------------------------------------------------
            CASE (NSCARC_DIMENSION_THREE)
         
               !< determine fine IW's, which must be merged to one coarse IW
               SELECT CASE (ABS(IOR0))
                  CASE (1)
                     IWF(1) = FF%IWG_MARKER + (2*IZ-2)*SLF%NY + 2*IY - 1
                     IWF(3) = FF%IWG_MARKER + (2*IZ-1)*SLF%NY + 2*IY - 1
                  CASE (2)
                     IWF(1) = FF%IWG_MARKER + (2*IZ-2)*SLF%NX + 2*IX - 1
                     IWF(3) = FF%IWG_MARKER + (2*IZ-1)*SLF%NX + 2*IX - 1
                  CASE (3)
                     IWF(1) = FF%IWG_MARKER + (2*IY-2)*SLF%NX + 2*IX - 1
                     IWF(3) = FF%IWG_MARKER + (2*IY-1)*SLF%NX + 2*IX - 1
               END SELECT
               IWF(2) = IWF(1)+1
               IWF(4) = IWF(3)+1
            
               !< set fine cell neighbors (they must be the same for all fine IW's)
               DO I=1,4
                  NOMF(I) = WF(IWF(I))%NOM
               ENDDO
                  
               IF (NOMF(1)/=NOMF(2) .OR. NOMF(1)/=NOMF(3) .OR. NOMF(1)/=NOMF(4)) THEN
                  WRITE(*,*) 'SCARC_SETUP_FACE: Inconsistent neighbors on IOR=', IOR0,' not allowed!'
                  STOP
               ENDIF
               WC(IWC)%NOM = NOMF(1) 
            
               !< set corresponding pressure_bc_index on coarser level
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
            
               !< in case of an internal boundary set WALL(10:15,IWC)
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
                           WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                           STOP
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
                           WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                           STOP
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
                           WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                           STOP
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
         
                  !< Allocate and specify ICN and ICE arrays for OSLC
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

OS%NIC_S = 0
FOUND = .FALSE.

SEARCH_LOOP: DO IW=1, OSL%NCG

   !< neighborship structure already known from finest level
   IF (SL%WALL(IW)%NOM/=NOM) CYCLE SEARCH_LOOP
   OS%NIC_S = OS%NIC_S + 1
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

OS%I_MIN_R = IMIN
OS%I_MAX_R = IMAX
OS%J_MIN_R = JMIN
OS%J_MAX_R = JMAX
OS%K_MIN_R = KMIN
OS%K_MAX_R = KMAX

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) '============ SETUP_EXCHANGE_DIMENSIONS FOR NL=',NL
   WRITE(SCARC_LU,*) 'IMIN=',IMIN
   WRITE(SCARC_LU,*) 'IMAX=',IMAX
   WRITE(SCARC_LU,*) 'JMIN=',JMIN
   WRITE(SCARC_LU,*) 'JMAX=',JMAX
   WRITE(SCARC_LU,*) 'KMIN=',KMIN
   WRITE(SCARC_LU,*) 'KMAX=',KMAX
ENDIF

END SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS

 
!> ----------------------------------------------------------------------------------------------------
!> Allocate several global structures for data exchange 
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBAL_DIMENSIONS
INTEGER :: NM, NL, IERR=0
INTEGER :: IREFINE

!> Allocate arrays which are used for (global) communciations
ALLOCATE(NC_GLOBAL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'NC_GLOBAL', IERR)
NC_GLOBAL = 0

ALLOCATE(NC_LOCAL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'NC_LOCAL', IERR)
NC_LOCAL = 0

ALLOCATE(SP_LOCAL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBAL', 'SP_LOCAL', IERR)
SP_LOCAL = 0.0_EB

!> Compute global number of cells for all levels
IREFINE=0
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

   MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      NC_LOCAL(NM)=SCARC(NM)%LEVEL(NL)%NC
   ENDDO MESHES_LOOP

   !! Determine global number of cells for all levels 
   IF (NMESHES>1) THEN
      IF (USE_MPI) THEN
         CALL MPI_ALLREDUCE(NC_LOCAL(MYID+1),NC_GLOBAL(NL),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
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

END SUBROUTINE SCARC_SETUP_GLOBAL_DIMENSIONS


!> ----------------------------------------------------------------------------------------------------
!> Allocate several global structures for data exchange 
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_NEIGHBORS
INTEGER :: NM, NOM, NL, IERR
TYPE (OSCARC_TYPE), POINTER :: OS, OSO


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

 
!> ----------------------------------------------------------------------------------------------------
!> Setup system of equation:
!> Define matrix stencils and initialize matrices and boundary conditions on all needed levels
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM
INTEGER :: NM, NL

SELECT CASE (TYPE_MULTIGRID)
   CASE (NSCARC_MULTIGRID_GEOMETRIC)
      DO NL=NLEVEL_MIN, NLEVEL_MAX
         CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NL)
      ENDDO
   CASE DEFAULT
      CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NLEVEL_MIN)
END SELECT


MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

WRITE(*,*) 'ALLOCATING MATRIX FOR MESH ',NM
   SELECT_SOLVER: SELECT CASE (TYPE_METHOD)
   
      !< ---------------------------------------------------------------------------------------------
      !< Krylov method (CG/BICG) as main solver, different preconditioners possible
      !< ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_KRYLOV)
   
         SELECT_PRECON: SELECT CASE (TYPE_PRECON)
   
            !< ---------------------------------------------------------------------------------------
            !< in case of multigrid as preconditioner:
            !< ---------------------------------------------------------------------------------------
            CASE (NSCARC_PRECON_MULTIGRID)
   
               SELECT_PRECON_MG: SELECT CASE (TYPE_MULTIGRID)
   
                  !< Geometric multigrid:
                  !<    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
                  CASE (NSCARC_MULTIGRID_GEOMETRIC)
                     
                     DO NL = NLEVEL_MIN, NLEVEL_MAX
                        CALL SCARC_SETUP_MATRIX  (NM, NL)
                        CALL SCARC_SETUP_BOUNDARY(NM, NL)
                     ENDDO 
   
                  !< Algebraic multigrid:
                  !<    -  use compact storage technique on all levels (no other choise possible!)
                  !<    -  assemble standard n-point-matrix only on finest level 
                  !<    -  construct all coarser levels by requested coarsening strategy
                  CASE (NSCARC_MULTIGRID_ALGEBRAIC)
   
                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
   
               END SELECT SELECT_PRECON_MG
   
            !< ---------------------------------------------------------------------------------------
            !< in case of one-level preconditioners (JACOBI/SSOR/FFT)
            !<    -  assemble standard n-point-matrix on finest level 
            !< ---------------------------------------------------------------------------------------
            CASE DEFAULT
   
               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
   
         END SELECT SELECT_PRECON
   
      !< ---------------------------------------------------------------------------------------------
      !< Multigrid as main solver
      !< ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_MULTIGRID)
   
         SELECT_MG: SELECT CASE (TYPE_MULTIGRID)
   
            !< ---------------------------------------------------------------------------------------
            !< Geometric multigrid:
            !<    -  assemble standard n-point-matrix hierarchy on all levels 
            !< ---------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_GEOMETRIC)
   
               DO NL = NLEVEL_MIN, NLEVEL_MAX
                  CALL SCARC_SETUP_MATRIX  (NM, NL)
                  CALL SCARC_SETUP_BOUNDARY(NM, NL)
               ENDDO 

   
            !< ---------------------------------------------------------------------------------------
            !< Algebraic multigrid:
            !<    -  use compact storage technique (no other choice possible!)
            !<    -  assemble standard n-point-matrix only on finest level
            !<    -  construct all coarser levels later by requested coarsening strategy
            !< ---------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_ALGEBRAIC)
   
               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         END SELECT SELECT_MG
   
   END SELECT SELECT_SOLVER
   
ENDDO MESHES_LOOP


!> ------------------------------------------------------------------------------------------------
!> Exchange matrix stencil information in case of AMG
!> ------------------------------------------------------------------------------------------------
IF (NMESHES>1 .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_ALGEBRAIC) THEN

   !< set sizes for matrix stencils on overlapping parts
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_STENCIL, NLEVEL_MIN)
   CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_STENCIL, NLEVEL_MIN)

   !< exchange matrix entries on overlapping parts
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SYSTEM, NLEVEL_MIN)
   CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_SYSTEM, NLEVEL_MIN)

ENDIF


!> ------------------------------------------------------------------------------------------------
!> If multigrid is used, setup global coarse grid matrix in case of a direct coarse grid solver
!> ------------------------------------------------------------------------------------------------
IF ((TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_PRECON_MULTIGRID) .AND. &
     TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
   WRITE(*,*) 'ACHTUNG: SETUP_COARSE_MATRIX muss noch an die neue WALL-Struktur angepasst werden!'
   STOP
   !CALL SCARC_SETUP_COARSE_MATRIX(NLEVEL_MAX)
ENDIF


!> ------------------------------------------------------------------------------------------------
!> Only temporarily: Show all relevant arrays and structures
!> ------------------------------------------------------------------------------------------------
SELECT CASE (TYPE_MULTIGRID)
   CASE (NSCARC_MULTIGRID_GEOMETRIC)
      DO NL=NLEVEL_MIN, NLEVEL_MAX
         CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NL, 'SETUP_SYSTEM', 'MATRIX000')
         CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL, 'SETUP_SYSTEM', 'WALLINFO')
      ENDDO
   CASE (NSCARC_MULTIGRID_ALGEBRAIC)
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NLEVEL_MIN, 'SETUP_SYSTEM', 'MATRIX000')
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NLEVEL_MIN, 'SETUP_SYSTEM', 'WALLINFO')
END SELECT
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_BCINDEX, NLEVEL_MIN, 'SETUP_SYSTEM', 'PRESSURE_BC_INDEX')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_ACELL  , NLEVEL_MIN, 'SETUP_SYSTEM', 'WALL_CELL')
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
INTEGER :: IX, IY, IZ
INTEGER :: IC, IP, IERR
INTEGER :: NA, NC
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

NA = SL%NAE
NC = SL%NCE              ! ONLY TEMPORARILY
NC = SL%NCO

ALLOCATE (SL%A_ROW(NC+10), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%A_ROW', IERR)
SL%A_ROW = 0

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) '------------------ CALLING SETUP_MATRIX FOR MESH ',NM
   WRITE(SCARC_LU,*) 'SETUP_MATRIX: NA=',SL%NAE
   WRITE(SCARC_LU,*) 'SETUP_MATRIX: NC=',SL%NCO
ENDIF

SELECT CASE (TYPE_DIMENSION)

   !! --------------------- 2D -----------------
   CASE (NSCARC_DIMENSION_TWO)

      ALLOCATE (SL%A(NA+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%A', IERR)
      SL%A = 0.0_EB
      
      ALLOCATE (SL%A_COL(NA+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%A_COL', IERR)
      SL%A_COL = 0

      IP  = 1
      IY  = 1
      DO IZ = 1, SL%NZ
         DO IX = 1, SL%NX
  
            IC = (IZ-1) * SL%NX + IX

            CALL SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)

            CALL SCARC_SETUP_MATRIX_SUBDIAG  ( 3, IC, IP, IZ, NM, NL)
            CALL SCARC_SETUP_MATRIX_SUBDIAG  ( 1, IC, IP, IX, NM, NL)
            CALL SCARC_SETUP_MATRIX_SUBDIAG  (-1, IC, IP, IX, NM, NL)
            CALL SCARC_SETUP_MATRIX_SUBDIAG  (-3, IC, IP, IZ, NM, NL)

         ENDDO
      ENDDO
   
   !! --------------------- 3D -----------------
   CASE (NSCARC_DIMENSION_THREE)
   
      ALLOCATE (SL%A(SL%NAE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%A', IERR)
      SL%A = 0.0_EB
      
      ALLOCATE (SL%A_COL(SL%NAE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SL%A_COL', IERR)
      SL%A_COL = 0.0_EB


      !< Compute single matrix entries and corresponding row and column pointers
      !< Along internal boundaries use placeholders for the neighboring matrix entries
      !< which will be communicated in a following step
      IP  = 1
      DO IZ = 1, SL%NZ
         DO IY = 1, SL%NY
            DO IX = 1, SL%NX
       
               IC = (IZ-1) * SL%NX * SL%NY + (IY-1) * SL%NX + IX
       
               CALL SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)
      
               CALL SCARC_SETUP_MATRIX_SUBDIAG  ( 3, IC, IP, IZ, NM, NL)
               CALL SCARC_SETUP_MATRIX_SUBDIAG  ( 2, IC, IP, IY, NM, NL)
               CALL SCARC_SETUP_MATRIX_SUBDIAG  ( 1, IC, IP, IX, NM, NL)
               CALL SCARC_SETUP_MATRIX_SUBDIAG  (-1, IC, IP, IX, NM, NL)
               CALL SCARC_SETUP_MATRIX_SUBDIAG  (-2, IC, IP, IY, NM, NL)
               CALL SCARC_SETUP_MATRIX_SUBDIAG  (-3, IC, IP, IZ, NM, NL)

            ENDDO
         ENDDO
      ENDDO
    
END SELECT
SL%A_ROW(SL%NC+1) = IP
SL%NA             = IP -1                                   !< set correct number of matrix entries


IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'SETUP_MATRIX: SL%NC=',SL%NC
   WRITE(SCARC_LU,*) 'SETUP_MATRIX: SL%NA=',SL%NA
   WRITE(SCARC_LU,*) 'SETUP_MATRIX: SIZE(SL%A)=',SIZE(SL%A)
   WRITE(SCARC_LU,*) 'SETUP_MATRIX: SIZE(SL%A_ROW)=',SIZE(SL%A_ROW)
   WRITE(SCARC_LU,*) 'SETUP_MATRIX: SIZE(SL%A_COL)=',SIZE(SL%A_COL)
ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX


!> ------------------------------------------------------------------------------------------------
!> Set main diagonal entry for matrix
!> In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil 
!> If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)
INTEGER, INTENT(IN)  :: IC, IX, IY, IZ
INTEGER, INTENT(IN)  :: NM, NL
INTEGER, INTENT(INOUT) :: IP
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

SL%A(IP) = SL%A(IP) - 2.0_EB/(SL%DXL(IX-1)*SL%DXL(IX))

IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) &
   SL%A(IP) = SL%A(IP) - 2.0_EB/(SL%DYL(IY-1)*SL%DYL(IY))

SL%A(IP) = SL%A(IP) - 2.0_EB/(SL%DZL(IZ-1)*SL%DZL(IZ))
            
SL%A_ROW(IC) = IP                                     
SL%A_COL(IP) = IC                                     
IP = IP + 1

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,'(a)') '------------- SETUP_MATRIX_MAINDIAG:'
   WRITE(SCARC_LU,'(5(a,i4))') 'IC=',IC,': IP=',IP-1,': IX=',IX,': IY=',IY,': IZ=',IZ
   WRITE(SCARC_LU,'(a, 10f12.6)') 'DXL=',SL%DXL
   WRITE(SCARC_LU,'(a, 10f12.6)') 'DYL=',SL%DYL
   WRITE(SCARC_LU,'(a, 10f12.6)') 'DZL=',SL%DZL
   WRITE(SCARC_LU,'(a,i4,a,f12.6)') 'A(',IP-1,')    =',SL%A(IP-1)
   WRITE(SCARC_LU,'(a,i4,a,i4)') 'A_ROW(',IC,')=',SL%A_ROW(IC)
   WRITE(SCARC_LU,'(a,i4,a,i4)') 'A_COL(',IP-1,')=',SL%A_COL(IP-1)
ENDIF
   
END SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG


!> ------------------------------------------------------------------------------------------------
!> Set subdiagonal entry for matrix
!> In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil 
!> If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG (IOR0, IC, IP, IPTR, NM, NL)
INTEGER, INTENT(IN)  :: IC, IOR0, IPTR
INTEGER, INTENT(IN)  :: NM, NL
INTEGER, INTENT(INOUT) :: IP
INTEGER  :: IW
REAL(EB) :: DSCAL, DH1, DH2
LOGICAL  :: NO_BDRY_CELL
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

DH1 = SL%FACE(IOR0)%DH(IPTR-1)
DH2 = SL%FACE(IOR0)%DH(IPTR)


!< set bounds and step sizes depending on orientation of face (lower or upper subdiagonal)
IF (IOR0 > 0) THEN
   NO_BDRY_CELL = IPTR > 1
   DSCAL = 2.0_EB/(DH1*(DH1+DH2))
ELSE
   NO_BDRY_CELL = IPTR < SL%FACE(IOR0)%NFC 
   DSCAL = 2.0_EB/(DH2*(DH1+DH2))
ENDIF

!< if IC is an internal cell of the mesh, compute normal matrix contribution for corresponding subdiagonal
IF (NO_BDRY_CELL) THEN

   SL%A(IP)     = SL%A(IP) + DSCAL
   SL%A_COL(IP) = IC + SL%FACE(IOR0)%IOFFSET_MATRIX
   IP = IP + 1

!< if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
ELSE IF (SL%FACE(IOR0)%NUM_NEIGHBORS /= 0) THEN

   IF (CELL_WITH_NEIGHBOR(SL%WALL, IC, IW, SL%FACE(IOR0)%IWG_MARKER, SL%FACE(IOR0)%NFW)) THEN
      SL%A(IP)     = SL%A(IP) + DSCAL
      !SL%A_COL(IP) = SL%WALL(IW)%ICE(1)
      SL%A_COL(IP) = SL%WALL(IW)%ICO
      IP = IP + 1
   ENDIF  

ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,'(a)') '------------- SETUP_MATRIX_SUBDIAG:'
   WRITE(SCARC_LU,'(4(a,i4))') 'IC=',IC,': IP=',IP-1,': IPTR=',IPTR
   WRITE(SCARC_LU,'(a, 10f12.6)') 'DH=',SL%FACE(IOR0)%DH 
   WRITE(SCARC_LU,*) 'NO_BDRY_CELL=',NO_BDRY_CELL
   WRITE(SCARC_LU,*) 'DSCAL = ', DSCAL
   WRITE(SCARC_LU,*) 'DH1 = ', DH1
   WRITE(SCARC_LU,*) 'DH2 = ', DH2
   WRITE(SCARC_LU,*) 'IOR0 = ', IOR0
   WRITE(SCARC_LU,*) 'SL%FACE(IOR0)%NUM_NEIGHBORS=',SL%FACE(IOR0)%NUM_NEIGHBORS
   WRITE(SCARC_LU,*) 'SL%FACE(IOR0)%IWG_MARKER=',SL%FACE(IOR0)%IWG_MARKER
   WRITE(SCARC_LU,*) 'SL%FACE(IOR0)%IOFFSET_MATRIX=',SL%FACE(IOR0)%IOFFSET_MATRIX
   WRITE(SCARC_LU,*) 'SL%FACE(IOR0)%NFW=',SL%FACE(IOR0)%NFW
   WRITE(SCARC_LU,*) 'SL%FACE(IOR0)%NFC=',SL%FACE(IOR0)%NFC
   WRITE(SCARC_LU,'(a,i4,a,f12.6)') 'A(',IP-1,')    =',SL%A(IP-1)
   WRITE(SCARC_LU,'(a,i4,a,i4)') 'A_ROW(',IC,')=',SL%A_ROW(IC)
   WRITE(SCARC_LU,'(a,i4,a,i4)') 'A_COL(',IP-1,')=',SL%A_COL(IP-1)
ENDIF


END SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG


!> ------------------------------------------------------------------------------------------------
!> Determine if cell IC has a neighbor and, if yes, save corresponding IW-value
!> ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION CELL_WITH_NEIGHBOR(WALL0, IC, IW, IW0, IL0)
TYPE (SCARC_WALL_TYPE), DIMENSION(:), INTENT(IN):: WALL0
INTEGER, INTENT(IN)    :: IC, IW0, IL0
INTEGER, INTENT(INOUT) :: IW

CELL_WITH_NEIGHBOR = .FALSE.

SEARCH_WALLCELLS_LOOP: DO IW = IW0, IW0+IL0-1
  IF (IC == WALL0(IW)%ICW .AND. WALL0(IW)%NOM /= 0) THEN
     CELL_WITH_NEIGHBOR = .TRUE.
     EXIT SEARCH_WALLCELLS_LOOP
  ENDIF
ENDDO SEARCH_WALLCELLS_LOOP

RETURN
END FUNCTION CELL_WITH_NEIGHBOR

!> ----------------------------------------------------------------------------------------------------
!> Extract overlapping matrix data
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OVERLAPS (NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM, LAST_NOM
INTEGER :: IERR, IW, IWC, IG, IC, JC, JCA, ICW, ICN, ICE, ICE0, ICE2, ICG, ICG2
INTEGER :: ICOL, ICCE, ICPL, NCOL, NCOL2
INTEGER :: ICC, ICOLG, ICOLE, IOR0, NCPL, NCE0
INTEGER, POINTER, DIMENSION(:,:) :: AUX
TYPE (SCARC_TYPE)         , POINTER :: S
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SLF, SLC
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSLF, OSLC

IERR = 0
SELECT_TYPE: SELECT CASE (NTYPE)

   !! --------------------------------------------------------------------------------------------
   !! set stencil sizes on overlapping parts
   !! --------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_STENCIL)

      STENCIL_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

         S => SCARC(NM)
         SLF => S%LEVEL(NL)

         STENCIL_EXTERNAL_CELL_LOOP: DO ICE = SLF%NC+1, SLF%NCE
      
           IW = ABS(SLF%ICG_TO_IWG(ICE))

           IF (IW > 0) THEN
              ICG = SLF%WALL(IW)%ICG(1)
           ELSE
              IW  = ABS(IW)
              ICG = SLF%WALL(ABS(IW))%ICG2(1)
           ENDIF

           NOM = SLF%WALL(IW)%NOM
           OSLF => S%OSCARC(NOM)%LEVEL(NL)

           SLF%A_ROW(ICE+1) = SLF%A_ROW(ICE) + OSLF%A_SIZE(ICG)
          
         ENDDO STENCIL_EXTERNAL_CELL_LOOP

      ENDDO STENCIL_MESHES_LOOP

   !! --------------------------------------------------------------------------------------------
   !! after exchange of system matrix
   !! --------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SYSTEM)

      SYSTEM_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
         
         S => SCARC(NM)
         SLF => S%LEVEL(NL)
         
         NCE0 = SLF%NCE+1
         DO IC = SLF%NC+1, SLF%NCE

            IW = SLF%ICG_TO_IWG(IC)
            IF (IW > 0) THEN
               ICN = SLF%WALL(IW)%ICN(1)
               ICE = SLF%WALL(IW)%ICE(1)
            ELSE
               IW = ABS(IW)
               ICN = SLF%WALL(IW)%ICN2(1)
               ICE = SLF%WALL(IW)%ICE2(1)
            ENDIF

            NOM  = SLF%WALL(IW)%NOM
            OSLF => S%OSCARC(NOM)%LEVEL(NL)

            ICOL = SLF%A_ROW(ICE)
            JC   = ABS(SLF%A_COL(ICOL))

            OSLF%ICN_TO_ICE(JC) = ICE

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'NOM=',NOM,': ICN_TO_ICE(',JC,')=',OSLF%ICN_TO_ICE(JC)

         ENDDO

         ALLOCATE (AUX(5,SLF%NC+1:SLF%NCE), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'AUX', IERR)
         AUX = 0

         DO IC = SLF%NC+1, SLF%NCE

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'SETUP_OVERLAPS FOR SYSTEM_MATRIX'
   WRITE(SCARC_LU,*) 'NC  = ', SLF%NC
   WRITE(SCARC_LU,*) 'NCE = ', SLF%NCE
   WRITE(SCARC_LU,*) '------------------ IC=',IC
ENDIF
            IW = SLF%ICG_TO_IWG(IC)
         WRITE(*,*) 'E: HIER WALL%ICN ändern auf OSLF!'
            IF (IW > 0) THEN
         !<   ICN = SLF%WALL(IW)%ICN(1)
            ELSE
              IW = ABS(IW)
         !<   ICN = SLF%WALL(IW)%ICN2(1)
            ENDIF

            NOM = SLF%WALL(IW)%NOM
 
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'IW =',IW
   WRITE(SCARC_LU,*) 'ICN=',ICN
ENDIF
            OSLF => S%OSCARC(NOM)%LEVEL(NL)
            SLF%CELL_MAP(1, IC) = NOM

          !< Caution: length must possibly be increased 

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,10i4)') '============= OVERLAPS 0: A_COL:',&
     (SLF%A_COL(ICOL), ICOL = SLF%A_ROW(IC), SLF%A_ROW(IC+1)-1)

            ICPL = 0
            DO ICOL = SLF%A_ROW(IC),SLF%A_ROW(IC+1)-1
               JC  = SLF%A_COL(ICOL)
               JCA = ABS(SLF%A_COL(ICOL))
               IF (JC < 0 .AND. JCA <= OSLF%NC) THEN
                  IF (OSLF%ICN_TO_ICE(JCA)==0) THEN
                     OSLF%ICN_TO_ICE(JCA) = NCE0
                     SLF%CELL_MAP(2, IC) = NCE0
                     NCE0 = NCE0 + 1
                  ELSE
                     SLF%CELL_MAP(2, IC) = OSLF%ICN_TO_ICE(JCA)
                  ENDIF
                  SLF%CELL_MAP(3, IC) = JC
                  SLF%A_COL(ICOL) = SLF%CELL_MAP(2,IC)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,3i3,a,i3)') &
    'OVERLAPS A: ICOL=',ICOL,': IC=',IC,': NOM =',NOM,': ICN_TO_ICE(',JCA,')=',&
     OSLF%ICN_TO_ICE(JCA),': CELL_MAP=',&
     SLF%CELL_MAP(1,IC),SLF%CELL_MAP(2,IC),SLF%CELL_MAP(2,IC), ': NCE0=',NCE0

               ELSE
                  SLF%A_COL(ICOL) = JC
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,2i3)') &
   'OVERLAPS B: ICOL=',ICOL,': IC=',IC,': NOM =',NOM,': A_COL(',ICOL,') =',SLF%A_COL(ICOL)

               ENDIF
            ENDDO
         
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,10i4)') '============= OVERLAPS C: A_COL:',&
     (SLF%A_COL(ICOL), ICOL = SLF%A_ROW(IC), SLF%A_ROW(IC+1)-1)
         ENDDO

         !SLF%NCE2 = NCE0-1              !< needed where ????? OUTDATED ???
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,*) 'OVERLAPS C: FINAL NCE = ', SLF%NCE



         DEALLOCATE (AUX)

      ENDDO SYSTEM_MESHES_LOOP

   !! --------------------------------------------------------------------------------------------
   !! after exchange of prolongation matrix
   !! --------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_PROLONGATION)
   
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL 000')

      PROLONGATION_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
         
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
            NCOL  = OSLF%P_ROW(ICG +1)-OSLF%P_ROW(ICG)         !< wozu wird das nochmal gebraucht ??
            !IF (TYPE_LAYER == NSCARC_LAYER_TWO) THEN
               !ICE2  = SLF%WALL(IW)%ICE2(1)
               !ICG2  = SLF%WALL(IW)%ICG2(1)
               NCOL2 = OSLF%P_ROW(ICG2+1)-OSLF%P_ROW(ICG2)     !< wozu wird das nochmal gebraucht ??
            !ENDIF

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) '============= SETUP_OVERLAPS: IW=',IW
   WRITE(SCARC_LU,*) 'ICW  =',ICW
   WRITE(SCARC_LU,*) 'ICC  =',ICC
   WRITE(SCARC_LU,*) 'ICE  =',ICE
   WRITE(SCARC_LU,*) 'ICG  =',ICG
   WRITE(SCARC_LU,*) 'NCOL =',NCOL
ENDIF
            IF (ICC >= NSCARC_CELL_TYPE_COARSE) THEN
               IWC = IWC + 1
               SLC%INTERNAL_BDRY_CELL(ICC) = NSCARC_LAYER_ONE
               SLC%WALL_INDEX(ICC, IOR0)   = IWC
               SLC%WALL(IWC)%NOM = NOM
               SLC%WALL(IWC)%ICW = ICC
               SLC%WALL(IWC)%IOR = IOR0
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'A: SLC%INTERNAL_BDRY_CELL(',ICC,')=',SLC%INTERNAL_BDRY_CELL(ICC)
   WRITE(SCARC_LU,*) 'A: SLC%WALL_INDEX(',ICC,',', IOR0,')=',SLC%WALL_INDEX(ICC, IOR0)
   WRITE(SCARC_LU,*) 'A: SLC%WALL(',IWC,')%NOM=',SLC%WALL(IWC)%NOM 
   WRITE(SCARC_LU,*) 'A: SLC%WALL(',IWC,')%ICW=',SLC%WALL(IWC)%ICW
   WRITE(SCARC_LU,*) 'A: SLC%WALL(',IWC,')%IOR=',SLC%WALL(IWC)%IOR 
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
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'B: ICCE=',ICCE
   WRITE(SCARC_LU,*) 'B: ICOL=',ICOL
   WRITE(SCARC_LU,*) 'B: OSLF%ICG0=',OSLF%ICG0
   WRITE(SCARC_LU,*) 'B: SLF%CELL_TYPE(',ICE,')=',SLF%CELL_TYPE(ICE)
   WRITE(SCARC_LU,*) 'B: SLC%EXT_PTR(',ICCE,')=',OSLF%CELL_TYPE(ICG)
   WRITE(SCARC_LU,*) 'B: OSLC%ICN_TO_ICE(',OSLF%P_COL(ICOL),')=',OSLC%ICN_TO_ICE(OSLF%P_COL(ICOL))
   !WRITE(SCARC_LU,*) 'OSLC%ICG_TO_ICE(',OSLF%ICG0,')=',OSLC%ICG_TO_ICE(OSLF%ICG0)
ENDIF
               !OSLC%ICG_TO_ICE(OSLF%ICG0) = ICCE
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

         IF (TYPE_DEBUG >= NSCARC_DEBUG_NONE) THEN
            DO IC = 1, SLC%NC
               WRITE(SCARC_LU,*) 'SLC%INTERNAL_BDRY_CELL(',IC,')=',SLC%INTERNAL_BDRY_CELL(IC)
            ENDDO
            DO IC = 1, SLC%NC
               WRITE(SCARC_LU,'(a,i3,a,8i4)') 'SLC%WALL_INDEX(',IC,', -3:3)=',SLC%WALL_INDEX(IC,-3:3)
            ENDDO
         ENDIF

         !< Replace negative cell numbers of neighboring cells for internal wall cells
         !< Note that this only holds true for fine grid cells 
         !< (weights of coarse cells don't reach in neighboring mesh)
         IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) THEN

         PROLONGATION_WALL_LOOP2: DO IW = 1, SLF%NW

            NOM  = SLF%WALL(IW)%NOM
            IF (NOM == 0) CYCLE PROLONGATION_WALL_LOOP2

            OSLF => S%OSCARC(NOM)%LEVEL(NL)

            ICW  = SLF%WALL(IW)%ICW
            ICC  = SLF%CELL_TYPE(ICW)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') &
      'B0: ============= : IW=',IW,': NOM=',NOM,': ICW=',ICW,': ICC=',ICC
ENDIF

            IF (ICC < NSCARC_CELL_TYPE_COARSE) THEN

               DO ICOL = SLF%P_ROW(ICW), SLF%P_ROW(ICW+1)-1
                  JC   = SLF%P_COL(ICOL)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') &
      'B1: ============= : ICC=',ICC,': JC=',JC
ENDIF
                  !< Additionally identify coarse cells from second layer 
                  !< adjacent to internal boundary 
                  IF (JC < 0) THEN
                     IC = OSLF%ICN_TO_ICE(ABS(JC))
                     SLF%P_COL(ICOL) = SLF%CELL_TYPE(IC)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
      'B2: ============= : IC=',IC,': SLF%P_COL(',ICOL,')=',SLF%P_COL(ICOL)
ENDIF
                  ELSE IF (JC <= SLF%NCC) THEN
                     IF (SLC%INTERNAL_BDRY_CELL(JC) > NSCARC_LAYER_ONE) THEN
                        IWC = IWC + 1
                        SLC%INTERNAL_BDRY_CELL(JC) = NSCARC_LAYER_TWO
                        SLC%WALL_INDEX(ICC, IOR0)  = -IWC
                        SLC%WALL(IWC)%NOM = NOM
                        SLC%WALL(IWC)%ICW = JC
                        SLC%WALL(IWC)%IOR = SLF%WALL(IW)%IOR

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'B3: IWC=',IWC
   WRITE(SCARC_LU,*) 'B3: SLC%INTERNAL_BDRY_CELL(',JC,')=',SLC%INTERNAL_BDRY_CELL(JC)
   WRITE(SCARC_LU,*) 'B3: SLC%WALL(',IWC,')%NOM=',NOM
   WRITE(SCARC_LU,*) 'B3: SLC%WALL(',IWC,')%ICW=',JC
   WRITE(SCARC_LU,*) 'B3: SLC%WALL(',IWC,')%IOR=',SLF%WALL(IWC)%IOR 
ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO PROLONGATION_WALL_LOOP2

   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL 001')
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION', 'FINAL CELL_TYPE 001')

         PROLONGATION_ECELL_LOOP: DO ICE0 = SLF%NC+1, SLF%NCE

            IW  = SLF%ICG_TO_IWG(ICE0)
            NOM = SLF%WALL(ABS(IW))%NOM
            OSLF => S%OSCARC(NOM)%LEVEL(NL)
            OSLC => S%OSCARC(NOM)%LEVEL(NL+1)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'C: =======================  IW=',IW, ': NOM=',NOM
   WRITE(SCARC_LU,'(a,5i4)') 'OSLF%P_ROW=',(OSLF%P_ROW(ICG), ICG=1,5)
   WRITE(SCARC_LU,'(a,10i4)') 'OSLF%P_COL=',(OSLF%P_COL(ICG), ICG=1,10)
ENDIF

            !< Replace positive cell numbers of neighboring cells in first layer
         WRITE(*,*) 'G: HIER WALL%ICN ändern auf OSLF!'
            IF (IW > 0) THEN
               !ICG = SLF%WALL(IW)%ICG(1)
               !ICE = SLF%WALL(IW)%ICE(1)
               ICOLE = SLF%P_ROW(ICE)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,'(a,i4,a,i4,a,i4)') 'C00: ICG=',ICG,': ICE=',ICE,': ICOLE=',ICOLE
ENDIF
               DO ICOLG = OSLF%P_ROW(ICG), OSLF%P_ROW(ICG+1)-1
                  SLF%P_COL(ICOLE) = OSLF%P_COL(ICOLG)
                  SLF%P_COL(ICOLE) = OSLF%ICN_TO_ICE(OSLF%P_COL(ICOLG))
                  SLF%P(ICOLE) = OSLF%P(ICOLG)
                  ICOLE = ICOLE + 1

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,'(a,i4,a,i4,a,i4,a,i4)') 'C01: ICOLG=',ICOLG,': JC=',JC,': OSLF%NCC=',OSLF%NCC
ENDIF
!>                 IF (JC > 0 .AND. JC <= OSLF%NCC) THEN
!>                    IC = OSLC%ICN_TO_ICE(JC)
!>                    !OSLF%P_COL(ICOLG) = IC
!>                    SLF%P_COL(ICOLE)  = IC
!>
!IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
!>  WRITE(SCARC_LU,'(a,i4,a,i4,a,i4,a,i4)') &
!>       'C10: OSLF%P_COL(',ICOLG,')=',IC,': SLF%P_COL(',ICOLE,')=',SLF%P_COL(ICOLE)
!ENDIF
!>                 ELSE
!>                    !OSLF%P_COL(ICOLG) = SLF%CELL_TYPE(ABS(JC))
!>                    SLF%P_COL(ICOLE) = SLF%CELL_TYPE(ABS(JC))
!>
!IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
!>  WRITE(SCARC_LU,'(a,i4,a,i4,a,i4,a,i4)') &
!>       'C20: OSLF%P_COL(',ICOLG,')=',IC,': SLF%P_COL(',ICOLE,')=',SLF%P_COL(ICOLE)
!ENDIF
!>                 ENDIF
!>                 SLF%P(ICOLE) = OSLF%P(ICOLG)
!>                 ICOLE = ICOLE + 1
               ENDDO
               SLF%P_ROW(ICE+1) = ICOLE
   
            !< Replace positive cell numbers of neighboring cells in second layer (if requested)
            ELSE
               IW = ABS(IW)
         WRITE(*,*) 'D: HIER WALL%ICN ändern auf OSLF!'
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
     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SLC%NW=',SLC%NW

         ENDIF     !GMG-ENDIF

      ENDDO PROLONGATION_MESHES_LOOP

 !! --------------------------------------------------------------------------------------------
 !! after exchange of GMG-matrix
 !! --------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_GMG)
   
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'CALLING SETUP_OVERLAP_GMG on level ', NL
      GMG_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
         
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '======================================='
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '============ NM=',NM
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '======================================='
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

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '============ IW=',IW,': NOM=',NOM, '; IG=',IG
            OSLF => S%OSCARC(NOM)%LEVEL(NL)
            OSLC => S%OSCARC(NOM)%LEVEL(NL+1)
    
            IOR0 = SLF%WALL(IW)%IOR
            ICW  = SLF%WALL(IW)%ICW
            ICC  = SLF%CELL_TYPE(ICW)

         WRITE(*,*) 'I: HIER WALL%ICN ändern auf OSLF!'
            !ICE   = SLF%WALL(IW)%ICE(1)
            !ICG   = SLF%WALL(IW)%ICG(1)
            NCOL  = OSLF%P_ROW(ICG +1)-OSLF%P_ROW(ICG)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'SETUP_OVERLAPS:'
   WRITE(SCARC_LU,*) 'ICW  =',ICW
   WRITE(SCARC_LU,*) 'ICC  =',ICC
   WRITE(SCARC_LU,*) 'ICE  =',ICE
   WRITE(SCARC_LU,*) 'ICG  =',ICG
   WRITE(SCARC_LU,*) 'NCOL =',NCOL
ENDIF
            IF (ICC >= NSCARC_CELL_TYPE_COARSE) THEN

               IWC  = IWC + 1
               ICCE = ICCE + 1
               NCPL = 1

               SLC%INTERNAL_BDRY_CELL(ICC) = NSCARC_LAYER_ONE
               SLC%WALL_INDEX(ICC, IOR0)   = IWC

         WRITE(*,*) 'J: HIER WALL%ICN ändern auf OSLF!'
               SLC%WALL(IWC)%NOM  = NOM
               !SLC%WALL(IWC)%ICW  = ICC
               SLC%WALL(IWC)%IOR  = IOR0
               !SLC%WALL(IWC)%NCPL = NCPL

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'SLC%WALL(',IWC,')%NOM=',SLC%WALL(IWC)%NOM 
   WRITE(SCARC_LU,*) 'SLC%WALL(',IWC,')%ICW=',SLC%WALL(IWC)%ICW
   WRITE(SCARC_LU,*) 'SLC%WALL(',IWC,')%IOR=',SLC%WALL(IWC)%IOR , SLC%INTERNAL_BDRY_CELL(ICC)
ENDIF

               !IF (SLF%CELL_TYPE(ICE) < NSCARC_CELL_TYPE_COARSE) THEN
             !<  WRITE(*,*) 'Wrong cominbation of coarse cells in OVERLAP_GMG, stop!'
             !<  STOP
               !NDIF

         WRITE(*,*) 'H: HIER WALL%ICN ändern auf OSLF!'
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
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'ICCE=',ICCE
   WRITE(SCARC_LU,*) 'SLF%CELL_TYPE(',ICE,')=',SLF%CELL_TYPE(ICE)
   WRITE(SCARC_LU,*) 'SLC%EXT_PTR(',ICCE,')=',SLC%EXT_PTR(ICCE)
   WRITE(SCARC_LU,*) 'OSLC%ICN_TO_ICE(',OSLF%P_COL(ICOL),')=',OSLC%ICN_TO_ICE(OSLF%P_COL(ICOL))
   WRITE(SCARC_LU,*) 'OSLC%ICG_TO_ICE(',IG,')=',OSLC%ICG_TO_ICE(IG)
   !WRITE(SCARC_LU,*) 'SLC%WALL(',IWC,')%ICE(1)=',SLC%WALL(IWC)%ICE(1)
   !WRITE(SCARC_LU,*) 'SLC%WALL(',IWC,')%ICG(1)=',SLC%WALL(IWC)%ICG(1)
   !WRITE(SCARC_LU,*) 'SLC%WALL(',IWC,')%ICN(1)=',SLC%WALL(IWC)%ICN(1)
ENDIF
            ENDIF

         ENDDO GMG_WALL_LOOP1

         IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
            DO IC = 1, SLC%NC
               WRITE(SCARC_LU,*) 'SLC%INTERNAL_BDRY_CELL(',IC,')=',SLC%INTERNAL_BDRY_CELL(IC)
            ENDDO
            DO IC = 1, SLC%NC
               WRITE(SCARC_LU,'(a,i3,a,8i4)') 'SLC%WALL_INDEX(',IC,', -3:3)=',SLC%WALL_INDEX(IC,-3:3)
            ENDDO
         ENDIF

      ENDDO GMG_MESHES_LOOP

END SELECT SELECT_TYPE

END SUBROUTINE SCARC_SETUP_OVERLAPS


!> ------------------------------------------------------------------------------------------------
!> Set pointer for different structures on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP
REAL(EB) :: DBC
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL


SL => SCARC(NM)%LEVEL(NL)

SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)

   !! ------------------ 2D -----------------------
   CASE (NSCARC_DIMENSION_TWO)

      WALLCELLS_LOOP2D: DO IW = 1, SL%NW
       
         IOR0 = SL%WALL(IW)%IOR
         IF (ABS(IOR0) == 2) CYCLE          !< 2D: in case of y-boundary cycle
   
         I    = SL%WALL(IW)%IXW
         K    = SL%WALL(IW)%IZW

         NOM  = SL%WALL(IW)%NOM
   
         SL%WALL(IW)%ICW = (K-1)*SL%NX + I

         SELECT CASE (IOR0)
            CASE (1)
               IC = (K-1) * SL%NX + I
               DBC= SL%DXI2
            CASE (-1)
               IC = K * SL%NX
               DBC= SL%DXI2
            CASE (3)
               IC = I
               DBC= SL%DZI2
            CASE (-3)
               IC = (SL%NZ-1) * SL%NX + I
               DBC= SL%DZI2
         END SELECT
   
         IP = SL%A_ROW(IC)
         SELECT CASE (SL%WALL(IW)%BTYPE)
            CASE (DIRICHLET)                      !< set Dirichlet BC's along open boundary cells
               SL%A(IP) = SL%A(IP) - DBC
            !CASE (INTERNAL)                      !< do nothing along internal boundaries (only debugging)
            CASE (NEUMANN)                        !< set Neumann BC's at all other nodes
               SL%A(IP) = SL%A(IP) + DBC
         END SELECT
   
      ENDDO WALLCELLS_LOOP2D
   

   !! ------------------ 2D -----------------------
   CASE (NSCARC_DIMENSION_THREE)

      WALLCELLS_LOOP3D: DO IW = 1, SL%NW
   
         IOR0 = SL%WALL(IW)%IOR

         I    = SL%WALL(IW)%IXW
         J    = SL%WALL(IW)%IYW
         K    = SL%WALL(IW)%IZW

         NOM  = SL%WALL(IW)%NOM
   
         SL%WALL(IW)%ICW = (K-1)*SL%NX*SL%NY + (J-1)*SL%NX + I

         SELECT CASE (IOR0)
            CASE (1)
               IC = (K-1) * SL%NX * SL%NY + (J-1) * SL%NX + I
               DBC= SL%DXI2
            CASE (-1)
               IC = (K-1) * SL%NX * SL%NY + J * SL%NX 
               DBC= SL%DXI2
            CASE (2)
               IC = (K-1) * SL%NX * SL%NY + I
               DBC= SL%DYI2
            CASE (-2)
               IC = (K-1) * SL%NX * SL%NY + (SL%NY-1) * SL%NX + I
               DBC= SL%DYI2
            CASE (3)
               IC = (J-1) * SL%NX + I
               DBC= SL%DZI2
            CASE (-3)
               IC = (SL%NZ-1) * SL%NX * SL%NY + (J-1) * SL%NX + I
               DBC= SL%DZI2
         END SELECT
   
         IP = SL%A_ROW(IC)
         SELECT CASE (SL%WALL(IW)%BTYPE)
            CASE (DIRICHLET)                      !< set Dirichlet BC's at open and null boundary cells
               SL%A(IP) = SL%A(IP) - DBC
            !CASE (INTERNAL)                      !< do nothing along internal boundaries (only debugging)
            CASE (NEUMANN)                        !< set Neumann BC's at all other cells
               SL%A(IP) = SL%A(IP) + DBC
         END SELECT
   
      ENDDO WALLCELLS_LOOP3D
     
END SELECT SELECT_DIMENSION


END SUBROUTINE SCARC_SETUP_BOUNDARY


REAL(EB) FUNCTION SET_EDGE_BOUNDARY(ITYPE, SCAL)
INTEGER, INTENT(IN):: ITYPE
REAL(EB), INTENT(IN):: SCAL
SELECT CASE (ITYPE)
   CASE (DIRICHLET) 
       SET_EDGE_BOUNDARY = 5.0_EB*SCAL
   CASE (NEUMANN) 
       SET_EDGE_BOUNDARY = 3.0_EB*SCAL
   CASE DEFAULT
      WRITE(*,*) 'CONFUSION IN SET_EDGE_BOUNDARY'
      stop
END SELECT

RETURN
END FUNCTION SET_EDGE_BOUNDARY


REAL(EB) FUNCTION SET_VERTEX_BOUNDARY(ITYPE1, ITYPE2, SCAL)
INTEGER, INTENT(IN):: ITYPE1, ITYPE2
REAL(EB), INTENT(IN):: SCAL
SELECT CASE (ITYPE1)
   CASE (DIRICHLET) 
      SELECT CASE(ITYPE2)
         CASE (DIRICHLET) 
            SET_VERTEX_BOUNDARY = 6.0_EB*SCAL
         CASE (NEUMANN) 
            SET_VERTEX_BOUNDARY = 4.0_EB*SCAL
         CASE (INTERNAL) 
            SET_VERTEX_BOUNDARY = 5.0_EB*SCAL
         CASE DEFAULT
            WRITE(*,*) 'CONFUSION IN SET_VERTEX_BOUNDARY, DIRICHLET'
            stop
      END SELECT
   CASE (NEUMANN) 
      SELECT CASE(ITYPE2)
         CASE (DIRICHLET) 
            SET_VERTEX_BOUNDARY = 4.0_EB*SCAL
         CASE (NEUMANN) 
            SET_VERTEX_BOUNDARY = 2.0_EB*SCAL
         CASE (INTERNAL) 
            SET_VERTEX_BOUNDARY = 3.0_EB*SCAL
         CASE DEFAULT
            WRITE(*,*) 'CONFUSION IN SET_VERTEX_BOUNDARY, NEUMANN'
            stop
      END SELECT
END SELECT
RETURN 
END FUNCTION SET_VERTEX_BOUNDARY



!> ----------------------------------------------------------------------------------------------------
!> Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_VECTORS
INTEGER :: IERR
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL
TYPE (SCARC_PRECON_TYPE) , POINTER :: SP
INTEGER :: NM, NL

IERR = 0

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

      SL => SCARC(NM)%LEVEL(NL)
   
      SELECT_METHOD: SELECT CASE (TRIM(SCARC_METHOD))
   
         !< working and auxiliary vectors for global CG/BICG-method
         CASE ('KRYLOV')
   
            ALLOCATE (SL%X(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'X', IERR)
            SL%X = 0.0_EB
       
            ALLOCATE (SL%F(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'B', IERR)
            SL%F = 0.0_EB
       
            ALLOCATE (SL%D(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'D', IERR)
            SL%D = 0.0_EB
   
            ALLOCATE (SL%W(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'W', IERR)
            SL%W = 0.0_EB
       
            ALLOCATE (SL%G(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'G', IERR)
            SL%G = 0.0_EB
       
            ALLOCATE (SL%Y(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'Y', IERR)
            SL%Y = 0.0_EB
   
            ALLOCATE (SL%Z(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'Z', IERR)
            SL%Z = 0.0_EB

            IF (TYPE_PRECON == NSCARC_PRECON_MULTIGRID) THEN

               ALLOCATE (SL%X2(SL%NCE), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'X2', IERR)
               SL%X2 = 0.0_EB
         
               ALLOCATE (SL%D2(SL%NCE), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'D2', IERR)
               SL%D2 = 0.0_EB
   
               IF (NL==NLEVEL_MAX) THEN
      
                  ALLOCATE (SL%W2(SL%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'W2', IERR)
                  SL%W2 = 0.0_EB
          
                  ALLOCATE (SL%G2(SL%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'G2', IERR)
                  SL%G2 = 0.0_EB
          
                  ALLOCATE (SL%Y2(SL%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Y2', IERR)
                  SL%Y2 = 0.0_EB
      
               ENDIF
             
            ENDIF
   
            IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
               SP => SCARC(NM)%PRECON(NL)
               ALLOCATE (SP%FFT(1:SL%NX+1, 1:SL%NY+1, 1:SL%NZ+1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'FFT', IERR)
               SP%FFT = 0.0_EB
            ENDIF
   
         !< working and auxiliary vectors for global GMG/AMG-method
         CASE ('MULTIGRID')
   
            ALLOCATE (SL%X(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'X', IERR)
            SL%X = 0.0_EB
       
            ALLOCATE (SL%F(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'F', IERR)
            SL%F = 0.0_EB
       
            ALLOCATE (SL%D(SL%NCE), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'D', IERR)
            SL%D = 0.0_EB
   
            IF (NL==NLEVEL_MAX) THEN
   
               ALLOCATE (SL%W(SL%NCE), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'W', IERR)
               SL%W = 0.0_EB
       
               ALLOCATE (SL%G(SL%NCE), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'G', IERR)
               SL%G = 0.0_EB
       
               ALLOCATE (SL%Y(SL%NCE), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'Y', IERR)
               SL%Y = 0.0_EB
   
            ENDIF
   
            IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
               SP => SCARC(NM)%PRECON(NL)
               ALLOCATE (SP%FFT(1:SL%NX+1, 1:SL%NY+1, 1:SL%NZ+1), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'FFT', IERR)
               SP%FFT = 0.0_EB
            ENDIF
   
      END SELECT SELECT_METHOD


   ENDDO LEVEL_LOOP
ENDDO  MESHES_LOOP

END SUBROUTINE SCARC_SETUP_VECTORS

!> ------------------------------------------------------------------------------------------------
!> Set sizes for transfer matrices
!> ------------------------------------------------------------------------------------------------
SUBROUTINE  SCARC_SETUP_SIZES(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: IW, IC, IG, ICP, ICE, IERR, ICCE
INTEGER :: NM, NOM
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: SLF, SLC
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSLF, OSLC


SELECT CASE (NTYPE)

   !! --------------------------------------------------------------------------------------------------
   !! Define sizes for system matrix A (including extended regions related to overlaps)
   !! --------------------------------------------------------------------------------------------------
   CASE (NSCARC_SIZE_MATRIX)

      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'STARTING SETUP_SIZES'

      LEVEL_SYSTEM_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SLF => SCARC(NM)%LEVEL(NL)
      
         !< Determine number of couplings in matrix stencil and matrix size in own part of mesh
         SELECT CASE (TYPE_DIMENSION)
            CASE (NSCARC_DIMENSION_TWO)
               SLF%NPOINTS = 5
            CASE (NSCARC_DIMENSION_THREE)
               SLF%NPOINTS = 7
         END SELECT
         SLF%NA = SLF%NC * SLF%NPOINTS

         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
            WRITE(SCARC_LU,*) 'NPOINTS=',SLF%NPOINTS
            WRITE(SCARC_LU,*) 'NA     =',SLF%NA
         ENDIF

         !< Determine sizes of overlapped parts for later communication with corresponding neighbors
         DO IW = 1, SLF%NW
            NOM = SLF%WALL(IW)%NOM
            IF (NOM /= 0) THEN
      
               OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
      
               IC = SLF%WALL(IW)%ICW
               OSLF%NA0 = OSLF%NA0 + SLF%NPOINTS
      
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
                  WRITE(SCARC_LU,'(4(a,i4))') 'IW=',IW,': NOM=',NOM,': IC=',IC,': NA0=',OSLF%NA0

               !IF (TYPE_LAYER == NSCARC_LAYER_TWO) THEN
               !   IC = SLF%WALL(IW)%ICW2
               !   OSLF%NA0 = OSLF%NA0 + SLF%NPOINTS
               !ENDIF

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'NA0(',NOM,')=', OSLF%NA0
            ENDIF
         ENDDO
      
      ENDDO LEVEL_SYSTEM_MESHES_LOOP
      
      !< Exchange sizes for AMG matrices and vectors
      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SIZE, NL)
      
      
      !< Determine extended sizes for extended prolongation and restriction matrix
      LEVEL_SYSTEM_MESHES_LOOP2: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SLF => SCARC(NM)%LEVEL(NL)
         SLF%NAE = SLF%NA
      
         ALLOCATE (SLF%CELL_MAP(3,SLF%NC+1:SLF%NCE), STAT=IERR)
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


      TRANSFER_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SLF => SCARC(NM)%LEVEL(NL)
         SLC => SCARC(NM)%LEVEL(NL+1)
         
         !< Determine dimensions of restriction and prolongation matrices in own mesh
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
    !<         SLF%CELL_TYPE(IC)  = SLF%NCC
            ELSE
               SLF%NCF = SLF%NCF + 1
               SLF%NP  = SLF%NP  + SLF%A_ROW(IC+1)-SLF%A_ROW(IC) - 1
               SLF%NR  = SLF%NP
            ENDIF
         ENDDO

    !<   SLF%NCCE = SLF%NCC
    !<   SLF%NCFE = SLF%NCF
    !<   DO IC = SLF%NC+1, SLF%NCE
    !<      IF (SLF%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
    !<         SLF%NCCE = SLF%NCCE + 1
    !<         SLF%CELL_TYPE(IC)  = SLF%NCCE
    !<      ELSE
    !<         SLF%NCFE = SLF%NCFE + 1
    !<      ENDIF
    !<   ENDDO

         SLF%NCCI = SLF%NCC
         SLF%NCCE = SLF%NCC
         !SLF%NCE0 = SLF%NC+1      !< really needed? OUTDATED ?
      
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'NCC =',SLF%NCC
   WRITE(SCARC_LU,*) 'NCCE=',SLF%NCCE
   WRITE(SCARC_LU,*) 'NP  =',SLF%NP 
   WRITE(SCARC_LU,*) 'NR  =',SLF%NR 
   WRITE(SCARC_LU,*) 'CELL_TYPE  ='
   WRITE(SCARC_LU,'(4i4)') (SLF%CELL_TYPE(IC), IC=1, SLF%NCE) 
ENDIF

         !< Determine number of coarse and fine cells and check correctness of computation
         IF (SLF%NCC + SLF%NCF /= SLF%NC) THEN
            WRITE(*,*) 'Error in AMG standard coarsening, N_CELLS_COARSE + N_CELLS_FINE = ', SLF%NCC + SLF%NCF, &
                       ' differs from N_CELLS = ', SLF%NC, ' on level ', NL
       !<   STOP
         ENDIF
         

         !< define variables for overlapping parts
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
      
         !< Determine sizes of overlapped parts for later communication with corresponding neighbors
         ICCE = SLF%NCC
         IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) SLF%NCW = 0

         DO IW = 1, SLF%NW
            NOM = SLF%WALL(IW)%NOM
            IF (NOM /= 0) THEN
      
               OSLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
      
         WRITE(*,*) 'K: HIER WALL%ICN ändern auf OSLF!'
               IC  = SLF%WALL(IW)%ICW
               !ICE = SLF%WALL(IW)%ICE(1)

               IF (TYPE_COARSENING<NSCARC_COARSENING_GMG.AND.SLF%CELL_TYPE(IC)>0) SLF%NCW = SLF%NCW + 1

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'CELL_TYPE(',IC,')=',SLF%CELL_TYPE(IC)

               IF (SLF%CELL_TYPE(ICE) >= NSCARC_CELL_TYPE_COARSE) THEN
                  OSLF%NP0  = OSLF%NP0  + 1
                  OSLF%ICG0  = OSLF%ICG0  + 1
                  OSLF%NR0  = OSLF%NP0
                  OSLF%NCC0 = OSLF%NCC0 + 1
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'OSLF(',NOM,')%ICG0 =',OSLF%ICG0
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'OSLF(',NOM,')%NCC0=',OSLF%NCC0
               ELSE
                  OSLF%NP0  = OSLF%NP0  + SLF%A_ROW(IC+1)-SLF%A_ROW(IC) - 1
                  OSLF%NR0  = OSLF%NP0  
                  OSLF%NCF0 = OSLF%NCF0 + 1
               ENDIF

         WRITE(*,*) 'L: HIER WALL%ICN ändern auf OSLF!, schauen wegen NLAYER2'
      
               !IF (NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) THEN
             !<   !IC   = SLF%WALL(IW)%ICW2
             !<   !ICE2 = SLF%WALL(IW)%ICE2(1)
             !<   IF (SLF%CELL_TYPE(ICE2) >= NSCARC_CELL_TYPE_COARSE) THEN
             !<      OSLF%NCC0 = OSLF%NCC0 + 1
             !<      OSLF%NCG  = OSLF%NCG  + 1
             !<      OSLF%NP0  = OSLF%NP0  + 1
             !<      OSLF%NR0  = OSLF%NP0
             !<   ELSE
             !<      OSLF%NCF0 = OSLF%NCF0 + 1
             !<      OSLF%NP0  = OSLF%NP0  + SLF%A_ROW(IC+1)-SLF%A_ROW(IC) - 1
             !<      OSLF%NR0  = OSLF%NP0
             !<   ENDIF
               !ENDIF
               OSLF%NCC=OSLF%NCC0
               OSLF%NCF=OSLF%NCF0
      
            ENDIF
         ENDDO

         !!!
         !< Determine new numbering for coarse cells in interior of mesh
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
      
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'NCC =',SLF%NCC
   WRITE(SCARC_LU,*) 'NCCE=',SLF%NCCE
   WRITE(SCARC_LU,*) 'NP  =',SLF%NP 
   WRITE(SCARC_LU,*) 'NR  =',SLF%NR 
   WRITE(SCARC_LU,*) 'CELL_TYPE  ='
   WRITE(SCARC_LU,'(4i4)') (SLF%CELL_TYPE(IC), IC=1, SLF%NCE) 
ENDIF

      ENDDO TRANSFER_MESHES_LOOP
      
      !< Exchange sizes for AMG matrices and vectors
      IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_TRANSFER_SIZE, NL)
      
      !< Determine extended sizes for extended prolongation and restriction matrix
      TRANSFER_MESHES_LOOP2: DO NM = NMESHES_MIN, NMESHES_MAX
      
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'NM =',NM

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
      
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '==============  NM=',NM,': NOM=',NOM
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'OSLF%NCC =',OSLF%NCC,': OSLF%NC=',OSLF%NC
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SLF%NPE =',SLF%NPE,': SLF%NP=',SLF%NP
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SLF%NRE =',SLF%NRE,': SLF%NR=',SLF%NR
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SLF%NCCE=',SLF%NCCE,': SLF%NCCE=',SLF%NCCE
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SLF%NCFE=',SLF%NCFE,': SLF%NCFE=',SLF%NCFE
         ENDDO TRANSFER_OTHER_MESHES_LOOP2

      ENDDO TRANSFER_MESHES_LOOP2

      !< Determine extended sizes for extended prolongation and restriction matrix
      TRANSFER_MESHES_LOOP3: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SLF => SCARC(NM)%LEVEL(NL)
         SLC => SCARC(NM)%LEVEL(NL+1)

         SLC%NC  = SLF%NCC
         SLC%NCE = SLF%NCCE
         SLC%NW  = SLF%NCW

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'SLC%NC =',SLC%NC
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'SLC%NCE=',SLC%NCE
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'SLC%NW =',SLC%NW
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'ALLOCATING SLC%WALL in length ', SLC%NW, ' to ', SLC%NW
         ALLOCATE (SLC%WALL(1:SLC%NW), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'SLC%WALL', IERR)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'ALLOCATING SLC%CELL_TYPE in length ', 1, ' to ', SLC%NCE
         ALLOCATE (SLC%CELL_TYPE(1:SLC%NCE), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'SLC%CELL_TYPE', IERR)
         SLC%CELL_TYPE = NSCARC_CELL_TYPE_NONE

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'ALLOCATING SLC%MEASURE in length ', 1, ' to ', SLC%NCE
         ALLOCATE (SLC%MEASURE(1:SLC%NCE), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'SLC%MEASURE', IERR)
         SLC%MEASURE = NSCARC_MEASURE_NONE

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'ALLOCATING SLC%INTERNAL_BDRY_CELL in length ', 1, ' to ', SLC%NC
         ALLOCATE(SLC%INTERNAL_BDRY_CELL(SLC%NC), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_WALLINFO','INTERNAL_BDRY_CELL',IERR)
         SLC%INTERNAL_BDRY_CELL = NSCARC_LAYER_NONE

         ALLOCATE(SLC%WALL_INDEX(SLC%NC, -3:3), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_WALLINFO','WALL_INDEX',IERR)
         SLC%WALL_INDEX = 0

         IF (SLC%NCE > SLC%NC) THEN
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'ALLOCATING SLC%ICG_TO_IWG in length ', SLC%NC+1, ' to ', SLC%NCE
            ALLOCATE(SLC%ICG_TO_IWG(SLC%NC+1:SLC%NCE), STAT=IERR)
            CALL ChkMemErr('SCARC_SETUP_WALLINFO','ICG_TO_IWG',IERR)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'ALLOCATING SLC%EXT_PTR in length ', SLC%NC+1, ' to ', SLC%NCE
            ALLOCATE(SLC%EXT_PTR(SLC%NC+1:SLC%NCE), STAT=IERR)
            CALL ChkMemErr('SCARC_SETUP_WALLINFO','EXT_PTR',IERR)
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

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_ACELL , NLEVEL_MIN,'SETUP_COARSENING','WALL_CELL')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NLEVEL_MIN,'SETUP_COARSENING', 'MATRIX')

!> Determine number of multigrid levels
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX-1
!LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MIN
   
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '======== STARTING COARENING LEVEL ', NL
   !! ---------------------------------------------------------------------------------------------
   !! Determine coarser meshes corresponding to requested coarsening strategy
   !!  --- allocate necessary arrays
   !!  --- setup measures of single cells
   !!  --- setup CELL_TYPEs of single cells
   !!  --- setup sizes of transformation matrices (prolongation/restriction)
   !! ---------------------------------------------------------------------------------------------
   !IF (NL == NLEVEL_MIN) THEN

   DO NM = NMESHES_MIN, NMESHES_MAX

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

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'ALLOCATING ST_COL in length ', SLF%NAE+1

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

   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPES', 'CELL_TYPE AAA')
         
   !! Set dimensions for coarser mesh and define sizes of prolongation and restriction matrices
   CALL SCARC_SETUP_SIZES (NSCARC_SIZE_TRANSFER, NL)

   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPES', 'CELL_TYPE BBB')

   !CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE,  NL)

   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPES', 'CELL_TYPE CCC')

   !! --------------------------------------------------------------------------------------------
   !! Allocate and define grid transfer matrices
   !!  NCE  : number of extended cells in fine grid 
   !!  NCCE : number of extended cells in coarse grid 
   !!  NCW  : number of wall cells in coarse grid
   !! --------------------------------------------------------------------------------------------
   DO NM = NMESHES_MIN, NMESHES_MAX

      SLF => SCARC(NM)%LEVEL(NL)            

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'ALLOCATING P in size ', SLF%NPE
   WRITE(SCARC_LU,*) 'ALLOCATING R in size ', SLF%NRE
   WRITE(SCARC_LU,*) 'CURRENT FINE   LEVEL ', NL
   WRITE(SCARC_LU,*) 'CURRENT COARSE LEVEL ', NL+1
ENDIF

      !< ---------------- allocate prolongation matrix including row and column pointers
      ALLOCATE (SLF%P(SLF%NPE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P', IERR)
      SLF%P = 0.0_EB
   
      ALLOCATE (SLF%P_COL(SLF%NPE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_COL', IERR)
      SLF%P_COL = 0.0_EB
   
      ALLOCATE (SLF%P_ROW(SLF%NCE+10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_ROW', IERR)
      SLF%P_ROW = 0.0_EB
   
      !< ---------------- allocate auxiliary arrays to mark relevant positions in A and P
      ALLOCATE (SLF%A_TAG(SLF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'A_TAG', IERR)
      SLF%A_TAG = 0
   
      ALLOCATE (SLF%P_TAG(SLF%NCCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_TAG', IERR)
      SLF%P_TAG = 0
   

      !< ---------------- allocate restriction matrix including row and column pointers
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
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R_COL', IERR)
         OSLF%P_COL = 0.0_EB
   
         ALLOCATE (OSLF%R_ROW(OSLF%NC+10), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R_ROW', IERR)
         OSLF%P_ROW = 0.0_EB
   
      ENDDO OTHER_MESH_LOOP2
   ENDDO

   !! determine prolongation and restriction matrix for the corresponding level
   IF (TYPE_COARSENING >= NSCARC_COARSENING_GMG) THEN
      CALL SCARC_SETUP_RESTRICTION(NL)
      CALL SCARC_SETUP_PROLONGATION(NL)
      CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_GMG, NL)
   ELSE

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'PART TWO - 1'

      CALL SCARC_SETUP_PROLONGATION(NL)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'PART TWO - 2'
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL PART TWO -2')

      CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_PROLONGATION, NL)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'PART TWO - 3'
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL PART TWO -3')

      CALL SCARC_MATRIX_TRANSPOSE(SLF%P, SLF%P_ROW, SLF%P_COL, SLF%R, SLF%R_ROW, SLF%R_COL, &
                                   SLF%NC, SLF%NCE, SLF%NCC, SLF%NCCE )

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'PART TWO - 4'
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL PART TWO -4')
   ENDIF

   !!!
   !! ONLY DEBUGGING - START
   !!!
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL 002')
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_RESTRICTION, NL, 'SETUP_RESTRICTION', 'FINAL REST 002')
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL+1, 'SETUP_SYSTEM', 'WALLINFO 002')

   IF (TYPE_LATEX == NSCARC_LATEX_ALL) THEN
      CALL SCARC_LATEX_GRID(NSCARC_LATEX_STAGGERED,NL)
      CALL SCARC_LATEX_GRID(NSCARC_LATEX_EQUAL,NL)
      CALL SCARC_LATEX_GRID(NSCARC_LATEX_NUMBER,NL)
   ENDIF
   !IF (TYPE_DEBUG > NSCARC_DEBUG_LESS .AND. NL>=1) THEN

   DO NM = NMESHES_MIN, NMESHES_MAX
      SLF => SCARC(NM)%LEVEL(NL)            

      CALL SCARC_LATEX_MATRIX2(SLF%A, SLF%A_ROW, SLF%A_COL, SLF%NC  , SLF%NC  , NM, NL, 'A')
      CALL SCARC_LATEX_MATRIX2(SLF%R, SLF%R_ROW, SLF%R_COL, SLF%NCCE, SLF%NCE , NM, NL, 'R')
      CALL SCARC_LATEX_MATRIX2(SLF%P, SLF%P_ROW, SLF%P_COL, SLF%NCE , SLF%NCCE, NM, NL, 'P')

      CALL SCARC_MATLAB_MATRIX(SLF%A, SLF%A_ROW, SLF%A_COL, SLF%NC  , SLF%NC  , NM, NL, 'A')
      CALL SCARC_MATLAB_MATRIX(SLF%R, SLF%R_ROW, SLF%R_COL, SLF%NCCE, SLF%NCE , NM, NL, 'R')
      CALL SCARC_MATLAB_MATRIX(SLF%P, SLF%P_ROW, SLF%P_COL, SLF%NCE , SLF%NCCE, NM, NL, 'P')
   ENDDO
   !ENDIF
   !!!
   !! ONLY DEBUGGING - END
   !!!

   !! -----------------------------------------------------------------------------------------
   !! Allocate coarse grid matrix including pointer arrays
   !! Note: number of cells on coarse level corresponds to number of c-points on fine level
   !! Compute coarse grid matrix by multiplication with restriction and prolongation matrix:
   !!  A_coarse := R * A_fine * P
   !! -----------------------------------------------------------------------------------------
   DO NM = NMESHES_MIN, NMESHES_MAX

      SLF => SCARC(NM)%LEVEL(NL)               !< Pointer to fine level
      SLC => SCARC(NM)%LEVEL(NL+1)             !< Pointer to coarse level

      SLC%NC  = SLF%NCC
      SLC%NCE = SLF%NCCE
      SLC%NW  = SLF%NCW
   
      !< Allocate WALL-information for coarser grids (own and neighboring grids)
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
    !<   ICC = SLF%CELL_TYPE(IC)
      !ENDDO


   ENDDO
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROL A')

!>  IF (NMESHES>1) CALL SCARC_SETUP_WALLINFO_AMG (NL)          !HIER NOCHMAL CHECKEN ----

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
!>  DO NM = NMESHES_MIN, NMESHES_MAX
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
SUBROUTINE PRINT_MATRIX(A, A_ROW, A_COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , DIMENSION(:), INTENT(IN) :: A_COL
INTEGER , INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(1), INTENT(IN):: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60):: MATRIX
REAL(EB):: MATRIX_LINE(1000)

WRITE (MATRIX, '(A,A1,A,i2.2,A,i2.2,A)') 'matrix/',CNAME,'_mesh',NM,'_level',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=MATRIX)
      
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= A_ROW(IC), A_ROW(IC+1)-1       
         IF (A_COL(ICOL)==JC) MATRIX_LINE(JC)=A(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX,1000) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO

CLOSE(MMATRIX)

1000 FORMAT(F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',&
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,';')
END SUBROUTINE PRINT_MATRIX

!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE PRINT_MATRIX2(A, A_ROW, A_COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , DIMENSION(:), INTENT(IN) :: A_COL
INTEGER , INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(1), INTENT(IN):: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60):: MATRIX
REAL(EB):: MATRIX_LINE(1000)

WRITE (MATRIX, '(A,A1,A,i2.2,A,i2.2,A)') 'matrix2/',CNAME,'_mesh',NM,'_level',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=MATRIX)
      
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= A_ROW(IC), A_ROW(IC+1)-1       
         IF (A_COL(ICOL)==JC) MATRIX_LINE(JC)=A(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX,1000) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO

CLOSE(MMATRIX)

1000 FORMAT(40f8.2)
END SUBROUTINE PRINT_MATRIX2


!> ----------------------------------------------------------------------------------------------------
!> Store subdivision information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COORDINATES_AMG(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER IXC, IYC, IZC
TYPE (SCARC_LEVEL_TYPE), POINTER :: SLC, SLF

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) '======= CALLING SETUP_COORDINATES_AMG'

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

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,'(a,8f8.2)') 'SLC%XCOR:', (SLC%XCOR(IXC), IXC=0,SLC%NX)
   WRITE(SCARC_LU,'(a,8f8.2)') 'SLC%ZCOR:', (SLC%ZCOR(IZC), IZC=0,SLC%NZ)
ENDIF

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

IOR_LAST    = 0
NEIGHBORS   = 0
         
MESHES_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX
   
   SLC => SCARC(NM)%LEVEL(NL+1)

   SLC%SUBDIVISION = 0

   WALLCELLS_LOOP: DO IW = 1, SLC%NW
   
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
               WRITE(*,*) 'More than 20 neighbors at one face not allowed yet!'
               STOP
            ENDIF
         ENDDO NEIGHBOR_LOOP
      ENDIF
      
      IOR_LAST = IOR0
   
   ENDDO WALLCELLS_LOOP

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
     WRITE(SCARC_LU,*) 'SETTING SUBDIVISION for mesh ', NM, ' on level ', NL
     WRITE(SCARC_LU,'(a,7i8)') 'SUBDIVISION (1,.): ', (SLC%SUBDIVISION(1, IOR0), IOR0=-3,3)
     WRITE(SCARC_LU,'(a,7i8)') 'SUBDIVISION (2,.): ', (SLC%SUBDIVISION(2, IOR0), IOR0=-3,3)
     WRITE(SCARC_LU,'(a,7i8)') 'SUBDIVISION (3,.): ', (SLC%SUBDIVISION(3, IOR0), IOR0=-3,3)
   ENDIF
      
ENDDO MESHES_LOOP1
         
END SUBROUTINE SCARC_SETUP_SUBDIVISION_AMG



!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLINFO_AMG(NL)
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
MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SLC => SCARC(NM)%LEVEL(NL+1)

   DO IW = 1, SLC%NW
   
      NOM = SLC%WALL(IW)%NOM
      OSLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)
   
      ICW = SLC%WALL(IW)%ICW

      IF (TYPE_COARSENING == NSCARC_COARSENING_GMG) NCPL = 1

      WRITE(*,*) 'L: HIER WALL%ICN ändern auf OSLF!'
!>     SLC%WALL(IW)%NCPL = NCPL

!>     ALLOCATE(SLC%WALL(IW)%ICN(NCPL), STAT=IERR)
!>     CALL ChkMemErr('SCARC_SETUP_WALLINFO_AMG','ICN',IERR)

!>     ALLOCATE(SLC%WALL(IW)%ICE(NCPL), STAT=IERR)
!>     CALL ChkMemErr('SCARC_SETUP_WALLINFO_AMG','ICE',IERR)

!>     ALLOCATE(SLC%WALL(IW)%ICG(NCPL), STAT=IERR)
!>     CALL ChkMemErr('SCARC_SETUP_WALLINFO_AMG','ICG',IERR)

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

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL+1, 'SETUP_WALLINFO_AMG', 'WALL pos 0 ')
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MESHINFO, NL+1)

MESHES_LOOP2: DO NM = NMESHES_MIN, NMESHES_MAX

   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
      OS  => SCARC(NM)%OSCARC(NOM)      

      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE OTHER_MESHES_LOOP
      OSLC  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)      

      ALLOCATE(OSLC%WALL(1:OSLC%NW), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_WALLINFO_AMG','OSLC%WALL',IERR)

   ENDDO OTHER_MESHES_LOOP
ENDDO MESHES_LOOP2

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL+1, 'SETUP_WALLINFO_AMG', 'WALL pos 1 ')
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WALLINFO, NL+1)
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NL+1, 'SETUP_WALLINFO_AMG', 'WALL pos 2')

END SUBROUTINE SCARC_SETUP_WALLINFO_AMG


!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_LATEX_MATRIX(A, A_ROW, A_COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , DIMENSION(:), INTENT(IN) :: A_COL
INTEGER , INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(1), INTENT(IN):: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60):: MATRIX
REAL(EB):: MATRIX_LINE(1000)

WRITE (MATRIX, '(A,A1,A,i2.2,A,i2.2,A)') 'matrix/',CNAME,'_mesh',NM,'_level',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=MATRIX)
      
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= A_ROW(IC), A_ROW(IC+1)-1       
         IF (A_COL(ICOL)==JC) MATRIX_LINE(JC)=A(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX,1000) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO

CLOSE(MMATRIX)

1000 FORMAT(F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',&
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,';')
END SUBROUTINE SCARC_LATEX_MATRIX

!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_LATEX_MATRIX2(A, A_ROW, A_COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , DIMENSION(:), INTENT(IN) :: A_COL
INTEGER , INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(1), INTENT(IN):: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60):: MATRIX
REAL(EB):: MATRIX_LINE(1000)

WRITE (MATRIX, '(A,A1,A,i2.2,A,i2.2,A)') 'matrix2/',CNAME,'_mesh',NM,'_level',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=MATRIX)
      
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= A_ROW(IC), A_ROW(IC+1)-1       
         IF (A_COL(ICOL)==JC) MATRIX_LINE(JC)=A(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX,1000) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO

CLOSE(MMATRIX)

1000 FORMAT(40f8.2)
END SUBROUTINE SCARC_LATEX_MATRIX2

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

WRITE (MATRIX, '(A,A1,A,i2.2,A,i2.2,A)') 'matlab/',CNAME,'_mesh',NM,'_level',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=MATRIX)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'PRINTING MATLAB INFORMATION FOR LEVEL ', NL
   WRITE(SCARC_LU,*) 'NC1  =',NC1
   WRITE(SCARC_LU,*) 'NC2  =',NC2
   WRITE(SCARC_LU,*) 'CNAME=',CNAME
ENDIF
      
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
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,1104) IC,(MATRIX_LINE(JC),JC=1,NC2)
      CASE (8)
         WRITE(MMATRIX,1008) (MATRIX_LINE(JC),JC=1,NC2)
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,1108) IC,(MATRIX_LINE(JC),JC=1,NC2)
      CASE (12)
         WRITE(MMATRIX,1012) (MATRIX_LINE(JC),JC=1,NC2)
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,1112) IC,(MATRIX_LINE(JC),JC=1,NC2)
      CASE (16)
         WRITE(MMATRIX,1016) (MATRIX_LINE(JC),JC=1,NC2)
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,1116) IC,(MATRIX_LINE(JC),JC=1,NC2)
      CASE (24)
         WRITE(MMATRIX,1024) (MATRIX_LINE(JC),JC=1,NC2)
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,1124) IC,(MATRIX_LINE(JC),JC=1,NC2)
      CASE (32)
         WRITE(MMATRIX,1032) (MATRIX_LINE(JC),JC=1,NC2)
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,1132) IC,(MATRIX_LINE(JC),JC=1,NC2)
      CASE (64)
         WRITE(MMATRIX,1064) (MATRIX_LINE(JC),JC=1,NC2)
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,1164) IC,(MATRIX_LINE(JC),JC=1,NC2)
      CASE DEFAULT
         WRITE(*,*) 'WRONG VALUE FOR NC2 in MATLAB_MATRIX: ', NC2
         STOP
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
SUBROUTINE SCARC_LATEX_WALL(MTABLE, NM, NL)
INTEGER, INTENT(IN):: NM, NL, MTABLE
INTEGER :: NOM, NOM0, NCPL
INTEGER :: IW, ICPL, ICE
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

SL => SCARC(NM)%LEVEL(NL)

!WRITE (CTABLE, '(A,A,A,i2.2,A,i2.2,A)') 'tables/',CNAME,'_mesh',NM,'_level',NL,'.tex'
!MTABLE=GET_FILE_NUMBER()
!OPEN(MTABLE,FILE=CTABLE)

WRITE(MTABLE,100) NM
WRITE(MTABLE,500) NM
WRITE(MTABLE,1000)
WRITE(MTABLE,1001) NM
WRITE(MTABLE,2000)
!WRITE(MTABLE,2001) 'IW','IC\_WALL','IC\_EXT','IC\_GHOST','IC\_GHOST'
WRITE(MTABLE,2001) 'IW','ICW','ICE','ICG','ICN'
NOM0 = 0
DO IW = 1, SL%NW
   NOM = SL%WALL(IW)%NOM
   NCPL = SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)%NCPL
   IF (NOM == 0) CYCLE
   IF (NOM /= NOM0) WRITE(MTABLE,2000)
   SELECT CASE(NCPL)
   CASE (1)
       WRITE(MTABLE,3001) IW, SL%WALL(IW)%ICW,&
                          (SL%WALL(IW)%ICE(ICPL), ICPL=1,NCPL), &
                          (SL%WALL(IW)%ICG(ICPL), ICPL=1,NCPL), &
                          (SL%WALL(IW)%ICN(ICPL), ICPL=1,NCPL)
    CASE (2)
       WRITE(MTABLE,3002) IW, SL%WALL(IW)%ICW,&
                          (SL%WALL(IW)%ICE(ICPL), ICPL=1,NCPL), &
                          (SL%WALL(IW)%ICG(ICPL), ICPL=1,NCPL), &
                          (SL%WALL(IW)%ICN(ICPL), ICPL=1,NCPL)
    CASE (3)
       WRITE(MTABLE,3003) IW, SL%WALL(IW)%ICW,&
                          (SL%WALL(IW)%ICE(ICPL), ICPL=1,NCPL), &
                          (SL%WALL(IW)%ICG(ICPL), ICPL=1,NCPL), &
                          (SL%WALL(IW)%ICN(ICPL), ICPL=1,NCPL)
    CASE (4)
       WRITE(MTABLE,3004) IW, SL%WALL(IW)%ICW,&
                          (SL%WALL(IW)%ICE(ICPL), ICPL=1,NCPL), &
                          (SL%WALL(IW)%ICG(ICPL), ICPL=1,NCPL), &
                          (SL%WALL(IW)%ICN(ICPL), ICPL=1,NCPL)
   END SELECT
   !WRITE(MTABLE,2000)
   NOM0 = NOM
ENDDO
WRITE(MTABLE,2000)
WRITE(MTABLE,4000)
WRITE(MTABLE,5000)


!> ICG_TO_IWG
WRITE(MTABLE,500) NM
WRITE(MTABLE,1004)
WRITE(MTABLE,1002) NM
WRITE(MTABLE,2000)
WRITE(MTABLE,2002) 'ICE', 'IW'
WRITE(MTABLE,2000)
DO ICE = SL%NC+1, SL%NCE
   WRITE(MTABLE,3000) ICE, SL%ICG_TO_IWG(ICE)
   !WRITE(MTABLE,2000)
ENDDO
WRITE(MTABLE,2000)
WRITE(MTABLE,4000)
WRITE(MTABLE,5000)

100  FORMAT( /,'\vspace{1cm}',/,'{\bf MESH ',i3,'} \\',/)
500  FORMAT( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/, &
             '%%%    Mesh(',i2,')',/,&
             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
1000 FORMAT( '\begin{minipage}{8cm}')
1004 FORMAT( '\begin{minipage}{5cm}')
1001 FORMAT( '\flushleft{\bf {\footnotesize Mesh(',i3,'):WALL:}}\\[2ex]',/,&
             '\begin{tabular}{|c|c|c|c|c|}')
1002 FORMAT( '\flushleft{\bf {\footnotesize Mesh(',i3,'):WALL\_PTR:}}\\[2ex]',/,&
             '\begin{tabular}{|c|c|}')
2000 FORMAT( '\hline')
2001 FORMAT ('{\footnotesize ',a,'} & {\footnotesize ',a,'} & {\footnotesize ',a,'} & {\footnotesize ',a,&
             '} & {\footnotesize ',a,'} \\')
2002 FORMAT ('{\footnotesize ',a,'} & {\footnotesize ',a,'} \\')
3000 FORMAT (i4,' &',i4, '\\')
3001 FORMAT (i4,' &',i4, ' &',i4,' &',i4,' &',i4,'\\')
3002 FORMAT (i4,' &',i4, ' &',i4,',',i4,' &',i4,',',i4,' &',i4,',',i4,'\\')
3003 FORMAT (i4,' &',i4, ' &',i4,',',i4,',',i4,' &',i4,',',i4,',',i4,' &',i4,',',i4,',',i4,'\\')
3004 FORMAT (i4,' &',i4, ' &',i4,',',i4,',',i4,',',i4,' &',i4,',',i4,',',i4,',',i4,' &',i4,',',i4,',',i4,',',i4,'\\')
4000 FORMAT( '\end{tabular}')
5000 FORMAT( '\end{minipage}',/,'\vspace{0.5cm}')
END SUBROUTINE SCARC_LATEX_WALL


!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_LATEX_NOM(MTABLE, NM, NOM, NL)
INTEGER, INTENT(IN):: NM, NOM, NL, MTABLE
INTEGER :: IG
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSL

OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

!WRITE (CTABLE, '(A,A,A,i2.2,A,I2.2,A,i2.2,A)') 'tables/',CNAME,'_mesh',NM,'_nom',nom,'_level',NL,'.tex'
!MTABLE=GET_FILE_NUMBER()
!OPEN(MTABLE,FILE=CTABLE)

!WRITE(MTABLE,100) NM, NOM

!> NOM: ICG_TO_IWG
WRITE(MTABLE, 500) NM,':WALL\_PTR(',NOM,'), GHOST\_PTR(',NOM,'), NOM\_PTR(',NOM,')'
WRITE(SCARC_LU, 500) NM,':WALL\_PTR(',NOM,'), GHOST\_PTR(',NOM,'), NOM\_PTR(',NOM,')'
WRITE(MTABLE,1000)
WRITE(MTABLE,1500) NM,NOM
WRITE(MTABLE,2000)
WRITE(MTABLE,3000) 'IG', 'IW','ICE','IW'
WRITE(SCARC_LU,3000) 'IG', 'IW','ICE','IW'
WRITE(MTABLE,2000)
DO IG = 1, OSL%NCG
   WRITE(MTABLE,3001) IG, OSL%ICG_TO_IWG(IG), OSL%ICG_TO_ICE(IG), OSL%ICN_TO_ICE(IG)
   WRITE(SCARC_LU,3001) IG, OSL%ICG_TO_IWG(IG), OSL%ICG_TO_ICE(IG), OSL%ICN_TO_ICE(IG)
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
INTEGER  :: NM, IC, JC, ICOL, IDIAG, IP, IS , IPOS
REAL(EB) :: ADIAG, ACOL, ROW_SCALE, ROW_SUM, MAX_ROW_SUM = 0.9_EB, THRESHOLD = 0.25_EB
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

STRENGTH_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SETUP_STRENGTH_MATRIX: ============== '

   SL => SCARC(NM)%LEVEL(NL)            

   STRENGTH_CELL_LOOP: DO IC = 1, SL%NC

      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '==============   IC =', IC, '==== NTYPE=',NTYPE

      IDIAG = SL%A_ROW(IC)
      ADIAG = SL%A(IDIAG)

      !< get row sum and scaling factor
      ROW_SCALE = 0.0_EB
      ROW_SUM   = ADIAG
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,f12.6)') 'ADIAG    =',ADIAG    
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,f12.6)') 'ROW_SCALE=',ROW_SCALE  
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,f12.6)') 'ROW_SUM  =',ROW_SUM  

      DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
         JC = SL%A_COL(ICOL)
         IF (JC <= SL%NC) THEN
            ACOL = SL%A(ICOL)
            ROW_SCALE = MAX(ROW_SCALE, ACOL)
            ROW_SUM   = ROW_SUM + ACOL
            IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i4,a,f12.6,a,f12.6, a, f12.6)') &
                'ICOL=',ICOL,': ACOL=',ACOL,': ROW_SCALE=',ROW_SCALE, ': ROW_SUM=',ROW_SUM
         ENDIF
      ENDDO
      
      ROW_SUM = ABS(ROW_SUM/ADIAG)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,f12.6)') 'ROW_SUM=',ROW_SUM

      !< get row entries of strength matrix S 
      SL%S_COL(IDIAG) = -1
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i4,a,i4)') 'S_COL(',IDIAG,')=',SL%S_COL(IDIAG)

      IF ((ROW_SUM > MAX_ROW_SUM) .AND. (MAX_ROW_SUM < 1.0_EB)) THEN
         !< set all dependencies to be weak
         DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
            JC = SL%A_COL(ICOL)
            IF (JC <= SL%NC) THEN
               SL%S_COL(ICOL) = -1
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i4,a,i4)') 'A: S_COL(',ICOL,')=',SL%S_COL(ICOL)
            ENDIF
         ENDDO
      ELSE
         !< set dependencies to be weak related to threshold
         DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
            JC = SL%A_COL(ICOL)
            IF (JC <= SL%NC) THEN
               IF (SL%A(ICOL) <= THRESHOLD * ROW_SCALE) SL%S_COL(ICOL) = -1
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i4,a,f12.6,a,f12.6,a,i3,a,i4)') &
                 'ZZZZ2: A(',ICOL,')=',SL%A(ICOL), ': T*R_S=',THRESHOLD*ROW_SCALE,' S_COL(',ICOL,')=',SL%S_COL(ICOL)
            !ELSE
            !<   SL%S_COL(ICOL) = -1
            ENDIF
          ENDDO
      ENDIF

   ENDDO STRENGTH_CELL_LOOP

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      WRITE(SCARC_LU,*) 'A: ---------------------- S_ROW:', SL%NC
      WRITE(SCARC_LU,'(4i9)') (SL%S_ROW(IC), IC=1,SL%NC+1)
      WRITE(SCARC_LU,*) '---------------------- S_COL:'
      DO IC = 1, SL%NC
         WRITE(SCARC_LU,'(i5,a,20i9)') IC,':',(SL%S_COL(IP),IP=SL%S_ROW(IC),SL%S_ROW(IC+1)-1)
      ENDDO
      WRITE(SCARC_LU,*) '---------------------- S:'
      DO IC = 1, SL%NC
         WRITE(SCARC_LU,'(i5,a,20f9.2)') IC,':',(SL%S(IP),IP=SL%S_ROW(IC),SL%S_ROW(IC+1)-1)
      ENDDO
   ENDIF

   !! Compress strength matrix
   IS = 1
   STRENGTH_CELL_LOOP2: DO IC = 1, SL%NC
      SL%S_ROW(IC) = IS
      DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
         !IF (SL%S_COL(ICOL) > -1.AND.SL%S_COL(ICOL)<=SL%NC) THEN
         IF (SL%S_COL(ICOL) > -1) THEN
            SL%S_COL(IS) = SL%S_COL(ICOL)
            WRITE(SCARC_LU,*) 'S_COL(',IS,')=',SL%S_COL(ICOL)
            IS = IS + 1
         ENDIF
      ENDDO
   ENDDO STRENGTH_CELL_LOOP2
   SL%S_ROW(SL%NC+1) = IS
   
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      WRITE(SCARC_LU,*) 'B: ---------------------- S_ROW:', SL%NC
      WRITE(SCARC_LU,'(4i9)') (SL%S_ROW(IC), IC=1,SL%NC+1)
      WRITE(SCARC_LU,*) '---------------------- S_COL:'
      DO IC = 1, SL%NC
         WRITE(SCARC_LU,'(i5,a,20i9)') IC,':',(SL%S_COL(IP),IP=SL%S_ROW(IC),SL%S_ROW(IC+1)-1)
      ENDDO
      WRITE(SCARC_LU,*) 'SIZE(ST_ROW)=',SIZE(SL%ST_ROW)
      WRITE(SCARC_LU,*) 'SIZE(ST_COL)=',SIZE(SL%ST_COL)
   ENDIF

   DO IC = 1, SL%NCE+1
      SL%ST_ROW(IC) = 0
   ENDDO

   IS = SL%S_ROW(SL%NC+1)-1
   DO ICOL = 1, IS
      SL%ST_ROW(SL%S_COL(ICOL)+1) = SL%ST_ROW(SL%S_COL(ICOL)+1) + 1
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
         WRITE(SCARC_LU,*) 'A: ST_ROW(',SL%S_COL(ICOL)+1,')=',SL%ST_ROW(SL%S_COL(ICOL)+1)
   ENDDO 
   SL%ST_ROW(1) = 1

   DO IC = 1, SL%NCE
      SL%ST_ROW(IC+1)= SL%ST_ROW(IC+1) + SL%ST_ROW(IC) 
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
         WRITE(SCARC_LU,*) 'B: ST_ROW(',IC+1,')=',SL%ST_ROW(IC+1)
   ENDDO 
   DO IC = 1, SL%NC
      DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
         IPOS = SL%S_COL(ICOL)
         SL%ST_COL(SL%ST_ROW(IPOS)) = IC
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
         WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
           'C: IC=',IC,': ICOL=',ICOL,': IPOS=', IPOS,': ST_COL(',SL%ST_ROW(IPOS),')=',IC
         SL%ST_ROW(IPOS) = SL%ST_ROW(IPOS) + 1
      ENDDO
   ENDDO 
   DO IC = SL%NCE+1, 2, -1
      SL%ST_ROW(IC) = SL%ST_ROW(IC-1)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
         WRITE(SCARC_LU,*) 'D: ST_ROW(',IC,')=',SL%ST_ROW(IC)
   ENDDO
   SL%ST_ROW(1) = 1

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      WRITE(SCARC_LU,*) '---------------------- ST_ROW:', SL%NC
      WRITE(SCARC_LU,'(4i9)') (SL%ST_ROW(IC), IC=1,SL%NCE+1)
      WRITE(SCARC_LU,*) '---------------------- ST_COL:'
      DO IC = 1, SL%NCE
         WRITE(SCARC_LU,'(i5,a,20i9)') IC,':',(SL%ST_COL(IP),IP=SL%ST_ROW(IC),SL%ST_ROW(IC+1)-1)
      ENDDO
   ENDIF

ENDDO STRENGTH_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_STRENGTH_MATRIX


!> ------------------------------------------------------------------------------------------------
!> Determine measure of cells corresponding to requested coarsening type
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COLORING(NTYPE, NL)
INTEGER, PARAMETER  :: NCOUPLINGS = 20
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, IC, IC0, JC, KC, ICOL, JCOL, IRAND, REMAINING_CELLS, ICG
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

      RS3_MESH_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SL => SCARC(NM)%LEVEL(NL)            

         REMAINING_CELLS = 0

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      WRITE(SCARC_LU,*) 'SCARC_SETUP_COLORING: INIT: NL=',NL
      WRITE(SCARC_LU,*) 'SCARC_SETUP_COLORING: CELL_TYPE: ', SL%NC
      WRITE(SCARC_LU,*) (SL%CELL_TYPE(IC), IC=1, SL%NC)
   ENDIF

         !< Currently the measures are computed as row sums of ST (number of influences for IC)
         RS3_MEASURE_LOOP0: DO IC = 1, SL%NC
            SL%MEASURE(IC) = SL%ST_ROW(IC+1)-SL%ST_ROW(IC) 
         ENDDO RS3_MEASURE_LOOP0


         !< Subdivide in fine and coarse grid cells
         RS3_MEASURE_LOOP1: DO IC = 1, SL%NC
            
            IF (SL%S_ROW(IC+1)-SL%S_ROW(IC) == 0) THEN
               SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE0
               SL%MEASURE(IC) = 0.0_EB
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'A: CELL_TYPE(',IC,') =',SL%CELL_TYPE(IC),': MEASURE(',IC,')=',SL%MEASURE(IC),': Remaining cells=',REMAINING_CELLS
             !< IF (AGGRESSIVE2) SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE0
            ELSE
               SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
               REMAINING_CELLS   = REMAINING_CELLS + 1
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'B: CELL_TYPE(',IC,') =',SL%CELL_TYPE(IC),': MEASURE(',IC,')=',SL%MEASURE(IC),': Remaining cells=',REMAINING_CELLS
            ENDIF
         ENDDO RS3_MEASURE_LOOP1

         RS3_MEASURE_LOOP2: DO IC = 1, SL%NC
            MEASURE = SL%MEASURE(IC)
            IF (SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_FINE0 .AND. SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_COARSE0) THEN
               IF (SL%MEASURE(IC) > 0.0_EB) THEN
                  !SL%MEASURE(IC) = SL%MEASURE(IC) 
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'C: CELL_TYPE(',IC,') =',SL%CELL_TYPE(IC),': MEASURE(',IC,')=',SL%MEASURE(IC),': Remaining cells=',REMAINING_CELLS
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
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'D: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS
                        ELSE
                           NEW_MEASURE = SL%MEASURE(JC)+1
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'E: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS
                          
                        ENDIF
                     ENDIF
                  ENDDO
                  REMAINING_CELLS = REMAINING_CELLS - 1
               ENDIF
            ENDIF
         ENDDO RS3_MEASURE_LOOP2


CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_COLORING', 'MEASURE AFTER FIRST RS3 LOOP ')
      
         RS3_CYCLE_LOOP: DO WHILE (REMAINING_CELLS > 0)
         

          !< get maximum (remaining) measure for all cells
            MEASURE_MAX = MAXVAL(SL%MEASURE(1:SL%NC))
            IF (MEASURE_MAX <= EPS) EXIT RS3_CYCLE_LOOP
         
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'MEASURE_MAX=',MEASURE_MAX

            RS3_CELL_LOOP: DO IC = 1, SL%NC
         
             !< Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, SL%MEASURE(IC))) THEN
         
                  SL%MEASURE(IC)  = NSCARC_MEASURE_NONE
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
                  REMAINING_CELLS = REMAINING_CELLS - 1
         
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'F: CELL_TYPE(',IC,') =',SL%CELL_TYPE(IC),': MEASURE(',IC,')=',SL%MEASURE(IC),': Remaining cells=',REMAINING_CELLS

                !< Determine set of fine cells 
                  DO ICOL = SL%ST_ROW(IC), SL%ST_ROW(IC+1)-1

                   !< IF JC hasn't been marked yet, set it to be a fine cell which is no longer measured
                     JC = SL%ST_COL(ICOL)

                     IF (JC > SL%NC) CYCLE
                     IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) THEN

                        SL%MEASURE(JC)  = NSCARC_MEASURE_NONE
                        SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_FINE
                        REMAINING_CELLS = REMAINING_CELLS - 1
      
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'G: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS

                      !<  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                        DO JCOL = SL%S_ROW(JC), SL%S_ROW(JC+1)-1
                           KC = SL%S_COL(JCOL)
                           IF (SL%CELL_TYPE(KC)==NSCARC_CELL_TYPE_NONE) THEN
                              SL%MEASURE(KC) = SL%MEASURE(KC) + 1.0_EB
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'H: CELL_TYPE(',KC,') =',SL%CELL_TYPE(KC),': MEASURE(',KC,')=',SL%MEASURE(KC),': Remaining cells=',REMAINING_CELLS
   
                           ENDIF
                        ENDDO 
                     ENDIF
                  ENDDO 
WRITE(SCARC_LU,*) '====================== REMAINING_CELLS = ', REMAINING_CELLS
!CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELL_TYPE', 'RS3_MEASURE ')
!CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_CELL_TYPE', 'RS3_CELL_TYPE ')
!WRITE(SCARC_LU,*) '====================== REMAINING_CELLS = ', REMAINING_CELLS
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_COLORING', 'MEASURE AFTER SECOND RS3 LOOP ')
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE AFTER SECOND RS3 LOOP ')
                  EXIT RS3_CELL_LOOP
               ENDIF
            ENDDO RS3_CELL_LOOP

            DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
               JC = SL%S_COL(ICOL)
               IF (JC > SL%NC) CYCLE
               IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) THEN
                  MEASURE = SL%MEASURE(JC) - 1
                  SL%MEASURE(JC) = MEASURE 
                  IF (MEASURE > 0.0_EB) THEN
                     SL%MEASURE(JC) = MEASURE
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'I: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS
                  ELSE
                     SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_FINE
                     REMAINING_CELLS = REMAINING_CELLS - 1 
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'J: CELL_TYPE(',JC,') =',SL%CELL_TYPE(JC),': MEASURE(',JC,')=',SL%MEASURE(JC),': Remaining cells=',REMAINING_CELLS
                     DO JCOL = SL%S_ROW(JC), SL%S_ROW(JC+1)-1
                        KC = SL%S_COL(JCOL)
                        IF (SL%CELL_TYPE(KC)==NSCARC_CELL_TYPE_NONE) THEN
                           SL%MEASURE(KC) = SL%MEASURE(KC) + 1.0_EB
                           MEASURE_MAX = MAX(MEASURE_MAX, SL%MEASURE(KC))
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i6,a,i3,a,f12.6,a,i3)') &
         'K: CELL_TYPE(',KC,') =',SL%CELL_TYPE(KC),': MEASURE(',KC,')=',SL%MEASURE(KC),': Remaining cells=',REMAINING_CELLS
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

       !< exchange information with neighboring meshes
         IF (NMESHES > 1) THEN

            DO IC = 1, SL%NCE
               SL%GRAPH(IC) = -1
            ENDDO

            IC0 = 1
            DO IC = 1, SL%NC
        IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'GRAPH::: =============== IC =',IC 
        IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'GRAPH:::0  CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
               IF (ICT2 /= IC) ICT = -1
        IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'GRAPH:::0  ICT=',ICT
   
               IF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE) THEN

                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                     JC = SL%S_COL(ICOL)
                     IF (JC <= SL%NC) THEN
                     IF (SL%CELL_TYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                        SL%GRAPH(JC) = IC
        IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'GRAPH::: GRAPH(',JC,')=',SL%GRAPH(JC)
                     ENDIF
                     ENDIF
                  ENDDO
   
                !< Hier fehlt noch die Abfrage nach Nachbarmesh !!!
   
                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                     JC = SL%S_COL(ICOL)
                     IF (JC <= SL%NC) THEN
                     IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_FINE) THEN
        IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'GRAPH:::1  CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
                  !< DIESER PART WIRD ERST FÃœR ANDERE VERFEINERUNGEN AKTIV
                  !< ACHTUNG: DANN NOCHMAL ÃœBERPRÃœFEN
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
        IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'GRAPH:::2  CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
                                 IF (ICT > -1) THEN
                                    SL%CELL_TYPE(ICT) = NSCARC_CELL_TYPE_FINE
        IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'GRAPH:::3  CELL_TYPE(',ICT,')=',SL%CELL_TYPE(ICT)
                                    ICT = -1
                                 ENDIF
                                 !< Hier fehlt noch Nachbaranteil
                                 BNONEMPTY = .FALSE.
                                 EXIT
                              ELSE
                                 ICT  = JC
                                 ICT2 = IC
                                 SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_COARSE
        IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'GRAPH::: CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
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
               IC0 = IC0 + 1 !> Achtung, nochmal prÃ¼fen, eventuell IC0 als Index verwenden?
            ENDDO
         ENDIF


         DO IC = 1, SL%NC
     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SUSI: GRAPH(',IC,')=',SL%GRAPH(IC)
         ENDDO
         DO IC = 1, SL%NC
     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SUSI: CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
         ENDDO
      ENDDO RS3_MESH_LOOP

    !< Exchange CELL_TYPE-data along internal boundaries
      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)


      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE BEFORE LAST EXCHANGE')
      CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE AFTER LAST EXCHANGE')

    !< Third pass: identify fine cells along boundary and get their coarse neighbors
      RS3_MESH_LOOP2: DO NM = NMESHES_MIN, NMESHES_MAX
         DO IC = 1, SL%NC
            SL%GRAPH(IC) = -1
         ENDDO
      ENDDO RS3_MESH_LOOP2


!> ----------------------------------------------------------------------------------------
!> Falgout coarsening
!> ----------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_FALGOUT)

      FCELL = ZCELL
      FALGOUT_MESHES_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX

         SL => SCARC(NM)%LEVEL(NL)

         !< Reset the measures as row sums of ST (number of influences for IC) plus random number
         !< to make them unique
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

      !< Initial exchange of measure array
      TYPE_VECTOR = NSCARC_VECTOR_MEASURE
      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR , NL)

      FALGOUT_MESHES_LOOP2: DO NM = NMESHES_MIN, NMESHES_MAX

         !< reset CELL_TYPE array and eliminate former fine cells
         FALGOUT_INTERNAL_CELL_LOOP1: DO IC = 1, SL%NC
            BNEIGHBOR = .FALSE.
            !< Is it a cell with a neighboring cell in another mesh?
            DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
               JC = SL%S_COL(ICOL)
               IF (JC > SL%NC) THEN
                  BNEIGHBOR = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (BNEIGHBOR .OR. SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE) SL%CELL_TYPE(IC)=NSCARC_CELL_TYPE_NONE
         ENDDO FALGOUT_INTERNAL_CELL_LOOP1

         !< initialize GRAPH and CELL_TYPE on ghost cells
         FALGOUT_EXTENDED_CELL_LOOP: DO IC = SL%NC+1, SL%NCE
            SL%GRAPH(IC)    = IC
            SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
         ENDDO FALGOUT_EXTENDED_CELL_LOOP

         ICG = 1
         FALGOUT_INTERNAL_CELL_LOOP2: DO IC = 1, SL%NC
            IF (SL%INTERNAL_BDRY_CELL(IC).NE.0.OR.SL%CELL_TYPE(IC)<NSCARC_CELL_TYPE_NONE) &
               SL%CELL_TYPE(IC)=NSCARC_CELL_TYPE_NONE
            IF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_ZPNT) THEN
               IF (SL%MEASURE(IC) >= 1.0_EB .OR. (SL%S_ROW(IC+1)-SL%S_ROW(IC)) > 0) THEN
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
                  SL%GRAPH(ICG) = IC
                  ICG = ICG + 1
               ELSE
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
               ENDIF
            ELSE IF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_SFPNT) THEN
               SL%MEASURE(IC) = NSCARC_MEASURE_NONE
            ELSE
               SL%GRAPH(ICG) = IC
               ICG = ICG + 1
            ENDIF
         ENDDO FALGOUT_INTERNAL_CELL_LOOP2
         
         IGRAPH  = ICG-1
         IGRAPHE = SL%NCE-SL%NC

      ENDDO FALGOUT_MESHES_LOOP2


      FALGOUT_EXTERNAL_CELL_LOOP2: DO IC = SL%NC+1, SL%NCE
         SL%MEASURE(IC) = NSCARC_MEASURE_NONE
      ENDDO FALGOUT_EXTERNAL_CELL_LOOP2

     !> ------------------------------------------------------------
     !> ----> ONLY TEMPORARILY - SET MEASURES BY HAND
     !> ------------------------------------------------------------
      FALGOUT_MESHES_LOOP20: DO NM = NMESHES_MIN, NMESHES_MAX

      IF (NM == 1) THEN
    !< Mesh 1
      SCARC(1)%LEVEL(NL)%MEASURE(13) = 3.6618       
      SCARC(1)%LEVEL(NL)%MEASURE(14) = 4.9696       
      SCARC(1)%LEVEL(NL)%MEASURE(15) = 4.8591       
      SCARC(1)%LEVEL(NL)%MEASURE(16) = 4.0864
      SCARC(1)%LEVEL(NL)%MEASURE(9)  = 3.0272       
      SCARC(1)%LEVEL(NL)%MEASURE(10) = 4.6014       
      SCARC(1)%LEVEL(NL)%MEASURE(11) = 4.4802       
      SCARC(1)%LEVEL(NL)%MEASURE(12) = 4.8368
      SCARC(1)%LEVEL(NL)%MEASURE(5)  = 3.5116       
      SCARC(1)%LEVEL(NL)%MEASURE(6)  = 4.4809       
      SCARC(1)%LEVEL(NL)%MEASURE(7)  = 4.2316       
      SCARC(1)%LEVEL(NL)%MEASURE(8)  = 4.8414
      SCARC(1)%LEVEL(NL)%MEASURE(1)  = 2.0215       
      SCARC(1)%LEVEL(NL)%MEASURE(2)  = 3.3343       
      SCARC(1)%LEVEL(NL)%MEASURE(3)  = 3.6478       
      SCARC(1)%LEVEL(NL)%MEASURE(4)  = 3.9119

      ELSE IF (NM==2) THEN
    !< Mesh 2
      SCARC(2)%LEVEL(NL)%MEASURE(13) =  4.4928
      SCARC(2)%LEVEL(NL)%MEASURE(14) =  4.0042
      SCARC(2)%LEVEL(NL)%MEASURE(15) =  4.9126
      SCARC(2)%LEVEL(NL)%MEASURE(16) =  3.6161
      SCARC(2)%LEVEL(NL)%MEASURE(9)  =  4.7065
      SCARC(2)%LEVEL(NL)%MEASURE(10) =  4.5361
      SCARC(2)%LEVEL(NL)%MEASURE(11) =  4.8637
      SCARC(2)%LEVEL(NL)%MEASURE(12) =  3.3562
      SCARC(2)%LEVEL(NL)%MEASURE(5)  =  4.0444
      SCARC(2)%LEVEL(NL)%MEASURE(6)  =  4.6998
      SCARC(2)%LEVEL(NL)%MEASURE(7)  =  4.2786
      SCARC(2)%LEVEL(NL)%MEASURE(8)  =  3.5202
      SCARC(2)%LEVEL(NL)%MEASURE(1)  =  3.0215
      SCARC(2)%LEVEL(NL)%MEASURE(2)  =  3.4658
      SCARC(2)%LEVEL(NL)%MEASURE(3)  =  3.4034
      SCARC(2)%LEVEL(NL)%MEASURE(4)  =  2.3706

      ELSE IF (NM==3) THEN
    !< Mesh 3
      SCARC(3)%LEVEL(NL)%MEASURE(13) =  2.3237
      SCARC(3)%LEVEL(NL)%MEASURE(14) =  3.0387
      SCARC(3)%LEVEL(NL)%MEASURE(15) =  3.9660
      SCARC(3)%LEVEL(NL)%MEASURE(16) =  3.1458
      SCARC(3)%LEVEL(NL)%MEASURE(9)  =  3.3858
      SCARC(3)%LEVEL(NL)%MEASURE(10) =  4.4708
      SCARC(3)%LEVEL(NL)%MEASURE(11) =  4.2472
      SCARC(3)%LEVEL(NL)%MEASURE(12) =  4.8756
      SCARC(3)%LEVEL(NL)%MEASURE(5)  =  3.5771
      SCARC(3)%LEVEL(NL)%MEASURE(6)  =  4.9188
      SCARC(3)%LEVEL(NL)%MEASURE(7)  =  4.3256
      SCARC(3)%LEVEL(NL)%MEASURE(8)  =  4.1991
      SCARC(3)%LEVEL(NL)%MEASURE(1)  =  3.0215
      SCARC(3)%LEVEL(NL)%MEASURE(2)  =  4.5974
      SCARC(3)%LEVEL(NL)%MEASURE(3)  =  4.1590
      SCARC(3)%LEVEL(NL)%MEASURE(4)  =  4.8292

      ELSE IF (NM==4) THEN
    !< Mesh 4
      SCARC(4)%LEVEL(NL)%MEASURE(13) =  3.1547
      SCARC(4)%LEVEL(NL)%MEASURE(14) =  3.0733
      SCARC(4)%LEVEL(NL)%MEASURE(15) =  3.0195
      SCARC(4)%LEVEL(NL)%MEASURE(16) =  2.6755
      SCARC(4)%LEVEL(NL)%MEASURE(9)  =  4.0651
      SCARC(4)%LEVEL(NL)%MEASURE(10) =  4.4055
      SCARC(4)%LEVEL(NL)%MEASURE(11) =  4.6307
      SCARC(4)%LEVEL(NL)%MEASURE(12) =  3.3950
      SCARC(4)%LEVEL(NL)%MEASURE(5)  =  4.1099
      SCARC(4)%LEVEL(NL)%MEASURE(6)  =  4.1378
      SCARC(4)%LEVEL(NL)%MEASURE(7)  =  4.3727
      SCARC(4)%LEVEL(NL)%MEASURE(8)  =  3.8780
      SCARC(4)%LEVEL(NL)%MEASURE(1)  =  4.0215
      SCARC(4)%LEVEL(NL)%MEASURE(2)  =  4.7289
      SCARC(4)%LEVEL(NL)%MEASURE(3)  =  4.9146
      SCARC(4)%LEVEL(NL)%MEASURE(4)  =  2.2879

      ENDIF

      ENDDO FALGOUT_MESHES_LOOP20
      


      !< Coloring loop until all cells are coarse or fine
      FALGOUT_MESHES_LOOP3: DO NM = NMESHES_MIN, NMESHES_MAX

         ICYCLE = 0
         FALGOUT_GRAPH_LOOP3: DO WHILE (BCYCLE)

!IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'DUDU: CALLING MEASURE_ADD in CYCLE ', ICYCLE, NMESHES

            CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_COLORING', 'MEASURE0 IN FALGOUT WHILE')

            IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MEASURE_ADD, NL)

            !CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_COLORING', 'MEASURE0 IN FALGOUT WHILE')

            IF (ICYCLE > 0) THEN


               ICG = 1
               FALGOUT_GRAPH_LOOP31: DO WHILE (ICG <= IGRAPH)

                  IC = SL%GRAPH(ICG)

                  !< if cell isn't marked as coarse cell yet and has measure less than 1, mark it as fine cell
                  !< take care that all dependencies have been taken into account
                  IF (SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_COARSE .AND. SL%MEASURE(IC) < NSCARC_MEASURE_ONE) THEN
                     SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
                     DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                        JC = SL%S_COL(ICOL)
                        IF (JC < 0) CYCLE
                        SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
                     ENDDO
                  ENDIF

                  !< if cell is already marked as fine or coarse, set its measure to zero and extract it from
                  !< the graph (put it at the end of the graph array and decrease number of relevant graph entries)
                  IF (SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_NONE) THEN
                     SL%MEASURE(IC) = NSCARC_MEASURE_NONE
                     SL%GRAPH(ICG) = SL%GRAPH(IGRAPH)
                     SL%GRAPH(IGRAPH) = IC
                     IGRAPH = IGRAPH - 1
                     ICG = ICG - 1
                  ENDIF
                 
                  ICG = ICG + 1
               ENDDO FALGOUT_GRAPH_LOOP31
       
            ENDIF

            !< exchanges measure and CELL_TYPEs of neighboring cells
            TYPE_VECTOR = NSCARC_VECTOR_MEASURE
            IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)
            !IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)

            !< Get global number of relevant graph entries
            IF (NMESHES>1 .AND. USE_MPI) THEN
               CALL MPI_ALLREDUCE (IGRAPH, IGRAPH_GLOBAL, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
            ENDIF

          !< If all cells have been encountered for, leave loop
            IF (IGRAPH_GLOBAL == 0) EXIT FALGOUT_GRAPH_LOOP3

          !< search for an independent set of points with maximal measure 
            IF (ICYCLE > 0) THEN
              CALL SCARC_SETUP_GRAPHSET(IGRAPH, NM, NL)
            ENDIF
            ICYCLE = ICYCLE + 1 

          !< exchanges CELL_TYPEs of neighboring cells (with new information from graphset)
            IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)

            ICG = SL%NC+1
            FALGOUT_GRAPH_LOOP32: DO WHILE (ICG <= SL%NC+IGRAPHE)
               IC = SL%GRAPH(ICG)

               IF (IC < 0) CYCLE
               IF (SL%CELL_TYPE(IC) < NSCARC_CELL_TYPE_NONE) THEN
                  SL%GRAPH(ICG) = SL%GRAPH(SL%NC+IGRAPHE)
                  SL%GRAPH(SL%NC+IGRAPHE) = IC
                  IGRAPHE = IGRAPHE - 1
               ENDIF
               ICG = ICG + 1
            ENDDO FALGOUT_GRAPH_LOOP32

            SL%MEASURE(SL%NC+1:SL%NCE) = NSCARC_MEASURE_NONE

            FALGOUT_GRAPH_LOOP33: DO ICG = 1, IGRAPH

               IC = SL%GRAPH(ICG)
   
             !< Coarse cells don't interpolate from influencing neighbors
               IF (SL%CELL_TYPE(IC) > NSCARC_CELL_TYPE_NONE) THEN
   
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
   
                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                     JC = SL%S_COL(ICOL)

                   !< remove edge from S and decrement measures of unmarked neighbors correspondingly
                     IF (JC > -1) THEN
                        SL%S_COL(ICOL) = - SL%S_COL(ICOL) - 1
                        IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) SL%MEASURE(JC) = SL%MEASURE(JC) - 1
                     ENDIF
                  ENDDO
   
               ELSE
   
                !< dependencies which are already marked
                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1

                     JC = SL%S_COL(ICOL)
                     IF (JC < 0) JC = -JC - 1

                     IF (SL%CELL_TYPE(JC) > NSCARC_CELL_TYPE_NONE) THEN

                   !< remove edge from S and temporarily reset CELL_TYPE
                        IF (SL%S_COL(ICOL) > -1) SL%S_COL(ICOL) = -SL%S_COL(ICOL) - 1
                        SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_CPNT            
                 
                     ELSE IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_SFPNT) THEN    !< necessary ??
                        IF (SL%S_COL(ICOL) > -1) SL%S_COL(ICOL) = -SL%S_COL(ICOL) - 1
                     ENDIF
   
                  ENDDO
   
                !< dependencies which aren't marked yet
                  DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
    
                     JC = SL%S_COL(ICOL)
                     IF (JC > -1 .AND. JC<=SL%NC) THEN
                     
                        BREAK = .TRUE.

                      !< check if there are common C-points
                        DO JCOL = SL%S_ROW(JC), SL%S_ROW(JC+1)-1
                           KC = SL%S_COL(JCOL)
                           IF (KC < 0) KC = -KC - 1

                           IF (KC <= SL%NC) THEN
                            !< remove edge from S and update measure
                              IF (SL%CELL_TYPE(KC) == NSCARC_CELL_TYPE_CPNT) THEN
                                 SL%S_COL(ICOL) = - SL%S_COL(ICOL) - 1
                                 SL%MEASURE(JC) = SL%MEASURE(JC) - 1 
                                 BREAK = .FALSE.
                                 EXIT
                              ENDIF
                           ENDIF
                        ENDDO
   
                        IF (BREAK) THEN
                           DO JCOL = SL%S_ROW(JC), SL%S_ROW(JC+1)-1
                              KC = SL%S_COL(JCOL)
                              IF (KC < 0) KC = -KC - 1
                              IF (KC > SL%NC) THEN
                               !< remove edge from S and update measure
                                 IF (SL%CELL_TYPE(KC) == NSCARC_CELL_TYPE_CPNT) THEN
                                    SL%S_COL(ICOL) = - SL%S_COL(ICOL) - 1
                                    SL%MEASURE(JC) = SL%MEASURE(JC) - 1 
                                    BREAK = .FALSE.
                                    EXIT
                                 ENDIF
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF

             !< reset CELL_TYPEs
               DO ICOL = SL%S_ROW(IC), SL%S_ROW(IC+1)-1
                  JC = SL%S_COL(ICOL)
                  IF (JC < 1) JC = -JC - 1
                  IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_CPNT) SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_COARSE
               ENDDO

            ENDDO FALGOUT_GRAPH_LOOP33
            CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_COLORING', 'MEASURE1 IN FALGOUT WHILE')
            CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_COLORING', 'CELL_TYPE1 IN FALGOUT WHILE')
            CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_GRAPH   , NL, 'SETUP_COLORING', 'GRAPH1 IN FALGOUT WHILE')

         ENDDO FALGOUT_GRAPH_LOOP3

       !< reset S-matrix
         DO ICOL = 1, SL%S_ROW(SL%NC+1)
            IF (SL%S_COL(ICOL) < 0) SL%S_COL(ICOL) = -SL%S_COL(ICOL)-1
         ENDDO
                  
      ENDDO FALGOUT_MESHES_LOOP3

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

!>
!> First mark every cell from the 'internal' graphset with measure bigger than one as a coarse cell
!>
DO IG = 1, IGRAPH
   IC = SL%GRAPH(IG)
   IF (SL%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
   ENDIF
ENDDO

!> Do the same with cells on the overlapping areas
DO IG = SL%NC+1, SL%NCE
   IC = SL%GRAPH(IG)
   IF (IC < 0) CYCLE
   IF (SL%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
   ENDIF
ENDDO


!>
!> remove nodes from the initial independent set depending on their measure
!> For each cell consider every connection in stencil and set cell with smaller measure to Zero-Type
!>
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
            ELSE IF (SL%MEASURE(JC) > SL%MEASURE(IC)) THEN
               SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
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
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i4,a,i4,a,i4)') 'PC(',INX,',',INZ,')=',PC(INX,INZ)
         CELL_TYPE(PC(INX,INZ)) = NSCARC_CELL_TYPE_SFINE
      ENDDO
   ENDDO
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'INX0=',INX0,': INZ0=',INZ0, NY1, NY2, INCRY
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
      RS3_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SL => SCARC(NM)%LEVEL(NL)
         REMAINING_CELLS = SL%NC

         RS3_CYCLE_LOOP: DO WHILE (REMAINING_CELLS > 0)
         

            !< get maximum (remaining) measure for all cells
            MEASURE_MAX = MAXVAL(SL%MEASURE(1:SL%NC))
            IF (MEASURE_MAX <= EPS) EXIT RS3_CYCLE_LOOP
         
            RS3_CELL_LOOP: DO IC = 1, SL%NC
         
               !< Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, SL%MEASURE(IC))) THEN
         
                  SL%MEASURE(IC)  = NSCARC_MEASURE_NONE
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
                  REMAINING_CELLS = REMAINING_CELLS - 1
         
                  !< Determine set of fine cells 
                  DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1

                     IF (STRONGLY_COUPLED(SL%A, SL%A_ROW, IC, ICOL)) THEN

                        !< IF JC hasn't been marked yet, set it to be a fine cell which is no longer measured
                        JC = SL%A_COL(ICOL)
                        IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) THEN

                           SL%MEASURE(JC)  = NSCARC_MEASURE_NONE
                           SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_FINE
                           REMAINING_CELLS = REMAINING_CELLS - 1
         
                           !<  increase measures of cells KC adjacent to fine cells JC based on strong couplings
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
WRITE(SCARC_LU,*) '====================== REMAINING_CELLS = ', REMAINING_CELLS
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

      FALGOUT_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

WRITE(SCARC_LU,*) 'A'
         !< copy column pointers for matrix
         DO IA = 1, SL%NA
            SL%S_COL(IA) = SL%A_COL(IA)
         ENDDO

WRITE(SCARC_LU,*) 'B'
         !< this part may be changed for other coarsening strategies (with respect to cfmarker != 1)
         ICOUNT = 1
         FALGOUT_CELL_LOOP: DO IC = 1, SL%NC
WRITE(SCARC_LU,*) 'PUH1:  IC=',IC
            IF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE) THEN
               IF (SL%MEASURE(IC) >= 1.0_EB .AND. SL%A_ROW(IC+1)-SL%A_ROW(IC)+1 > 0) THEN
                  SL%CELL_TYPE(IC)  = NSCARC_CELL_TYPE_NONE
WRITE(SCARC_LU,*) 'PUH2:  SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
                  SL%GRAPH(ICOUNT) = IC
WRITE(SCARC_LU,*) 'C:  SL%GRAPH(',ICOUNT,')=',SL%GRAPH(ICOUNT)
                  ICOUNT = ICOUNT + 1
                  IF (ICOUNT > SL%NC) THEN
                     WRITE(*,*) 'WARNING: GRAPH_SIZE MUST BE INCREASED!!'
                     STOP
                  ENDIF
               ELSE
                  SL%CELL_TYPE(IC)  = NSCARC_CELL_TYPE_FINE
WRITE(SCARC_LU,*) 'D: SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
               ENDIF
            ELSEIF (SL%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE0) THEN
               SL%MEASURE(IC) = NSCARC_MEASURE_NONE
            ELSE
               SL%GRAPH(ICOUNT) = IC
WRITE(SCARC_LU,*) 'F:  SL%GRAPH(',ICOUNT,')=',SL%GRAPH(ICOUNT)
               ICOUNT = ICOUNT + 1
            ENDIF
         ENDDO FALGOUT_CELL_LOOP

WRITE(SCARC_LU,*) '-----------------------------------------------'
DO IC = 1, SL%NC
   WRITE(SCARC_LU,*) 'E: SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
ENDDO
WRITE(SCARC_LU,*) '-----------------------------------------------'
DO IC = 1, SL%NC
   WRITE(SCARC_LU,*) 'E: SL%MEASURE(',IC,')=',SL%MEASURE(IC)
ENDDO
WRITE(SCARC_LU,*) '-----------------------------------------------'
DO IC = 1, SL%NC
   WRITE(SCARC_LU,*) 'E: SL%GRAPH(',IC,')=',SL%GRAPH(IC)
ENDDO
WRITE(SCARC_LU,*) '-----------------------------------------------'
DO IC = 1, SL%NC
   WRITE(SCARC_LU,*) 'E: SL%A_ROW(',IC,')=',SL%A_ROW(IC)
ENDDO
WRITE(SCARC_LU,*) '-----------------------------------------------'
DO IC = 1, SL%NA
   WRITE(SCARC_LU,*) 'E: SL%S_COL(',IC,')=',SL%S_COL(IC)
ENDDO
WRITE(SCARC_LU,*) '-----------------------------------------------'

         IGRAPH = ICOUNT
         ILOOP  = 0
         FALGOUT_COARSENING_LOOP: DO

WRITE(SCARC_LU,*) 'iter =',ILOOP

            !< Care: Is graph_size always SL%NC ?
            IF (ILOOP > 0) THEN

WRITE(SCARC_LU,*) 'HALLO2 ILOOP=',ILOOP

               DO IG = 1, IGRAPH

                  IC = SL%GRAPH(IG)
                  JC = IG

                  !< make IC a fine cell and look for all dependencies
                  IF ((SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_COARSE) .AND. (SL%MEASURE(IC) < 1)) THEN
                     SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE          
WRITE(SCARC_LU,*) 'Setting %CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
                     DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1       
                        IF (SL%S_COL(IROW) > -1) SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_NONE
WRITE(SCARC_LU,*) ' jS: ',IROW, ' Setting SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
                     ENDDO
                  ENDIF

                  !< remove cell from graph
                  IF (SL%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_NONE) THEN
                     SL%MEASURE(JC) = NSCARC_MEASURE_NONE
                     IGRAPH = IGRAPH - 1
                     SL%GRAPH(JC) = SL%GRAPH(IGRAPH)
                     SL%GRAPH(IGRAPH) = IC
WRITE(SCARC_LU,*) 'graph_size=',IGRAPH,' graph_array(',JC,')=',SL%GRAPH(JC),': graph_array(',IGRAPH,')=',SL%GRAPH(IGRAPH)
                     JC = JC - 1               !< correct ???
                  ENDIF

               ENDDO
   
            ENDIF
            ILOOP = ILOOP + 1

            FALGOUT_GRAPH_LOOP: DO IG = 1, IGRAPH

               IC = SL%GRAPH(IG)

WRITE(SCARC_LU,*) 'SUSI ===== IC = ',IC,' ===================================================================='
DO JC = 1, SL%NC
   WRITE(SCARC_LU,*) 'E: SL%GRAPH(',JC,')=',SL%GRAPH(JC)
ENDDO
WRITE(SCARC_LU,*) '------------------------------------------------------------------------------'
WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,10i3)') 'before i=',IC,&
                                 ': CF_marker(',IC,')=',SL%CELL_TYPE(IC),&
                                 ': graph_array(',ic,')=',SL%GRAPH(IC), &
                                 ': S_diag_i(',IC,')=',SL%A_ROW(IC)
WRITE(SCARC_LU,*) '------------------------------------------------------------------------------'

               !< C-points are not influenced by neighbors (no interpolation)
               IF (SL%CELL_TYPE(IC) > 0) THEN

WRITE(SCARC_LU,*) '================== COARSE'
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE            !< define it as C-point

WRITE(SCARC_LU,*) 'A0: SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)

                  DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1       
                     JC = SL%S_COL(IROW)
WRITE(SCARC_LU,*) 'A1: --------- IROW=',IROW,': JC=', JC 
                     IF (JC > -1) THEN
                        SL%S_COL(IROW) = -SL%S_COL(IROW) - 1  !< remove edge 
WRITE(SCARC_LU,*) 'A2: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
WRITE(SCARC_LU,*) 'A3: SL%CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
                        IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_NONE) THEN
WRITE(SCARC_LU,*) 'A4: SL%MEASURE (',JC,')=',SL%MEASURE(JC)
                           SL%MEASURE(JC) = SL%MEASURE(JC) - 1.0_EB   !< decrement measure of unmarked neigh.
WRITE(SCARC_LU,*) 'A5: SL%MEASURE (',JC,')=',SL%MEASURE(JC)
                        ENDIF
                     ENDIF
                  ENDDO

               ELSE

WRITE(SCARC_LU,*) '================== FINE'
                  !< dependencies which have already been marked
                  DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1       
                     JC = SL%S_COL(IROW)
WRITE(SCARC_LU,*) 'C1: --------- IROW=',IROW,': JC=', JC 
                     IF (JC < 0) JC = -JC-1
WRITE(SCARC_LU,*) 'C2: --------- IROW=',IROW,': JC=', JC 
WRITE(SCARC_LU,*) 'C3: SL%CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
                     IF (SL%CELL_TYPE(JC) > 0) THEN
WRITE(SCARC_LU,*) 'C4: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                        IF (SL%S_COL(IROW) > -1) THEN
                           SL%S_COL(IROW) = -SL%S_COL(IROW) -1   !< remove edge
WRITE(SCARC_LU,*) 'C5: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                           SL%CELL_TYPE(JC)    = NSCARC_CELL_TYPE_COMMON   !< temporarily set CELL_TYPE to common
WRITE(SCARC_LU,*) 'C6: SL%CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
                        ENDIF
                     ELSEIF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_FINE0) THEN
                        IF (SL%S_COL(IROW) > -1) THEN
                           SL%S_COL(IROW) = -SL%S_COL(IROW) -1   !< remove edge
WRITE(SCARC_LU,*) 'C7: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                        ENDIF
                     ENDIF
                  ENDDO

                  !< dependencies which haven't been marked yet
                  DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1       
                     IF (SL%S_COL(IROW) > -1) THEN
WRITE(SCARC_LU,*) 'E: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                        BREMOVE = .TRUE.
                        JC = SL%S_COL(IROW)
                        DO JROW = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1 
                           KC = SL%S_COL(JROW)
WRITE(SCARC_LU,*) 'E1: SL%S_COL(',JROW,')=',SL%S_COL(JROW)
                           IF (KC < 0) KC = -KC-1                        !< check for all dependencies !!
                           IF (SL%CELL_TYPE(KC) == NSCARC_CELL_TYPE_COMMON) THEN
                              SL%S_COL(IROW) = -SL%S_COL(IROW)-1
WRITE(SCARC_LU,*) 'E2: SL%S_COL(',IROW,')=',SL%S_COL(IROW)
                              SL%MEASURE(JC) = SL%MEASURE(JC) - 1.0_EB
WRITE(SCARC_LU,*) 'E3: SL%MEASURE (',JC,')=',SL%MEASURE(JC)
                              BREMOVE = .FALSE.
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO

               ENDIF

               !< reset CELL_TYPES
               DO IROW = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1       
                  JC = SL%S_COL(IROW)
WRITE(SCARC_LU,*) 'H1: jS = ', IROW,' j=',JC
                  IF (JC < 0) JC = -JC-1
WRITE(SCARC_LU,*) 'H2: jS = ', IROW,' j=',JC
                  IF (SL%CELL_TYPE(JC) == NSCARC_CELL_TYPE_COMMON) SL%CELL_TYPE(JC)=NSCARC_CELL_TYPE_COARSE
WRITE(SCARC_LU,*) 'H3: SL%CELL_TYPE(',JC,')=',SL%CELL_TYPE(JC)
               ENDDO

WRITE(SCARC_LU,*) '------------------------------------------------------------------------------'
WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,10i3)') 'after  i=',IC,&
                                 ': CF_marker(',IC,')=',SL%CELL_TYPE(IC),&
                                 ': graph_array(',ic,')=',SL%GRAPH(IC), &
                                 ': S_diag_i(',IC,')=',SL%A_ROW(IC)
            ENDDO FALGOUT_GRAPH_LOOP

WRITE(SCARC_LU,*) 'AFTER RESET'
          
WRITE(SCARC_LU,*) '-----------------------------------------------'
DO IC = 1, SL%NC
   WRITE(SCARC_LU,*) 'Z: SL%CELL_TYPE(',IC,')=',SL%CELL_TYPE(IC)
ENDDO
WRITE(SCARC_LU,*) '-----------------------------------------------'
DO IC = 1, SL%NC
   WRITE(SCARC_LU,*) 'Z: SL%MEASURE(',IC,')=',SL%MEASURE(IC)
ENDDO

         ENDDO FALGOUT_COARSENING_LOOP
      
      ENDDO FALGOUT_LOOP

   !! ---------------------------------------------------------------------------------------------
   !! set CELL_TYPEs for GMG-like interpolation
   !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_GMG)

      SELECT_GMG_DIMENSION: SELECT CASE (TYPE_DIMENSION)

         CASE (NSCARC_DIMENSION_TWO)

            GMG_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

               SLF => SCARC(NM)%LEVEL(NL)
               SLC => SCARC(NM)%LEVEL(NL+1)

               !< Analyze grid sizes
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

            GMG3_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

               SLF => SCARC(NM)%LEVEL(NL)
               SLC => SCARC(NM)%LEVEL(NL+1)

               !< Analyze grid sizes
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
                  WRITE(*,*) 'CELL NUNMBERS NOT DIVISABLE BY 2 or 3, STOP!'
                  STOP
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

      BDRY_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SL => SCARC(NM)%LEVEL(NL)
      
         !< First define coarse cells on boundary layer(adjacent to neighboring meshes)
         BDRY_CYCLE_LOOP: DO
         
            !!!
            !< get maximum (remaining) measure for all cells
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

               !< Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, SL%MEASURE(IC))) THEN
         
                  SL%MEASURE(IC)  = 0.0_EB
                  SL%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
         
                  !< Determine set of fine cells 
                  DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
         
                   !< JC is set to be a fine cell which is no longer measured
                     JC = SL%A_COL(ICOL)
         
                     SL%MEASURE(JC)  = 0.0_EB
                     SL%CELL_TYPE(JC) = NSCARC_CELL_TYPE_SFINE
         
                     !<  increase measures of cells KC adjacent to fine cells JC based on strong couplings
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


SELECT_INTERPOLATION: SELECT CASE (TYPE_INTERPOL)

   !! ------------------------------------------------------------------------------------------
   !! Standard interpolation
   !! ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_STANDARD)

      MESHES_LOOP_STANDARD: DO NM = NMESHES_MIN, NMESHES_MAX
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

      MESHES_LOOP_GMG: DO NM = NMESHES_MIN, NMESHES_MAX
      
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

      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NL, 'SETUP_PROLONGATION', 'Matrix structure ')
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION','CELL_TYPE INIT')

      MESHES_LOOP_CLASSICAL: DO NM = NMESHES_MIN, NMESHES_MAX

         SL => SCARC(NM)%LEVEL(NL)
         
         IP = 1
         INTERNAL_CELL_LOOP_CLASSICAL: DO IC = 1, SL%NC
         
            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0
         
            !!!
            !< If IC is a coarse cell, its value is taken
            !!!
            IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
         
               NEIGHBOR(1)= SL%CELL_TYPE(IC)
               WEIGHTS(1) = 1.0_EB
               NWEIGHTS   = 1

            !!!
            !< If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
            !!!
            ELSE 
         
             !< Get main diagonal entry a_ii for that fine cell
               IDIAG = SL%A_ROW(IC)
               SUM_DIAG = SL%A(IDIAG)
         
               !< First search for all neighboring coarse grid cells (store them in NEIGHBOR)
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

               !< Then search for the strongly and weakly coupled fine grid cells
               DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1

                  IW2 = 1
                  SUM_COARSE = 0.0_EB

                  JC = SL%A_COL(ICOL)

                  SELECT CASE (SL%CELL_TYPE(JC))

                     CASE (NSCARC_CELL_TYPE_SFINE)

                        !< search for couplings KC of the strongly coupled JC which belong to the 
                        !< coarse interpolatory set of IC
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

            !< 
            !< Define corresponding entry in interpolation matrix P by means of the upper weights
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

      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NL, 'SETUP_PROLONGATION', 'TYPE2: Matrix structure ')
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION','TYPE2: CELL_TYPE INIT')

      MESHES_LOOP_CLASSICAL2: DO NM = NMESHES_MIN, NMESHES_MAX

         SL => SCARC(NM)%LEVEL(NL)
         
         !< Build Interpolation
         ICOL0 = 1
         INTERNAL_CELL_LOOP2_CLASSICAL2: DO IC = 1, SL%NC
     
          !< For a C-point IC, use identity and set mapping
            IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               SL%P_ROW(IC) = ICOL0
               SL%P_COL(ICOL0) = SL%CELL_TYPE(IC)
               SL%P(ICOL0) = 1.0_EB

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)') &
      'C: IC=', IC,': P_ROW(',IC,')=', SL%P_ROW(IC),': P_COL(',ICOL0,')=',SL%P_COL(ICOL0), &
      ': P(',ICOL0,')==',SL%P(ICOL0)

               ICOL0 = ICOL0 + 1

            !< For a F-point IC, build interpolation
            ELSE

               SL%P_ROW(IC)= ICOL0               !< diagonal part of P
               ICOL_FIRST = ICOL0

               DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                  JC = SL%A_COL(ICOL)

                  !< If JC is a C-point, initialize interpolation to zero and set column number in P_COL
                  IF (SL%CELL_TYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     SL%CELL_MARKER(JC) = ICOL0
                     SL%P_COL(ICOL0) = SL%CELL_TYPE(JC)
                     SL%P(ICOL0) = 0.0_EB
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)') &
      'D: IC=', IC,': CELL_MARKER(',JC,')=', SL%CELL_MARKER(JC),': P_COL(',ICOL0,')=',SL%P_COL(ICOL0), &
      ': P(',ICOL0,')==',SL%P(ICOL0)
                     ICOL0 = ICOL0 + 1
                  !< If JC is a F-point, set it to be a strong F-point with relevant connections
                  ELSE IF (SL%CELL_TYPE(JC) /= NSCARC_CELL_TYPE_WFINE) THEN
                     SL%CELL_MARKER(JC) = NSCARC_CELL_TYPE_SFINE
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)') &
      'E: IC=', IC,': CELL_MARKER(',JC,')=', SL%CELL_MARKER(JC)
                  ENDIF

               ENDDO
               ICOL_LAST = ICOL0 - 1

               !< consider the i-th row of A, start with diagonal part
               IDIAG = SL%A_ROW(IC)
               DATA_DIAG = SL%A(IDIAG)
   
               DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                  JC = SL%A_COL(ICOL)
                  
                  !< ---------------------------------------------------------------------------------------
                  !< If JC is a strongly coupled C-point to IC, the value a_(IC, JC) must be incorporated to the 
                  !< interpolation weight
                  IF (SL%CELL_MARKER(JC) >= ICOL_FIRST) THEN
                     JCOL = SL%CELL_MARKER(JC)
                     SL%P(JCOL) = SL%P(JCOL) + SL%A(ICOL)
   
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,f12.6,a,i3,a,i3,a,i3)') &
      'G: IC=', IC,': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,': P(',JCOL,')=',SL%P(JCOL)

                  !< ---------------------------------------------------------------------------------------
                  !< If JC is a strongly coupled F-point to IC, the value a_(IC, JC) must be distributed to C-points
                  !< which are strongly coupled to IC (no distribution to the diagonal part)
                  ELSE IF (SL%CELL_MARKER(JC) >= NSCARC_CELL_TYPE_SFINE) THEN
                     DATA_SUM  = 0.0_EB
                     SIGN0 = 1
                     IF (SL%A(IDIAG) > 0) SIGN0=-1
            
                     !< search JC in the row of A and get sum of the couplings to C-points strongly coupled to IC
                     DO JCOL = SL%A_ROW(JC), SL%A_ROW(JC+1)-1
                        KC = SL%A_COL(JCOL)
                        IF (SL%CELL_MARKER(KC) > ICOL_FIRST .AND. SIGN0*SL%A(JCOL) >0) THEN
                           DATA_SUM = DATA_SUM + SL%A(JCOL)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,f12.6,a,i3,a,i3,a,i3)') &
      'H: IC=', IC,': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,': KC=',KC,': DATA_SUM=',DATA_SUM
                        ENDIF
                     ENDDO
   
                     IF (DATA_SUM .NE. 0) THEN
      
                        DATA_INTERPOL = SL%A(ICOL)/DATA_SUM
                     
                        !< loop over row of A for JC and spread data
                        DO JCOL = SL%A_ROW(JC), SL%A_ROW(JC+1)-1
                           KC = SL%A_COL(JCOL)
                           IF (SL%CELL_MARKER(KC) >= ICOL_FIRST .AND. SIGN0*SL%A(JCOL) > 0) THEN
                              KCOL = SL%CELL_MARKER(KC)
                              SL%P(KCOL) = SL%P(KCOL) + DATA_INTERPOL * SL%A(JCOL)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3, a, i3,a,i3,a,i3,a,i3,a,i3,a,f12.6,a,i3,a,i3)') &
      'I: IC=', IC,': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,': KC=',KC,': P(',KCOL,')=',SL%P(KCOL)
                           ENDIF
                        ENDDO
      
      
                     ENDIF
   
                  !< ---------------------------------------------------------------------------------------
                  !< If JC is a weakly coupled F-point to IC, the value a_(IC, JC) must be distributed to the diagonal
                  ELSE IF (SL%CELL_TYPE(JC) .NE. NSCARC_CELL_TYPE_WFINE) THEN
                     DATA_DIAG = DATA_DIAG + SL%A(ICOL)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3, a, i3,a,i3,a,i3,a,i3,a,f12.6,a,i3,a,i3,a,i3)') &
      'J: IC=', IC,': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,': KC=',KC,': DATA_DIAG=',DATA_DIAG
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
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6,a,i3,a,i3)') &
      'K: IC=', IC,': P(',ICOL,')=',SL%P(ICOL)
                  ENDDO
               ENDIF
            ENDIF
            
         ENDDO INTERNAL_CELL_LOOP2_CLASSICAL2
         SL%P_ROW(IC) = ICOL0
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)') &
      'F: IC=', IC,': P_ROW(',IC,')=', SL%P_ROW(IC)

 
         
      ENDDO MESHES_LOOP_CLASSICAL2

      WRITE(SCARC_LU,*) 'SIZE(P)=',SIZE(SL%P)
      WRITE(SCARC_LU,*) 'SIZE(P_ROW)=',SIZE(SL%P_ROW)
      WRITE(SCARC_LU,*) 'SIZE(P_COL)=',SIZE(SL%P_COL)
      WRITE(SCARC_LU,*) 'SIZE(R)=',SIZE(SL%R)
      WRITE(SCARC_LU,*) 'SIZE(R_ROW)=',SIZE(SL%R_ROW)
      WRITE(SCARC_LU,*) 'SIZE(R_COL)=',SIZE(SL%R_COL)
      WRITE(SCARC_LU,*) 'P:'
      WRITE(SCARC_LU,'(f12.6)') (SL%P(IC), IC=1, SIZE(SL%P))
      WRITE(SCARC_LU,*) 'P_ROW:'
      WRITE(SCARC_LU,'(i4)') (SL%P_ROW(IC), IC=1, SIZE(SL%P_ROW))
      WRITE(SCARC_LU,*) 'P_COL:'
      WRITE(SCARC_LU,'(i4)') (SL%P_COL(IC), IC=1, SIZE(SL%P_COL))



      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'TYPE2: final prolongation ')
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_RESTRICTION, NL, 'SETUP_PROLONGATION', 'TYPE2: final restriction ')
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION','TYPE2: final CELL_TYPE INIT')

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)

   !! ---------------------------------------------------------------------------------------------
   !! Direct interpolation
   !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_DIRECT)

      MESHES_LOOP_DIRECT: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SL => SCARC(NM)%LEVEL(NL)
         CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION','CELL_TYPE INIT')
         
         IP = 1
         INTERNAL_CELL_LOOP: DO IC = 1, SL%NC
         
            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0
         
            !!!
            !< If IC is a coarse cell, its value is taken
            !!!
            IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
         
               NEIGHBOR(1)= SL%CELL_TYPE(IC)
               NEIGHBOR(1)= IC
               WEIGHTS(1) = 1.0_EB
               NWEIGHTS   = 1
         
            !!!
            !< If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
            !!!
            ELSE 
         
               !< Get main diagonal entry a_ii for that fine cell
               IDIAG = SL%A_ROW(IC)
         
               !< Select type of fine cell (weakly/strongly coupled)
               SELECT_FPOINT_TYPE: SELECT CASE(SL%CELL_TYPE(IC))

                  !!!
                  !< Strongly coupled fine cell IC
                  !< approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                  !< Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                  !!!
                  CASE (NSCARC_CELL_TYPE_SFINE)
         
                     !< Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                     SUM_COUPLED = 0.0_EB
                     SUM_CPOINTS = 0.0_EB
         
                     DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                        JC = SL%A_COL(ICOL)
                        SUM_COUPLED = SUM_COUPLED + SL%A(ICOL)
                        IF (SL%CELL_TYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + SL%A(ICOL)
                     ENDDO
         
                     SCAL = - SUM_COUPLED/SUM_CPOINTS
         
                     !< for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
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
                  !< Weakly coupled fine cell IC:
                  !< Determine strongly coupled fine cells JC surrounding IC and, in turn, replace 
                  !< each of them by a mean value of their surrounding strongly coupled coarse cells
                  !!!
                  CASE (NSCARC_CELL_TYPE_WFINE)
         
                     IW = 1                                             !< weights counter
                     DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1         !< loop over couplings of weakly coupled fine cell
         
                        !< regard subdiagonal matrix entries a_ij of weakly coupled fine cell JC to IC
                        JC = SL%A_COL(ICOL)
         
                        !< Find all surrounding (coupled) points of JC and compute scaling factor 
                        !< compute scaling factor SCAL = a_ij/a_ii * [sum_(k in N_j) a_jk]/[sum_(l in C_j) ajl] 
                        SUM_COUPLED = 0.0_EB
                        SUM_CPOINTS = 0.0_EB
         
                        DO JCOL = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                           KC = SL%A_COL(JCOL)
                           SUM_COUPLED = SUM_COUPLED + SL%A(JCOL)
                           IF (SL%CELL_TYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + SL%A(JCOL)
                        ENDDO
         
                        IF (SUM_CPOINTS == 0.0_EB) THEN
                           WRITE(*,*) 'Error in SCARC_INTERPOLATION on mesh ',NM,' for cell ',IC,' with neighbor ', JC, &
                                      'on level ', NL,' stopping program!'
                           WRITE(*,*) '  SUM_COUPLED =',SUM_COUPLED
                           WRITE(*,*) '  SUM_CPOINTS =',SUM_CPOINTS
                           STOP
                        ENDIF
                        SCAL =  SL%A(ICOL)/SL%A(IDIAG) * SUM_COUPLED/SUM_CPOINTS
                            
                        !< Get diagonal matrix a_jj for point JC
                        JDIAG = SL%A_ROW(JC)
         
                        !< Compute interpolation weights for all strong coarse cell couplings KC of JC
                        !< note that a coarse cell KC may be considered several times for different JC's
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
         
                     !< make weights unique (add weights for multiple coarse cells)
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
         
         
            !< 
            !< Define corresponding entry in interpolation matrix P by means of the upper weights
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
         !SL%NROW0 = IP                !< really neede - OUTDATED ?
         
         !< Determine new coarse number and store it in CELL_TYPE
         INTERNAL_CELL_LOOP2: DO IC = 1, SL%NC
            DO ICOL = SL%P_ROW(IC), SL%P_ROW(IC+1)-1
              JC = SL%P_COL(ICOL)
              IF (JC <= SL%NC) THEN                 
                 SL%P_COL(ICOL) = SL%CELL_TYPE(JC)         
              ENDIF
            ENDDO
         ENDDO INTERNAL_CELL_LOOP2


         !< In the multi-mesh case determine the new coarse cell numbers for neighboring cells
         IF (NMESHES > 1) THEN
            GHOST_CELL_LOOP: DO IW = 1, SL%NW
               IF (SL%WALL(IW)%NOM==0) CYCLE GHOST_CELL_LOOP
               ICG = SL%WALL(IW)%ICG(1)
               ICA = SL%WALL(IW)%ICW
               DO ICOL = SL%P_ROW(ICA), SL%P_ROW(ICA+1)-1
                 JC = SL%P_COL(ICOL)
                 IF (JC > SL%NC) THEN
                    JCO = SL%WALL(IW)%ICN(1)
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

      MESHES_LOOP_DIRECT_BDRY: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SL => SCARC(NM)%LEVEL(NL)
         CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELL_TYPE, NL, 'SETUP_PROLONGATION','CELL_TYPE INIT')
         
         IP = 1
         INTERNAL_CELL_BDRY_LOOP: DO IC = 1, SL%NC
         
            DIRECT_BDRY_IF: IF (SL%INTERNAL_BDRY_CELL(IC)/=0) THEN

               VALUES   = 0.0_EB
               WEIGHTS  = 0.0_EB
               NEIGHBOR = 0
         
               !!!
               !< If IC is a coarse cell, its value is taken
               !!!
               IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
         
                  NEIGHBOR(1)= SL%CELL_TYPE(IC)
                  NEIGHBOR(1)= IC
                  WEIGHTS(1) = 1.0_EB
                  NWEIGHTS   = 1
         

               !!!
               !< If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
               !!!
               ELSE 


                  !< Get main diagonal entry a_ii for that fine cell
                  IDIAG = SL%A_ROW(IC)
         
                  !< Strongly coupled fine cell IC
                  !< approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                  !< Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                  !< Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                  SUM_COUPLED = 0.0_EB
                  SUM_CPOINTS = 0.0_EB
         
                  DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                     JC = SL%A_COL(ICOL)
                     SUM_COUPLED = SUM_COUPLED + SL%A(ICOL)
                     IF (JC > SL%NC) CYCLE
                     IF (SL%INTERNAL_BDRY_CELL(JC)/=0.AND.SL%CELL_TYPE(JC)>0) SUM_CPOINTS = SUM_CPOINTS + SL%A(ICOL)
                  ENDDO
         
                  SCAL = - SUM_COUPLED/SUM_CPOINTS
         
                  !< for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
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
               !< If IC is a coarse cell, its value is taken
               !!!
               IF (SL%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
         
                  NEIGHBOR(1)= SL%CELL_TYPE(IC)
                  NEIGHBOR(1)= IC
                  WEIGHTS(1) = 1.0_EB
                  NWEIGHTS   = 1
         
               !!!
               !< If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
               !!!
               ELSE 
         
                  !< Get main diagonal entry a_ii for that fine cell
                  IDIAG = SL%A_ROW(IC)
         
                  !< Select type of fine cell (weakly/strongly coupled)
                  SELECT_FPOINT_BDRY_TYPE: SELECT CASE(SL%CELL_TYPE(IC))

                     !!!
                     !< Strongly coupled fine cell IC
                     !< approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                     !< Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                     !!!
                     CASE (NSCARC_CELL_TYPE_SFINE)
         
                        !< Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                        SUM_COUPLED = 0.0_EB
                        SUM_CPOINTS = 0.0_EB
         
                        DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                           JC = SL%A_COL(ICOL)
                           SUM_COUPLED = SUM_COUPLED + SL%A(ICOL)
                           IF (SL%CELL_TYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + SL%A(ICOL)
                        ENDDO
         
                        SCAL = - SUM_COUPLED/SUM_CPOINTS
         
                        !< for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
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
                     !< Weakly coupled fine cell IC:
                     !< Determine strongly coupled fine cells JC surrounding IC and, in turn, replace 
                     !< each of them by a mean value of their surrounding strongly coupled coarse cells
                     !!!
                     CASE (NSCARC_CELL_TYPE_WFINE)
            
                        IW = 1                                             !< weights counter
                        DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1         !< loop over couplings of weakly coupled fine cell
         
                           !< regard subdiagonal matrix entries a_ij of weakly coupled fine cell JC to IC
                           JC = SL%A_COL(ICOL)
            
                           !< Find all surrounding (coupled) points of JC and compute scaling factor 
                           !< compute scaling factor SCAL = a_ij/a_ii * [sum_(k in N_j) a_jk]/[sum_(l in C_j) ajl] 
                           SUM_COUPLED = 0.0_EB
                           SUM_CPOINTS = 0.0_EB
            
                           DO JCOL = SL%A_ROW(JC)+1, SL%A_ROW(JC+1)-1
                              KC = SL%A_COL(JCOL)
                              SUM_COUPLED = SUM_COUPLED + SL%A(JCOL)
                              IF (SL%CELL_TYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + SL%A(JCOL)
                           ENDDO
         
                           IF (SUM_CPOINTS == 0.0_EB) THEN
                              WRITE(*,*) 'Error in SCARC_INTERPOLATION on mesh ',NM,' for cell ',IC,' with neighbor ', JC, &
                                         'on level ', NL,' stopping program!'
                              WRITE(*,*) '  SUM_COUPLED =',SUM_COUPLED
                              WRITE(*,*) '  SUM_CPOINTS =',SUM_CPOINTS
                              STOP
                           ENDIF
                           SCAL =  SL%A(ICOL)/SL%A(IDIAG) * SUM_COUPLED/SUM_CPOINTS
                            
                           !< Get diagonal matrix a_jj for point JC
                           JDIAG = SL%A_ROW(JC)
         
                           !< Compute interpolation weights for all strong coarse cell couplings KC of JC
                           !< note that a coarse cell KC may be considered several times for different JC's
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
         
                        !< make weights unique (add weights for multiple coarse cells)
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
            !< 
            !< Define corresponding entry in interpolation matrix P by means of the upper weights
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
         !SL%NROW0 = IP                !< really neede - OUTDATED ?
         
         !< Determine new coarse number and store it in CELL_TYPE
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

      MESHES_LOOP_MULTIPASS: DO NM = NMESHES_MIN, NMESHES_MAX
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
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i4,a,7i4)') 'IZ=',IZ,': PZ=',PZ(1:INCRZ), NY1, NY2, INCRY

DO IX = NX1, NX2, INCRX

   DO INX = 1, INCRX
      PX(INX) = IX + INX - 1
   ENDDO

   R_ROW(ICC) = IP

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i4,a,4i4)') 'IX=',IX,': PX=',PX(1:INCRX)
   DO INZ = 1, INCRZ
      DO INX = 1, INCRX
         PC(INX,INZ) = (PZ(INZ)-1)*NX + PX(INX)
         IF (CELL_TYPE(PC(INX,INZ))== NSCARC_CELL_TYPE_COARSE) CELL_TYPE(PC(INX,INZ))=ICC
         R_COL(IP)   = PC(INX,INZ)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i4,a,i4,a,i4)') 'PC(',INX,',',INZ,')=',PC(INX,INZ)
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
SUBROUTINE SCARC_MATRIX_TRANSPOSE(M, M_ROW, M_COL, MT, MT_ROW, MT_COL, NC, NCE, NCC, NCCE)
INTEGER,  DIMENSION(:)  , INTENT(IN)  :: M_ROW , M_COL
INTEGER,  DIMENSION(:)  , INTENT(OUT) :: MT_ROW, MT_COL
REAL(EB), DIMENSION(:)  , INTENT(IN)  :: M
REAL(EB), DIMENSION(:)  , INTENT(OUT) :: MT
INTEGER, INTENT(IN)  :: NC, NCE, NCC, NCCE
INTEGER :: IC, JC, ICOL, JCOL

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,*) 'STARTING MATRIX_TRANSPOSE'
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'NC=',NC
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'NCE=',NCE
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'NCC=',NCC
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'NCCE=',NCCE
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'M_ROW'
   WRITE(SCARC_LU,'(8i5)') (M_ROW(IC), IC=1, NCE+1)
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'M_COL'
   WRITE(SCARC_LU,'(8i5)') (M_COL(IC), IC=1, NCE)
ENDIF

!> identify the number of non-zero entries in every column of M (corresponds to a row in MT) 
!> And store it in the MT_ROW-array (caution: starts from the second position)
MT_ROW(1) = 1
DO ICOL = 1, M_ROW(NCE+1)-1
   IC = M_COL(ICOL)
   MT_ROW(IC+1) = MT_ROW(IC+1) + 1
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,i3)') 'ICOL=',ICOL, ': MT_ROW(',IC,')=',MT_ROW(IC)
ENDDO

IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'MT_ROW 2:'
   WRITE(SCARC_LU,'(8i5)') (MT_ROW(IC), IC=1, NCE+1)
ENDIF

!> shift it to the first position while summing up the entries for M_ROW
DO IC = 2, NCCE+1
   MT_ROW(IC) = MT_ROW(IC) + MT_ROW(IC-1)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,i3,i3)') 'MT_ROW(',IC,')=',MT_ROW(IC), MT_ROW(IC-1)
ENDDO

IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'MT_ROW 3:'
   WRITE(SCARC_LU,'(8i5)') (MT_ROW(IC), IC=1, NCCE+1)
ENDIF

!> set correct entries for MT
DO IC = 1, NCE
   DO ICOL = M_ROW(IC), M_ROW(IC+1)-1
      JC   = M_COL(ICOL)
      JCOL = MT_ROW(JC)
!IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(8i4,a,i3,a,i3)') (MT_ROW(I),I=1,8),'JCOL=',JCOL,': JC=',JC
      MT_COL(JCOL) = IC
      MT(JCOL)     = M(ICOL) 
      MT_ROW(JC)   = MT_ROW(JC) + 1
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3, a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)') &
          'IC=',IC, ': ICOL=',ICOL,': JC=',JC,': JCOL=',JCOL,&
          ': MT_COL(',JCOL,')=', MT_COL(JCOL), ': MT_ROW(',JCOL,')=',  MT_ROW(JC), ': MT(',JCOL,')=',MT(JCOL)
   ENDDO
ENDDO

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
   WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'MT_ROW 4:'
   WRITE(SCARC_LU,'(8i5)') (MT_ROW(IC), IC=1, NCCE+1)
ENDIF

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'IC=',NCCE+1, ': MT_ROW(',NCCE+1,')=',MT_ROW(NCCE+1)
DO IC = NCCE, 2, -1
   MT_ROW(IC) = MT_ROW(IC-1)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'IC=',IC, ': MT_ROW(',IC,')=',MT_ROW(IC)
ENDDO
MT_ROW(1)=1
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'IC=',1, ': MT_ROW(',1,')=',MT_ROW(1)


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

MESHES_LOOP_SYSTEM_AMG: DO NM = NMESHES_MIN, NMESHES_MAX

   SLF => SCARC(NM)%LEVEL(NL)
   SLC => SCARC(NM)%LEVEL(NL+1)

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SCARC_SETUP_SYSTEM_AMG: NC  =', SLF%NC 
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SCARC_SETUP_SYSTEM_AMG: NCC =', SLF%NCC

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
!>      !< loop over all entries in row IC0 of R
!>        LOOP_NBR_R_ENTRIES: DO ICOL1 = SLF%R_ROW(IC0), SLF%R_ROW(IC0+1)-1
!>  
!>           IC1 = SLF%R_COL(ICOL1) 
!>           IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '--------- ICOL1 =', ICOL1,'   : IC1=',IC1
!>           
!>         !< loop over all entries in row IC1 of A
!>           LOOP_NBR_A_ENTRIES: DO ICOL2 = SLF%A_ROW(IC1), SLF%A_ROW(IC1+1)-1
!>  
!>              IC2 = SLF%A_COL(ICOL2)
!>              IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '   ------ ICOL2 =', ICOL2,'   : IC2=',IC2
!>     
!>              IF (SLF%A_TAG(IC2) /= IC0) THEN
!>     
!>               !< mark IC2 as already considered
!>                 SLF%A_TAG(IC2) = IC0
!>                 IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '          A_TAG(',IC2,') =', SLF%A_TAG(IC2)
!>       
!>               !< loop over all entries in row IC2 of P
!>                 LOOP_NBR_P_ENTRIES: DO ICOL3 = SLF%P_ROW(IC2), SLF%P_ROW(IC2+1)-1
!>     
!>                    IC3 = SLF%P_COL(ICOL3)
!>                    IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '    !> ---- ICOL3 =', ICOL3,'   : IC3=',IC3
!>                    !IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
!>                  !<   '          P_TAG(',IC3,') =', SLF%P_TAG(IC3),': IROW_INIT=',IROW_INIT,' : IROW=',IROW
!>     
!>                  !< verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
!>                    IF (SLF%P_TAG(IC3) < IROW_INIT) THEN
!>                       SLF%P_TAG(IC3) = IROW
!>                       IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '          P_TAG(',IC3,') =', SLF%P_TAG(IC3)
!>                       IROW = IROW + 1
!>                       IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '          IROW =', IROW
!>                    ENDIF
!>                 ENDDO LOOP_NBR_P_ENTRIES
!>     
!>              ENDIF
!>  
!>           ENDDO LOOP_NBR_A_ENTRIES
!>        ENDDO LOOP_NBR_R_ENTRIES
!>  
!>      !< Store counters
!>        IROW_SAVE = IROW
!>  
!>
!>     ENDDO LOOP_NBR_CELLS
!>
!>  ENDIF


   !! Determine size of matrix RAP
   !! loop over interior c-cells
   LOOP1_OWN_CELLS: DO IC0 = 1, SLF%NCCE

      IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) THEN
         WRITE(SCARC_LU,*) '================ MARKMARK1:  IC0 =', IC0-1
         WRITE(SCARC_LU,*) 'IROW     = ', IROW
         WRITE(SCARC_LU,*) 'IROW_INIT= ', IROW_INIT
         WRITE(SCARC_LU,*) '================ A_MARKER'
         DO IC1 = 1, SLF%NCE
            WRITE(SCARC_LU,'(a,i4,a,i4)') 'A_MARKER(',IC1,')=', SLF%A_TAG(IC1)-1
         ENDDO
         WRITE(SCARC_LU,*) '================ P_MARKER'
         DO IC1 = 1, SLF%NCCE
            WRITE(SCARC_LU,'(a,i4,a,i4)') 'P_MARKER(',IC1,')=', SLF%P_TAG(IC1)-1
         ENDDO
         WRITE(SCARC_LU,*) '================ A_ROW'
         DO IC1 = 1, 12
            WRITE(SCARC_LU,'(a,i4,a,i4)') 'A_ROW(',IC1,')=', SLF%A_ROW(IC1)
         ENDDO
         WRITE(SCARC_LU,*) '================ A_COL'
         DO IC1 = 1, 12
            WRITE(SCARC_LU,'(a,i4,a,i4)') 'A_COL(',IC1,')=', SLF%A_COL(IC1)
         ENDDO
         WRITE(SCARC_LU,*) '================ R_ROW'
         DO IC1 = 1, 12
            WRITE(SCARC_LU,'(a,i4,a,i4)') 'R_ROW(',IC1,')=', SLF%R_ROW(IC1)
         ENDDO
         WRITE(SCARC_LU,*) '================ R_COL'
         DO IC1 = 1, 12
            WRITE(SCARC_LU,'(a,i4,a,i4)') 'R_COL(',IC1,')=', SLF%R_COL(IC1)
         ENDDO
      ENDIF

      IROW_INIT = IROW

      IF (BSQUARE) THEN
         SLF%P_TAG(IC0) = IROW
         IROW = IROW + 1
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SQUARE P_TAG(',IC0,') =', SLF%P_TAG(IC0),'   : IROW=',IROW
      ENDIF

      !< loop over all entries in row IC0 of R
      LOOP1_R_ENTRIES: DO ICOL1 = SLF%R_ROW(IC0), SLF%R_ROW(IC0+1)-1

         IC1 = SLF%R_COL(ICOL1) 
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '--------- ICOL1 =', ICOL1,'   : IC1=',IC1
         
         !< loop over all entries in row IC1 of A
         LOOP1_A_ENTRIES: DO ICOL2 = SLF%A_ROW(IC1), SLF%A_ROW(IC1+1)-1

            IC2 = SLF%A_COL(ICOL2)
            IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '   ------ ICOL2 =', ICOL2,'   : IC2=',IC2
   
            IF (SLF%A_TAG(IC2) /= IC0) THEN
   
               !< mark IC2 as already considered
               SLF%A_TAG(IC2) = IC0
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '          A_TAG(',IC2,') =', SLF%A_TAG(IC2)
     
               !< loop over all entries in row IC2 of P
               LOOP1_P_ENTRIES: DO ICOL3 = SLF%P_ROW(IC2), SLF%P_ROW(IC2+1)-1
   
                  IC3 = SLF%P_COL(ICOL3)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '    !> ---- ICOL3 =', ICOL3,'   : IC3=',IC3
                  !IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
                !<   '          P_TAG(',IC3,') =', SLF%P_TAG(IC3),': IROW_INIT=',IROW_INIT,' : IROW=',IROW
   
                !< verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
                  IF (IC3 <= 0) WRITE (SCARC_LU,*) '------------!!!> CAUTION1: IC3 <=0 :' , IC3
                  IF (IC3 >0) THEN   !< ONLY TEMPORARILY
                  IF (SLF%P_TAG(IC3) < IROW_INIT) THEN
                     SLF%P_TAG(IC3) = IROW
                     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '          P_TAG(',IC3,') =', SLF%P_TAG(IC3)
                     IROW = IROW + 1
                     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '          IROW =', IROW
                  ENDIF
                  ENDIF   !< ONLY TEMPORARILY
               ENDDO LOOP1_P_ENTRIES
   
            ENDIF

         ENDDO LOOP1_A_ENTRIES
      ENDDO LOOP1_R_ENTRIES

      !< Store counters
      IROW_SAVE = IROW

   ENDDO LOOP1_OWN_CELLS

   !! Determine size of matrix RAP  == SLC%A 
   SLC%NA = IROW
   SLC%NAE = SLC%NA + 100  !< ONLY TEMPORARILY

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'AllOCATION RAP in length ', SLC%NA

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

      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,*) '================ MARKMARK2:  IC0 =', IC0-1
         WRITE(SCARC_LU,*) 'IROW     = ', IROW
         WRITE(SCARC_LU,*) 'IROW_INIT= ', IROW_INIT
         WRITE(SCARC_LU,*) '================ A_MARKER'
         DO IC1 = 1, SLF%NCE
            WRITE(SCARC_LU,'(a,i4,a,i4)') 'A_MARKER(',IC1,')=', SLF%A_TAG(IC1)-1
         ENDDO
         WRITE(SCARC_LU,*) '================ P_MARKER'
         DO IC1 = 1, SLF%NCCE
            WRITE(SCARC_LU,'(a,i4,a,i4)') 'P_MARKER(',IC1,')=', SLF%P_TAG(IC1)-1
         ENDDO
      ENDIF

      IROW_INIT = IROW

      SLC%A_ROW(IC0) = IROW_INIT

      IF (BSQUARE) THEN
         SLF%P_TAG(IC0) = IROW
         SLC%A_COL(IROW) = IC0
         SLC%A(IROW) = 0.0_EB
         IROW = IROW + 1
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      WRITE(SCARC_LU,*) 'SQUARE P_TAG(',IC0,')         =', SLF%P_TAG(IC0)
      WRITE(SCARC_LU,*) 'SQUARE RAP_diag_data(',IC0,') =', SLC%A(IROW)
      WRITE(SCARC_LU,*) 'SQUARE RAP_diag_j(',IC0,')    =', SLC%A_COL(IROW)
      WRITE(SCARC_LU,*) 'SQUARE IROW=',IROW
      WRITE(SCARC_LU,*) 'NC =',SLF%NC
      WRITE(SCARC_LU,*) 'NCC=',SLF%NCC
   ENDIF
      ENDIF

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      WRITE(SCARC_LU,*) 'MYDEBUG: 0: ======== RAP 10: A_diag_i  ======================'
      DO IC1 = 1, SLF%NC+1
         !WRITE(SCARC_LU,*) 'A_ROW(',IC1,')=',SLF%A_ROW(IC1)
         WRITE(SCARC_LU,'(a,i5)') 'MYDEBUG: 0: ', SLF%A_ROW(IC1)
      ENDDO
      WRITE(SCARC_LU,*) 'MYDEBUG: 0: ======== RAP 10: A_diag_i  ======================'
      DO IC1 = 1, SLF%NA
         !WRITE(SCARC_LU,*) 'A_COL(',IC1,')=',SLF%A_COL(IC1)
         WRITE(SCARC_LU,'(a,i5)') 'MYDEBUG: 0: ', SLF%A_COL(IC1)
      ENDDO
      WRITE(SCARC_LU,*) 'MYDEBUG: 0: ======== RAP 10: A_diag_i  ======================'
      DO IC1 = 1, SLF%NA
         !WRITE(SCARC_LU,*) 'A(',IC1,')=',SLF%A(IC1)
         WRITE(SCARC_LU,'(a,f12.6)') 'MYDEBUG: 0: ', SLF%A(IC1)
      ENDDO
      WRITE(SCARC_LU,*) '==============================================='
      WRITE(SCARC_LU,*) 'MATRIX R_ROW:'
      DO IC1 = 1, SLF%NCC+1
         WRITE(SCARC_LU,*) 'R_ROW(',IC1,')=',SLF%R_ROW(IC1)
      ENDDO
      WRITE(SCARC_LU,*) 'MATRIX R_COL:'
      DO IC1 = 1, SLF%NR
         WRITE(SCARC_LU,*) 'R_COL(',IC1,')=',SLF%R_COL(IC1)
      ENDDO
      WRITE(SCARC_LU,*) 'MATRIX R:'
      DO IC1 = 1, SLF%NR
         WRITE(SCARC_LU,*) 'R(',IC1,')=',SLF%R(IC1)
      ENDDO
   ENDIF

      !< loop over all entries in row IC0 of R
      LOOP2_R_ENTRIES: DO ICOL1 = SLF%R_ROW(IC0), SLF%R_ROW(IC0+1)-1

         IC1 = SLF%R_COL(ICOL1) 
         R_VALUE = SLF%R(ICOL1)

         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6)') &
                     '--------- ICOL1 =', ICOL1,'   : IC1 =',IC1,': R_VALUE=',R_VALUE
         
         !< loop over all entries in row IC1 of A
         LOOP2_A_ENTRIES: DO ICOL2 = SLF%A_ROW(IC1), SLF%A_ROW(IC1+1)-1

            IC2 = SLF%A_COL(ICOL2)
            RA_VALUE = R_VALUE * SLF%A(ICOL2)

         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6,a,i3,a,i3)') &
                     '--------- ICOL2 =', ICOL2,'   : IC2 =',IC2,': RA_VALUE=',RA_VALUE, ': A_TAG(',IC2,')=',SLF%A_TAG(IC2)
   
            !< Hasn't cell IC2 been considered before? (new values for RAP can only be done for unmarked cells)
            IF (SLF%A_TAG(IC2) /= IC0) THEN
   
               !< mark IC2 as already considered
               SLF%A_TAG(IC2) = IC0
     
               !< loop over all entries in row IC2 of P
               LOOP2_P_ENTRIES1: DO ICOL3 = SLF%P_ROW(IC2), SLF%P_ROW(IC2+1)-1
   
                  IC3 = SLF%P_COL(ICOL3)
                  RAP_VALUE = RA_VALUE * SLF%P(ICOL3)

         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6)') &
                     '---A----- ICOL3 =', ICOL3,'   : IC3 =',IC3,': RAP_VALUE=',RAP_VALUE
   
                    !< verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
                  IF (IC3 <= 0) WRITE (SCARC_LU,*) '------------!!!> CAUTION2: IC3 <=0 : ' , IC3
                  IF (IC3 >0) THEN     !< ONLY TEMPORARILY

                  IF (SLF%P_TAG(IC3) < IROW_INIT) THEN
                     SLF%P_TAG(IC3)  = IROW
                     SLC%A_COL(IROW) = SLF%P_COL(ICOL3)
                     SLC%A(IROW)     = RAP_VALUE
                     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '   A1     P_TAG(',IC3,')     =', SLF%P_TAG(IC3)
                     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '   A1     SLC%A_COL(',IROW,') =', SLC%A_COL(IROW)
                     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '   A1     SLC%A(',IROW,')     =', SLC%A(IROW)
                     IROW = IROW + 1
                     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '   A1     IROW =', IROW
                  ELSE
                     SLC%A(SLF%P_TAG(IC3)) = SLC%A(SLF%P_TAG(IC3)) + RAP_VALUE
                     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '   A2     SLC%A(',SLF%P_TAG(IC3),&
                                                         ')     =', SLC%A(SLF%P_TAG(IC3))
                  ENDIF
                  ENDIF  !< ONLY TEMPORARILY
               ENDDO LOOP2_P_ENTRIES1
   
            !< or has it been already considered
            ELSE

               !< loop over all entries in row IC2 of P
               LOOP2_P_ENTRIES2: DO ICOL3 = SLF%P_ROW(IC2), SLF%P_ROW(IC2+1)-1
   
                  IC3 = SLF%P_COL(ICOL3)
                  RAP_VALUE = RA_VALUE * SLF%P(ICOL3)
                  SLC%A(SLF%P_TAG(IC3)) = SLC%A(SLF%P_TAG(IC3)) + RAP_VALUE

         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6,a,i3,a,f12.6)') &
                     '---B----- ICOL3 =', ICOL3,'   : IC3 =',IC3,': RAP_VALUE=',RAP_VALUE,&
                     ': A(',SLF%P_TAG(IC3),')=',SLC%A(SLF%P_TAG(IC3))

                  IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) '    !> ---- ICOL3 =', ICOL3,'   : IC3=',IC3
                  IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
                     '          P_TAG(',IC3,') =', SLF%P_TAG(IC3),': IROW_INIT=',IROW_INIT,' : IROW=',IROW
   
               ENDDO LOOP2_P_ENTRIES2
            ENDIF

         ENDDO LOOP2_A_ENTRIES
      ENDDO LOOP2_R_ENTRIES

      !< Store counters
      IROW_SAVE = IROW

   ENDDO LOOP2_C_CELLS
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
      WRITE(SCARC_LU,*) 'MYDEBUG: 0: ======== RAP 11: RAP_diag_i  ================= NCC:', SLF%NCC
      DO IC1 = 1, SLF%NCC+1
         !WRITE(SCARC_LU,*) 'A_ROW(',IC1,')=',SLF%A_ROW(IC1)
         WRITE(SCARC_LU,'(a,i5)') 'MYDEBUG: 0: ', SLC%A_ROW(IC1)-1
      ENDDO
      WRITE(SCARC_LU,*) 'MYDEBUG: 0: ======== RAP 11: RAP_diag_i  ================= NA :', SLF%NA
      DO IC1 = 1, SLC%NA
         !WRITE(SCARC_LU,*) 'A_COL(',IC1,')=',SLF%A_COL(IC1)
         WRITE(SCARC_LU,'(a,i5)') 'MYDEBUG: 0: ', SLC%A_COL(IC1)-1
      ENDDO
      WRITE(SCARC_LU,*) 'MYDEBUG: 0: ======== RAP 11: RAP_diag_i  ================= NA :=', SLF%NA
      DO IC1 = 1, SLC%NA
         !WRITE(SCARC_LU,*) 'A(',IC1,')=',SLF%A(IC1)
         WRITE(SCARC_LU,'(a,f12.6)') 'MYDEBUG: 0: ', -SLC%A(IC1)
      ENDDO
   ENDIF


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

      MESHES_LOOP_GMG: DO NM = NMESHES_MIN, NMESHES_MAX
         SL => SCARC(NM)%LEVEL(NL)
         IP = 1
         ICC = 1

         !< Analyze grid sizes
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

      MESHES_LOOP_GMG3: DO NM = NMESHES_MIN, NMESHES_MAX
         SL => SCARC(NM)%LEVEL(NL)
         IP = 1
         ICC = 1

         !< Analyze grid sizes
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
            WRITE(*,*) 'CELL NUNMBERS NOT DIVISABLE BY 2 or 3, STOP!'
            STOP
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

      MESHES_LOOP_CLASSICAL2: DO NM = NMESHES_MIN, NMESHES_MAX
         SL => SCARC(NM)%LEVEL(NL)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'CALLING MATRIX_TRANSPOSE HERE'
         CALL SCARC_MATRIX_TRANSPOSE(SL%P, SL%P_ROW, SL%P_COL, SL%R, SL%R_ROW, SL%R_COL, &
                                      SL%NC, SL%NCE, SL%NCC, SL%NCCE)
      ENDDO MESHES_LOOP_CLASSICAL2


 !! ------------------------------------------------------------------------------------------
 !! DEFAULT
 !! ------------------------------------------------------------------------------------------
   CASE DEFAULT
  
      MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
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

      MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
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
               
!WRITE(SCARC_LU,*) '============= NM=',NM,': NOM=',NOM,': IG0=',IG0
                  OTHER_ROW_LOOP: DO ICP2 = OSL%P_ROW(IG),OSL%P_ROW(IG+1)-1
!WRITE(SCARC_LU,*) 'OSL%P_COL(',ICP2,')=',OSL%P_COL(ICP2)
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

      MESHES_GMG: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SLF => SCARC(NM)%LEVEL(NL)             !< pointer to fine level
         SLC => SCARC(NM)%LEVEL(NL+1)           !< pointer to coarse level
         
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
!WRITE(SCARC_LU,*) 'MATRIX(',IC,',',ICP2,')=',MATRIX(IC,ICP2)
                  
                  AUX2 = AUX2 + R_WEIGHT(IC, ICP1, NM, NL) * AUX1
      
               ENDDO
               VAL(ICP2) = AUX2
      
            ENDDO ICP2_GMG_LOOP1
         
!WRITE(SCARC_LU,'(4f12.6)') ((MATRIX(IC0,ICP0),ICP0=1,4),IC0=1,16)

            !< analyze new matrix line and store it corresponding to compact storage technique:
            !< (diagonal entry first)
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

!WRITE(SCARC_LU,*) 'ENDING MATRIX'

 !! ---------------------------------------------------------------------------------
 !! Default case
 !! ---------------------------------------------------------------------------------
   CASE DEFAULT

      MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SLF => SCARC(NM)%LEVEL(NL)             !< pointer to fine level
         SLC => SCARC(NM)%LEVEL(NL+1)           !< pointer to coarse level
         
         IF (NMESHES==1) THEN
            ALLOCATE (VAL(SLF%NC), STAT=IERR)
         ELSE
            ALLOCATE (VAL(SLF%NCE), STAT=IERR)
         ENDIF
         CALL CHKMEMERR ('SCARC_TRANSFER_MATRIX', 'VAL', IERR)
         VAL = 0.0_EB
         
         IP = 1
         ICP1_LOOP: DO ICP1 = 1, SLF%NCC
      
      
            !< ---------------------------------------------------------------------------------
            !< First: loop over the internal coarse cells
            !< ---------------------------------------------------------------------------------
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
         
            !< analyze new matrix line and store it corresponding to compact storage technique:
            !< (diagonal entry first)
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
REAL (EB) :: TNOW_SOLVER

TNOW_SOLVER = SECOND()

!WRITE(*,*) 'STARTING SCARC_SOLVER'

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

 !! Krylov method (CG/BICG)
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV: SELECT CASE (TYPE_KRYLOV)
         CASE (NSCARC_KRYLOV_CG)                           
!WRITE(*,*) 'STARTING SCARC_METHOD_CG'
            CALL SCARC_METHOD_CG  (NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F, TYPE_PRECON)
         CASE (NSCARC_KRYLOV_BICG)                         
!WRITE(*,*) 'STARTING SCARC_METHOD_BICG'
            CALL SCARC_METHOD_BICG(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F, TYPE_PRECON)
      END SELECT SELECT_KRYLOV

 !! Multigrid method
   CASE (NSCARC_METHOD_MULTIGRID)

!WRITE(*,*) 'STARTING SCARC_METHOD_MULTIGRID'
      CALL SCARC_METHOD_MULTIGRID(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F, TYPE_SMOOTH)

END SELECT SELECT_METHOD

!WRITE(*,*) 'LEAVING SCARC_SOLVER'

TUSED_SCARC(NSCARC_TIME_SOLVER  ,:)=TUSED_SCARC(NSCARC_TIME_SOLVER  ,:)+SECOND()-TNOW_SOLVER
TUSED_SCARC(NSCARC_TIME_TOTAL,:)=TUSED_SCARC(NSCARC_TIME_TOTAL,:)+SECOND()-TNOW_SOLVER
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
FUNCTION POINT_TO_VECTOR(NVECTOR, NM, NL)
REAL(EB), POINTER, DIMENSION(:) :: POINT_TO_VECTOR
INTEGER, INTENT(IN):: NVECTOR, NM, NL

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_F)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%F
   CASE (NSCARC_VECTOR_X)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%X
   CASE (NSCARC_VECTOR_D)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%D
   CASE (NSCARC_VECTOR_Y)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%Y
   CASE (NSCARC_VECTOR_G)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%G
   CASE (NSCARC_VECTOR_W)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%W
   CASE (NSCARC_VECTOR_Z)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%Z
   CASE (NSCARC_VECTOR_X2)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%X2
   CASE (NSCARC_VECTOR_D2)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%D2
   CASE (NSCARC_VECTOR_Y2)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%Y2
   CASE (NSCARC_VECTOR_G2)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%G2
   CASE (NSCARC_VECTOR_W2)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%W2
   CASE (NSCARC_VECTOR_Z2)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%Z2
   CASE (NSCARC_VECTOR_MEASURE)
      POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%MEASURE
END SELECT

RETURN
END FUNCTION POINT_TO_VECTOR


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

DO NM = NMESHES_MIN, NMESHES_MAX
   SL => SCARC(NM)%LEVEL(NL)
   V  => POINT_TO_VECTOR (NVECTOR, NM, NL)
   DO IC = 1, SL%NC
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
REAL(EB) :: VCO
INTEGER , POINTER, DIMENSION(:)     :: AC_ROW, AC_COL
INTEGER , POINTER :: NC
INTEGER :: NM, IC, JC, ICOL

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
DO NM = NMESHES_MIN, NMESHES_MAX

   V1 => POINT_TO_VECTOR (NVECTOR1, NM, NL)
   V2 => POINT_TO_VECTOR (NVECTOR2, NM, NL)
   
   NC => SCARC(NM)%LEVEL(NL)%NC                                       !< number of cells

   AC     => SCARC(NM)%LEVEL(NL)%A                                    !< system matrix
   AC_ROW => SCARC(NM)%LEVEL(NL)%A_ROW                                !< row pointer
   AC_COL => SCARC(NM)%LEVEL(NL)%A_COL                                !< column pointer

!WRITE(SCARC_LU,'(a, 14f10.5)') 'V1:', (V1(IC),IC=1,12)
!WRITE(SCARC_LU,'(a, 14f10.5)') 'V2:', (V2(IC),IC=1,12)

   DO IC = 1, NC
           
      !< diagonal entry
      ICOL = AC_ROW(IC)                             
      JC   = AC_COL(ICOL)

      V2 (IC) = AC(ICOL)* V1(JC)
      VCO = V2(IC)
!WRITE(SCARC_LU,'(a,i4,a,5f14.6,i4)') 'MATVEC 1: V2(',IC,')=',V2(IC), 0.0_EB, V1(JC), AC(ICOL), AC(ICOL)*V1(JC), JC

      !< subdiagonal entries
      DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1          
         JC = AC_COL(ICOL)
         V2(IC) =  V2(IC) + AC(ICOL)* V1(JC)
!WRITE(SCARC_LU,'(a,i4,a,5f14.6,i4)') 'MATVEC 2: V2(',IC,')=',V2(IC), VCO, V1(JC), AC(ICOL), AC(ICOL)*V1(JC), JC
         VCO = V2(IC)
      ENDDO

   ENDDO

ENDDO 

CALL SCARC_DEBUG_LEVEL (NVECTOR1, 'SCARC_MATVEC_PRODUCT', 'V1 FINAL', NL)
CALL SCARC_DEBUG_LEVEL (NVECTOR2, 'SCARC_MATVEC_PRODUCT', 'V2 FINAL', NL)

END SUBROUTINE SCARC_MATVEC_PRODUCT


!> ------------------------------------------------------------------------------------------------
!> Compute global scalarproductt (including global data exchange) 
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_SCALAR_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER, INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), DIMENSION(:)    , POINTER ::  V1, V2
INTEGER , POINTER :: NC
INTEGER  :: NM, IERR, NL0, IC
REAL(EB) :: SP_GLOBAL
!>MKL code
!REAL(EB) :: DDOT
!EXTERNAL :: DDOT

DO NM = NMESHES_MIN, NMESHES_MAX

   V1 => POINT_TO_VECTOR (NVECTOR1, NM, NL)
   V2 => POINT_TO_VECTOR (NVECTOR2, NM, NL)
   
   NC  => SCARC(NM)%LEVEL(NL)%NC

      
!WRITE(SCARC_LU,'(a, 4f25.15)') 'SCALAR: V1:', (V1(IC),IC=1,NC)
!WRITE(SCARC_LU,'(a, 4f25.15)') 'SCALAR: V2:', (V2(IC),IC=1,NC)
   !IF (SCARC_MKL) THEN
 !!   SP_LOCAL(NM) = DDOT(NC, V1, 1, V2, 1)
   !ELSE
      SP_LOCAL(NM) = 0.0_EB
      DO IC = 1, NC
         SP_LOCAL(NM) = SP_LOCAL(NM) + V1(IC) * V2(IC)
      ENDDO
   !ENDIF

ENDDO


!> ------------------------------------------------------------------------------------------------
!> Get global scalarproduct by a global summation of the local values
!> ------------------------------------------------------------------------------------------------
IERR = 0
NL0  = NL
SP_GLOBAL   = 0.0_EB
IF (NMESHES>1 .AND. USE_MPI) THEN
   CALL MPI_ALLREDUCE (SP_LOCAL(MYID+1), SP_GLOBAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   DO NM=1,NMESHES
      SP_GLOBAL = SP_GLOBAL + SP_LOCAL(NM)
   ENDDO
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
INTEGER  :: NM, IERR, IC
REAL(EB) :: SP_GLOBAL  
REAL(EB) :: DDOT
EXTERNAL :: DDOT

DO NM = NMESHES_MIN, NMESHES_MAX

   V  => POINT_TO_VECTOR (NVECTOR1, NM, NL)
   NC => SCARC(NM)%LEVEL(NL)%NC

   !IF (SCARC_MKL) THEN
   !   SP_LOCAL(NM) = DDOT(NC, VC, 1, VC, 1)
   !ELSE

      SP_LOCAL(NM) = 0.0_EB
      DO IC = 1, NC
         SP_LOCAL(NM) = SP_LOCAL(NM) + V(IC) * V(IC)
      ENDDO

   !ENDIF

ENDDO


!> ------------------------------------------------------------------------------------------------
!> Get global scalarproduct by a global summation of the local values
!> ------------------------------------------------------------------------------------------------
IERR = 0
SP_GLOBAL = 0.0_EB
IF (NMESHES>1 .AND. USE_MPI) THEN
   CALL MPI_ALLREDUCE (SP_LOCAL(MYID+1), SP_GLOBAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   DO NM=1,NMESHES
      SP_GLOBAL = SP_GLOBAL + SP_LOCAL(NM)
   ENDDO
ENDIF
SP_GLOBAL = SQRT (SP_GLOBAL/REAL(NC_GLOBAL(NL), EB))   

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
!>MKL code
!INTEGER , POINTER :: NC
!EXTERNAL :: DAXPBY

DO NM = NMESHES_MIN, NMESHES_MAX

   V1 => POINT_TO_VECTOR(NVECTOR1, NM, NL)
   V2 => POINT_TO_VECTOR(NVECTOR2, NM, NL)

   !IF (SCARC_MKL) THEN
 !!   NC  => SCARC(NM)%LEVEL(NL)%NC
 !!   CALL DAXPBY(NC, SCAL1, V1, 1, SCAL2, V2, 1)
   !ELSE
     V2 = SCAL1 * V1 + SCAL2 * V2
   !ENDIF

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
!> MKL code
!INTEGER , POINTER :: NC
!EXTERNAL :: DCOPY, DSCAL

DO NM = NMESHES_MIN, NMESHES_MAX

   V1 => POINT_TO_VECTOR(NVECTOR1, NM, NL)
   V2 => POINT_TO_VECTOR(NVECTOR2, NM, NL)

   !IF (SCARC_MKL) THEN
 !!   NC  => SCARC(NM)%LEVEL(NL)%NC
 !!   CALL DCOPY(NC, V1, 1, V2, 1)
 !!   CALL DSCAL(NC, SCAL1, V2, 1)
   !ELSE
      V2 = SCAL1 * V1
   !ENDIF

ENDDO

END SUBROUTINE SCARC_VECTOR_COPY


!> ------------------------------------------------------------------------------------------------
!> Clear vector
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_CLEAR(NVECTOR, NL)
INTEGER , INTENT(IN):: NVECTOR, NL
REAL(EB), DIMENSION(:)    , POINTER ::  VC
INTEGER  :: NM

DO NM = NMESHES_MIN, NMESHES_MAX
   VC => POINT_TO_VECTOR(NVECTOR, NM, NL)
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

DO NM = NMESHES_MIN, NMESHES_MAX
   VC => POINT_TO_VECTOR(NVECTOR, NM, NL)
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
INTEGER  :: NPRECON, NM, I, J, K, IC, ICOL
REAL(EB) :: AUX, OMEGA=1.5_EB
TYPE (MESH_TYPE), POINTER :: M

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
 !! Multigrid preconditioner
 !! ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID (NSCARC_SCOPE_PRECON, NVECTOR2, TYPE_SMOOTH)


 !! ----------------------------------------------------------------------------------------
 !! Jacobi preconditioner
 !! ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_JACOBI)

      DO NM = NMESHES_MIN, NMESHES_MAX

         V2 => POINT_TO_VECTOR(NVECTOR2, NM, NL)
   
         NC     => SCARC(NM)%LEVEL(NL)%NC
         AC     => SCARC(NM)%LEVEL(NL)%A
         AC_ROW => SCARC(NM)%LEVEL(NL)%A_ROW

         DO IC = 1, NC
            V2 (IC) = V2 (IC) / AC (AC_ROW(IC))
!WRITE(SCARC_LU,'(a,i4,a,2f16.6,i4)') 'JACOBI: V2(',IC,')=',V2(IC),AC(AC_ROW(IC)), AC_ROW(IC)
         ENDDO

      ENDDO
   
 !! ----------------------------------------------------------------------------------------
 !! SSOR preconditioner
 !! ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_SSOR)
   
      DO NM = NMESHES_MIN, NMESHES_MAX
      
         V2 => POINT_TO_VECTOR(NVECTOR2, NM, NL)
   
         NC     => SCARC(NM)%LEVEL(NL)%NC
         AC     => SCARC(NM)%LEVEL(NL)%A
         AC_ROW => SCARC(NM)%LEVEL(NL)%A_ROW
         AC_COL => SCARC(NM)%LEVEL(NL)%A_COL

         !< use only matrix superdiagonals
         FORWARD_CELL_LOOP: DO IC = 1, NC

            AUX = 0.0_EB
            LOWER_DIAG_LOOP: DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1
               IF (AC_COL(ICOL) >= IC) CYCLE LOWER_DIAG_LOOP
               IF (AC_COL(ICOL) <= NC) AUX = AUX + AC(ICOL) * V2(AC_COL(ICOL))
            ENDDO LOWER_DIAG_LOOP
            V2(IC) = (V2(IC) - AUX * OMEGA) / AC(AC_ROW(IC))
!WRITE(SCARC_LU,'(a,i4,a,2f16.6,i4)') 'SSOR1: V2(',IC,')=',V2(IC),AC(AC_ROW(IC)), AC_ROW(IC)
         
         ENDDO FORWARD_CELL_LOOP
         
         !< use only matrix subdiagonals
         BACKWARD_CELL_LOOP: DO IC = NC-1, 1, -1
            
            AUX = 0.0_EB
            UPPER_DIAG_LOOP: DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1
               IF (AC_COL(ICOL) <= IC) CYCLE UPPER_DIAG_LOOP
               IF (AC_COL(ICOL) <= NC) AUX = AUX + AC(ICOL) * V2(AC_COL(ICOL))
            ENDDO UPPER_DIAG_LOOP
            V2(IC) = V2(IC) - AUX * OMEGA / AC(AC_ROW(IC))
!WRITE(SCARC_LU,'(a,i4,a,2f16.6,i4)') 'SSOR2: V2(',IC,')=',V2(IC),AC(AC_ROW(IC)), AC_ROW(IC)

         ENDDO BACKWARD_CELL_LOOP
   
      ENDDO
      
 !! ----------------------------------------------------------------------------------------
 !! FFT preconditioner
 !! ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_FFT)
      
      DO NM = NMESHES_MIN, NMESHES_MAX
      
         V1 => POINT_TO_VECTOR(NVECTOR1, NM, NL)
         V2 => POINT_TO_VECTOR(NVECTOR2, NM, NL)
   
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


END SUBROUTINE SCARC_PRECONDITIONING 

!> ------------------------------------------------------------------------------------------------
!> Save and reset settings of CALLing parent-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SAVE_PARENT(PARENT)
TYPE (SCARC_PARENT_TYPE), INTENT(OUT):: PARENT
PARENT%TYPE_METHOD   = TYPE_METHOD
PARENT%TYPE_SCOPE    = TYPE_SCOPE 
PARENT%TYPE_PRECON   = TYPE_PRECON
PARENT%TYPE_SMOOTH   = TYPE_SMOOTH
PARENT%TYPE_CYCLE    = TYPE_CYCLE 
PARENT%TYPE_ACCURACY = TYPE_ACCURACY
END SUBROUTINE SCARC_SAVE_PARENT


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
END SUBROUTINE SCARC_RESET_PARENT


!> ------------------------------------------------------------------------------------------------
!> Perform global CG-method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CG(NSCOPE, NVECTOR, NPRECON)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR, NPRECON
INTEGER   :: NL
INTEGER   :: ISTATE, ITE
REAL (EB) :: SIGMA0, SIGMA1, ALPHA0, GAMMA0
REAL (EB) :: TNOW_KRYLOV
TYPE (SCARC_SCOPE_TYPE) :: CG
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW_KRYLOV = SECOND()


!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Set environment variables and define working level
!>   - Initialize solution, right hand side vector and auxiliary vectors
!>   - Define iterations parameters
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(CG, NSCOPE, NPRECON, NVECTOR, NL)
CALL SCARC_SETUP_SOLVER(NL)

!CALL SCARC_SETUP_VALUE( 1.0_EB, CG%X, NL)
!CALL SCARC_SETUP_VALUE(10.0_EB, CG%F, NL)

CALL SCARC_DEBUG_LEVEL (CG%X, 'SCARC_METHOD_CG', 'X INIT', NL)
CALL SCARC_DEBUG_LEVEL (CG%F, 'SCARC_METHOD_CG', 'F INIT', NL)

!> ------------------------------------------------------------------------------------------------
!> Compute initial residual and perform initial preconditioning
!> ------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (CG%X, CG%W, NL)                                !<  W := A*X
CALL SCARC_VECTOR_SUM     (CG%F, CG%W, -1.0_EB, 1.0_EB, NL)               !<  W := W - F

CALL SCARC_DEBUG_LEVEL (CG%W, 'SCARC_METHOD_CG', 'W INIT', NL)

CG%RESIN = SCARC_L2NORM (CG%W, NL)                                        !<  RESIN := ||W||

CALL SCARC_CONVERGENCE_INFO(CG%RESIN, 0, NL, CG%CROUTINE)
CALL SCARC_VECTOR_COPY     (CG%W, CG%G, 1.0_EB, NL)                       !<  G := W
CALL SCARC_PRECONDITIONING (CG%W, CG%G, NL, NSCOPE)                       !<  G := PRECON(W)

SIGMA0 = SCARC_SCALAR_PRODUCT(CG%W, CG%G, NL)                             !<  SIGMA0 := (W,G)

CALL SCARC_DEBUG_LEVEL (CG%G, 'SCARC_METHOD_CG', 'G PRECON INIT', NL)

CALL SCARC_VECTOR_COPY (CG%G, CG%D, -1.0_EB, NL)                          !<  D := -G

!> ------------------------------------------------------------------------------------------------
!> Perform conjugate gradient looping
!> ------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, CG%NIT
 
   CALL SCARC_MATVEC_PRODUCT (CG%D, CG%Y, NL)                             !<  Y := A*D

CALL SCARC_DEBUG_LEVEL (CG%Y, 'SCARC_METHOD_CG', 'Y LOOP', NL)

   ALPHA0 = SCARC_SCALAR_PRODUCT (CG%D, CG%Y, NL)                         !<  ALPHA0 := (D,Y)
   ALPHA0 = SIGMA0/ALPHA0                                                 

!WRITE(SCARC_LU,*) 'ALPHA0=',ALPHA0,': SIGMA0=',SIGMA0

   CALL SCARC_VECTOR_SUM (CG%D, CG%X, ALPHA0, 1.0_EB, NL)                 !<  X := ALPHA0*D + X
   CALL SCARC_VECTOR_SUM (CG%Y, CG%W, ALPHA0, 1.0_EB, NL)                 !<  W := ALPHA0*Y + W

CALL SCARC_DEBUG_LEVEL (CG%D, 'SCARC_METHOD_CG', 'D LOOP', NL)
CALL SCARC_DEBUG_LEVEL (CG%X, 'SCARC_METHOD_CG', 'X LOOP', NL)
CALL SCARC_DEBUG_LEVEL (CG%W, 'SCARC_METHOD_CG', 'W LOOP', NL)

   CG%RES = SCARC_L2NORM (CG%W, NL)                                       !<  RES := ||W||
   ISTATE = SCARC_CONVERGENCE_STATE (CG, ITE, NL)                         !<  RES < TOL ??
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP
IF (TYPE_DEBUG >NSCARC_DEBUG_LESS.AND.MYID==0) &
   WRITE(0,'(a,i3,a,e14.5)') ' CG  -Iteration  =',ITE,': Residuum=',SCARC_RESIDUAL

   CALL SCARC_VECTOR_COPY     (CG%W, CG%G, 1.0_EB, NL)                    !<  G := W
   CALL SCARC_PRECONDITIONING (CG%W, CG%G, NL, NSCOPE)                    !<  G := PRECON(W)

CALL SCARC_DEBUG_LEVEL (CG%G, 'SCARC_METHOD_CG', 'G PRECON ', NL)
CALL SCARC_DEBUG_LEVEL (CG%W, 'SCARC_METHOD_CG', 'W PRECON ', NL)

   SIGMA1 = SCARC_SCALAR_PRODUCT (CG%W, CG%G, NL)                         !<  SIGMA1 := (W,G)
   GAMMA0 = SIGMA1/SIGMA0                                                   
   SIGMA0 = SIGMA1                                                         

   CALL SCARC_VECTOR_SUM (CG%G, CG%D, -1.0_EB, GAMMA0, NL)                !<  D := -G + GAMMA0*D

ENDDO CG_LOOP
 
!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of CG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector 
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries 
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(CG, ITE, ISTATE)
!IF (TYPE_DEBUG >=NSCARC_DEBUG_NONE.AND.MYID==0) &
!>  WRITE(0,'(a,e14.5)') '                                       !> ---->  Konvergenzrate=',SCARC_CAPPA


IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_FINISH_SOLVER  (NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS (NLEVEL_MIN)
!ELSE IF (TYPE_SCOPE == NSCARC_SCOPE_COARSE) THEN
!>  CALL SCARC_VECTOR_COPY (CG%X, CG%F, 1.0_EB, NL)                         
ENDIF

CALL SCARC_DEBUG_LEVEL (CG%X, 'SCARC_METHOD_CG', 'X FINAL', NL)
CALL SCARC_DEBUG_LEVEL (CG%F, 'SCARC_METHOD_CG', 'F FINAL', NL)

CALL SCARC_RESET_PARENT(PARENT)

!WRITE(*,*) 'FINAL: SCARC(1)%X:'
!WRITE(*,'(8f15.6)') (SCARC(1)%LEVEL(NL)%X(IC), IC=1, SCARC(1)%LEVEL(NL)%NC)

TUSED_SCARC(NSCARC_TIME_KRYLOV,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_METHOD_CG


!> ------------------------------------------------------------------------------------------------
!> Perform global BICGstab-method based on global possion-matrix 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_BICG(NSCOPE, NVECTOR, NPRECON)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR, NPRECON
INTEGER   :: NL = NSCARC_LEVEL_NONE
INTEGER   :: ISTATE, ITE
REAL (EB) :: ALPHA0, ALPHA1, ALPHA2, RHO0, RHO1, DTHETA, DBETA
REAL (EB) :: TNOW_KRYLOV
TYPE (SCARC_SCOPE_TYPE)  :: BICG
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW_KRYLOV = SECOND()

!> ------------------------------------------------------------------------------------------------
!> Initialization
!>   - Set environment variables and define working level
!>   - Initialize solution, right hand side vector and auxiliary vectors
!>   - Define iterations parameters
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(BICG, NSCOPE, NPRECON, NVECTOR, NL)
CALL SCARC_SETUP_SOLVER(NL)
   
!> ------------------------------------------------------------------------------------------------
!> Compute initial defect and perform (double) initial preconditioning
!> ------------------------------------------------------------------------------------------------
ALPHA0 = 1.0_EB
RHO0   = 1.0_EB
RHO1   = 0.0_EB
DBETA  = 0.0_EB
DTHETA = 1.0_EB
 
!CALL SCARC_VECTOR_COPY    (BICG%F, BICG%W, 1.0_EB, NL)                         !<  W := F
!CALL SCARC_PRECONDITIONING (BICG%W, BICG%W, NL, NSCOPE)                        !<  W := PRECON(W)
CALL SCARC_MATVEC_PRODUCT  (BICG%X, BICG%W, NL)                                 !<  W := A*X
CALL SCARC_VECTOR_SUM      (BICG%F, BICG%W, 1.0_EB, -1.0_EB, NL)                !<  W := F - W
CALL SCARC_PRECONDITIONING (BICG%W, BICG%W, NL, NSCOPE)                         !<  W := PRECON(W)

BICG%RESIN = SCARC_L2NORM (BICG%W, NL)                                          !<  RESIN := ||W||
CALL SCARC_CONVERGENCE_INFO (BICG%RESIN, 0, NL, BICG%CROUTINE)

CALL SCARC_VECTOR_COPY (BICG%W, BICG%G, 1.0_EB, NL)                             !<  G := W
   
!> ------------------------------------------------------------------------------------------------
!> Perform bi-conjugate gradient looping:
!> ------------------------------------------------------------------------------------------------
BICG_LOOP: DO ITE = 1, BICG%NIT

   RHO1  = SCARC_SCALAR_PRODUCT (BICG%G, BICG%W, NL)                            !< RHO1 := (G,W)
   DBETA = (RHO1*DTHETA)/(RHO0*ALPHA0)                                     
   RHO0  = RHO1

   CALL SCARC_VECTOR_SUM      (BICG%W, BICG%Z, 1.0_EB       , DBETA , NL)       !< Z := W + DBETA*Z
   CALL SCARC_VECTOR_SUM      (BICG%Y, BICG%Z, -DBETA*ALPHA0, 1.0_EB, NL)       !< Z := -DBETA*ALPHA0*Y + Z
   CALL SCARC_MATVEC_PRODUCT  (BICG%Z, BICG%Y, NL)                              !< Y := A*Z
   CALL SCARC_PRECONDITIONING (BICG%Y, BICG%Y, NL, NSCOPE)                      !< Z := PRECON(Z)

   DTHETA = SCARC_SCALAR_PRODUCT (BICG%G, BICG%Y, NL)                           !< DTHETA := (G,Y)
   DTHETA = RHO1/DTHETA

   CALL SCARC_VECTOR_SUM      (BICG%Y, BICG%W, -DTHETA, 1.0_EB, NL)             !< W := -DTHETA*Y + W
   CALL SCARC_MATVEC_PRODUCT  (BICG%W, BICG%D, NL)                              !< D := A*W
   CALL SCARC_PRECONDITIONING (BICG%D, BICG%D, NL, NSCOPE)                      !< D := PRECON(D)
   
   ALPHA1 = SCARC_SCALAR_PRODUCT (BICG%D, BICG%W, NL)                           !< ALPHA1 := (D,W)
   ALPHA2 = SCARC_SCALAR_PRODUCT (BICG%D, BICG%D, NL)                           !< ALPHA2 := (D,D)
   ALPHA0 = ALPHA1/ALPHA2

   CALL SCARC_VECTOR_SUM (BICG%Z, BICG%X,  DTHETA, 1.0_EB, NL)                  !< X :=  DTHETA*Z + X
   CALL SCARC_VECTOR_SUM (BICG%W, BICG%X,  ALPHA0, 1.0_EB, NL)                  !< X :=  ALPHA0*W + X
   CALL SCARC_VECTOR_SUM (BICG%D, BICG%W, -ALPHA0, 1.0_EB, NL)                  !< W := -ALPHA0*D + W

   BICG%RES = SCARC_L2NORM (BICG%W, NL)                                         !< RES := ||W||

   ISTATE = SCARC_CONVERGENCE_STATE(BICG, ITE, NL)                              !< RES < TOL ???
IF (TYPE_DEBUG >NSCARC_DEBUG_INFO1.AND.MYID==0) &
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

TUSED_SCARC(NSCARC_TIME_KRYLOV,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_METHOD_BICG


!> ------------------------------------------------------------------------------------------------
!> Perform geometric multigrid method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MULTIGRID(NSCOPE, NVECTOR, NPRECON)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR, NPRECON
INTEGER   :: NL = NSCARC_LEVEL_NONE
INTEGER   :: ISTATE, ICYCLE, ITE
REAL (EB) :: TNOW_MULTIGRID
TYPE (SCARC_SCOPE_TYPE)  :: MG
TYPE (SCARC_PARENT_TYPE) :: PARENT

TNOW_MULTIGRID = SECOND()

!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Set environment variables and define working level
!>   - Initialize solution, right hand side vector 
!>   - Define iterations parameters (NL is set to finest level)
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(MG, NSCOPE, NPRECON, NVECTOR, NL)
CALL SCARC_SETUP_SOLVER(NL)

CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MG', 'F init ', NL)
CALL SCARC_DEBUG_LEVEL (MG%X, 'SCARC_METHOD_MG', 'X init ', NL)


!> ------------------------------------------------------------------------------------------------
!> Compute initial defect:  RESIN := || F - A*X ||
!>   - Initialize cycle counts for MG-iteration
!>   - Perform initial matrix-vector product on finest level
!>   - calculate norm of initial residual on finest level
!> ------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (MG%X, MG%D, NL)                                       !<  D := A*X
CALL SCARC_VECTOR_SUM     (MG%F, MG%D, 1.0_EB, -1.0_EB, NL)                      !<  D := F - D 

CALL SCARC_DEBUG_LEVEL (MG%D, 'SCARC_METHOD_MG', 'D init ', NL)

ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_SETUP, NL)
MG%RESIN = SCARC_L2NORM (MG%D, NL)                                               !<  RESIN := ||D||

CALL SCARC_CONVERGENCE_INFO(MG%RESIN, 0, NL, MG%CROUTINE)

MG%NIT=10    !!! ONLY TEMPORARILY

!> ------------------------------------------------------------------------------------------------

!> Perform multigrid-looping (start each iteration on finest level) 
!> ------------------------------------------------------------------------------------------------
MULTIGRID_LOOP: DO ITE = 1, MG%NIT
 
   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLE_EXIT)

      !< presmoothing  (smoothing/restriction till coarsest level is reached)
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)
         CALL SCARC_SMOOTHING    (MG, NSCARC_CYCLE_PRESMOOTH, NL, NSCOPE)        !< D_fine   := smooth(defect)
         CALL SCARC_RESTRICTION  (MG%D, MG%F, NL, NL+1)                          !< F_coarse := rest(D_fine)
CALL SCARC_DEBUG_LEVEL (MG%D, 'SCARC_METHOD_MG', 'D after restrict', NL)
CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MG', 'F after restrict', NL+1)
         CALL SCARC_VECTOR_CLEAR (MG%X, NL+1)                                    !< X_coarse := 0.0
         NL = NL + 1                                                             !< set coarser level
      ENDDO PRESMOOTHING_LOOP

      !< coarse grid solver
      SELECT CASE (TYPE_COARSE)
         CASE (NSCARC_COARSE_ITERATIVE)
            CALL SCARC_METHOD_CG (NSCARC_SCOPE_COARSE, MG%F, NSCARC_PRECON_SSOR) ! X_coarse := exact_sol(.)
         CASE (NSCARC_COARSE_DIRECT)
            CALL SCARC_METHOD_GE (MG%X, MG%F)
      END SELECT

CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MG', 'F after coarse', NL)
CALL SCARC_DEBUG_LEVEL (MG%X, 'SCARC_METHOD_MG', 'X after coarse', NL)

      !< postsmoothing (smoothing/restriction till finest level is reached again)
      POSTSMOOTHING_LOOP: DO WHILE (NL > NLEVEL_MIN) 
         NL=NL-1
         CALL SCARC_PROLONGATION (MG%X, MG%D, NL+1, NL)                          !< D_fine := prol(X_coarse)
CALL SCARC_DEBUG_LEVEL (MG%D, 'SCARC_METHOD_MG', 'D after prol', NL)
         CALL SCARC_VECTOR_SUM   (MG%D, MG%X, 1.0_EB, 1.0_EB, NL)                !< X_fine := D_fine + X_fine
         CALL SCARC_SMOOTHING    (MG, NSCARC_CYCLE_POSTSMOOTH, NL, NSCOPE)       !< D_fine := smooth(defect)
         ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_PROCEED, NL)                  !< perform requested cycle
         IF (ICYCLE /= NSCARC_CYCLE_POSTSMOOTH) CYCLE CYCLE_LOOP
      ENDDO POSTSMOOTHING_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) THEN
      WRITE(*,*) 'ERROR in SCARC_MULTIGRID, wrong level ', NL
      STOP
   ENDIF

 !! ---------------------------------------------------------------------------------------------
 !! Compute norm of new residual on finest level and  leave loop correspondingly
 !! ---------------------------------------------------------------------------------------------
   CALL SCARC_MATVEC_PRODUCT (MG%X, MG%D, NL)                                    !< D := A*X
   CALL SCARC_VECTOR_SUM     (MG%F, MG%D, 1.0_EB, -1.0_EB, NL)                   !< D := F - D

   MG%RES = SCARC_L2NORM (MG%D, NL)                                              !< RES := ||D||
   ISTATE = SCARC_CONVERGENCE_STATE(MG, ITE, NL)                                 !< convergence ?
IF (TYPE_DEBUG >NSCARC_DEBUG_NONE.AND.MYID==0) &
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

!CALL SCARC_DEBUG_LEVEL (MG%X, 'SCARC_METHOD_MULTIGRID', 'X FINAL', NLEVEL_MIN)
!CALL SCARC_DEBUG_LEVEL (MG%F, 'SCARC_METHOD_MULTIGRID', 'F FINAL', NLEVEL_MIN)
CALL SCARC_RESET_PARENT(PARENT)

TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)=TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)+SECOND()-TNOW_MULTIGRID
TUSED_SCARC(NSCARC_TIME_TOTAL    ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL    ,:)+SECOND()-TNOW_MULTIGRID
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
      
 !! ---------------------------------------------------------------------------------------------
 !! reset cycle counts at beginning of each new multigrid iteration
 !! ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_RESET)

      DO NM = NMESHES_MIN, NMESHES_MAX
         DO NL0 = NLEVEL_MIN, NLEVEL_MAX
            SCARC(NM)%CYCLE_COUNT(1,NL0)=SCARC(NM)%CYCLE_COUNT(2,NL0)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLE_NEXT

 !! ---------------------------------------------------------------------------------------------
 !! determine where to proceed with cycling
 !! ---------------------------------------------------------------------------------------------
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


!> ------------------------------------------------------------------------------------------------
!> Perform smoothing 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SMOOTHING(MG, NTYPE, NL, NSCOPE)
INTEGER , INTENT(IN) :: NTYPE, NL, NSCOPE
TYPE (SCARC_SCOPE_TYPE), INTENT(INOUT) :: MG
INTEGER :: ITE, ISTATE
REAL(EB):: TNOW_SMOOTH
LOGICAL :: BMATVEC, BL2NORM
TYPE (SCARC_SCOPE_TYPE) :: SM

TNOW_SMOOTH = SECOND()

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
   CALL SCARC_MATVEC_PRODUCT (SM%X, SM%D, NL)                                  !<  D := A*X
   CALL SCARC_VECTOR_SUM     (SM%F, SM%D, 1.0_EB, -1.0_EB, NL)                 !<  D := F - D
ENDIF

CALL SCARC_DEBUG_LEVEL (SM%D, 'SCARC_SMOOTHER', 'D init ', NL)

IF (BL2NORM.AND.BMATVEC) THEN
   SM%RESIN = SCARC_L2NORM (SM%D, NL)                                          !<  RESIN := ||D||
ELSE
   SM%RESIN = SCARC_RESIDUAL
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Smoothing loop
!> ------------------------------------------------------------------------------------------------
SMOOTH_LOOP: DO ITE=1, SM%NIT
 
   CALL SCARC_PRECONDITIONING (SM%D, SM%D, NL, NSCOPE)                         !<  D := PRECON (D)
CALL SCARC_DEBUG_LEVEL (SM%D, 'SCARC_SMOOTHER', 'D precon ', NL)
   CALL SCARC_VECTOR_SUM      (SM%D, SM%X, SM%OMEGA, 1.0_EB, NL)               !<  X := OMEGA*D + X
   CALL SCARC_MATVEC_PRODUCT  (SM%X, SM%D, NL)                                 !<  D := A*X
CALL SCARC_DEBUG_LEVEL (SM%D, 'SCARC_SMOOTHER', 'D matvec ', NL)
   CALL SCARC_VECTOR_SUM      (SM%F, SM%D, 1.0_EB, -1.0_EB, NL)                !<  D := F - D
CALL SCARC_DEBUG_LEVEL (SM%D, 'SCARC_SMOOTHER', 'D sum ', NL)

   IF (BL2NORM.OR.ITE==SM%NIT) THEN
      SM%RES = SCARC_L2NORM (SM%D, NL)                                         !<  RES := ||D||
!WRITE(SCARC_LU,*) 'SMOOTHER RESIDUUM =',SM%RES
      ISTATE = SCARC_CONVERGENCE_STATE(SM, ITE, NL)
      IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP                     !<  RES < TOL ?
   ENDIF
!WRITE(SCARC_LU,*) 'SMOOTHER LOOP ', ITE,': RESIDUAL =',SM%RES, BL2NORM

ENDDO SMOOTH_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> ------------------------------------------------------------------------------------------------
IF (BL2NORM) CALL SCARC_CONVERGENCE_RATE(SM, ITE, NL)

TUSED_SCARC(NSCARC_TIME_SMOOTH,:)=TUSED_SCARC(NSCARC_TIME_SMOOTH,:)+SECOND()-TNOW_SMOOTH
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_SMOOTH
END SUBROUTINE SCARC_SMOOTHING


!> ------------------------------------------------------------------------------------------------
!> Gaussian elimination for coarse grid solution
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_GE(NVECTORX, NVECTORF)
INTEGER, INTENT(IN) :: NVECTORX, NVECTORF
INTEGER :: NM, IC, IOFFSET, IERR
REAL(EB), POINTER, DIMENSION(:)     :: VCX, VCF
REAL (EB) :: TNOW_COARSE
TYPE (SCARC_TYPE), POINTER :: SM
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL

TNOW_COARSE = SECOND()

IERR=0
SM => SCARC(NMASTER)

!> ------------------------------------------------------------------------------------------------
!> Parallel version
!> ------------------------------------------------------------------------------------------------
IF (USE_MPI) THEN

   DO NM = 1, NMESHES
      SM%COUNTS1(NM-1) = NC_COARSE(NM)
      SM%DISPLS1(NM-1) = SM%OFFSET(NM)
   ENDDO

   DO NM = NMESHES_MIN, NMESHES_MAX
   
      SL => SCARC(NM)%LEVEL(NLEVEL_MAX)
      IOFFSET = SM%OFFSET(NM)

      DO IC = 1, SL%NZ
         SM%X_COARSE (IC) = VCF(IC)
      ENDDO
 
      WRITE(*,*) 'SCARC_METHOD_GE, LEVEL, NOT FINISHED YET'
      STOP
   
   ENDDO
   
   WRITE(*,*) 'STILL MKL VERSION !!'
   !IF (MYID+1 == NMASTER) THEN
 !!   CALL DGETRS('N', NC_COARSE0, 1, SM%A_COARSE, NC_COARSE0, SM%PIVOT, &
 !!               SM%X_COARSE, NC_COARSE0, IERR)
   !ENDIF
   
   DO NM = NMESHES_MIN, NMESHES_MAX
   
      SL => SCARC(NM)%LEVEL(NLEVEL_MAX)
      IOFFSET = SM%OFFSET(NM)
  
      WRITE(*,*) 'ACHTUNG; WAS IST HIER MIT BVECTOR???', NVECTORX, NVECTORF
      !VBX => POINT_TO_HVECTOR (NVECTORX, NM)
 
      DO IC = 1, SL%NC
         VCX(IC) = SM%X_COARSE (IC) 
      ENDDO

   ENDDO

!> ------------------------------------------------------------------------------------------------
!> Serial version
!> ------------------------------------------------------------------------------------------------
ELSE

   DO NM = NMESHES_MIN, NMESHES_MAX
   
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
   
   DO NM = NMESHES_MIN, NMESHES_MAX
   
      SL => SCARC(NM)%LEVEL(NLEVEL_MAX)
      IOFFSET = SM%OFFSET(NM)
  
      WRITE(*,*) 'ACHTUNG; WAS IST HIER MIT BVECTOR???'
      !VBX => POINT_TO_HVECTOR (NVECTORX, NM)
 
      DO IC = 1, SL%NC
         VCX(IC) = SM%X_COARSE (IC) 
      ENDDO

   ENDDO

ENDIF

TUSED_SCARC(NSCARC_TIME_COARSE,:)=TUSED_SCARC(NSCARC_TIME_COARSE,:)+SECOND()-TNOW_COARSE
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_COARSE
END SUBROUTINE SCARC_METHOD_GE


 
!> ------------------------------------------------------------------------------------------------
!> Setup environement in every solver CALL (i.e. set pointers to used vectors)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SCOPE(SCOPE, NSCOPE, NPRECON, NRHS, NL)
INTEGER, INTENT(IN)  :: NSCOPE, NPRECON, NRHS
INTEGER, INTENT(OUT) :: NL
TYPE (SCARC_SCOPE_TYPE), INTENT(OUT):: SCOPE

TYPE_SCOPE  = NSCOPE

SELECT CASE (NSCOPE)
   CASE (NSCARC_SCOPE_COARSE)
      TYPE_METHOD = NSCARC_METHOD_KRYLOV
   CASE (NSCARC_SCOPE_PRECON)
      TYPE_METHOD = NSCARC_METHOD_MULTIGRID
   CASE DEFAULT
      TYPE_METHOD = TYPE_METHOD
END SELECT

SCOPE%RESIN = 1.0_EB
SCOPE%RES   = 0.0_EB

SELECT CASE (TYPE_METHOD)
   
 !! ---------------------------------------------------------------------------------------------------
 !! Krylov method
 !! ---------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      TYPE_PRECON = NPRECON

      SELECT CASE (TYPE_KRYLOV)
      
         !< CG-method
         CASE (NSCARC_KRYLOV_CG)
      
            SCOPE%F = NRHS                                             !< set correct right hand side vector

            SELECT CASE (NSCOPE)
               CASE (NSCARC_SCOPE_MAIN)

                  NL = NLEVEL_MIN

                  SCOPE%CROUTINE = 'SCARC_GLOBAL_CG'
                  SCOPE%OMEGA    =  SCARC_PRECON_OMEGA
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
                  SCOPE%OMEGA    =  SCARC_COARSE_OMEGA
                  SCOPE%EPS      =  SCARC_COARSE_ACCURACY
                  SCOPE%NIT      =  SCARC_COARSE_ITERATIONS

                  IF (TYPE_METHOD0 == NSCARC_METHOD_MULTIGRID) THEN              !<  pure MG-method
                     SCOPE%X = NSCARC_VECTOR_X                                   !<  take first-stage MG-vectors
                     SCOPE%D = NSCARC_VECTOR_D 
                     SCOPE%G = NSCARC_VECTOR_G 
                     SCOPE%Y = NSCARC_VECTOR_Y 
                     SCOPE%W = NSCARC_VECTOR_W 
                  ELSE                                                           !< CG-MG-method
                     SCOPE%X = NSCARC_VECTOR_X2                                  !< take second-stage MG-vectors
                     SCOPE%D = NSCARC_VECTOR_D2
                     SCOPE%G = NSCARC_VECTOR_G2
                     SCOPE%Y = NSCARC_VECTOR_Y2
                     SCOPE%W = NSCARC_VECTOR_W2
                  ENDIF

                  !TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
            END SELECT
      
         !< BICG-method
         CASE (NSCARC_KRYLOV_BICG)
      
            SCOPE%F = NRHS                                             !< set correct right hand side vector

            SELECT CASE (NSCOPE)
               CASE (NSCARC_SCOPE_MAIN)

                  NL = NLEVEL_MIN

                  SCOPE%CROUTINE = 'SCARC_GLOBAL_BICG'
                  SCOPE%OMEGA    =  SCARC_PRECON_OMEGA
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
                  SCOPE%OMEGA    =  SCARC_COARSE_OMEGA
                  SCOPE%EPS      =  SCARC_COARSE_ACCURACY
                  SCOPE%NIT      =  SCARC_COARSE_ITERATIONS

                  IF (TYPE_METHOD0 == NSCARC_METHOD_MULTIGRID) THEN              !<  pure MG-method
                     SCOPE%X = NSCARC_VECTOR_X                                   !<  take first-stage MG-vectors
                     SCOPE%D = NSCARC_VECTOR_D 
                     SCOPE%G = NSCARC_VECTOR_G 
                     SCOPE%Y = NSCARC_VECTOR_Y 
                     SCOPE%W = NSCARC_VECTOR_W 
                     SCOPE%Z = NSCARC_VECTOR_Z 
                  ELSE                                                           !< CG-MG-method
                     SCOPE%X = NSCARC_VECTOR_X2                                  !< take second-stage MG-vectors
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

      TYPE_SMOOTH = NPRECON
      TYPE_PRECON = NPRECON

      SCOPE%EPS   = SCARC_MULTIGRID_ACCURACY
      SCOPE%NIT   = SCARC_MULTIGRID_ITERATIONS
      SCOPE%OMEGA = SCARC_SMOOTH_OMEGA

      NL = NLEVEL_MIN
      SCOPE%F = NRHS                                                   !< set correct right hand side vector

      !< select scope (multigrid as main solver or preconditioner)
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

END SELECT
   
END SUBROUTINE SCARC_SETUP_SCOPE


!> ----------------------------------------------------------------------------------------------------
!> Set initial solution corresponding to boundary data in BXS, BXF, ...
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SOLVER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, I, J, K, IC
REAL(EB):: F_OLD!, LINE(10)
TYPE (MESH_TYPE), POINTER ::  M
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL


SELECT_METHOD: SELECT CASE (TYPE_METHOD)

 !! In case of a Krylov method initialize auxiliary arrays
   CASE (NSCARC_METHOD_KRYLOV)

      DO NM = NMESHES_MIN, NMESHES_MAX
         SL => SCARC(NM)%LEVEL(NL)
         DO IC = SL%NC+1, SL%NCE
            SL%D(IC) = 0.0_EB
            SL%G(IC) = 0.0_EB
            SL%Y(IC) = 0.0_EB
            SL%W(IC) = 0.0_EB
         ENDDO
      ENDDO

      IF (TYPE_KRYLOV == NSCARC_KRYLOV_BICG) THEN
         DO NM = NMESHES_MIN, NMESHES_MAX
            SL => SCARC(NM)%LEVEL(NL)
            DO IC = SL%NC+1, SL%NCE
               SL%Z(IC) = 0.0_EB
            ENDDO
         ENDDO
      ENDIF
      
 !! In case of a multigrid method with coarse grid solution by CG, clear CG-vectors on max level
   CASE (NSCARC_METHOD_MULTIGRID)
      
      DO NM = NMESHES_MIN, NMESHES_MAX
         SL => SCARC(NM)%LEVEL(NL)
         DO IC = SL%NC+1, SL%NCE
            SL%D = 0.0_EB
         ENDDO
         SL => SCARC(NM)%LEVEL(NLEVEL_MAX)
         DO IC = SL%NC+1, SL%NCE
            SL%G = 0.0_EB
            SL%W = 0.0_EB
         ENDDO
         
      ENDDO

END SELECT SELECT_METHOD

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   
 !! Initialize solution and right hand side vector corresponding to boundary conditions
   SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
   
      !< -------------- 2D --------------
      CASE (NSCARC_DIMENSION_TWO)
   
         DO NM = NMESHES_MIN, NMESHES_MAX
   
            M  => MESHES(NM)
            SL => SCARC(NM)%LEVEL(NL)
         
            !< get right hand side and initial vector
            DO K = 1, M%KBAR
               DO I = 1, M%IBAR
                  IC = (K-1)*M%IBAR + I
                  SL%F(IC) = M%PRHS (I, 1, K)
                  IF (PREDICTOR) THEN
                     SL%X(IC) = M%H (I, 1, K)
                  ELSE
                     SL%X(IC) = M%HS(I, 1, K)
                  ENDIF
               ENDDO
            ENDDO
      
      
            !< set boundary conditions at exterior boundaries (corresponding to pois.f90)
            LEVEL_WALLCELLS_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
         
               I    = SL%WALL(IW)%IXW 
               J    = SL%WALL(IW)%IYW 
               K    = SL%WALL(IW)%IZW 
               IC   = SL%WALL(IW)%ICW 
               IOR0 = SL%WALL(IW)%IOR
      
               IF (J /= 1) THEN
                  WRITE(*,*) 'Wrong index for J =',J,' in SCARC_SETUP_SOLVER !!!'
                  STOP
               ENDIF

               SELECT CASE (IOR0)
                  CASE (1)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DXI2 * M%BXS(1,K)         !< Dirichlet
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) + SL%DXI * M%BXS(1,K)                   !< Neumann
                     ENDIF
                  CASE (-1)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DXI2 *M%BXF(1,K)          !< Dirichlet
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DXI *M%BXF(1,K)                    !< Neumann
                     ENDIF
                  CASE (3)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DZI2 * M%BZS(I,1)         !< Dirichlet
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) + SL%DZI * M%BZS(I,1)                   !< Neumann
                     ENDIF
                  CASE (-3)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DZI2 * M%BZF(I,1)         !< Dirichlet
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DZI  * M%BZF(I,1)                  !< Neumann
                     ENDIF
               END SELECT
            
            ENDDO LEVEL_WALLCELLS_LOOP2D
   
         ENDDO
      !< -------------- 3D --------------
      CASE (NSCARC_DIMENSION_THREE)
      
         DO NM = NMESHES_MIN, NMESHES_MAX
   
            M  => MESHES(NM)
            SL => SCARC(NM)%LEVEL(NL)
      
            !< get right hand side and initial vector
            DO K = 1, M%KBAR
               DO J = 1, M%JBAR
                  DO I = 1, M%IBAR
                     IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                     SL%F(IC) = M%PRHS (I, J, K)
                     IF (PREDICTOR) THEN
                        SL%X(IC) = M%H (I, J, K)
                     ELSE
                        SL%X(IC) = M%HS(I, J, K)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
      
            LEVEL_WALLCELLS_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
            
               I    = SL%WALL(IW)%IXW 
               J    = SL%WALL(IW)%IYW 
               K    = SL%WALL(IW)%IZW 
               IC   = SL%WALL(IW)%ICW 
               IOR0 = SL%WALL(IW)%IOR
            
               SELECT CASE (IOR0)
                  CASE (1)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        F_OLD = SL%F(IC)
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DXI2 * M%BXS(J,K)         !< Dirichlet
!IF (J==18.AND.(K==18)) &
!WRITE(*,'(a,i3,a,i3,a,i3,a,3f18.10)') 'A:F(',I,',',J,',',K,')=',SL%F(IC), F_OLD, M%BXS(J,K)
!IF (J==18.AND.(K==19)) &
!WRITE(*,'(a,i3,a,i3,a,i3,a,3f18.10)') 'A:F(',I,',',J,',',K,')=',SL%F(IC), F_OLD, M%BXS(J,K)
!IF (J==19.AND.(K==18)) &
!WRITE(*,'(a,i3,a,i3,a,i3,a,3f18.10)') 'A:F(',I,',',J,',',K,')=',SL%F(IC), F_OLD, M%BXS(J,K)
!IF (J==19.AND.(K==19)) &
!WRITE(*,'(a,i3,a,i3,a,i3,a,3f18.10)') 'A:F(',I,',',J,',',K,')=',SL%F(IC), F_OLD, M%BXS(J,K)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        F_OLD = SL%F(IC)
                        SL%F(IC) = SL%F(IC) + SL%DXI * M%BXS(J,K)                   !< Neumann
!IF (J==18.AND.(K==18)) &
!WRITE(*,'(a,i3,a,i3,a,i3,a,3f18.10)') 'A:F(',I,',',J,',',K,')=',SL%F(IC), F_OLD, M%BXS(J,K)
!IF (J==18.AND.(K==19)) &
!WRITE(*,'(a,i3,a,i3,a,i3,a,3f18.10)') 'A:F(',I,',',J,',',K,')=',SL%F(IC), F_OLD, M%BXS(J,K)
!IF (J==19.AND.(K==18)) &
!WRITE(*,'(a,i3,a,i3,a,i3,a,3f18.10)') 'A:F(',I,',',J,',',K,')=',SL%F(IC), F_OLD, M%BXS(J,K)
!IF (J==19.AND.(K==19)) &
!WRITE(*,'(a,i3,a,i3,a,i3,a,3f18.10)') 'A:F(',I,',',J,',',K,')=',SL%F(IC), F_OLD, M%BXS(J,K)
                     ENDIF
                  CASE (-1)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DXI2 *M%BXF(J,K)          !< Dirichlet
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'C:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DXI *M%BXF(J,K)                    !< Neumann
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'D:F(',I,',',J,',',K,')=',SL%F(IC)
                     ENDIF
                  CASE (2)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DYI2 * M%BYS(I,K)         !< Dirichlet
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'E:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) + SL%DYI * M%BYS(I,K)                   !< Neumann
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'F:F(',I,',',J,',',K,')=',SL%F(IC)
                     ENDIF
                  CASE (-2)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DYI2 *M%BYF(I,K)          !< Dirichlet
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'G:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DYI *M%BYF(I,K)                    !< Neumann
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'H:F(',I,',',J,',',K,')=',SL%F(IC)
                     ENDIF
                  CASE (3)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DZI2 * M%BZS(I,J)         !< Dirichlet
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'I:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) + SL%DZI * M%BZS(I,J)                   !< Neumann
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'J:F(',I,',',J,',',K,')=',SL%F(IC)
                     ENDIF
                  CASE (-3)
                     IF (SL%WALL(IW)%BTYPE==DIRICHLET) THEN
                        SL%F(IC) = SL%F(IC) - 2.0_EB * SL%DZI2 * M%BZF(I,J)         !< Dirichlet
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'K:F(',I,',',J,',',K,')=',SL%F(IC)
                     ELSE IF (SL%WALL(IW)%BTYPE==NEUMANN) THEN
                        SL%F(IC) = SL%F(IC) - SL%DZI  * M%BZF(I,J)                  !< Neumann
!WRITE(*,'(a,i3,a,i3,a,i3,a,f18.10)') 'L:F(',I,',',J,',',K,')=',SL%F(IC)
                     ENDIF
               END SELECT
            
            ENDDO LEVEL_WALLCELLS_LOOP3D
   
         ENDDO
   
!WRITE(*,'(a10,10f12.6)') 'F(,5,5)', LINE(1:10)

   END SELECT SELECT_DIMENSION

ELSE IF (TYPE_SCOPE == NSCARC_SCOPE_COARSE) THEN
   DO NM = NMESHES_MIN, NMESHES_MAX
      SCARC(NM)%LEVEL(NLEVEL_MAX)%X = 0.0_EB
   ENDDO
ENDIF


END SUBROUTINE SCARC_SETUP_SOLVER

!> ------------------------------------------------------------------------------------------------
!> Print out residual information for loop ITE
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_INFO(RES, ITE, NL, CROUTINE)
INTEGER, INTENT(IN) :: ITE, NL
REAL(EB), INTENT(IN) :: RES
CHARACTER(*), INTENT(IN) :: CROUTINE
INTEGER:: NM

DO NM = NMESHES_MIN, NMESHES_MAX
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(SCARC_LU,1000) TRIM(CROUTINE), NM, NL, ITE,  RES
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS.AND.MYID==0) write(*,1000) TRIM(CROUTINE), NM, NL, ITE,  RES
ENDDO

1000 FORMAT (5X,A30,': mesh', i4,': level=',i4,': #ite= ',i4,': res =',e25.16)
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

IF (TYPE_DEBUG >= NSCARC_DEBUG_INFO2) WRITE(SCARC_LU,1000) TRIM(SCOPE%CROUTINE), NL, ITE, SCOPE%RES
IF (TYPE_DEBUG >= NSCARC_DEBUG_INFO2.AND.MYID==0) WRITE(*,1000) TRIM(SCOPE%CROUTINE), NL, ITE, SCOPE%RES

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

IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH.AND.TYPE_METHOD==TYPE_METHOD0.AND.TYPE_SCOPE==NSCARC_SCOPE_MAIN) THEN
   WRITE(SCARC_LU,2000) SCOPE%CROUTINE, ITERATIONS, CAPPA 
   IF (MYID==0) WRITE(*,2000) SCOPE%CROUTINE, ITERATIONS, CAPPA 
ELSE IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(SCARC_LU,2000) SCOPE%CROUTINE, ITERATIONS, CAPPA 
   IF (MYID==0) WRITE(*,2000) SCOPE%CROUTINE, ITERATIONS, CAPPA 
ENDIF
IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) THEN
   WRITE(SCARC_LU,3000)
   IF (MYID==0) WRITE(*,3000)
ENDIF
IF (TYPE_DEBUG == NSCARC_DEBUG_MUCH.AND.TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   WRITE(SCARC_LU,3000)
   IF (MYID==0) WRITE(*,3000)
ENDIF

2000 FORMAT (5X,A30,': iterations: ',i6,':  convergence rate =',e14.6)
3000 FORMAT ('==========================================================================')
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

TYPE_VECTOR = NVECTOR_FI
IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC.AND.TYPE_COARSENING < NSCARC_COARSENING_GMG) &
   CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_VECTOR, NL_FI)

SELECT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

 !! ------------------------- Geometric multigrid -------------------------------------------
   CASE (NSCARC_MULTIGRID_GEOMETRIC)

      DO NM = NMESHES_MIN, NMESHES_MAX
   
         NX_CO => SCARC(NM)%LEVEL(NL_CO)%NX
         NY_CO => SCARC(NM)%LEVEL(NL_CO)%NY
         NZ_CO => SCARC(NM)%LEVEL(NL_CO)%NZ
   
         NX_FI = 2*NX_CO
         NY_FI = 2*NY_CO
         NZ_FI = 2*NZ_CO
   
         DC_FI  => POINT_TO_VECTOR(NVECTOR_FI, NM, NL_FI)
         FC_CO  => POINT_TO_VECTOR(NVECTOR_CO, NM, NL_CO)

         SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !< ----------------- 2D ---------------
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
         
IF (TYPE_DEBUG>= NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'RESTRICT1: FC_CO(',IC_CO,')=',FC_CO(IC_CO)
                  ENDDO
               ENDDO
            
            !< ----------------- 3D ---------------
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

      DO NM = NMESHES_MIN, NMESHES_MAX

         DC_FI => POINT_TO_VECTOR(NVECTOR_FI, NM, NL_FI)
         FC_CO => POINT_TO_VECTOR(NVECTOR_CO, NM, NL_CO)

         NC_CO => SCARC(NM)%LEVEL(NL_CO)%NC

         R     => SCARC(NM)%LEVEL(NL_FI)%R
         R_ROW => SCARC(NM)%LEVEL(NL_FI)%R_ROW
         R_COL => SCARC(NM)%LEVEL(NL_FI)%R_COL
   
         DO IC_CO = 1, NC_CO
            AUX = 0.0_EB
            DO ICOL = R_ROW(IC_CO), R_ROW(IC_CO+1)-1
               IC = R_COL(ICOL)
               AUX0 = AUX
               AUX = AUX + DC_FI(IC) * R(ICOL)
            ENDDO
            FC_CO(IC_CO) = AUX
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'RESTRICT2: FC_CO(',IC_CO,')=',FC_CO(IC_CO)
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


SELECT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

 !! ------------------------- Geometric multigrid -------------------------------------------
   CASE (NSCARC_MULTIGRID_GEOMETRIC)

      DO NM = NMESHES_MIN, NMESHES_MAX
   
         NX_CO => SCARC(NM)%LEVEL(NL_CO)%NX
         NY_CO => SCARC(NM)%LEVEL(NL_CO)%NY
         NZ_CO => SCARC(NM)%LEVEL(NL_CO)%NZ
   
         NX_FI = 2*NX_CO
         NY_FI = 2*NY_CO
         NZ_FI = 2*NZ_CO
   
         XC_CO  => POINT_TO_VECTOR(NVECTOR_CO, NM, NL_CO)
         DC_FI  => POINT_TO_VECTOR(NVECTOR_FI, NM, NL_FI)
   
         SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !< ----------------- 2D ---------------
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
IF (TYPE_DEBUG>= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'PROL1:: DC_FI(',IC_FI(I),')=',DC_FI(IC_FI(I))
                     ENDDO
         
                  ENDDO
               ENDDO
         
                        
            !< ----------------- 3D ---------------
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
   
         END SELECT SELECT_DIMENSION

      ENDDO

 !! ------------------------- Algebraic multigrid -------------------------------------------
   CASE (NSCARC_MULTIGRID_ALGEBRAIC)

SCAL = 1.0_EB
IF (TYPE_COARSENING >= NSCARC_COARSENING_GMG) SCAL=2.0_EB

      DO NM = NMESHES_MIN, NMESHES_MAX

         XC_CO  => POINT_TO_VECTOR(NVECTOR_CO, NM, NL_CO)
         DC_FI  => POINT_TO_VECTOR(NVECTOR_FI, NM, NL_FI)

         NC_FI  => SCARC(NM)%LEVEL(NL_FI)%NC
         NCE_FI => SCARC(NM)%LEVEL(NL_FI)%NCE

         P      => SCARC(NM)%LEVEL(NL_FI)%P
         P_ROW  => SCARC(NM)%LEVEL(NL_FI)%P_ROW
         P_COL  => SCARC(NM)%LEVEL(NL_FI)%P_COL
         P_PTR  => SCARC(NM)%LEVEL(NL_FI)%P_PTR
   
         DO IC = 1, NC_FI
            AUX = 0.0_EB
            DO ICOL = P_ROW(IC), P_ROW(IC+1)-1
               IC_CO = P_COL(ICOL)
               AUX = XC_CO(IC_CO) * P(ICOL)
            ENDDO
            DC_FI(IC) = SCAL*AUX
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'PROL2:: DC_FI(',IC,')=',DC_FI(IC)
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

DO NM = NMESHES_MIN, NMESHES_MAX

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
            IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
            HP(I, J, K) = SL%X(IC)
         ENDDO
      ENDDO
   ENDDO

ENDDO

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


DO NM = NMESHES_MIN, NMESHES_MAX

 !! point to correct pressure vector on mesh 'NM'
   M  => MESHES(NM)
   SL => SCARC(NM)%LEVEL(NL)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF
   
 !! compute ghost cell values
   WALLCELLS_LOOP: DO IW = 1, SL%NW
   
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
   ENDDO WALLCELLS_LOOP
   
ENDDO


!> -----------------------------------------------------------------------------------------------
!> Perform data exchange to achieve consistency of ghost values along internal boundaries 
!> -----------------------------------------------------------------------------------------------
!WRITE(*,*) 'SUSI: BEFORE EXCHANGE_PRESURE'
IF (NMESHES > 1) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_PRESSURE, NL)
!WRITE(*,*) 'SUSI: AFTER  EXCHANGE_PRESURE'

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
RECEIVE_MESH_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   RECEIVE_OMESH_LOOP: DO NOM = 1, NMESHES
    
      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_OMESH_LOOP

      S   => SCARC(NM)                           
      OS  => SCARC(NM)%OSCARC(NOM)               

      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE RECEIVE_OMESH_LOOP

      SELECT_EXCHANGE_TYPE: SELECT CASE (TYPE_EXCHANGE)

      
         !< ---------------------------------------------------------------------------------------
         !< Exchange information about neighboring step size along internal boundary
         !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER0(1),1,MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !< ---------------------------------------------------------------------------------------
         !< Exchange information about neighboring wall data
         !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALLINFO)

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'RECEIVING IN SIZE ', SIZE(OS%RECV_INTEGER)
            OSL => OS%LEVEL(NL)
            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),SIZE(OS%RECV_INTEGER),MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
   
         !< ---------------------------------------------------------------------------------------
         !< Exchange information about neighboring grid dimensions
         !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MESHINFO)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),5,MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !< ---------------------------------------------------------------------------------------
         !< Exchange information about neighboring step size along internal boundary
         !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WIDTHINFO)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_REAL0(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !< ---------------------------------------------------------------------------------------
         !< Exchange number of neighboring cells for AMG method (compact type only)
         !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_TRANSFER_SIZE)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),9,MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !< ---------------------------------------------------------------------------------------
         !< Exchange number of neighboring cells for AMG method (compact type only)
         !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),1,MPI_INTEGER,SNODE, &
                           TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

         !< ---------------------------------------------------------------------------------------
         !< Exchange neighboring celltypes
         !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_TYPE)

            NREQ_SCARC = NREQ_SCARC+1
            CALL MPI_IRECV(OS%RECV_INTEGER(1),SIZE(OS%RECV_INTEGER),MPI_INTEGER,&
                           SNODE,TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)


         !< ---------------------------------------------------------------------------------------
         !< Perform exchanges for 
         !<    - internal values for matrix-vector multiplication 
         !<    - internal boundariy values
         !<    - internal subdiagonal matrix values
         !<    - internal subdiagonal or ghost matrix values
         !<    - internal measure/CELL_TYPE values
         !< ---------------------------------------------------------------------------------------
         CASE DEFAULT

            NREQ_SCARC = NREQ_SCARC+1
!WRITE(SCARC_LU,*) 'RECEIVING IN LENGTH ', SIZE(OS%RECV_REAL)
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
INTEGER  :: NM, NOM, NCPL, NCOL, NROW
INTEGER  :: IERR, IW, IWL, IWG, IWW, IW0, IROW, ICOL, IPTR, ICPL, ICELL_TYPE
INTEGER  :: IOR0, IC, JC, ICC, JCC, ICE, IC1, JC0, ICO, ICN, ICG, ICW
INTEGER  :: I, J, K, II, JJ, KK, LL
REAL(EB) :: ZSUM, ZSUM1
INTEGER , POINTER, DIMENSION(:)     ::  RECV_INTEGER
REAL(EB), POINTER, DIMENSION(:)     ::  RECV_REAL
REAL(EB), POINTER, DIMENSION(:)     ::  VECTOR
REAL(EB), POINTER, DIMENSION(:,:,:) :: HVECTOR
TYPE (SCARC_TYPE)      , POINTER :: S 
TYPE (OSCARC_TYPE)     , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, OSL, SLO

IERR = 0

!WRITE(SCARC_LU,*) 'STARTING EXCHANGE_SEND, TYPE_EXCHANGE=',TYPE_EXCHANGE

!> ------------------------------------------------------------------------------------------------
!> Collect data for sending corresponding to requested exchange type
!> ------------------------------------------------------------------------------------------------
MESH_PACK_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

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

           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_BASIC: pack neighboring sizes for wallinfo exchange
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC) 

IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(a,i4,a,i4,a,3i4)') 'BASIC: NM=',NM,': NOM=',NOM,': SENDING ',OS%LEVEL(NL)%NWL

            IF (RNODE /= SNODE) THEN
               OS%SEND_INTEGER0(1)=OS%LEVEL(NL)%NWL
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER0(1),1,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF

           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_WIDTHINFO: pack neighboring grid resolution information
           !< ---------------------------------------------------------------------------------------
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

            IF (RNODE /= SNODE) THEN
               OS%SEND_REAL0(1) = OSL%DH
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(a,i4,a,i4,a,f12.6)') 'RESOLUTION: NM=',NM,': NOM=',NOM,': SENDING ',OS%SEND_REAL0(1)
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL0(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF


           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_WALLINFO: pack neighboring wallinfo information
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALLINFO)

            IF (RNODE /= SNODE) THEN

               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'NWL=',OSL%NWL
               IPTR=1
               DO IWL = 1, OSL%NWL
                  IWG = OSL%IWL_TO_IWG(IWL)
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
                 WRITE(SCARC_LU,'(a,i4,a,i4,a,i4)') 'IWL=',IWL,': IWG=',IWG,': NCPLS=',OSL%NCPLS
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
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(i4,i4)') (I, OS%SEND_INTEGER(I), I=IPTR,IPTR+11)
                  IPTR = IPTR + 12
                  DO ICPL=1,OSL%NCPLS
                     OS%SEND_INTEGER(IPTR)=SL%WALL(IWG)%ICE(ICPL)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(i4,i4)') (I, OS%SEND_INTEGER(I), I=IPTR,IPTR)
                     IPTR = IPTR + 1
                  ENDDO
                  DO ICPL=1,OSL%NCPLS
                     OS%SEND_INTEGER(IPTR)=SL%WALL(IWG)%ICN(ICPL)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(i4,i4)') (I, OS%SEND_INTEGER(I), I=IPTR,IPTR)
                     IPTR = IPTR + 1
                  ENDDO
               ENDDO
               NREQ_SCARC = NREQ_SCARC+1
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'SENDING IN SIZE ', SIZE(OS%SEND_INTEGER)
               CALL MPI_ISEND(OS%SEND_INTEGER(1),SIZE(OS%SEND_INTEGER),MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)

            ENDIF


           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_MESHINFO: pack neighboring mesh information
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MESHINFO) 

            IF (RNODE /= SNODE) THEN
               SL  =>  S%LEVEL(NL)
               OS%SEND_INTEGER(1)=SL%NX
               OS%SEND_INTEGER(2)=SL%NY
               OS%SEND_INTEGER(3)=SL%NZ
               OS%SEND_INTEGER(4)=SL%NC
               OS%SEND_INTEGER(5)=SL%NW
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(a,i4,a,i4,a,5i4)') 'MESHINFO: NM=',NM,': NOM=',NOM,': SENDING ',OS%SEND_INTEGER(1:5)
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER(1),5,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF



           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_PRESSURE: pack overlapping parts of a H or HS
         ! (pressure vectors for PREDICTOR and CORRECTOR
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_PRESSURE)

!WRITE(*,*) 'STARTING TO PACK PRESSURE ', OSL%NCG
            IF (PREDICTOR) THEN
               HVECTOR => POINT_TO_HVECTOR (NSCARC_VECTOR_H, NM)
            ELSE
               HVECTOR => POINT_TO_HVECTOR (NSCARC_VECTOR_H, NM)
            ENDIF
            
            LL  = 1
            PACK_PRESSURE: DO ICG=1, OSL%NCG
               DO KK=OSL%WALL(ICG)%IZN(1),OSL%WALL(ICG)%IZN(2)
                  DO JJ=OSL%WALL(ICG)%IYN(1),OSL%WALL(ICG)%IYN(2)
                     DO II=OSL%WALL(ICG)%IXN(1),OSL%WALL(ICG)%IXN(2)
                        OS%SEND_REAL(LL) = HVECTOR(II,JJ,KK)
!WRITE(SCARC_LU,*) 'PACKING H(',II,JJ,KK,')=',HVECTOR(II,JJ,KK)
                        LL = LL+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO PACK_PRESSURE

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_VECTOR: pack overlapping parts of a given vector
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_VECTOR)

            VECTOR => POINT_TO_VECTOR(TYPE_VECTOR, NM, NL)

            IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
               WRITE(SCARC_LU,*) 'VECTOR TO BE SENT:', SL%NC, TYPE_VECTOR, OSL%NCG, OSL%NCPLR
               WRITE(SCARC_LU,'(4f18.12)') (VECTOR(ICG), ICG=1,SL%NC)
            ENDIF

            LL = 1
            PACK_VECTOR: DO ICG = 1, OSL%NCG
               DO ICPL=1,OSL%NCPLR
                  ICN = OSL%WALL(ICG)%ICN(ICPL) 
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'ICG=',ICG,': ICN=', ICN, OSL%NCPLR
                  OS%SEND_REAL(LL) = VECTOR(ICN)
               IF (TYPE_DEBUG > NSCARC_DEBUG_LESS)  &
                  WRITE(SCARC_LU,'(a,i6,a,f12.6,2i6)') 'SEND_REAL(',LL,')=',OS%SEND_REAL(ICG), ICG, ICN
                  LL = LL + 1
               ENDDO
            ENDDO PACK_VECTOR

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


         !! ---------------------------------------------------------------------------------------
         !! EXCHANGE_MEASURE_ADD: pack neighboring measure information (AMG only) with adding
         !! ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MEASURE_ADD)

            LL  = 1
            PACK_MEASURE_ADD: DO IWL=1, OSL%NWL
               ICW = OSL%ICG_TO_ICW(IWL)
               OS%SEND_REAL(LL) = SL%MEASURE(ICW)
               IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,'(a,i4,a,f12.6,3i3)') &
                  'PACKA VECTOR(',IC,')=',VECTOR(ICW), IWL, OSL%WALL(IWL)%ICE(1)
               LL = LL+1
            ENDDO PACK_MEASURE_ADD
                  
            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


           !< ---------------------------------------------------------------------------------------
           !< Send data along internal boundaries corresponding to requested exchange type
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_TYPE)

            WRITE(*,*) 'HALLO, HIER NOCHMAL CHECKEN !!!'

            LL  = 1
            PACK_CELL_TYPE: DO IWL=1, OSL%NWL
               ICW = OSL%IWL_TO_ICW(IWL)
               OS%SEND_INTEGER(LL) = SL%CELL_TYPE(ICW)
            IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) WRITE(SCARC_LU,*) 'SENDING VECTOR(',IC,')=',VECTOR(IC)
               LL = LL+1
            ENDDO PACK_CELL_TYPE

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER(1), SIZE(OS%SEND_INTEGER), MPI_INTEGER, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_MATRIX_SIZE: Send sizes of overlapping matrix areas
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE) 

            IF (RNODE /= SNODE) THEN
               OSL => OS%LEVEL(NL)
               OS%SEND_INTEGER(1) = OSL%NA0
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,*) 'SIZE_MATRIX: SENDING IN SIZE 1: ', OS%SEND_INTEGER(1)
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_INTEGER(1),1,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF


           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_MATRIX_STENCIL: pack subdiagonal entries of system matrix on overlapping parts
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG)

            LL = 1
            IWW = 0
        
            WRITE(*,*) 'EXCHANGE_MATRIX_SUBDIAG, PACK: Achtung, hier nochmal checken für AMG!!'

            PACK_MATRIX_SUBDIAG1: DO IWL=1,OSL%NCG
               IOR0 = ABS(OSL%WALL(IWL)%IOR)
               LL = LL + 1
               DO ICPL = 1, OSL%NCPL
                  IC = OSL%WALL(IWL)%ICN(ICPL)
                  OS%SEND_REAL(LL) = SL%DI2(IOR0)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) &
                     'PACK_MATRIX_SUBDIAG1: IC=',IC,': OS%SEND_REAL(',LL,')=', OS%SEND_REAL(LL)
                  LL = LL+1
               ENDDO
            ENDDO PACK_MATRIX_SUBDIAG1
         
            IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
               WRITE(SCARC_LU,*) 'NSCARC_EXCHANGE_MATRIX_SUBDIAG: SEND_REAL'
               WRITE(SCARC_LU,'(8f12.6)') (OS%SEND_REAL(IC), IC=1, SIZE(OS%SEND_REAL))
            ENDIF

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF



           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_MATRIX_STENCIL: pack matrix stencil information
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_STENCIL)

             !< Pack first cell layer
            LL = 1
            OS%SEND_REAL = 0.0_EB

            PACK_MATRIX_STENCIL: DO IWL=1, OSL%NCG
               IF (OSL%WALL(IWL)%NOM/=NM) CYCLE PACK_MATRIX_STENCIL
         
               NCPL = OSL%WALL(IWL)%NCPL
      
               DO ICPL = 1, NCPL
                  IC = OSL%WALL(IWL)%ICN(ICPL)
      
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'ICPL=',ICPL,': IC=',IC, ': IWL=',IWL
      
                  OS%SEND_REAL(LL)   = REAL(IWL,EB)
                  NCOL = 0
                  DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
                     JC = SL%A_COL(ICOL)
                     IF (JC > SL%NC) THEN
                        IW0 = SL%ICG_TO_IWG(JC)
                        !IF (SL%OSL%WALL(IW0)%NOM == NOM) NCOL = NCOL + 1
                        NCOL = NCOL + 1
                     ELSE
                        NCOL = NCOL + 1
                     ENDIF
                  ENDDO
                  OS%SEND_REAL(LL+1) = REAL(NCOL)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6)') &
                                         'PACK_MATRIX_STENCIL: IC=',IC,': OS%SEND_REAL(',LL+1,')=', OS%SEND_REAL(LL+1)
      
                  LL = LL + 2
               ENDDO 
      
                !< Pack second cell layer if requested
               IF (NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) THEN
         
                  DO ICPL = 1, NCPL
                     IC = OSL%WALL(IWL)%ICN2(ICPL)
         
                     OS%SEND_REAL(LL) = REAL(IWL,EB)
                     NCOL = 0
                     DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
                        JC = SL%A_COL(ICOL)
                        IF (JC > SL%NC) THEN
                           IW0 = SL%ICG_TO_IWG(JC)
                           IF (OSL%WALL(IW0)%NOM == NOM) NCOL = NCOL + 1
                        ELSE
                           NCOL = NCOL + 1
                        ENDIF
                     ENDDO
                     OS%SEND_REAL(LL+1) = REAL(NCOL)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6)') &
                                          'PACK_MATRIX_STENCIL: IC=',IC,': OS%SEND_REAL(',LL+1,')=', OS%SEND_REAL(LL+1)
         
                     LL = LL + 2
                  ENDDO
               ENDIF
            ENDDO PACK_MATRIX_STENCIL
      
            OS%SEND_REAL(LL) = -999.0_EB
            LL = LL + 1

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_MATRIX_SYSTEM: pack overlapping parts of system matrix
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SYSTEM)


            LL = 1
            OS%SEND_REAL = 0.0_EB
      
            !< Pack first cell layer
            PACK_MATRIX_SYSTEM: DO IWL=1, OSL%NCG
               IF (OSL%WALL(IWL)%NOM/=NM) CYCLE PACK_MATRIX_SYSTEM
         
               NCPL = OSL%WALL(IWL)%NCPL
      
               DO ICPL = 1, NCPL
                  IC = OSL%WALL(IWL)%ICN(ICPL)
      
                  OS%SEND_REAL(LL)   = REAL(IWL,EB)
                  LL = LL + 1
      
                  MATRIX_COLUMN_LOOP: DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
                     JC = SL%A_COL(ICOL)
                     IF (JC > SL%NC) THEN
                        IW0 = SL%ICG_TO_IWG(JC)
                        !IF (SL%OSL%WALL(IW0)%NOM /= NOM) CYCLE MATRIX_COLUMN_LOOP
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6)') &
                                         '==========================    ICG_TO_IWG(',JC,')=',SL%ICG_TO_IWG(JC) 
                        OS%SEND_REAL(LL) = - REAL(SL%ICG_TO_IWG(JC),EB)
                        OS%SEND_REAL(LL) = - REAL(JC)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6)') &
                                         'PACK_MATRIX_SYSTEM: IC=',IC,': OS%SEND_REAL(',LL,')=', OS%SEND_REAL(LL)
                     ELSE
                        OS%SEND_REAL(LL) =   REAL(JC,EB)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6)') &
                                         'PACK_MATRIX_SYSTEM: IC=',IC,': OS%SEND_REAL(',LL,')=', OS%SEND_REAL(LL)
                     ENDIF
                     OS%SEND_REAL(LL+1) = SL%A(ICOL)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,f12.6)') &
                                         'PACK_MATRIX_SYSTEM: IC=',IC,': OS%SEND_REAL(',LL+1,')=', OS%SEND_REAL(LL+1)
                     LL = LL + 2
                  ENDDO MATRIX_COLUMN_LOOP
      
               ENDDO 
      
               !< Pack second cell layer if requested
               IF (NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) THEN
         
                  DO ICPL = 1, NCPL
                     IC = OSL%WALL(IWL)%ICN2(ICPL)
                     MATRIX_COLUMN_LOOP2: DO ICOL = SL%A_ROW(IC), SL%A_ROW(IC+1)-1
                        JC = SL%A_COL(ICOL)
                        IF (JC > SL%NC) THEN
                           IW0 = SL%ICG_TO_IWG(JC)
                           !IF (SL%OSL%WALL(IW0)%NOM /= NOM) CYCLE MATRIX_COLUMN_LOOP2
                           OS%SEND_REAL(LL) = - REAL(SL%ICG_TO_IWG(JC),EB)
                        ELSE
                           OS%SEND_REAL(LL) =   REAL(JC,EB)
                        ENDIF
                        OS%SEND_REAL(LL+1) = SL%A(ICOL)
                        LL = LL + 2
                     ENDDO MATRIX_COLUMN_LOOP2
         
                  ENDDO 
               ENDIF
            ENDDO PACK_MATRIX_SYSTEM
      
            OS%SEND_REAL(LL) = -999.0_EB
      
            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF

      
         !! ---------------------------------------------------------------------------------------
         !! EXCHANGE_MATRIX_PROL: pack overlapping parts of prolongation matrix
         !! ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_PROL)
      
            LL = 1
      
            PACK_MATRIX_PROL: DO IWL=1, OSL%NCG
               IF (OSL%WALL(IWL)%NOM/=NM) CYCLE PACK_MATRIX_PROL
               NCPL = OSL%WALL(IWL)%NCPL
      
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,'(a,i4,a,i4,a,2f12.6)') '========================== IWL=',IWL,': NCPL=',NCPL
      ENDIF
               OS%SEND_REAL(LL)   = REAL(IWL,EB)
               LL = LL + 1
               PROL_FIRST_LAYER_LOOP: DO IC=OSL%WALL(IWL)%ICN(1),OSL%WALL(IWL)%ICN(NCPL)
      
                  OS%SEND_REAL(LL)   = REAL(SL%CELL_TYPE(IC),EB)
                  OS%SEND_REAL(LL+1) = REAL(SL%P_ROW(IC+1)-SL%P_ROW(IC),EB)
      
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,'(a,i4,a,i4,a,2f12.6)') 'IC=',IC,': IC=',IC,': SEND ', OS%SEND_REAL(LL:LL+1)
      ENDIF
                  LL = LL + 2
      
                  DO ICOL = SL%P_ROW(IC), SL%P_ROW(IC+1)-1
                     JC = SL%P_COL(ICOL)
                     IF (JC >= SL%NC) THEN
                        OS%SEND_REAL(LL) = - REAL(JC,EB)
                     ELSE
                        OS%SEND_REAL(LL) =   REAL(JC,EB)
                     ENDIF
                     OS%SEND_REAL(LL+1) = SL%P(ICOL)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,'(a,i4,a,i4,a,2f12.6)') 'IC=',IC,': JC=',JC,': SEND ', OS%SEND_REAL(LL:LL+1)
      ENDIF
                     LL = LL + 2
                  ENDDO 
      
               ENDDO PROL_FIRST_LAYER_LOOP
      
               IF (NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) THEN
         
                  PROL_SECOND_LAYER_LOOP: DO IC=OSL%WALL(IWL)%ICN2(1),OSL%WALL(IWL)%ICN2(NCPL)
         
                     OS%SEND_REAL(LL)   = REAL(SL%CELL_TYPE(IC),EB)
                     OS%SEND_REAL(LL+1) = REAL(SL%P_ROW(IC+1)-SL%P_ROW(IC),EB)
         
                     LL = LL + 2
         
                     DO ICOL = SL%P_ROW(IC), SL%P_ROW(IC+1)-1
                        JC = SL%P_COL(ICOL)
                        IF (JC >= SL%NC) THEN
                           OS%SEND_REAL(LL) = - REAL(JC,EB)
                        ELSE
                           OS%SEND_REAL(LL) =   REAL(JC,EB)
                        ENDIF
                        OS%SEND_REAL(LL+1) = SL%P(ICOL)
                        LL = LL + 2
                     ENDDO 
         
                  ENDDO PROL_SECOND_LAYER_LOOP
      
               ENDIF
      
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
      
            LL = 1
            PACK_MATRIX_REST: DO IWL=1, OSL%NCG
               IF (OSL%WALL(IWL)%NOM/=NM) CYCLE PACK_MATRIX_REST
               NCPL = OSL%WALL(IWL)%NCPL
      
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,'(a,i4,a,i4,a,2f12.6)') '========================== IWL=',IWL,': NCPL=',NCPL
      ENDIF
               OS%SEND_REAL(LL)   = REAL(IWL,EB)
               LL = LL + 1
               REST_FIRST_LAYER_LOOP: DO IC=OSL%WALL(IWL)%ICN(1),OSL%WALL(IWL)%ICN(NCPL)
      
                  ICC = SL%CELL_TYPE(IC)
                  IF (ICC > 0) THEN
                     OS%SEND_REAL(LL)   = REAL(ICC)
                     OS%SEND_REAL(LL+1) = REAL(SL%R_ROW(ICC+1)-SL%R_ROW(ICC),EB)
      
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,'(a,i4,a,i4,a,2f12.6)') 'IC=',IC,': ICC=',ICC,': SEND ', OS%SEND_REAL(LL:LL+1)
      ENDIF
                     LL = LL + 2
      
                     DO ICOL = SL%R_ROW(ICC), SL%R_ROW(ICC+1)-1
                        JCC = SL%R_COL(ICOL)
                        IF (JCC >= SL%NCC) THEN
                           OS%SEND_REAL(LL) = - REAL(JCC,EB)
                        ELSE
                           OS%SEND_REAL(LL) =   REAL(JCC,EB)
                        ENDIF
                        OS%SEND_REAL(LL+1) = SL%R(ICOL)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,'(a,i4,a,i4,a,2f12.6)') 'IC=',IC,': ICC=',ICC,': SEND ', OS%SEND_REAL(LL:LL+1)
      ENDIF
                        LL = LL + 2
                     ENDDO 
      
                  ENDIF
      
               ENDDO REST_FIRST_LAYER_LOOP
      
            ENDDO PACK_MATRIX_REST
            OS%SEND_REAL(LL) = -999.0_EB

            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF


           !< ---------------------------------------------------------------------------------------
           !< EXCHANGE_TRANSER_SIZE: pack sizes of transfer matrices (AMG only)
           !< ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_TRANSFER_SIZE) 

            IF (RNODE /= SNODE) THEN
               SL  =>  S%LEVEL(NL)
               OSL => OS%LEVEL(NL)
               NREQ_SCARC = NREQ_SCARC+1
               OS%SEND_INTEGER(1) = SL%NC
               OS%SEND_INTEGER(2) = SL%NW
               OS%SEND_INTEGER(3) = SL%NCE
               OS%SEND_INTEGER(4) = SL%NW
               OS%SEND_INTEGER(5) = SL%NCCI
               OS%SEND_INTEGER(6) = OSL%NP0
               OS%SEND_INTEGER(7) = OSL%NR0
               OS%SEND_INTEGER(8) = OSL%NCC0
               OS%SEND_INTEGER(9) = OSL%NCF0
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,*) 'SIZE_MATRIX: NM=',NM,': NOM=',NOM,': SENDING ',OS%SEND_INTEGER(1:9)
               CALL MPI_ISEND(OS%SEND_INTEGER(1),9,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
            ENDIF


      END SELECT OMESH_PACK_SELECT
   ENDDO OMESH_PACK_LOOP
ENDDO MESH_PACK_LOOP


!> ------------------------------------------------------------------------------------------------
!> Information from Mesh NM is received by Mesh NOM  (NOM receiver, NM sender)
!> ------------------------------------------------------------------------------------------------
IF (TYPE_DEBUG == NSCARC_DEBUG_INFO2) &
   WRITE(SCARC_LU,*) 'BEFORE WAITALL ', NREQ_SCARC,': TYPE_EXCHANGE=',TYPE_EXCHANGE

IF (USE_MPI.AND.NREQ_SCARC/=0) CALL MPI_WAITALL(NREQ_SCARC,REQ_SCARC(1:NREQ_SCARC),MPI_STATUSES_IGNORE,IERR)

IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
   WRITE(SCARC_LU,*) 'AFTER WAITALL', IERR


!> ------------------------------------------------------------------------------------------------
!> Extract communication data from corresponding RECEIVE-buffers
!> ------------------------------------------------------------------------------------------------
MESH_UNPACK_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   OMESH_UNPACK_LOOP: DO NOM=1,NMESHES
    
      SNODE  = PROCESS(NM)
      RNODE  = PROCESS(NOM)

      OS => SCARC(NM)%OSCARC(NOM)

      OMESH_UNPACK_IF: IF (OS%NICMAX_S/=0 .AND. OS%NICMAX_R/=0) THEN

         SL   => SCARC(NM)%LEVEL(NL)
         OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

         IF (RNODE/=SNODE) THEN
            RECV_INTEGER  => SCARC(NM)%OSCARC(NOM)%RECV_INTEGER
            RECV_REAL     => SCARC(NM)%OSCARC(NOM)%RECV_REAL
         ELSE
            RECV_INTEGER  => SCARC(NOM)%OSCARC(NM)%SEND_INTEGER
            RECV_REAL     => SCARC(NOM)%OSCARC(NM)%SEND_REAL
         ENDIF
            
         OMESH_UNPACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

               !< ------------------------------------------------------------------------------------
              !< EXCHANGE_BASIC: unpack basic neighboring sizes for wallinfo exchange
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_BASIC) 

               IF (RNODE/=SNODE) THEN
                  OSL%NCG  = OS%RECV_INTEGER0(1)
               ELSE
                  OSL%NCG  = SCARC(NOM)%OSCARC(NM)%LEVEL(NL)%NWL
               ENDIF

               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
                  WRITE(SCARC_LU,*) '======== NSCARC_EXCHANGE_BASIC:'
                  WRITE(SCARC_LU,*) 'NM=',NM,': NOM=',NOM,' : OSL%NCG =',OSL%NCG
               ENDIF

              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_WIDTHINFO: unpack neighboring step size information
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_WIDTHINFO) 

               IF (RNODE/=SNODE) THEN
                  OSL%DH = OS%RECV_REAL0(1)
               ELSE
                  OSL%DH = SCARC(NOM)%OSCARC(NM)%LEVEL(NL)%DH
               ENDIF

               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
                 WRITE(SCARC_LU,'(a,i4,a,2f12.6)') 'WIDTH1: FACE(', OSL%IOR,')%DH =',OSL%DH
                 WRITE(SCARC_LU,'(a,10f8.4)') 'WIDTH: DXL:',(SL%DXL(I), I=0,SL%NX)
                 WRITE(SCARC_LU,'(a,10f8.4)') 'WIDTH: DYL:',(SL%DYL(J), J=0,SL%NY)
                 WRITE(SCARC_LU,'(a,10f8.4)') 'WIDTH: DZL:',(SL%DZL(K), K=0,SL%NZ)
               ENDIF

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

               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
                 WRITE(SCARC_LU,'(a,i4,a,2f12.6)') 'WIDTH2: FACE(', OSL%IOR,')%DH =',OSL%DH
                 WRITE(SCARC_LU,'(a,10f8.4)') 'WIDTH: DXL:',(SL%DXL(I), I=0,SL%NX)
                 WRITE(SCARC_LU,'(a,10f8.4)') 'WIDTH: DYL:',(SL%DYL(J), J=0,SL%NY)
                 WRITE(SCARC_LU,'(a,10f8.4)') 'WIDTH: DZL:',(SL%DZL(K), K=0,SL%NZ)
               ENDIF


              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_WALLINFO: unpack neighboring wall information
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_WALLINFO)
               IF (RNODE/=SNODE) THEN
                  IPTR=1
                  DO ICG = 1, OSL%NCG
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(a,i4,a,i4,a,i4)') 'RECEIVING ICG=',ICG,': NCPLR=',OSL%NCPLR
                     OSL%WALL(ICG)%IXG    = OS%RECV_INTEGER(IPTR    )
                     OSL%WALL(ICG)%IYG    = OS%RECV_INTEGER(IPTR + 1)
                     OSL%WALL(ICG)%IZG    = OS%RECV_INTEGER(IPTR + 2)
                     OSL%WALL(ICG)%IXW    = OS%RECV_INTEGER(IPTR + 3)
                     OSL%WALL(ICG)%IYW    = OS%RECV_INTEGER(IPTR + 4)
                     OSL%WALL(ICG)%IZW    = OS%RECV_INTEGER(IPTR + 5)
                     OSL%WALL(ICG)%IXN(1) = OS%RECV_INTEGER(IPTR + 6)
                     OSL%WALL(ICG)%IXN(2) = OS%RECV_INTEGER(IPTR + 7)
                     OSL%WALL(ICG)%IYN(1) = OS%RECV_INTEGER(IPTR + 8)
                     OSL%WALL(ICG)%IYN(2) = OS%RECV_INTEGER(IPTR + 9)
                     OSL%WALL(ICG)%IZN(1) = OS%RECV_INTEGER(IPTR +10)
                     OSL%WALL(ICG)%IZN(2) = OS%RECV_INTEGER(IPTR +11)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(i4,i4)') (I, OS%RECV_INTEGER(I), I=IPTR,IPTR+11)
                     IPTR = IPTR + 12
                     ALLOCATE (OSL%WALL(ICG)%ICE(OSL%NCPLR))
                     DO ICPL=1,OSL%NCPLR
                        OSL%WALL(ICG)%ICE(ICPL) = OS%RECV_INTEGER(IPTR)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(i4,i4)') (I, OS%RECV_INTEGER(I), I=IPTR,IPTR)
                        IPTR = IPTR + 1
                     ENDDO
                     ALLOCATE (OSL%WALL(ICG)%ICN(OSL%NCPLR))
                     DO ICPL=1,OSL%NCPLR
                        OSL%WALL(ICG)%ICN(ICPL) = OS%RECV_INTEGER(IPTR)
IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
   WRITE(SCARC_LU,'(i4,i4)') (I, OS%RECV_INTEGER(I), I=IPTR,IPTR)
                        IPTR = IPTR + 1
                     ENDDO
                  ENDDO
               ELSE
                  SLO => SCARC(NOM)%LEVEL(NL)
                  DO ICG = 1, OSL%NCG
                     IWG = SCARC(NOM)%OSCARC(NM)%LEVEL(NL)%IWL_TO_IWG(ICG)
                     OSL%WALL(ICG)%IXG    = SLO%WALL(IWG)%IXG
                     OSL%WALL(ICG)%IYG    = SLO%WALL(IWG)%IYG
                     OSL%WALL(ICG)%IZG    = SLO%WALL(IWG)%IZG
                     OSL%WALL(ICG)%IXW    = SLO%WALL(IWG)%IXW
                     OSL%WALL(ICG)%IYW    = SLO%WALL(IWG)%IYW
                     OSL%WALL(ICG)%IZW    = SLO%WALL(IWG)%IZW
                     OSL%WALL(ICG)%IXN(1) = SLO%WALL(IWG)%IXN(1)
                     OSL%WALL(ICG)%IXN(2) = SLO%WALL(IWG)%IXN(2)
                     OSL%WALL(ICG)%IYN(1) = SLO%WALL(IWG)%IYN(1)
                     OSL%WALL(ICG)%IYN(2) = SLO%WALL(IWG)%IYN(2)
                     OSL%WALL(ICG)%IZN(1) = SLO%WALL(IWG)%IZN(1)
                     OSL%WALL(ICG)%IZN(2) = SLO%WALL(IWG)%IZN(2)
                     ALLOCATE (OSL%WALL(ICG)%ICE(OSL%NCPLR))
                     DO ICPL=1,OSL%NCPLR
                        OSL%WALL(ICG)%ICE(ICPL) = SLO%WALL(IWG)%ICE(ICPL)
                     ENDDO
                     ALLOCATE (OSL%WALL(ICG)%ICN(OSL%NCPLR))
                     DO ICPL=1,OSL%NCPLR
                        OSL%WALL(ICG)%ICN(ICPL) = SLO%WALL(IWG)%ICN(ICPL)
                     ENDDO
                  ENDDO
               ENDIF


              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_WALLINFO: unpack neighboring mesh information
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MESHINFO) 
               IF (RNODE /= SNODE) THEN
                  OSL%NX = OS%RECV_INTEGER(1)
                  OSL%NY = OS%RECV_INTEGER(2)
                  OSL%NZ = OS%RECV_INTEGER(3)
                  OSL%NC = OS%RECV_INTEGER(4)
                  OSL%NW = OS%RECV_INTEGER(5)
               ELSE
                  SLO => SCARC(NOM)%LEVEL(NL)
                  OSL%NX =  SLO%NX
                  OSL%NY =  SLO%NY
                  OSL%NZ =  SLO%NZ
                  OSL%NC =  SLO%NC
                  OSL%NW =  SLO%NW
               ENDIF


              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_WALLINFO: unpack neighboring matrix size information
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SIZE) 
               IF (RNODE /= SNODE) THEN
                  OSL%NA = OS%RECV_INTEGER(1)
               ELSE
                  OSL%NA = SCARC(NOM)%OSCARC(NM)%LEVEL(NL)%NA0
               ENDIF


              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_WALLINFO: unpack neighboring transfer size information
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_TRANSFER_SIZE) 
               IF (RNODE /= SNODE) THEN
                  OSL%NC  =  OS%RECV_INTEGER(1)
                  OSL%NW  =  OS%RECV_INTEGER(2)
                  OSL%NCE =  OS%RECV_INTEGER(3)
                  OSL%NW  =  OS%RECV_INTEGER(4)
                  OSL%NCCI=  OS%RECV_INTEGER(5)
                  OSL%NP  =  OS%RECV_INTEGER(6)
                  OSL%NR  =  OS%RECV_INTEGER(7)
                  OSL%NCC =  OS%RECV_INTEGER(8)
                  OSL%NCF =  OS%RECV_INTEGER(9)
               ELSE
                  SL  => SCARC(NOM)%LEVEL(NL)
                  SLO => SCARC(NOM)%OSCARC(NM)%LEVEL(NL)
                  OSL%NC  =  SL%NC
                  OSL%NW  =  SL%NW
                  OSL%NCE =  SL%NCE
                  OSL%NCCI=  SL%NCCI
                  OSL%NP  =  SLO%NP0
                  OSL%NR  =  SLO%NR0
                  OSL%NCC =  SLO%NCC0
                  OSL%NCF =  SLO%NCF0
               ENDIF

              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_PRESSURE: unpack overlapping parts of a H or HS
            ! (pressure vectors for PREDICTOR and CORRECTOR
              !< ------------------------------------------------------------------------------------
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
!WRITE(SCARC_LU,'(a,3i4,a,f12.6)') 'UNPACKING H(', I, J, K,')=',HVECTOR(I, J, K)
               
               ENDDO UNPACK_PRESSURE
            
              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_VECTOR: unpack overlapping parts of a given vector
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_VECTOR) 

               VECTOR => POINT_TO_VECTOR (TYPE_VECTOR, NM, NL)
            
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,'(a,16f12.6)') 'RECV_REAL:',(RECV_REAL(ICG),ICG=1,OSL%NCG)

               LL = 1
               UNPACK_VECTOR: DO IWL = 1, OSL%NWL
                  ZSUM = 0.0_EB
                  IWG = OSL%IWL_TO_IWG(IWL)
                  ICO = OSL%IWL_TO_ICO(IWL)
                  ICO = SL%WALL(IWG)%ICO
                  DO ICPL = 1, OSL%NCPL
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'IWL=',IWL,': ICO=',ICO,': ICPL=',ICPL,': ICG=',ICG
                     ZSUM = ZSUM + RECV_REAL(LL)
                     LL = LL + 1
                  ENDDO
                  VECTOR(ICO) = ZSUM/REAL(OSL%NCPL,EB)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,'(a,i4,a,i4,a,f12.6)') &
                     'IWL=',IWL,':  UNPACK VECTOR(',ICO,')=',VECTOR(ICO)
               ENDDO UNPACK_VECTOR
               
            
              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_MEASURE_ADD: unpack neighboring measure information (AMG only) with adding
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MEASURE_ADD) 
            
               LL = 1
               UNPACK_MEASURE_ADD: DO
               
                  IW = NINT(RECV_REAL(LL))
                  IF (IW==-999) EXIT UNPACK_MEASURE_ADD
                  LL = LL + 1
               
                  ICW = SL%WALL(IW)%ICW
                  SL%MEASURE(ICW) = SL%MEASURE(ICW)  + RECV_REAL(LL)
                  IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i4,a,3f12.6)') &
                     'UNPACKA:  SL%MEASURE(',ICW,')=', SL%MEASURE(ICW),  RECV_REAL(LL)
               
                  LL = LL + 1
               
               ENDDO UNPACK_MEASURE_ADD
            
              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_CELL_TYPE: unpack neighboring CELL_TYPE information
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_CELL_TYPE) 
            
            LL = 0
            UNPACK_CELL_TYPE: DO
            
               IW = NINT(RECV_REAL(LL+1))
               IF (IW==-999) EXIT UNPACK_CELL_TYPE
            
               ZSUM1=0.0_EB
               NCPL = SL%WALL(IW)%NCPL
               DO ICPL = 1, NCPL
                  IC1 = SL%WALL(IW)%ICN(ICPL)
                  ZSUM1=ZSUM1+RECV_REAL(LL+2)
                  LL = LL+2
               ENDDO
            
               IC1 = SL%WALL(IW)%ICG(1)
               OSL%CELL_TYPE(IC1) = NINT(ZSUM1/REAL(SL%WALL(IW)%NCPL,EB))
            
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'RECEIVING OSL%CELL_TYPE(',IC1,')=',OSL%CELL_TYPE(IC1)
            
               IC1 = SL%WALL(IW)%ICE(1)
               OSL%CELL_TYPE(IC1) = NINT(ZSUM1/REAL(SL%WALL(IW)%NCPL,EB))
            
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'RECEIVING  OSL%CELL_TYPE(',IC1,')=', OSL%CELL_TYPE(IC1)
            
            ENDDO UNPACK_CELL_TYPE
            
            
              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_MATRIX_SUBDIAG: unpack subdiagonal information from !neighboring system matrix
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG) 
            
               WRITE(*,*) 'EXCHANGE_MATRIX_SUBDIAG, UNPACK: Achtung, hier nochmal checken für AMG!!'
               IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
                  WRITE(SCARC_LU,*) 'NSCARC_EXCHANGE_MATRIX_SUBDIAG: RECV_REAL'
                  WRITE(SCARC_LU,'(8f12.6)') (RECV_REAL(IC), IC=1, SIZE(RECV_REAL))
               ENDIF

               LL = 0
               UNPACK_MATRIX_SUBDIAG: DO
               
                  IW = NINT(RECV_REAL(LL+1))
                  IF (IW==-999) EXIT UNPACK_MATRIX_SUBDIAG
                  ZSUM=0.0_EB
                  NCPL = SL%WALL(IW)%NCPL

         
                  LL = LL + 1
                  DO IC=SL%WALL(IW)%ICN(1),SL%WALL(IW)%ICN(NCPL)
                     ZSUM=ZSUM+RECV_REAL(LL+1)
                     LL = LL+1
                  ENDDO
               
                  IC=SL%WALL(IW)%ICW
               
                  IROW = SL%A_ROW(IC)
                  DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                     IF (IW == -SL%A_COL(ICOL)) THEN
                        SL%A(ICOL)     = ZSUM/REAL(SL%WALL(IW)%NCPL,EB)
                        SL%A_COL(ICOL) = SL%WALL(IW)%ICE(1)
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)') &
                 'UNPACK_MATRIX_SUBDIAG1: NM=',NM,': NOM=',NOM,': A_COL(',ICOL,')=', SL%A_COL(ICOL),' A(',ICOL,')=',SL%A(ICOL)
                     ENDIF
                  ENDDO
               
                  IF (NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) THEN
         
                     IW = NINT(RECV_REAL(LL+1))
                     IF (IW==-999) EXIT UNPACK_MATRIX_SUBDIAG
                     ZSUM=0.0_EB
                     NCPL = SL%WALL(IW)%NCPL
                     LL = LL + 1
                     DO IC=SL%WALL(IW)%ICN2(1),SL%WALL(IW)%ICN2(NCPL)
                        ZSUM=ZSUM+RECV_REAL(LL+1)
                        LL = LL+1
                     ENDDO
                  
                     IC=SL%WALL(IW)%ICW
                  
                     IROW = SL%A_ROW(IC)
                     DO ICOL = SL%A_ROW(IC)+1, SL%A_ROW(IC+1)-1
                        IF (IW == -SL%A_COL(ICOL)) THEN
                           SL%A(ICOL)     = ZSUM/REAL(SL%WALL(IW)%NCPL,EB)
                           SL%A_COL(ICOL) = SL%WALL(IW)%ICE(1)
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)')& 
                 'UNPACK_MATRIX_SUBDIAG1: NM=',NM,': NOM=',NOM,': A_COL(',ICOL,')=', SL%A_COL(ICOL),' A(',ICOL,')=',SL%A(ICOL)
                        ENDIF
                     ENDDO
                  
                  ENDIF
         
               ENDDO UNPACK_MATRIX_SUBDIAG
            
         
              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_MATRIX_STENCIL: unpack stencil information from neighboring system matrix
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_STENCIL) 
            
               UNPACK_MATRIX_STENCIL: DO
            
                  IW = NINT(RECV_REAL(LL))
                  IF (IW==-888.OR.IW==-999) EXIT UNPACK_MATRIX_STENCIL
            
                  !< ------------------------  First layer  -----------------------------------------
                  ICG = SL%WALL(IW)%ICG(1)
                  OSL%A_SIZE(ICG) = NINT(RECV_REAL(LL+1))
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') &
                'UNPACK_MATRIX_STENCIL1: NM=',NM,': NOM=',NOM,': A_SIZE(',ICG,')=',OSL%A_SIZE(ICG)
         
                  LL = LL + 2
         
                  IF (NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) THEN
                      IW = NINT(RECV_REAL(LL))
                      IF (IW==-999) EXIT UNPACK_MATRIX_STENCIL
            
                      !< Unpack second cell layer if requested
                      ICG = SL%WALL(IW)%ICG2(1)
                      OSL%A_SIZE(ICG) = NINT(RECV_REAL(LL+1))
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') &
                'UNPACK_MATRIX_STENCIL2: NM=',NM,': NOM=',NOM,': A_SIZE(',ICG,')=',OSL%A_SIZE(ICG)
                      LL = LL + 2
                   ENDIF
         
                ENDDO UNPACK_MATRIX_STENCIL
         
            
            
              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_MATRIX_SYSTEM: extract neighboring system matrix information
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SYSTEM) 
            
               LL = 1
               UNPACK_MATRIX_SYSTEM: DO
            
                  IW = NINT(RECV_REAL(LL))
                  IF (IW==-888.OR.IW==-999) THEN
                     OSL%NLEN_MATRIX_SYSTEM = LL
                     IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) &
                           'NLEN_MATRIX_SYSTEM(',NOM,')=',OSL%NLEN_MATRIX_SYSTEM, SIZE(RECV_REAL)
                     EXIT UNPACK_MATRIX_SYSTEM
                  ENDIF
                  LL = LL + 1
            
                  ICE = SL%WALL(IW)%ICE(1)
                  ICN = SL%WALL(IW)%ICN(1)
                  ICG = SL%WALL(IW)%ICG(1)
                  ICW = SL%WALL(IW)%ICW
         
                  DO ICOL = SL%A_ROW(ICE), SL%A_ROW(ICE+1)-1
                     JC0 = NINT(RECV_REAL(LL))
                     !IF (JC0 > 0) THEN
                     !<   JC = - JC0
                     !ELSE IF (JC0 < 0) THEN            !< CAUTION: correct ?
                     !<   JC = SL%WALL(IW)%ICW
                     !ENDIF
                     !SL%A_COL(ICOL)= JC
                     IF (JC0 > 0) THEN
                        SL%A_COL(ICOL)= - JC0
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(7(a,i3),a,f12.6)') &
              'A: UNPACK_MATRIX_SYSTEM1: NM=',NM,': NOM=',NOM,': JC0=',JC0,': ICE=',ICE,&
               ': A_COL(',ICOL,')=',SL%A_COL(ICOL),' A(',ICOL,')=',SL%A(ICOL)
                     ELSE
                        !SL%A_COL(ICOL)= ABS(JC0)       !< CAUTION: still not sure which one is correct
                        SL%A_COL(ICOL)= ICW
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(7(a,i3),a,f12.6)') &
              'B: UNPACK_MATRIX_SYSTEM1: NM=',NM,': NOM=',NOM,': JC0=',JC0,': ICE=',ICE,&
               ': A_COL(',ICOL,')=',SL%A_COL(ICOL),' A(',ICOL,')=',SL%A(ICOL)
                     ENDIF
                     SL%A(ICOL)    = RECV_REAL(LL+1)
                     LL = LL + 2
                  ENDDO 
         
                   !< Unpack second cell layer if requested
                  IF (NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) THEN
                     ICE = SL%WALL(IW)%ICE2(1)
                     DO ICOL = SL%A_ROW(ICE), SL%A_ROW(ICE+1)-1
                        SL%A_COL(ICOL)= - NINT(RECV_REAL(LL))
                        SL%A(ICOL)    = RECV_REAL(LL+1)
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.6)') &
              'UNPACK_MATRIX_SYSTEM2: NM=',NM,': NOM=',NOM,': ICE=',ICE,': A_COL(',ICOL,')=',JC,' A(',ICOL,')=',SL%A(ICOL)
                        LL = LL + 2
                     ENDDO 
                   ENDIF
                ENDDO UNPACK_MATRIX_SYSTEM
         
         
              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_MATRIX_PROL: Extract neighboring prolongation matrix information
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_PROL) 
            
               LL = 1
               NROW=1
               ICC = 1
               UNPACK_MATRIX_PROL: DO
            
                  IW = NINT(RECV_REAL(LL))
                  LL = LL + 1
                  IF (IW==-888.OR.IW==-999) EXIT UNPACK_MATRIX_PROL
            
                  ICG = SL%WALL(IW)%ICG(1)
                  ICE = SL%WALL(IW)%ICE(1)
                  ICN = SL%WALL(IW)%ICN(1)
                  ICELL_TYPE = NINT(RECV_REAL(LL))
                  NCOL      = NINT(RECV_REAL(LL+1))
                  
         IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
            WRITE(SCARC_LU,'(a,i4,a,f12.6,a,i4,a,f12.6)') &
                'RECV_REAL(',LL,')=',RECV_REAL(LL), ': RECV_REAL(', LL+1,')=', RECV_REAL(LL+1)
            WRITE(SCARC_LU,'(a,i4,a,i4,a,i4,a,i4,a,i4)') &
                '========================== IW=',IW,': NCOL=',NCOL,': ICG=',ICG,': ICE=',ICE,': ICN=',ICN
         ENDIF
                  OSL%CELL_TYPE(ICG) = ICELL_TYPE
                  OSL%P_ROW(ICG) = NROW
         IF (TYPE_DEBUG> NSCARC_DEBUG_NONE) THEN
            WRITE(SCARC_LU,*) 'UNPACK PROLONGATION: OSL%CELL_TYPE(',ICG,')=',OSL%CELL_TYPE(ICG)
            WRITE(SCARC_LU,*) 'UNPACK PROLONGATION: OSL%P_ROW(',ICG,')   =',OSL%P_ROW(ICG)
         ENDIF
         
                  LL = LL + 2
         
                  DO ICOL = NROW, NROW+NCOL-1
                     OSL%P_COL(ICOL)= NINT(RECV_REAL(LL))
                     OSL%P(ICOL)    = RECV_REAL(LL+1)
         IF (TYPE_DEBUG> NSCARC_DEBUG_NONE) THEN
            WRITE(SCARC_LU,*) 'UNPACK PROLONGATION: OSL%P_COL(',ICOL,')   =',OSL%P_COL(ICOL)
            WRITE(SCARC_LU,*) 'UNPACK PROLONGATION: OSL%P(',ICOL,')   =',OSL%P(ICOL)
         ENDIF
                     LL = LL + 2
                  ENDDO 
                  NROW = NROW + NCOL
                  OSL%P_ROW(ICG+1) = NROW
         
               ENDDO UNPACK_MATRIX_PROL
            

              !< ------------------------------------------------------------------------------------
              !< EXCHANGE_MATRIX_REST: exchange neighboring restriction matrix information
              !< ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_REST) 
            
               LL = 1
               NROW=1
               ICC = 1
               UNPACK_MATRIX_REST: DO
            
                  IW = NINT(RECV_REAL(LL))
                  LL = LL + 1
                  IF (IW==-888.OR.IW==-999) EXIT UNPACK_MATRIX_REST
            
                  !ICG = SL%WALL(IW)%ICG(1)
                  !ICE = SL%WALL(IW)%ICE(1)
                  !ICN = SL%WALL(IW)%ICN(1)
                  !ICELL_TYPE = NINT(RECV_REAL(LL))
                  !NCOL      = NINT(RECV_REAL(LL+1))
                 !< 
                 !< OSL%CELL_TYPE(ICG) = ICELL_TYPE
                 !< OSL%P_ROW(ICG) = NROW
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
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, SC1, SC2, SC3, SC4
TYPE (SCARC_LEVEL_TYPE), POINTER :: OSL, OSC1, OSC2, OSC12, OSC13, OSC21, OSC24, OSC31, OSC34, OSC42, OSC43

IF (TYPE_DEBUG < NSCARC_DEBUG_LESS) RETURN

SELECT CASE (NTYPE)

 !! ------------------------------------------------------------------------------------------------
 !! Debug system matrix A (corresponding to system type)
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIX)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(SCARC_LU,*) 'NA =',SL%NA
         WRITE(SCARC_LU,*) 'NC =',SL%NC
         WRITE(SCARC_LU,*) 'SIZE(A) =',SIZE(SL%A)
         WRITE(SCARC_LU,*) 'SIZE(A_COL) =',SIZE(SL%A_COL)
         WRITE(SCARC_LU,*) 'SIZE(A_ROW) =',SIZE(SL%A_ROW)
         WRITE(SCARC_LU,*) '---------------------- A_ROW:', SL%NC
         WRITE(SCARC_LU,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NC+1)
         WRITE(SCARC_LU,*) '---------------------- A_COL:'
         DO IC = 1, SL%NC
            WRITE(SCARC_LU,'(i5,a,20i9)') IC,':',(SL%A_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
         ENDDO
         WRITE(SCARC_LU,*) '---------------------- A:', SL%A_ROW(IC+1)- SL%A_ROW(IC)
         DO IC = 1, SL%NC
            WRITE(SCARC_LU,'(i5,a,20f9.2)') IC,':',(SL%A(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
         ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug extended system matrix A (corresponding to system type)
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIXE)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(SCARC_LU,*) 'PRINTING EXTENDED MATRIX, NCE=',SL%NCE
         WRITE(SCARC_LU,*) '---------------------- A_ROW:', SL%NCE
         WRITE(SCARC_LU,'(4i9)') (SL%A_ROW(IC), IC=1,SL%NCE+1)
         WRITE(SCARC_LU,*) '---------------------- A_COL:'
         DO IC = 1, SL%NCE
            WRITE(SCARC_LU,'(i5,a,20i9)') IC,':',(SL%A_COL(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
         ENDDO
         WRITE(SCARC_LU,*) '---------------------- A:'
         DO IC = 1, SL%NCE
            WRITE(SCARC_LU,'(i5,a,20f9.2)') IC,':',(SL%A(IP),IP=SL%A_ROW(IC),SL%A_ROW(IC+1)-1)
         ENDDO
         WRITE(SCARC_LU,*) 'SIZE(A) =',SIZE(SL%A)
         WRITE(SCARC_LU,*) 'SIZE(A_COL) =',SIZE(SL%A_COL)
         WRITE(SCARC_LU,*) 'SIZE(A_ROW) =',SIZE(SL%A_ROW)
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug prolongation matrix P (only for compact system)
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_PROLONGATION)

      DO NM = NMESHES_MIN, NMESHES_MAX
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         WRITE(SCARC_LU,*) '============= P_ROW:', NM
         WRITE(SCARC_LU,'(4i9)') (SL%P_ROW(IC), IC=1, SL%NCE+1)
         WRITE(SCARC_LU,*) '============= P_COL:', NM, SL%NCE
         DO IC=1,SL%NCE
            WRITE(SCARC_LU,'(i5,a,10i9)') IC,':',(SL%P_COL(IP), IP=SL%P_ROW(IC),SL%P_ROW(IC+1)-1)
         ENDDO
         WRITE(SCARC_LU,*) '============= P:', NM, SL%NCE
         DO IC=1,SL%NCE
            WRITE(SCARC_LU,'(i5,a,10f9.2)') IC,':',(SL%P(IP), IP=SL%P_ROW(IC),SL%P_ROW(IC+1)-1)
         ENDDO
         DO NOM = 1, NMESHES
            IF (NOM == NM) CYCLE 
            OS  => SCARC(NM)%OSCARC(NOM)      
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE 
            OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)      
            WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P_ROW:', OSL%NCG
            WRITE(SCARC_LU,'(4i9)') (OSL%P_ROW(IC), IC=1, OSL%NCG+1)
            WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P_COL:', OSL%NCG
            DO IC=1,OSL%NCG
               WRITE(SCARC_LU,'(i5,a,10i9)') IC,':',(OSL%P_COL(IP), IP=OSL%P_ROW(IC),OSL%P_ROW(IC+1)-1)
            ENDDO
            WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P:'
            DO IC=1,OSL%NCG
               WRITE(SCARC_LU,'(i5,a,10f9.2)') IC,':',(OSL%P(IP), IP=OSL%P_ROW(IC),OSL%P_ROW(IC+1)-1)
            ENDDO
         ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug restriction matrix R (only for compact system)
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_RESTRICTION)

      DO NM = NMESHES_MIN, NMESHES_MAX
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         WRITE(SCARC_LU,*) '============= R_ROW:', NM, SL%NCCE
         WRITE(SCARC_LU,'(4i9)') (SL%R_ROW(IC), IC=1, SL%NCCE)
         WRITE(SCARC_LU,*) '============= R_COL:', NM
         DO IC=1,SL%NCCE
            WRITE(SCARC_LU,'(i5,a,10i9)') IC,':',(SL%R_COL(IP), IP=SL%R_ROW(IC),SL%R_ROW(IC+1)-1)
         ENDDO
         WRITE(SCARC_LU,*) '============= R:', NM
         DO IC=1,SL%NCCE
            WRITE(SCARC_LU,'(i5,a,10f9.2)') IC,':',(SL%R(IP), IP=SL%R_ROW(IC),SL%R_ROW(IC+1)-1)
         ENDDO
         !DO NOM = 1, NMESHES
         DO NOM = 1, NMESHES
            IF (NOM == NM) CYCLE 
            OS  => SCARC(NM)%OSCARC(NOM)      
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE 
            OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)      
            WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%R_ROW:', OSL%NCC
            WRITE(SCARC_LU,'(4i9)') (OSL%R_ROW(IC), IC=1, OSL%NCC+1)
            WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%R_COL:', OSL%NCC
            DO IC=1,OSL%NCC
               WRITE(SCARC_LU,'(i5,a,10i9)') IC,':',(OSL%R_COL(IP), IP=OSL%R_ROW(IC),OSL%R_ROW(IC+1)-1)
            ENDDO
            WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P:', OSL%NCC
            DO IC=1,OSL%NCC
               WRITE(SCARC_LU,'(i5,a,10f9.2)') IC,':',(OSL%R(IP), IP=OSL%R_ROW(IC),OSL%R_ROW(IC+1)-1)
            ENDDO
         ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug FACEINFO
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_FACEINFO)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         M  => MESHES(NM)
         DO IOR0 = -3, 3
            IF (IOR0 == 0) CYCLE
            WRITE(SCARC_LU,*) '========================================='
            WRITE(SCARC_LU,*) '============= DEBUGGING FACE(',IOR0,'):'
            WRITE(SCARC_LU,*) '========================================='
            WRITE(SCARC_LU,*) 'FACE(.)%NFM:', SL%FACE(IOR0)%NFM
            WRITE(SCARC_LU,*) 'FACE(.)%NFX:', SL%FACE(IOR0)%NFX
            WRITE(SCARC_LU,*) 'FACE(.)%NFY:', SL%FACE(IOR0)%NFY
            WRITE(SCARC_LU,*) 'FACE(.)%NFZ:', SL%FACE(IOR0)%NFZ
            WRITE(SCARC_LU,*) 'FACE(.)%NWL:', SL%FACE(IOR0)%NFW
            WRITE(SCARC_LU,*) 'FACE(.)%IWG_MARKER:', SL%FACE(IOR0)%IWG_MARKER
            WRITE(SCARC_LU,*) '----------------------------------------------'
            WRITE(SCARC_LU,*) 'FACE(.)%BXY:', SL%FACE(IOR0)%BXY
            WRITE(SCARC_LU,*) 'FACE(.)%BXZ:', SL%FACE(IOR0)%BXZ
            WRITE(SCARC_LU,*) 'FACE(.)%BYZ:', SL%FACE(IOR0)%BYZ
            WRITE(SCARC_LU,*)
            !WRITE(SCARC_LU,*) 'FACE(.)%DH :', SL%FACE(IOR0)%DH 
            WRITE(SCARC_LU,*) 'FACE(.)%IOFFSET_MATRI :', SL%FACE(IOR0)%IOFFSET_MATRIX 
            WRITE(SCARC_LU,*)
            IF (SL%FACE(IOR0)%NUM_NEIGHBORS /= 0) THEN
               DO INBR=1,SL%FACE(IOR0)%NUM_NEIGHBORS
                  NOM = SL%FACE(IOR0)%NEIGHBORS(INBR)
                  OSL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  WRITE(SCARC_LU,'(a,i4,a,2f12.6)') '---OSL(',NOM,')%DH :',OSL%DH
                  WRITE(SCARC_LU,'(a,i4,a,2i4)') '---OSL(',NOM,')%IOR   :',OSL%IOR
                  WRITE(SCARC_LU,'(a,i4,a,2i4)') '---OSL(',NOM,')%NWL(.):',OSL%NWL
                  WRITE(SCARC_LU,'(a,i4,a,2i4)') '---OSL(',NOM,')%NCG(.):',OSL%NCG
                  WRITE(SCARC_LU,*)
                  WRITE(SCARC_LU,'(a,i4,a)') '------OSL(',NOM,')%ICG_TO_IWG:'
                  DO IW = 1, OSL%NCG
                     WRITE(SCARC_LU,'(16i4)') OSL%ICG_TO_IWG(IW)
                  ENDDO
                  WRITE(SCARC_LU,*)
                  WRITE(SCARC_LU,'(a,i4,a)') '------OSL(',NOM,')%ICG_TO_ICW:'
                  DO IW = 1, OSL%NCG
                     WRITE(SCARC_LU,'(16i4)') OSL%ICG_TO_ICW(IW)
                  ENDDO
                  WRITE(SCARC_LU,*)
                  WRITE(SCARC_LU,'(a,i4,a)') '------OSL(',NOM,')%ICG_TO_ICE:'
                  DO IW = 1, OSL%NCG
                     WRITE(SCARC_LU,'(16i4)') OSL%ICG_TO_ICE(IW)
                  ENDDO
                  WRITE(SCARC_LU,*)
                  WRITE(SCARC_LU,'(a,i4,a)') '------OSL(',NOM,')%ICG_TO_ICN:'
                  DO IW = 1, OSL%NCG
                     WRITE(SCARC_LU,'(16i4)') OSL%ICG_TO_ICN(IW)
                  ENDDO
                  WRITE(SCARC_LU,*)
                  WRITE(SCARC_LU,'(a,i4,a)') '------OSL(',NOM,')%IWL_TO_IWG:'
                  DO IW = 1, OSL%NWL
                     WRITE(SCARC_LU,'(16i4)') OSL%IWL_TO_IWG(IW)
                  ENDDO
                  WRITE(SCARC_LU,*)
                  WRITE(SCARC_LU,'(a,i4,a)') '------OSL(',NOM,')%IWL_TO_ICW:'
                  DO IW = 1, OSL%NWL
                     WRITE(SCARC_LU,'(16i4)') OSL%IWL_TO_ICW(IW)
                  ENDDO
                  WRITE(SCARC_LU,*)
                  WRITE(SCARC_LU,'(a,i4,a)') '------OSL(',NOM,')%IWL_TO_ICO:'
                  DO IW = 1, OSL%NWL
                     WRITE(SCARC_LU,'(16i4)') OSL%IWL_TO_ICO(IW)
                  ENDDO
                  WRITE(SCARC_LU,*)
                  WRITE(SCARC_LU,'(a,i4,a)') '------OSL(',NOM,')%IWL_TO_ICN:'
                  DO IW = 1, OSL%NWL
                     WRITE(SCARC_LU,'(16i4)') OSL%IWL_TO_ICN(IW, 1:OSL%NCPL)
                  ENDDO
                  WRITE(SCARC_LU,*)
                  WRITE(SCARC_LU,'(a,i4,a)') '------OSL(',NOM,')%IWL_TO_ICG:'
                  DO IW = 1, OSL%NCG
                     WRITE(SCARC_LU,'(16i4)') OSL%IWL_TO_ICG(IW, 1:OSL%NCPL)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug WALLINFO
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_WALLINFO)

!WRITE(SCARC_LU,*) ': NL=',NL, ': CNAME=',CNAME,': CROUTINE=',CROUTINE
      !IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC) THEN
      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         M  => MESHES(NM)

            WRITE(SCARC_LU,*) 'NCE =',SL%NCE
            WRITE(SCARC_LU,*) 'NG =',SL%NCG
            WRITE(SCARC_LU,*)
            WRITE(SCARC_LU,'(a,i4,a)') '------SL%ICE_TO_IWG:'
            DO IW = SL%NC+1, SL%NCE
               WRITE(SCARC_LU,'(16i4)') SL%ICE_TO_IWG(IW)
            ENDDO
            WRITE(SCARC_LU,'(a,i4,a)') '------SL%ICE_TO_IWL:'
            DO IW = SL%NC+1, SL%NCE
               WRITE(SCARC_LU,'(16i4)') SL%ICE_TO_IWL(IW)
            ENDDO
            WRITE(SCARC_LU,*)
            WRITE(SCARC_LU,'(a,i4,a)') '------SL%ICE_TO_ICG:'
            DO IW = SL%NC+1, SL%NCE
               WRITE(SCARC_LU,'(16i4)') SL%ICE_TO_ICG(IW)
            ENDDO
            WRITE(SCARC_LU,*)
            WRITE(SCARC_LU,'(a,i4,a)') '------SL%ICE_TO_ICW:'
            DO IW = SL%NC+1, SL%NCE
               WRITE(SCARC_LU,'(16i4)') SL%ICE_TO_ICW(IW)
            ENDDO
         IF (NL == 1) THEN
         WRITE(SCARC_LU,*) '============= WALL(.)%IXG:', NM
         WRITE(SCARC_LU,'(16i5)') (SL%WALL(IW)%IXG, IW=1,SL%NW)
         WRITE(SCARC_LU,*) '============= WALL(.)%IYG:', NM
         WRITE(SCARC_LU,'(16i5)') (SL%WALL(IW)%IYG, IW=1,SL%NW)
         WRITE(SCARC_LU,*) '============= WALL(.)%IZG:', NM
         WRITE(SCARC_LU,'(16i5)') (SL%WALL(IW)%IZG, IW=1,SL%NW)
         WRITE(SCARC_LU,*)
         WRITE(SCARC_LU,*) '============= WALL(.)%IXW:', NM
         WRITE(SCARC_LU,'(16i5)') (SL%WALL(IW)%IXW, IW=1,SL%NW)
         WRITE(SCARC_LU,*) '============= WALL(.)%IYW:', NM
         WRITE(SCARC_LU,'(16i5)') (SL%WALL(IW)%IYW, IW=1,SL%NW)
         WRITE(SCARC_LU,*) '============= WALL(.)%IZW:', NM
         WRITE(SCARC_LU,'(16i5)') (SL%WALL(IW)%IZW, IW=1,SL%NW)
         ENDIF
         WRITE(SCARC_LU,*)
         WRITE(SCARC_LU,*) '============= WALL(.)%IOR:', NM
         WRITE(SCARC_LU,'(16i5)') (SL%WALL(IW)%IOR, IW=1,SL%NW)
         WRITE(SCARC_LU,*)
         WRITE(SCARC_LU,*) '============= WALL(.)%NOM:', NM
         WRITE(SCARC_LU,'(16i5)') (SL%WALL(IW)%NOM, IW=1,SL%NW)
         WRITE(SCARC_LU,*)
         IF (NL == 1) THEN
         WRITE(SCARC_LU,*) '============= WALL(.)%NCPL:', NM
         WRITE(SCARC_LU,'(16i5)') (SL%WALL(IW)%NCPL, IW=1,SL%NW)
         ENDIF
         WRITE(SCARC_LU,*)
         WRITE(SCARC_LU,*) '============= WALL(.)%ICW:', NM
         DO IW=1,SL%NW
            WRITE(SCARC_LU,'(a,i3, a,16i5)') 'IW=',IW,':',SL%WALL(IW)%ICW
         ENDDO
         WRITE(SCARC_LU,*)
         WRITE(SCARC_LU,*) '============= WALL(.)%ICO:', NM
         DO IW=1,SL%NW
            IF (SL%WALL(IW)%NOM /=0) &
               WRITE(SCARC_LU,'(a,i3, a,16i5)') 'IW=',IW,':',SL%WALL(IW)%ICO
         ENDDO
         WRITE(SCARC_LU,*)
         WRITE(SCARC_LU,*) '============= WALL(.)%ICE:', NM
         DO IW=1,SL%NW
            IF (SL%WALL(IW)%NOM /=0) &
               WRITE(SCARC_LU,'(a,i3, a,16i5)') 'IW=',IW,':',SL%WALL(IW)%ICE
         ENDDO
         WRITE(SCARC_LU,*)
         WRITE(SCARC_LU,*) '============= WALL(.)%ICN:', NM
         DO IW=1,SL%NW
            IF (SL%WALL(IW)%NOM /=0) &
               WRITE(SCARC_LU,'(a,i3, a,16i5)') 'IW=',IW,':',SL%WALL(IW)%ICN
         ENDDO
         WRITE(SCARC_LU,*)
         WRITE(SCARC_LU,*) '============= WALL(.)%ICG:', NM
         DO IW=1,SL%NW
            IF (SL%WALL(IW)%NOM /=0) &
               WRITE(SCARC_LU,'(a,i3, a,16i5)') 'IW=',IW,':',SL%WALL(IW)%ICG
         ENDDO
         WRITE(SCARC_LU,*) '====================================================='
         WRITE(SCARC_LU,*) ' Plotting out M%WALL-structure'
         WRITE(SCARC_LU,*) '====================================================='
         WRITE(SCARC_LU,*) ' M%WALL(.)%WALL_INDEX'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%WALL_INDEX, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%SURF_INDEX'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%SURF_INDEX, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%BACK_INDEX'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%BACK_INDEX, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%BOUNDARY_TYPE'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%BOUNDARY_TYPE, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%OBST_INDEX'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%OBST_INDEX, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%PRESSURE_BC_INDEX'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%PRESSURE_BC_INDEX, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%PRESSURE_ZONE'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%PRESSURE_ZONE, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%VENT_INDEX'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%VENT_INDEX, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%NOM'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%NOM, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%NOM_IB(1)'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%NOM_IB(1), IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%NOM_IB(2)'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%NOM_IB(2), IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%NOM_IB(3)'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%NOM_IB(3), IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%NOM_IB(4)'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%NOM_IB(4), IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%NOM_IB(5)'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%NOM_IB(5), IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%NOM_IB(6)'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%NOM_IB(6), IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%ONE_D%II'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%ONE_D%II, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%ONE_D%JJ'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%ONE_D%JJ, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%ONE_D%KK'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%ONE_D%KK, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%ONE_D%IIG'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%ONE_D%IIG, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%ONE_D%JJG'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%ONE_D%JJG, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%ONE_D%KKG'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%ONE_D%KKG, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%ONE_D%N_LAYER_CELLS'
         WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%ONE_D%IOR, IW=1,SL%NW)
         WRITE(SCARC_LU,*) ' M%WALL(.)%ONE_D%IOR'
         !WRITE(SCARC_LU,'(16i5)') (M%WALL(IW)%ONE_D%N_LAYER_CELLS, IW=1,SL%NW)



      ENDDO
      !ENDIF


 !! ------------------------------------------------------------------------------------------------
 !! Debug PRESSURE_BC_INDEX
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_BCINDEX)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SL => SCARC(NM)%LEVEL(NL)
         WRITE(SCARC_LU, '(8i4)') (SL%WALL(J)%BTYPE, J=1,SL%NW)
      ENDDO
      
 !! ------------------------------------------------------------------------------------------------
 !! Debug WALL_CELL
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_ACELL)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(SCARC_LU, '(8i8)') (SL%WALL(IW)%ICW, IW=1,SL%NW)
      ENDDO
      
 !! ------------------------------------------------------------------------------------------------
 !! Debug GHOST_CELL
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_GCELL)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(SCARC_LU,*) '-------- GHOST,  NGE=',SL%NCG
         DO IW=1,SL%NW
            !IF (SL%WALL(IW)%NOM/=0) &
            !<   WRITE(SCARC_LU, '(8i8)') (SL%WALL(IW)%ICG(I), I=1,SL%WALL(IW)%NCPL)
         ENDDO
         !DO IW=1,SL%NW
         !<   IF (SL%WALL(IW)%NOM/=0.AND.NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) &
         !<      WRITE(SCARC_LU, '(8i8)') (SL%WALL(IW)%ICG2(I), I=1,SL%WALL(IW)%NCPL)
         !ENDDO
      ENDDO

 !! ------------------------------------------------------------------------------------------------
 !! Debug NOM_CELL
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_NCELL)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(SCARC_LU,*) '-------- NOM, NCE= ', SL%NCE
         DO IW=1,SL%NW
            !IF (SL%WALL(IW)%NOM/=0) &
            !<   WRITE(SCARC_LU, '(8i8)') (SL%WALL(IW)%ICN(I), I=1,SL%WALL(IW)%NCPL)
         ENDDO
         !DO IW=1,SL%NW
         !<   IF (SL%WALL(IW)%NOM/=0.AND.NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) &
         !<      WRITE(SCARC_LU, '(8i8)') (SL%WALL(IW)%ICN2(I), I=1,SL%WALL(IW)%NCPL)
         !ENDDO
      ENDDO
      
 !! ------------------------------------------------------------------------------------------------
 !! Debug SUBDIVISION
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_SUBDIVISION)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SL  => SCARC(NM)%LEVEL(NL)
         WRITE(SCARC_LU,*) 'SUBDIVISION IOR= 1 '
         WRITE(SCARC_LU,'(3i4)') (SL%SUBDIVISION(I, 1), I=1,3)
         WRITE(SCARC_LU,*) 'SUBDIVISION IOR=-1 '
         WRITE(SCARC_LU,'(3i4)') (SL%SUBDIVISION(I,-1), I=1,3)
         WRITE(SCARC_LU,*) 'SUBDIVISION IOR= 2 '
         WRITE(SCARC_LU,'(3i4)') (SL%SUBDIVISION(I, 2), I=1,3)
         WRITE(SCARC_LU,*) 'SUBDIVISION IOR=-2 '
         WRITE(SCARC_LU,'(3i4)') (SL%SUBDIVISION(I,-2), I=1,3)
         WRITE(SCARC_LU,*) 'SUBDIVISION IOR= 3 '
         WRITE(SCARC_LU,'(3i4)') (SL%SUBDIVISION(I, 3), I=1,3)
         WRITE(SCARC_LU,*) 'SUBDIVISION IOR=-3 '
         WRITE(SCARC_LU,'(3i4)') (SL%SUBDIVISION(I,-3), I=1,3)
      ENDDO
   
 !! ------------------------------------------------------------------------------------------------
 !! Debug MEASURE
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MEASURE)

      WRITE(SCARC_LU,1000) CROUTINE, CNAME, 1, NL
      IF (NMESHES == 4 .OR. NMESHES==1) THEN                !< only temporarily
      DO NM = NMESHES_MIN, NMESHES_MAX
         SL  => SCARC(NM)%LEVEL(NL)
         IF (NL == 1) THEN
            DO K = SL%NZ, 1, -1
               WRITE(SCARC_LU, '(10f11.5)') (SL%MEASURE((K-1)*SL%NX*SL%NY+I), I=1, SL%NX)
            ENDDO
            WRITE(SCARC_LU,*) '--------------------'
            WRITE(SCARC_LU, '(4f11.5)') (SL%MEASURE(IC), IC=SL%NC+1, SL%NCE)
         ELSE
            WRITE(SCARC_LU, '(4f11.5)') (SL%MEASURE(IC), IC=1, SL%NC)
         ENDIF
      ENDDO
      ELSE IF (NMESHES==17.AND.SCARC(1)%LEVEL(NL)%NX==4) THEN
      SC1 => SCARC(1)%LEVEL(NL)
      SC2 => SCARC(2)%LEVEL(NL)
      SC3 => SCARC(3)%LEVEL(NL)
      SC4 => SCARC(4)%LEVEL(NL)
      OSC12 => SCARC(1)%OSCARC(2)%LEVEL(NL)
      OSC13 => SCARC(1)%OSCARC(3)%LEVEL(NL)
      OSC21 => SCARC(2)%OSCARC(1)%LEVEL(NL)
      OSC24 => SCARC(2)%OSCARC(4)%LEVEL(NL)
      OSC31 => SCARC(3)%OSCARC(1)%LEVEL(NL)
      OSC34 => SCARC(3)%OSCARC(4)%LEVEL(NL)
      OSC42 => SCARC(4)%OSCARC(2)%LEVEL(NL)
      OSC43 => SCARC(4)%OSCARC(3)%LEVEL(NL)
      WRITE(SCARC_LU,'(5f6.2,6X,5f6.2)') SC3%MEASURE(13:16),OSC34%MEASURE(7), OSC43%MEASURE(7),SC4%MEASURE(13:16)
      WRITE(SCARC_LU,'(5f6.2,6X,5f6.2)') SC3%MEASURE( 9:12),OSC34%MEASURE(5), OSC43%MEASURE(5),SC4%MEASURE( 9:12)
      WRITE(SCARC_LU,'(5f6.2,6X,5f6.2)') SC3%MEASURE( 5: 8),OSC34%MEASURE(3), OSC43%MEASURE(3),SC4%MEASURE( 5: 8)
      WRITE(SCARC_LU,'(5f6.2,6X,5f6.2)') SC3%MEASURE( 1: 4),OSC34%MEASURE(1), OSC43%MEASURE(1),SC4%MEASURE( 1: 4)
      WRITE(SCARC_LU,'(4F6.2,18X,4F6.2)') OSC31%MEASURE( 1: 7:2),OSC42%MEASURE(1:7:2)
      WRITE(SCARC_LU,*)
      WRITE(SCARC_LU,'(4F6.2,18X,4F6.2)') OSC13%MEASURE( 1: 7:2),OSC24%MEASURE(1:7:2)
      WRITE(SCARC_LU,'(5f6.2,6X,5f6.2)') SC1%MEASURE(13:16),OSC12%MEASURE(7), OSC21%MEASURE(7),SC2%MEASURE(13:16)
      WRITE(SCARC_LU,'(5f6.2,6X,5f6.2)') SC1%MEASURE( 9:12),OSC12%MEASURE(5), OSC21%MEASURE(5),SC2%MEASURE( 9:12)
      WRITE(SCARC_LU,'(5f6.2,6X,5f6.2)') SC1%MEASURE( 5: 8),OSC12%MEASURE(3), OSC21%MEASURE(3),SC2%MEASURE( 5: 8)
      WRITE(SCARC_LU,'(5f6.2,6X,5f6.2)') SC1%MEASURE( 1: 4),OSC12%MEASURE(1), OSC21%MEASURE(1),SC2%MEASURE( 1: 4)
      ELSE IF (NMESHES==4) THEN
      SC1 => SCARC(1)%LEVEL(NL)
      SC2 => SCARC(2)%LEVEL(NL)
      SC3 => SCARC(3)%LEVEL(NL)
      SC4 => SCARC(4)%LEVEL(NL)
      OSC12 => SCARC(1)%OSCARC(2)%LEVEL(NL)
      OSC13 => SCARC(1)%OSCARC(3)%LEVEL(NL)
      OSC21 => SCARC(2)%OSCARC(1)%LEVEL(NL)
      OSC24 => SCARC(2)%OSCARC(4)%LEVEL(NL)
      OSC31 => SCARC(3)%OSCARC(1)%LEVEL(NL)
      OSC34 => SCARC(3)%OSCARC(4)%LEVEL(NL)
      OSC42 => SCARC(4)%OSCARC(2)%LEVEL(NL)
      OSC43 => SCARC(4)%OSCARC(3)%LEVEL(NL)
      WRITE(SCARC_LU,'(6f6.2,6X,6f6.2)') SC3%MEASURE(16:20),OSC34%MEASURE(7), OSC43%MEASURE(7),SC4%MEASURE(13:16)
      WRITE(SCARC_LU,'(6f6.2,6X,6f6.2)') SC3%MEASURE(11:15),OSC34%MEASURE(5), OSC43%MEASURE(5),SC4%MEASURE(9:12)
      WRITE(SCARC_LU,'(6f6.2,6X,6f6.2)') SC3%MEASURE( 6:10),OSC34%MEASURE(3), OSC43%MEASURE(3),SC4%MEASURE(5:8)
      WRITE(SCARC_LU,'(6f6.2,6X,6f6.2)') SC3%MEASURE( 1: 5),OSC34%MEASURE(1), OSC43%MEASURE(1),SC4%MEASURE(1: 4)
      WRITE(SCARC_LU,'(5F6.2,18X,5F6.2)') OSC31%MEASURE( 1: 9:2),OSC42%MEASURE(1:7:2)
      WRITE(SCARC_LU,*)
      WRITE(SCARC_LU,'(5F6.2,18X,5F6.2)') OSC13%MEASURE( 1: 9:2),OSC24%MEASURE(1:7:2)
      WRITE(SCARC_LU,'(6f6.2,6X,6f6.2)') SC1%MEASURE(16:20),OSC12%MEASURE(7), OSC21%MEASURE(7),SC2%MEASURE(13:16)
      WRITE(SCARC_LU,'(6f6.2,6X,6f6.2)') SC1%MEASURE(11:15),OSC12%MEASURE(5), OSC21%MEASURE(5),SC2%MEASURE(9:12)
      WRITE(SCARC_LU,'(6f6.2,6X,6f6.2)') SC1%MEASURE( 6:10),OSC12%MEASURE(3), OSC21%MEASURE(3),SC2%MEASURE( 5:8)
      WRITE(SCARC_LU,'(6f6.2,6X,6f6.2)') SC1%MEASURE( 1: 5),OSC12%MEASURE(1), OSC21%MEASURE(1),SC2%MEASURE( 1: 4)
      ELSE IF (NMESHES==2 .AND. SCARC(1)%LEVEL(1)%NX == 5) THEN
      SC1 => SCARC(1)%LEVEL(NL)
      SC2 => SCARC(2)%LEVEL(NL)
      OSC1 => SCARC(1)%OSCARC(2)%LEVEL(NL)
      OSC2 => SCARC(2)%OSCARC(1)%LEVEL(NL)
      WRITE(SCARC_LU,'(5F6.2,a,1F6.2,6X,1F6.2,a,4F6.2)') &
                               SC1%MEASURE(16:20),'|',OSC1%MEASURE(4), OSC2%MEASURE(4),'|',SC2%MEASURE(13:16)
      WRITE(SCARC_LU,'(5F6.2,a,1F6.2,6X,1F6.2,a,4F6.2)') &
                               SC1%MEASURE(11:15),'|',OSC1%MEASURE(3), OSC2%MEASURE(3),'|',SC2%MEASURE(9:12)
      WRITE(SCARC_LU,'(5F6.2,a,1F6.2,6X,1F6.2,a,4F6.2)') &
                               SC1%MEASURE(6:10) ,'|',OSC1%MEASURE(2), OSC2%MEASURE(2),'|',SC2%MEASURE(5:8)
      WRITE(SCARC_LU,'(5F6.2,a,1F6.2,6X,1F6.2,a,4F6.2)') &
                               SC1%MEASURE(1:5)  ,'|',OSC1%MEASURE(1), OSC2%MEASURE(1),'|',SC2%MEASURE(1:4)
      ELSE IF (SCARC(1)%LEVEL(1)%NX == 8) THEN
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(57:64),SC1%MEASURE(72),SC2%MEASURE(72),SC2%MEASURE(57:64)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(49:56),SC1%MEASURE(71),SC2%MEASURE(71),SC2%MEASURE(49:56)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(41:48),SC1%MEASURE(70),SC2%MEASURE(70),SC2%MEASURE(41:48)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(33:40),SC1%MEASURE(69),SC2%MEASURE(69),SC2%MEASURE(33:40)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(25:32),SC1%MEASURE(68),SC2%MEASURE(68),SC2%MEASURE(25:32)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(17:24),SC1%MEASURE(67),SC2%MEASURE(67),SC2%MEASURE(17:24)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(9:16) ,SC1%MEASURE(66),SC2%MEASURE(66),SC2%MEASURE(9:16)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(1:8)  ,SC1%MEASURE(65),SC2%MEASURE(65),SC2%MEASURE(1:8)
      ELSE IF (SCARC(1)%LEVEL(1)%NX == 9) THEN
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(64:72),SC1%MEASURE(80),SC2%MEASURE(72),SC2%MEASURE(57:64)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(55:63),SC1%MEASURE(79),SC2%MEASURE(71),SC2%MEASURE(49:56)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(46:54),SC1%MEASURE(78),SC2%MEASURE(70),SC2%MEASURE(41:48)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(37:45),SC1%MEASURE(77),SC2%MEASURE(69),SC2%MEASURE(33:40)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(28:36),SC1%MEASURE(76),SC2%MEASURE(68),SC2%MEASURE(25:32)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(19:27),SC1%MEASURE(75),SC2%MEASURE(67),SC2%MEASURE(17:24)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(10:18),SC1%MEASURE(74),SC2%MEASURE(66),SC2%MEASURE(9:16)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(1:9)  ,SC1%MEASURE(73),SC2%MEASURE(65),SC2%MEASURE(1:8)
      ENDIF

 !! ------------------------------------------------------------------------------------------------
 !! Debug CELL_TYPE
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_CELL_TYPE)

      IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
       IF (NMESHES == 4 .OR. NMESHES==1) THEN                !< only temporarily
         DO NM = NMESHES_MIN, NMESHES_MAX
           WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
           SL  => SCARC(NM)%LEVEL(NL)
           IF (NL == 1) THEN
              DO K = SL%NZ, 1, -1
                 WRITE(SCARC_LU, '(8i6)') (SL%CELL_TYPE((K-1)*SL%NX*SL%NY+I), I=1, SL%NX)
              ENDDO
              WRITE(SCARC_LU,*) '--------------------'
              WRITE(SCARC_LU, '(4i6)') (SL%CELL_TYPE(IC), IC=SL%NC+1, SL%NCE)
           ELSE
              WRITE(SCARC_LU, '(4i6)') (SL%CELL_TYPE(IC), IC=1, SL%NC)
           ENDIF
         ENDDO
       ELSEIF (NMESHES == 4.AND.NL==1) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(SCARC_LU,'(10i5)') SL%CELL_TYPE(1:SL%NC)
         SL  => SCARC(2)%LEVEL(NL)
         WRITE(SCARC_LU,'(10i5)') SL%CELL_TYPE(1:SL%NC)
         SL  => SCARC(3)%LEVEL(NL)
         WRITE(SCARC_LU,'(9i5)') SL%CELL_TYPE(1:SL%NC)
         SL  => SCARC(4)%LEVEL(NL)
         WRITE(SCARC_LU,'(9i5)') SL%CELL_TYPE(1:SL%NC)
       ELSEIF (NMESHES == 4.AND.NL==2) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(SCARC_LU,'(5i5)') SL%CELL_TYPE(1:SL%NC)
         SL  => SCARC(2)%LEVEL(NL)
         WRITE(SCARC_LU,'(5i5)') SL%CELL_TYPE(1:SL%NC)
         SL  => SCARC(3)%LEVEL(NL)
         WRITE(SCARC_LU,'(4i5)') SL%CELL_TYPE(1:SL%NC)
         SL  => SCARC(4)%LEVEL(NL)
         WRITE(SCARC_LU,'(4i5)') SL%CELL_TYPE(1:SL%NC)
       ELSEIF (NMESHES==1.AND.NL==1) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(SCARC_LU,'(9i5)') SL%CELL_TYPE(1:SL%NC)
       ELSEIF (NMESHES==1.AND.NL==2) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(SCARC_LU,'(3i5)') SL%CELL_TYPE(1:SL%NC)
       ENDIF
      ENDIF
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
      WRITE(SCARC_LU,1000) CROUTINE, CNAME, 1, NL
      IF (NMESHES == 1.OR.NL>100) THEN                !< only temporarily
      DO NM = NMESHES_MIN, NMESHES_MAX
         SL  => SCARC(NM)%LEVEL(NL)
         IF (NL == 1) THEN
            DO K = SL%NZ, 1, -1
               WRITE(SCARC_LU, '(10i5)') (SL%CELL_TYPE((K-1)*SL%NX*SL%NY+I), I=1, SL%NX)
            ENDDO
         ELSE
            WRITE(SCARC_LU, '(4i5)') (SL%CELL_TYPE(IC), IC=1, SL%NC)
         ENDIF
      ENDDO
      ELSE IF (NMESHES==4.AND.SCARC(1)%LEVEL(NL)%NC==4) THEN
      SC1 => SCARC(1)%LEVEL(NL)
      SC2 => SCARC(2)%LEVEL(NL)
      SC3 => SCARC(3)%LEVEL(NL)
      SC4 => SCARC(4)%LEVEL(NL)
      OSC12 => SCARC(1)%OSCARC(2)%LEVEL(NL)
      OSC13 => SCARC(1)%OSCARC(3)%LEVEL(NL)
      OSC21 => SCARC(2)%OSCARC(1)%LEVEL(NL)
      OSC24 => SCARC(2)%OSCARC(4)%LEVEL(NL)
      OSC31 => SCARC(3)%OSCARC(1)%LEVEL(NL)
      OSC34 => SCARC(3)%OSCARC(4)%LEVEL(NL)
      OSC42 => SCARC(4)%OSCARC(2)%LEVEL(NL)
      OSC43 => SCARC(4)%OSCARC(3)%LEVEL(NL)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC3%CELL_TYPE(13:16),OSC34%CELL_TYPE(7), OSC43%CELL_TYPE(7),SC4%CELL_TYPE(13:16)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC3%CELL_TYPE( 9:12),OSC34%CELL_TYPE(5), OSC43%CELL_TYPE(5),SC4%CELL_TYPE(9:12)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC3%CELL_TYPE( 5:8),OSC34%CELL_TYPE(3), OSC43%CELL_TYPE(3),SC4%CELL_TYPE( 5:8)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC3%CELL_TYPE( 1: 4),OSC34%CELL_TYPE(1), OSC43%CELL_TYPE(1),SC4%CELL_TYPE( 1: 4)
      WRITE(SCARC_LU,'(4I6,18X,4I6)') OSC31%CELL_TYPE( 1: 7:2),OSC42%CELL_TYPE(1:7:2)
      WRITE(SCARC_LU,*)
      WRITE(SCARC_LU,'(4I6,18X,4I6)') OSC13%CELL_TYPE( 1: 7:2),OSC24%CELL_TYPE(1:7:2)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC1%CELL_TYPE(13:16),OSC12%CELL_TYPE(7), OSC21%CELL_TYPE(7),SC2%CELL_TYPE(13:16)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC1%CELL_TYPE(9:12),OSC12%CELL_TYPE(5), OSC21%CELL_TYPE(5),SC2%CELL_TYPE(9:12)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC1%CELL_TYPE( 5:8),OSC12%CELL_TYPE(3), OSC21%CELL_TYPE(3),SC2%CELL_TYPE(5:8)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC1%CELL_TYPE( 1: 4),OSC12%CELL_TYPE(1), OSC21%CELL_TYPE(1),SC2%CELL_TYPE(1: 4)
      ELSE IF (NMESHES==4) THEN
      SC1 => SCARC(1)%LEVEL(NL)
      SC2 => SCARC(2)%LEVEL(NL)
      SC3 => SCARC(3)%LEVEL(NL)
      SC4 => SCARC(4)%LEVEL(NL)
      OSC12 => SCARC(1)%OSCARC(2)%LEVEL(NL)
      OSC13 => SCARC(1)%OSCARC(3)%LEVEL(NL)
      OSC21 => SCARC(2)%OSCARC(1)%LEVEL(NL)
      OSC24 => SCARC(2)%OSCARC(4)%LEVEL(NL)
      OSC31 => SCARC(3)%OSCARC(1)%LEVEL(NL)
      OSC34 => SCARC(3)%OSCARC(4)%LEVEL(NL)
      OSC42 => SCARC(4)%OSCARC(2)%LEVEL(NL)
      OSC43 => SCARC(4)%OSCARC(3)%LEVEL(NL)
      WRITE(SCARC_LU,'(6I6,6X,6I6)') SC3%CELL_TYPE(16:20),OSC34%CELL_TYPE(7), OSC43%CELL_TYPE(7),SC4%CELL_TYPE(13:16)
      WRITE(SCARC_LU,'(6I6,6X,6I6)') SC3%CELL_TYPE(11:15),OSC34%CELL_TYPE(5), OSC43%CELL_TYPE(5),SC4%CELL_TYPE(9:12)
      WRITE(SCARC_LU,'(6I6,6X,6I6)') SC3%CELL_TYPE( 6:10),OSC34%CELL_TYPE(3), OSC43%CELL_TYPE(3),SC4%CELL_TYPE( 5:8)
      WRITE(SCARC_LU,'(6I6,6X,6I6)') SC3%CELL_TYPE( 1: 5),OSC34%CELL_TYPE(1), OSC43%CELL_TYPE(1),SC4%CELL_TYPE( 1: 4)
      WRITE(SCARC_LU,'(5I6,18X,5I6)') OSC31%CELL_TYPE( 1: 9:2),OSC42%CELL_TYPE(1:7:2)
      WRITE(SCARC_LU,*)
      WRITE(SCARC_LU,'(5I6,18X,5I6)') OSC13%CELL_TYPE( 1: 9:2),OSC24%CELL_TYPE(1:7:2)
      WRITE(SCARC_LU,'(6I6,6X,6I6)') SC1%CELL_TYPE(16:20),OSC12%CELL_TYPE(7), OSC21%CELL_TYPE(7),SC2%CELL_TYPE(13:16)
      WRITE(SCARC_LU,'(6I6,6X,6I6)') SC1%CELL_TYPE(11:15),OSC12%CELL_TYPE(5), OSC21%CELL_TYPE(5),SC2%CELL_TYPE(9:12)
      WRITE(SCARC_LU,'(6I6,6X,6I6)') SC1%CELL_TYPE( 6:10),OSC12%CELL_TYPE(3), OSC21%CELL_TYPE(3),SC2%CELL_TYPE(5:8)
      WRITE(SCARC_LU,'(6I6,6X,6I6)') SC1%CELL_TYPE( 1: 5),OSC12%CELL_TYPE(1), OSC21%CELL_TYPE(1),SC2%CELL_TYPE(1: 4)
      ELSE IF (NMESHES==2 .AND. SCARC(1)%LEVEL(1)%NX == 5) THEN
      SC1 => SCARC(1)%LEVEL(NL)
      SC2 => SCARC(2)%LEVEL(NL)
      OSC1 => SCARC(1)%OSCARC(2)%LEVEL(NL)
      OSC2 => SCARC(2)%OSCARC(1)%LEVEL(NL)
      WRITE(SCARC_LU,'(7I6,6X,6I6)') SC1%CELL_TYPE(16:20),OSC1%CELL_TYPE(4),OSC1%CELL_TYPE(8), &
                                     OSC2%CELL_TYPE(8),OSC2%CELL_TYPE(4),SC2%CELL_TYPE(13:16)
      WRITE(SCARC_LU,'(7I6,6X,6I6)') SC1%CELL_TYPE(11:15),OSC1%CELL_TYPE(3),OSC1%CELL_TYPE(7), &
                                     OSC2%CELL_TYPE(7),OSC2%CELL_TYPE(3),SC2%CELL_TYPE(9:12)
      WRITE(SCARC_LU,'(7I6,6X,6I6)') SC1%CELL_TYPE(6:10) ,OSC1%CELL_TYPE(2),OSC1%CELL_TYPE(6), &
                                     OSC2%CELL_TYPE(6),OSC2%CELL_TYPE(2),SC2%CELL_TYPE(5:8)
      WRITE(SCARC_LU,'(7I6,6X,6I6)') SC1%CELL_TYPE(1:5)  ,OSC1%CELL_TYPE(1),OSC1%CELL_TYPE(5), &
                                     OSC2%CELL_TYPE(5),OSC2%CELL_TYPE(1),SC2%CELL_TYPE(1:4)
      ELSE IF (SCARC(1)%LEVEL(1)%NX == 8) THEN
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELL_TYPE(57:64),SC1%CELL_TYPE(72),SC2%CELL_TYPE(72),SC2%CELL_TYPE(57:64)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELL_TYPE(49:56),SC1%CELL_TYPE(71),SC2%CELL_TYPE(71),SC2%CELL_TYPE(49:56)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELL_TYPE(41:48),SC1%CELL_TYPE(70),SC2%CELL_TYPE(70),SC2%CELL_TYPE(41:48)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELL_TYPE(33:40),SC1%CELL_TYPE(69),SC2%CELL_TYPE(69),SC2%CELL_TYPE(33:40)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELL_TYPE(25:32),SC1%CELL_TYPE(68),SC2%CELL_TYPE(68),SC2%CELL_TYPE(25:32)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELL_TYPE(17:24),SC1%CELL_TYPE(67),SC2%CELL_TYPE(67),SC2%CELL_TYPE(17:24)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELL_TYPE(9:16) ,SC1%CELL_TYPE(66),SC2%CELL_TYPE(66),SC2%CELL_TYPE(9:16)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELL_TYPE(1:8)  ,SC1%CELL_TYPE(65),SC2%CELL_TYPE(65),SC2%CELL_TYPE(1:8)
      ELSE IF (SCARC(1)%LEVEL(1)%NX == 9) THEN
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELL_TYPE(64:72),SC1%CELL_TYPE(80),SC2%CELL_TYPE(72),SC2%CELL_TYPE(57:64)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELL_TYPE(55:63),SC1%CELL_TYPE(79),SC2%CELL_TYPE(71),SC2%CELL_TYPE(49:56)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELL_TYPE(46:54),SC1%CELL_TYPE(78),SC2%CELL_TYPE(70),SC2%CELL_TYPE(41:48)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELL_TYPE(37:45),SC1%CELL_TYPE(77),SC2%CELL_TYPE(69),SC2%CELL_TYPE(33:40)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELL_TYPE(28:36),SC1%CELL_TYPE(76),SC2%CELL_TYPE(68),SC2%CELL_TYPE(25:32)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELL_TYPE(19:27),SC1%CELL_TYPE(75),SC2%CELL_TYPE(67),SC2%CELL_TYPE(17:24)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELL_TYPE(10:18),SC1%CELL_TYPE(74),SC2%CELL_TYPE(66),SC2%CELL_TYPE(9:16)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELL_TYPE(1:9)  ,SC1%CELL_TYPE(73),SC2%CELL_TYPE(65),SC2%CELL_TYPE(1:8)
      ENDIF
      IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
      DO NM = NMESHES_MIN, NMESHES_MAX
         SL  => SCARC(NM)%LEVEL(NL)
         DO NOM = 1, NMESHES
            IF (NOM == NM) CYCLE 
            OS  => SCARC(NM)%OSCARC(NOM)      
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE 
            OSL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)      
            WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%CELL_TYPE:', OSL%NCG
            DO IC=1,OSL%NCG
               WRITE(SCARC_LU,'(i5,a,I5)') IC,':',OSL%CELL_TYPE(IC)
            ENDDO
         ENDDO
      ENDDO
      ENDIF
      ENDIF

 !! ------------------------------------------------------------------------------------------------
 !! Debug GRAPH
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_GRAPH)

      IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
       IF (NMESHES == 4 .OR. NMESHES==1) THEN                !< only temporarily
         DO NM = NMESHES_MIN, NMESHES_MAX
           WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
           SL  => SCARC(NM)%LEVEL(NL)
           IF (NL == 1) THEN
              DO K = SL%NZ, 1, -1
                 WRITE(SCARC_LU, '(8i6)') (SL%GRAPH((K-1)*SL%NX*SL%NY+I), I=1, SL%NX)
              ENDDO
              WRITE(SCARC_LU,*) '--------------------'
              WRITE(SCARC_LU, '(4i6)') (SL%GRAPH(IC), IC=SL%NC+1, SL%NCE)
           ELSE
              WRITE(SCARC_LU, '(4i6)') (SL%GRAPH(IC), IC=1, SL%NC)
           ENDIF
         ENDDO
       ELSEIF (NMESHES == 4.AND.NL==1) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(SCARC_LU,'(10i5)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(2)%LEVEL(NL)
         WRITE(SCARC_LU,'(10i5)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(3)%LEVEL(NL)
         WRITE(SCARC_LU,'(9i5)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(4)%LEVEL(NL)
         WRITE(SCARC_LU,'(9i5)') SL%GRAPH(1:SL%NC)
       ELSEIF (NMESHES == 4.AND.NL==2) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(SCARC_LU,'(5i5)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(2)%LEVEL(NL)
         WRITE(SCARC_LU,'(5i5)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(3)%LEVEL(NL)
         WRITE(SCARC_LU,'(4i5)') SL%GRAPH(1:SL%NC)
         SL  => SCARC(4)%LEVEL(NL)
         WRITE(SCARC_LU,'(4i5)') SL%GRAPH(1:SL%NC)
       ELSEIF (NMESHES==1.AND.NL==1) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(SCARC_LU,'(9i5)') SL%GRAPH(1:SL%NC)
       ELSEIF (NMESHES==1.AND.NL==2) THEN
         SL  => SCARC(1)%LEVEL(NL)
         WRITE(SCARC_LU,'(3i5)') SL%GRAPH(1:SL%NC)
       ENDIF
      ENDIF
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
      WRITE(SCARC_LU,1000) CROUTINE, CNAME, 1, NL
      IF (NMESHES == 1.OR.NL>100) THEN                !< only temporarily
      DO NM = NMESHES_MIN, NMESHES_MAX
         SL  => SCARC(NM)%LEVEL(NL)
         IF (NL == 1) THEN
            DO K = SL%NZ, 1, -1
               WRITE(SCARC_LU, '(10i5)') (SL%GRAPH((K-1)*SL%NX*SL%NY+I), I=1, SL%NX)
            ENDDO
         ELSE
            WRITE(SCARC_LU, '(4i5)') (SL%GRAPH(IC), IC=1, SL%NC)
         ENDIF
      ENDDO
      ENDIF
      ENDIF

 !! ------------------------------------------------------------------------------------------------
 !! Debug COARSE
 !! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_COARSE)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
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
            WRITE(SCARC_LU,'(12A4)') (LINE(I:I), I=1,SL%NX)
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
 
IF (TYPE_DEBUG < NSCARC_DEBUG_LESS) RETURN

!DO NM = NMESHES_MIN, NMESHES_MAX
DO NM = NMESHES_MIN, NMESHES_MIN
   
   M  => MESHES(NM)
   VC => POINT_TO_VECTOR (NVECTOR, NM, NLEVEL_MIN)
   
   NX8=MIN(8,M%IBAR)
   NY8=MIN(8,M%JBAR)
   NZ8=MIN(8,M%KBAR)
   
   WRITE(SCARC_LU,*) '======================================================================'
   WRITE(SCARC_LU,2000) CROUTINE, CNAME, NM, NLEVEL_MIN
   WRITE(SCARC_LU,*) '======================================================================'
   SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
   
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
         !WRITE(SCARC_LU,*)  '------------------------------------------------',&
         !<                   '---------------------------------------------------'
   
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
   END SELECT SELECT_DIMENSION

ENDDO


2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4) 
END SUBROUTINE SCARC_DEBUG_VECTOR


!> ------------------------------------------------------------------------------------------------
!> Only for debugging reasons: print out vector information on level NL 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_LEVEL (NVECTOR, CROUTINE, CNAME, NL)
INTEGER, INTENT(IN):: NVECTOR, NL
REAL (EB), POINTER, DIMENSION(:)     :: VC
REAL (EB):: VALUES(0:10)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
TYPE (SCARC_LEVEL_TYPE), POINTER :: SL
 
IF (TYPE_DEBUG < NSCARC_DEBUG_LESS) RETURN

         
DO NM = NMESHES_MIN, NMESHES_MAX

   SL => SCARC(NM)%LEVEL(NL)
   VC => POINT_TO_VECTOR (NVECTOR, NM, NL)
   NX8=MIN(4,SL%NX)
   NY8=MIN(4,SL%NY)
   NZ8=MIN(16,SL%NZ)

   WRITE(SCARC_LU,*) '=========================================================='
   WRITE(SCARC_LU,2001) CROUTINE, CNAME, NM, NL, SL%NC, SL%NCE, NX8, NY8, NZ8, NVECTOR
   WRITE(SCARC_LU,*) '=========================================================='
   !IF (NL == NLEVEL_MIN) THEN
         DO KK = NZ8, 1, - 1
            DO JJ = NY8, 1, - 1
               DO II=1,NX8
                  IC = (KK-1)*SL%NX*SL%NY + (JJ-1)*SL%NY + II
                  IF (ABS(VC(IC))<1.0E-14_EB) THEN
                     VALUES(II)=0.0_EB
                  ELSE
                     VALUES(II)=VC(IC)
                  ENDIF
               ENDDO
               WRITE(SCARC_LU, '(10e18.10)') (VALUES(II), II=1, NX8)
            ENDDO
         ENDDO
   !ENDIF

ENDDO

2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4, ': NX, NY, NZ=',3i4,': NVECTOR=',I4)
2001 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4, ': NC, NCE=',2i4, ': NX, NY, NZ=',3i4,': NVECTOR=',I4)
END SUBROUTINE SCARC_DEBUG_LEVEL


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

DO NM = NMESHES_MIN, NMESHES_MAX

   SL => SCARC(NM)%LEVEL(NL)
   VC => POINT_TO_VECTOR (NVECTOR, NM, NL)

   IF (NL == NLEVEL_MIN) THEN
         DO KK = SL%NX, 1, - 1
            DO JJ = SL%NY, 1, - 1
               WRITE(SCARC_LU, '(16e25.16)') (VALUES((KK-1)*SL%NX*SL%NY+(JJ-1)*SL%NY+II), II=1,SL%NX)
            ENDDO
         ENDDO
   ENDIF

ENDDO

END SUBROUTINE SCARC_DUMP_LEVEL

!> ------------------------------------------------------------------------------------------------
!> Produce latex information about coarsening of grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_LATEX_GRID(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NL0, MLATEX, I, K, IC, IC1, IC2, IC3, IC4, IC5, IC6
CHARACTER (60) :: CLATEX


SELECT CASE(NTYPE)

 !! ---------------------------------------------------------------------------------------------
 !! Produce Latex information about grid coarsening with different node sizes for different levels
 !! ---------------------------------------------------------------------------------------------
   CASE(NSCARC_LATEX_STAGGERED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
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
WRITE(SCARC_LU,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL0)%NCC
                   IF (SCARC(NM)%LEVEL(NL0)%NCC>=1) THEN
                   IF (NL==1) CYCLE
                   IC2=SCARC(NM)%LEVEL(NL0+1)%CELL_TYPE(IC1)
                   IF (IC2>0) THEN
                      WRITE(MLATEX,3002) 'L3','L3',I,K,I,K
WRITE(SCARC_LU,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL0+1)%NCC
                      IF (SCARC(NM)%LEVEL(NL0+1)%NCC>=1) THEN
                      IF (NL==2) CYCLE
                      IC3=SCARC(NM)%LEVEL(NL0+2)%CELL_TYPE(IC2)
                      IF (IC3>0) THEN
                         WRITE(MLATEX,3002) 'L4','L4',I,K,I,K
WRITE(SCARC_LU,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL0+2)%NCC
                         IF (SCARC(NM)%LEVEL(NL0+2)%NCC>=1) THEN
                         IF (NL==3) CYCLE
                         IC4=SCARC(NM)%LEVEL(NL0+3)%CELL_TYPE(IC3)
                         IF (IC4>0) THEN
                            WRITE(MLATEX,3002) 'L5','L5',I,K,I,K
WRITE(SCARC_LU,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL0+3)%NCC
                            IF (SCARC(NM)%LEVEL(NL0+3)%NCC>=1) THEN
                            IF (NL==4) CYCLE
                            IC5=SCARC(NM)%LEVEL(NL0+4)%CELL_TYPE(IC4)
                            IF (IC5>0) THEN
                               WRITE(MLATEX,3002) 'L6','L6',I,K,I,K
WRITE(SCARC_LU,*) 'IN LATEX_MATRIX for LEVEL ', NL,': NCC=',SCARC(NM)%LEVEL(NL+4)%NCC
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

      DO NM = NMESHES_MIN, NMESHES_MAX
      
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

      DO NM = NMESHES_MIN, NMESHES_MAX
      
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
!WRITE(SCARC_LU,*) 'IC=',IC,I,K
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
!> Bitte alles noch mal nachprÃƒÂ¼fen --------
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
 !! CD_VA_2D UrsprÂŸnglicher Testfall
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
     !< Traveling vortex mit glatter Verteilung, Gebiet [0,1]x[0,1], Wirbel im Mittelpunkt 0.5,0.5
      WRITE(*,*) '------------------- neuer Testfall --------------------------'
     CENTER_X    = 0.5_EB
          CENTER_Y    = 0.0_EB
       CENTER_Z    = 0.5_EB
     SCALE_RADIUS= 0.4_EB
          RHO_C       = 0.5_EB
          U_C         = 0.0_EB
          V_C         = 0.0_EB

     !< Dichte, zellzentriert
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

     !< u-Geschwindigkeit, flächenzentriert 
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

     !< v-Geschwindigkeit, flächenzentriert 
     DO K=1,KBAR
        DO J=0,JBAR
           DO I=1,IBAR
                   V(I,J,K) = 0.0_EB
           ENDDO
        ENDDO
     ENDDO

     !< w-Geschwindigkeit, flächenzentriert
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
!> Bitte alles noch mal nachprüfen
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
INTEGER:: NM, IERR, I
INTEGER , ALLOCATABLE, DIMENSION(:)   :: COUNTS_SCARC_TIMERS, DISPLS_SCARC_TIMERS
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: BUFFER
CHARACTER(40) :: NAME_SCARC(0:N_TIMERS_SCARC)
REAL(EB) :: TPCNT_SCARC(0:N_TIMERS_SCARC)

!> only temporarily - use routine only in debug mode
IF (TYPE_DEBUG == NSCARC_DEBUG_NONE) RETURN

IERR=0
 
ALLOCATE(COUNTS_SCARC_TIMERS(0:NUMPROCS-1))
ALLOCATE(DISPLS_SCARC_TIMERS(0:NUMPROCS-1))
ALLOCATE(BUFFER(N_TIMERS_SCARC, NMESHES))

COUNTS_SCARC_TIMERS = COUNTS_SCARC*N_TIMERS_SCARC
DISPLS_SCARC_TIMERS = DISPLS_SCARC*N_TIMERS_SCARC

BUFFER = TUSED_SCARC
IF (USE_MPI) CALL MPI_GATHERV(TUSED_SCARC(1,DISPLS_SCARC(MYID)+1),COUNTS_SCARC_TIMERS(MYID),&
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
   
   DO NM=1,NMESHES
      DO I=1,N_TIMERS_SCARC
         TPCNT_SCARC(I) = 100._EB*TUSED_SCARC(I,NM)/TUSED(1,NM)
      ENDDO
      WRITE(LU_OUTPUT,443) NM
      !WRITE(LU_OUTPUT,444) 'MAIN',TUSED(1,NM),100._EB
      WRITE(LU_OUTPUT,440) "Walltime for one solve ",SCARC_WALL_TIME
      WRITE(LU_OUTPUT,444) (NAME_SCARC(I),TUSED_SCARC(I,NM),TPCNT_SCARC(I),I=1,N_TIMERS_SCARC)
   ENDDO
ENDIF

   IF (TYPE_DEBUG /= NSCARC_DEBUG_NONE ) CLOSE(SCARC_LU)

443 FORMAT(//' ScaRC: CPU Time Usage, Mesh ',I3// &
         47X,' CPU (s)        %  '/ &
         7X,' -----------------------------------------------------------------')
440 FORMAT(7X,A40,1F15.6)
444 FORMAT(7X,A40,2F11.2)

END SUBROUTINE SCARC_TIMINGS


!> ------------------------------------------------------------------------------------------------
!> Stops the code gracefully after writing a message
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SHUTDOWN(CMESSAGE)
CHARACTER(*) CMESSAGE
LOGICAL :: EX

WRITE(LU_ERR,'(/A)') TRIM(CMESSAGE)
INQUIRE(FILE=FN_OUTPUT,EXIST=EX)
IF (.NOT.EX) OPEN(LU_OUTPUT,FILE=TRIM(CHID)//'.out',STATUS='REPLACE',FORM='FORMATTED')
WRITE(LU_OUTPUT,'(/A)') TRIM(CMESSAGE)
STOP

END SUBROUTINE SCARC_SHUTDOWN


!> ------------------------------------------------------------------------------------------------
!> Compute current revision number
!> ------------------------------------------------------------------------------------------------
SUBROUTINE GET_REV_SCRC(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') SCRCREV(INDEX(SCRCREV,':')+2:LEN_TRIM(SCRCREV)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') SCRCDATE

END SUBROUTINE GET_REV_SCRC

END MODULE SCRC
 

   !SELECT CASE (ABS(IOR0))
   !   CASE (1)
   !      N1 = 2
   !      N2 = 5 
   !      M1 = 3
   !      M2 = 6
   !   CASE (2)
   !      N1 = 1
   !      N2 = 5 
   !      M1 = 4
   !      M2 = 6
   !   CASE (3)
   !      N1 = 1
   !      N2 = 2 
   !      M1 = 4
   !      M2 = 5
   !END SELECT

!---------------------------------------------------------------------------------------------------
! Find the wall cell index corresponding to cell (II,JJ,KK) 
!---------------------------------------------------------------------------------------------------
!INTEGER FUNCTION SCARC_WALLCELLS_INDEX(IC, NM, NL)
!INTEGER, INTENT(IN) :: IC
!TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL
!
!SL => SCARC(NM)%LEVEL(NL)
!
!SCARC_WALLCELLS_INDEX = (KK-1)*SL%NX*SL%NY + (JJ-1)*SL%NX + II
!
!END FUNCTION SCARC_WALLCELLS_INDEX
!
!
!---------------------------------------------------------------------------------------------------
! Find the wall index corresponding to the -IOR face of cell (II,JJ,KK) 
!---------------------------------------------------------------------------------------------------
!SUBROUTINE SCARC_WALL_INDEX(II, JJ, KK, IOR, IW, NM, NL)
!INTEGER, INTENT(IN) :: NM, NL, IOR, II, JJ, KK
!INTEGER, INTENT(OUT) :: IW
!INTEGER :: IC
!TYPE (SCARC_LEVEL_TYPE) , POINTER :: SL
!
!SL => SCARC(NM)%LEVEL(NL)
!IC  = SL%CELL_INDEX(II,JJ,KK)
!
!IF (M%SOLID(IC)) THEN
!  SELECT CASE(IOR)
!     CASE(-1)
!        IF (II>0)      II = II-1
!     CASE( 1)
!        IF (II<M%IBP1) II = II+1
!     CASE(-2)
!        IF (JJ>0)      JJ = JJ-1
!     CASE( 2)
!        IF (JJ<M%JBP1) JJ = JJ+1
!     CASE(-3)
!        IF (KK>0)      KK = KK-1
!     CASE( 3)
!        IF (KK<M%KBP1) KK = KK+1
!  END SELECT
!ENDIF
!
!IC  = SL%CELL_INDEX(II,JJ,KK)
!IW  = SL%WALL_INDEX(IC,-IOR)
!
!IF (IW<=0) THEN
!  SELECT CASE(IOR)
!     CASE(-1)
!        IF (II>0)      IC = M%CELL_INDEX(II-1,JJ,KK)
!     CASE( 1)
!        IF (II<M%IBP1) IC = M%CELL_INDEX(II+1,JJ,KK)
!     CASE(-2)
!        IF (JJ>0)      IC = M%CELL_INDEX(II,JJ-1,KK)
!     CASE( 2)
!        IF (JJ<M%JBP1) IC = M%CELL_INDEX(II,JJ+1,KK)
!     CASE(-3)
!        IF (KK>0)      IC = M%CELL_INDEX(II,JJ,KK-1)
!     CASE( 3)
!        IF (KK<M%KBP1) IC = M%CELL_INDEX(II,JJ,KK+1)
!  END SELECT
!  IW = SL%WALL_INDEX(IC,-IOR)
!ENDIF
!
!END SUBROUTINE SCARC_WALL_INDEX
 
  
