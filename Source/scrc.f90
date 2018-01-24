#ifdef WITH_MKL
#define __MKL_PARDISO_F90
#define __MKL_CLUSTER_SPARSE_SOLVER_F90
!include "include/mkl_pardiso.f90"
!include "include/mkl_cluster_sparse_solver.f90"
#endif

MODULE SCRC

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
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
PUBLIC SCARC_SAVE_RHS                          !> Save right hand side computed in pres.f90 for later use

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

PUBLIC SCARC_FFT                               !> If .true. first FFT, then ScaRC
PUBLIC SCARC_PRESSURE_ITERATIONS               !> Internal counter for pressure iterations per time step
PUBLIC SCARC_TWOLEVEL                          !> Predefined ScaRC case

PUBLIC SCARC_KRYLOV                            !> Type of Krylov method
PUBLIC SCARC_KRYLOV_ITERATIONS                 !> Maximum number of iterations for Krylov method
PUBLIC SCARC_KRYLOV_ACCURACY                   !> Requested accuracy for Krylov method
PUBLIC SCARC_KRYLOV_INTERPOL                   !> Interpolation type for twolevel methods

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
PUBLIC SCARC_PRECON_OMEGA                      !> Relaxation parameter for preconditioning method

PUBLIC SCARC_COARSE                            !> Coarse grid solver for multigrid method
PUBLIC SCARC_COARSE_ITERATIONS                 !> Maximum number of iterations for coarse grid solver
PUBLIC SCARC_COARSE_ACCURACY                   !> Requested accuracy for coarse grid solver
PUBLIC SCARC_COARSE_MTYPE                      !> Type of direct coarse grid matrix (symmetric/nonsymmetric)

PUBLIC SCARC_CASE                              !> Predefined ScaRC case
PUBLIC SCARC_DEBUG                             !> Debugging parameter
PUBLIC SCARC_LAYER                             !> Number of layers for data exchange
PUBLIC TYPE_DEBUG                              !> Debug file unit

!> ------------------------------------------------------------------------------------------------
!> Miscellaneous declarations 
!> ------------------------------------------------------------------------------------------------
!> General definitions
CHARACTER(40) :: SCARC_METHOD    = 'KRYLOV'                 !> Requested solver method (KRYLOV/MULTIGRID)
CHARACTER(40) :: SCARC_INITIAL   = 'null'                   !> Initial solution (currently only default is used)

!> General iteration parameters
REAL (EB)     :: SCARC_RESIDUAL             =  0.0_EB       !> Residual of global selected solver
INTEGER       :: SCARC_ITERATIONS           =  0            !> Number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_CAPPA                =  0.0_EB       !> Convergence rate of selected ScarC solver
REAL (EB)     :: SCARC_ACCURACY_DIVERGENCE  =  1.E+6_EB     !> Divergence epsilon for all solvers
REAL (EB)     :: SCARC_ACCURACY_RELATIVE    =  1.E+0_EB     !> Minimum relative accuracy for all solvers
CHARACTER(40) :: SCARC_ACCURACY             = 'ABSOLUTE'    !> Accuracy type (ABSOLUTE/RELATIVE)

!> Parameter for DISCRET
CHARACTER(40) :: SCARC_DISCRETIZATION   = 'STRUCTURED'      !> Discretization of whole domain or with OBSTs 

!> Parameter for FFT-ScaRC combination
LOGICAL       :: SCARC_FFT = .FALSE.                        !> If .true. use combination of FFT and ScaRC
INTEGER       :: SCARC_PRESSURE_ITERATIONS = 0              !> internal counter for pressure iterations
CHARACTER(40) :: SCARC_TWOLEVEL = 'NONE'                    !> Type of two-level method

!> Parameter for MKL solver
CHARACTER(40) :: SCARC_MKL       = 'GLOBAL'                 !> Type of MKL solver (LOCAL->Pardiso/GLOBAL->Cluster_Sparse_solver)
CHARACTER(40) :: SCARC_MKL_MTYPE = 'SYMMETRIC'              !> Type of MKL solver (LOCAL->Pardiso/GLOBAL->Cluster_Sparse_solver)

!> Parameters for multigrid-type methods
CHARACTER(40) :: SCARC_MULTIGRID            = 'GEOMETRIC'   !> Type of MG-method (GEOMETRIC/ALGEBRAIC)
CHARACTER(40) :: SCARC_MULTIGRID_COARSENING = 'FALGOUT'     !> Coarsening strategy  (RS3/A1/A2/PMIS/FDS...)
CHARACTER(40) :: SCARC_MULTIGRID_INTERPOL   = 'DIRECT'      !> Interpolation strategy (DIRECT/RS/STANDARD)
CHARACTER(1)  :: SCARC_MULTIGRID_CYCLE      = 'V'           !> Cycling type  (F/V/W)
INTEGER       :: SCARC_MULTIGRID_LEVEL      = -1            !> User defined number of MG-levels (optionally)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS = 100           !> Max number of iterations
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY   = 1.E-8_EB      !> Requested accuracy for convergence

!> Parameters for Krylov-type methods
CHARACTER(40) :: SCARC_KRYLOV            = 'CG'             !> Type of Krylov-method (CG/BICG)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 500              !> Max number of iterations
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-8_EB         !> Requested accuracy for convergence
CHARACTER(40) :: SCARC_KRYLOV_INTERPOL   = 'NONE'           !> Interpolation type (only for two-level variants)

!> Parameters for smoothing method (used in multigrids-methods)
CHARACTER(40) :: SCARC_SMOOTH            = 'SSOR'           !> Smoother for MG (JACOBI/SSOR)
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 5                !> Max number of iterations
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-10_EB        !> Requested accuracy for convergence
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.80E+0_EB       !> Relaxation parameter

!> Parameters for preconditioning method (used in Krylov-methods)
CHARACTER(40) :: SCARC_PRECON            = 'FFT'            !> Preconditioner for CG/BICG (JACOBI/SSOR/FFT/PARDISO/MG)
INTEGER       :: SCARC_PRECON_ITERATIONS = 100              !> Max number of iterations
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-12_EB        !> Requested accuracy for convergence
REAL (EB)     :: SCARC_PRECON_OMEGA      = 0.80E+0_EB       !> Relaxation parameter

!> Parameters for coarse grid method
CHARACTER(40) :: SCARC_COARSE            = 'DIRECT'         !> Coarse grid solver (iterative CG-solver/direct MKL-solver)
CHARACTER(40) :: SCARC_COARSE_MTYPE      = 'SYMMETRIC'      !> Type of coarse grid matrix (nonsymmetric/symmetric)
INTEGER       :: SCARC_COARSE_ITERATIONS = 100              !> Max number of iterations for iterative solver
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-12_EB        !> Requested accuracy for iterative solver

!> Debugging parameters
CHARACTER(40) :: SCARC_DEBUG = 'NONE'                       !> Debugging level (NONE/LESS/MEDIUM/MUCH)
CHARACTER(40) :: SCARC_FN                                   !> File name for ScaRC debug messages
CHARACTER(40) :: SCARC_CASE = 'NONE'                        !> Debugging level (NONE/LESS/MEDIUM/MUCH)
CHARACTER(100):: SCARC_MESSAGE                              !> ScaRC messages
CHARACTER(40) :: SCARC_ROUTINE                              !> name of active SCARC-routine
INTEGER       :: LU_SCARC=0                                 !> Unit number for ScaRC debug file
INTEGER       :: SCARC_LAYER = 1                            !> Unit number for ScaRC debug file

!> order for the treatment of the single mesh faces
INTEGER :: FACE_ORDER_XYZ(6) = (/1,-1,2,-2,3,-3/)           !> Coordinate direction related order of mesh faces


!>PRIVATE
! Public declarations

!> ------------------------------------------------------------------------------------------------
!> Global constants
!> ------------------------------------------------------------------------------------------------
INTEGER, PARAMETER :: NSCARC_UNDEFINED = -1                    !> undefined value

INTEGER, PARAMETER :: NSCARC_DISCRET_STRUCTURED     =  1, &    !> structured discretization
                      NSCARC_DISCRET_UNSTRUCTURED   =  2       !> unstructured discretization

INTEGER, PARAMETER :: NSCARC_SCOPE_MAIN             =  1, &    !> method used as main solver
                      NSCARC_SCOPE_SMOOTH           =  2, &    !> method used as smoother
                      NSCARC_SCOPE_PRECON           =  3, &    !> method used as preconditiner
                      NSCARC_SCOPE_COARSE           =  4       !> method used as coarse grid solver

INTEGER, PARAMETER :: NSCARC_METHOD_KRYLOV          =  1, &    !> Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID       =  2, &    !> multigrid-method as global solver
                      NSCARC_METHOD_MKL             =  3       !> multigrid-method as global solver

INTEGER, PARAMETER :: NSCARC_KRYLOV_CG              =  1, &    !> CG   as Krylov solver
                      NSCARC_KRYLOV_BICG            =  2       !> BICG as Krylov solver

INTEGER, PARAMETER :: NSCARC_MULTIGRID_GEOMETRIC    =  1, &    !> geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC    =  2       !> algebraic multigrid

INTEGER, PARAMETER :: NSCARC_MKL_LOCAL              =  1, &    !> local use of MKL solver (Pardiso)
                      NSCARC_MKL_GLOBAL             =  2, &    !> global use of MKL solver (Cluster_Sparse_Solver)
                      NSCARC_MKL_COARSE             =  3       !> MKL solver on coarse grid level

INTEGER, PARAMETER :: NSCARC_EXCHANGE_BASIC         =  1, &    !> initialize wall information
                      NSCARC_EXCHANGE_WALLINFO      =  2, &    !> initialize wall information
                      NSCARC_EXCHANGE_VECTOR        =  3, &    !> matrix-vector communication
                      NSCARC_EXCHANGE_PRESSURE      =  4, &    !> vector values along internal boundaries
                      NSCARC_EXCHANGE_MEASURE       =  5, &    !> measure values along internal boundaries
                      NSCARC_EXCHANGE_MEASURE_ADD   =  6, &    !> measure values along internal boundaries
                      NSCARC_EXCHANGE_CELL_TYPE     =  7, &    !> cell types along internal boundaries
                      NSCARC_EXCHANGE_CELL_TYPE2    =  8, &    !> cell types II along internal boundaries
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
                      NSCARC_EXCHANGE_DISCRET       = 22       !> exchange DISCRET along internal boundaries

INTEGER, PARAMETER :: NSCARC_SMOOTH_JACOBI          =  1, &    !> smoothing by JACOBI-method
                      NSCARC_SMOOTH_SSOR            =  2, &    !> smoothing by SSOR-method
                      NSCARC_SMOOTH_FFT             =  3, &    !> smoothing by FFT-method
                      NSCARC_SMOOTH_PARDISO         =  4, &    !> smoothing by PARDISO-method
                      NSCARC_SMOOTH_CLUSTER         =  5       !> smoothing by PARDISO-method

INTEGER, PARAMETER :: NSCARC_PRECON_JACOBI          =  1, &    !> preconditioning by JACOBI-method
                      NSCARC_PRECON_SSOR            =  2, &    !> preconditioning by SSOR-method
                      NSCARC_PRECON_FFT             =  3, &    !> preconditioning by FFT-method
                      NSCARC_PRECON_PARDISO         =  4, &    !> preconditioning by PARDISO-method
                      NSCARC_PRECON_CLUSTER         =  5, &    !> preconditioning by CLUSTER-method
                      NSCARC_PRECON_MULTIGRID       =  6       !> preconditioning by MG-method

INTEGER, PARAMETER :: NSCARC_TWOLEVEL_NONE          =  0, &    !> no two levels, only one level
                      NSCARC_TWOLEVEL_ADD           =  1, &    !> additive 2-level method
                      NSCARC_TWOLEVEL_MUL           =  2, &    !> multiplicative 2-level method
                      NSCARC_TWOLEVEL_MUL2          =  3       !> multiplicative 2-level method

INTEGER, PARAMETER :: NSCARC_CYCLE_F                =  0, &    !> F-cycle for mg-method
                      NSCARC_CYCLE_V                =  1, &    !> V-cycle for mg-method
                      NSCARC_CYCLE_W                =  2, &    !> W-cycle for mg-method
                      NSCARC_CYCLE_SETUP            =  3, &    !> initialize cycle counts
                      NSCARC_CYCLE_RESET            =  4, &    !> reset cycle counts
                      NSCARC_CYCLE_PROCEED          =  5, &    !> proceed cycle counts
                      NSCARC_CYCLE_PRESMOOTH        =  6, &    !> presmoothing cycle
                      NSCARC_CYCLE_POSTSMOOTH       =  7, &    !> postsmoothing cycle
                      NSCARC_CYCLE_NEXT             =  8, &    !> perform next cycling loop
                      NSCARC_CYCLE_EXIT             =  9       !> exit cycling loop

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED          =  0, &    !> proceed loop
                      NSCARC_STATE_CONV0            =  1, &    !> convergence already for resin
                      NSCARC_STATE_CONV             =  2, &    !> convergence
                      NSCARC_STATE_DIVG             =  3       !> divergence

INTEGER, PARAMETER :: NSCARC_DEBUG_INFO0            =  0, &    !> info0  level of debugging requested
                      NSCARC_DEBUG_INFO1            =  1, &    !> info1  level of debugging requested
                      NSCARC_DEBUG_INFO2            =  2, &    !> info2  level of debugging requested
                      NSCARC_DEBUG_LESS             =  3, &    !> low    level of debugging requested
                      NSCARC_DEBUG_MEDIUM           =  4, &    !> medium level of debugging requested
                      NSCARC_DEBUG_MUCH             =  5, &    !> strong level of debugging requested
                      NSCARC_DEBUG_EXTREME          = 22, &    !> extreme level of debugging requested
                      NSCARC_DEBUG_MATLAB           = 23, &    !> extreme level of debugging requested
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
                      NSCARC_DEBUG_CELL_TYPE        = 18, &    !> show CELL_TYPE
                      NSCARC_DEBUG_GRAPH            = 19, &    !> show CELL_TYPE
                      NSCARC_DEBUG_COARSE           = 20, &    !> show coarse grid
                      NSCARC_DEBUG_PROLONGATION     = 21, &    !> show prolongation matrix
                      NSCARC_DEBUG_RESTRICTION      = 22       !> show restriction matrix

INTEGER, PARAMETER :: NSCARC_COARSENING_BASIC       =  1, &    !> basic coarsening
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

INTEGER, PARAMETER :: NSCARC_COARSE_ITERATIVE       =  1, &    !> iterative solution of coarse grid problem
                      NSCARC_COARSE_DIRECT          =  2       !> direct solution of coarse grid problem

INTEGER, PARAMETER :: NSCARC_SIZE_MATRIX            =  2, &    !> size of system matrix for compact system
                      NSCARC_SIZE_TRANSFER          =  3       !> size of transfer matrices for compact system

INTEGER, PARAMETER :: NSCARC_SOLVER_MAIN            =  1, &    !> main solver
                      NSCARC_SOLVER_PRECON          =  2, &    !> preconditioner
                      NSCARC_SOLVER_COARSE          =  3       !> coarse grid solver

INTEGER, PARAMETER :: NSCARC_VECTOR_H                  =  1, &    !> selection parameter for vector X
                      NSCARC_VECTOR_HS                 =  2, &    !> selection parameter for vector Y2
                      NSCARC_VECTOR_X                  =  3, &    !> selection parameter for vector Y2
                      NSCARC_VECTOR_F                  =  4, &    !> selection parameter for vector F
                      NSCARC_VECTOR_Y                  =  5, &    !> selection parameter for vector Y
                      NSCARC_VECTOR_G                  =  6, &    !> selection parameter for vector G
                      NSCARC_VECTOR_W                  =  7, &    !> selection parameter for vector R
                      NSCARC_VECTOR_D                  =  8, &    !> selection parameter for vector D
                      NSCARC_VECTOR_Z                  =  9, &    !> selection parameter for vector Z
                      NSCARC_VECTOR_MEASURE            = 10, &    !> selection parameter for vector MEASURE
                      NSCARC_VECTOR_CELL_TYPE          = 11, &    !> selection parameter for vector CELL_TYPE
                      NSCARC_VECTOR_GRAPH              = 12       !> selection parameter for vector GRAPH

INTEGER, PARAMETER :: NSCARC_MATRIX_SYSTEM          =  1, &    !> exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SUBDIAG         =  2, &    !> exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SUBDIAG_LOWER   =  3, &    !> exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SUBDIAG_UPPER   =  4, &    !> exchange subdiagonal matrix entries
                      NSCARC_MATRIX_TRANSFER        =  5, &    !> exchange prolongation matrix
                      NSCARC_MATRIX_STENCIL         =  6, &    !> exchange prolongation matrix
                      NSCARC_MATRIX_PROLONGATION    =  7, &    !> exchange prolongation matrix
                      NSCARC_MATRIX_RESTRICTION     =  8, &    !> exchange restriction matrix
                      NSCARC_MATRIX_GMG             =  9       !> exchange prolongation matrix

INTEGER, PARAMETER :: NSCARC_ACCURACY_ABSOLUTE      =  1, &    !> absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE      =  2       !> relative accuracy must be reached

REAL(EB), PARAMETER:: NSCARC_MEASURE_NONE           =  0.0_EB, &
                      NSCARC_MEASURE_ONE            =  1.0_EB, &  !> coarse-grid cell
                      NSCARC_MEASURE_COARSE         =  6.0_EB, &  !> coarse-grid cell
                      NSCARC_MEASURE_FINE           =  5.0_EB, &  !> fine-grid cell
                      NSCARC_MEASURE_SFINE          =  5.0_EB, &  !> strongly coupled fine-grid cell
                      NSCARC_MEASURE_WFINE          =  4.0_EB, &  !> weakly   coupled fine-grid cell
                      NSCARC_MEASURE_BDRY           =  1.0_EB     !> boundry weight

INTEGER, PARAMETER :: NSCARC_CELL_TYPE_COARSE       =  1, &    !> coarse-grid cell
                      NSCARC_CELL_TYPE_COARSE0      =  3, &    !> special coarse-grid cell
                      NSCARC_CELL_TYPE_COMMON       =  3, &    !> common cell
                      NSCARC_CELL_TYPE_FINE         = -1, &    !> fine-grid cell
                      NSCARC_CELL_TYPE_FINE0        = -3, &    !> special fine-grid cell
                      NSCARC_CELL_TYPE_SFINE        = -1, &    !> strongly coupled fine-grid cell
                      NSCARC_CELL_TYPE_WFINE        = -2, &    !> weakly   coupled fine-grid cell
                      NSCARC_CELL_TYPE_FPNT         = -1, &    !> special f-point
                      NSCARC_CELL_TYPE_ZPNT         = -2, &    !> special z-point
                      NSCARC_CELL_TYPE_SFPNT        = -3, &    !> special sf-point
                      NSCARC_CELL_TYPE_CPNT         =  2       !> special c-point

INTEGER, PARAMETER :: NSCARC_INTERPOL_STANDARD      =  1, &    !> standard interpolation
                      NSCARC_INTERPOL_CONSTANT      =  2, &    !> standard interpolation
                      NSCARC_INTERPOL_BILINEAR      =  3, &    !> standard interpolation
                      NSCARC_INTERPOL_CLASSICAL     =  4, &    !> classical interpolation
                      NSCARC_INTERPOL_CLASSICAL2    =  5, &    !> classical interpolation
                      NSCARC_INTERPOL_DIRECT        =  6, &    !> direct interpolation
                      NSCARC_INTERPOL_DIRECT_BDRY   =  7, &    !> direct interpolation with special boundary
                      NSCARC_INTERPOL_MULTIPASS     =  8, &    !> multipass interpolation
                      NSCARC_INTERPOL_GMG           =  9, &    !> GMG-like interpolation
                      NSCARC_INTERPOL_GMG3          = 10       !> GMG3-like interpolation

INTEGER, PARAMETER :: NSCARC_CASE_CD_NSA_2D         =  1, &    !> CD_NSA_2D-case
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

INTEGER, PARAMETER :: NSCARC_LEVEL_MIN              =  0, &    !> minimum multigrid level
                      NSCARC_LEVEL_MAX              = 15, &    !> maximum multigrid level
                      NSCARC_LEVEL_SINGLE           =  1, &    !> only one grid level needed
                      NSCARC_LEVEL_MULTI            =  2, &    !> multi grid levels needed
                      NSCARC_LEVEL_AMG              =  3       !> maximum number of grid levels 

INTEGER, PARAMETER :: NSCARC_LAYER_ONE              =  1, &    !> communication of one abutting layer
                      NSCARC_LAYER_TWO              =  2       !> communication of two abutting layers

INTEGER, PARAMETER :: NSCARC_DUMP_RHS               =  1, &    !> dump rhs
                      NSCARC_DUMP_PRES              =  2       !> dump pressure

INTEGER, PARAMETER :: NSCARC_STENCIL_CENTRAL        =  1, &    !> standard 5- or 7-point stencil
                      NSCARC_STENCIL_AMG            =  2       !> arbitrary AMG-stencil

INTEGER, PARAMETER :: NSCARC_DIAG_MAIN              =  1, &    !> standard 5- or 7-point stencil
                      NSCARC_DIAG_LOWER             =  2, &    !> standard 5- or 7-point stencil
                      NSCARC_DIAG_UPPER             =  3       !> standard 5- or 7-point stencil

INTEGER, PARAMETER :: NSCARC_COUPLING_MAX           = 10       !> maximum of possible couplings in stencil

INTEGER, PARAMETER :: NSCARC_NUM_FACES              =  6       !> number of faces per mesh
INTEGER, PARAMETER :: NSCARC_MAX_FACE_NEIGHBORS     = 10       !> max number neighbors per mesh face
INTEGER, PARAMETER :: NSCARC_MAX_STENCIL            =  7
INTEGER, PARAMETER :: NSCARC_MAX_MESH_NEIGHBORS     =  6*NSCARC_MAX_FACE_NEIGHBORS

INTEGER, PARAMETER  :: NSCARC_ZERO_INTEGER = 0                  !> zero initialization for integer variables
REAL(EB), PARAMETER :: NSCARC_ZERO_REAL    = 0.0_EB             !> zero initialization for real variables

INTEGER :: IERROR = 0
INTEGER :: ITE_CG, ITE_MG, ITE_SM

!> --------------------------------------------------------------------------------------------
!> Logical indicators for different methods 
!> --------------------------------------------------------------------------------------------
LOGICAL :: BCG = .FALSE.
LOGICAL :: BMG = .FALSE.
LOGICAL :: BCGGMG  = .FALSE.
LOGICAL :: BCGADD  = .FALSE.
LOGICAL :: BCGMUL  = .FALSE.
LOGICAL :: BGMG = .FALSE.
LOGICAL :: BAMG = .FALSE.
LOGICAL :: BTWOLEVEL   = .FALSE.
LOGICAL :: BMULTILEVEL = .FALSE.

!> --------------------------------------------------------------------------------------------
!> Unstructured grid stuff - based on Mvanella
!> --------------------------------------------------------------------------------------------
INTEGER,  PARAMETER :: IS_GASPHASE  = -1
INTEGER,  PARAMETER :: IS_CUTCFE    =  0
INTEGER,  PARAMETER :: IS_SOLID     =  1
INTEGER,  PARAMETER :: IS_UNDEFINED =-11

INTEGER,  PARAMETER :: IS_CGSC   = 1 ! Face media type: IS_GASPHASE, IS_SOLID or IS_CUTCFE.
INTEGER,  PARAMETER :: IS_UNKH   = 2 ! H unknown number.
INTEGER,  PARAMETER :: IS_NCVARS = 2 ! Number of face variables 

! Cartesian Cell centered variables, actual case initialized as CC_IBM=.FALSE.:
INTEGER :: CGSC=IS_CGSC, UNKH=IS_UNKH, NCVARS=IS_NCVARS

!> --------------------------------------------------------------------------------------------
!> Global variables
!> --------------------------------------------------------------------------------------------
!> use integer types for the user defined input data (based on SCARC_TYPE_... variables)
INTEGER :: TYPE_DISCRET     = NSCARC_DISCRET_STRUCTURED    !> Type of discretization (structured/unstructured)
INTEGER :: TYPE_SCOPE       = NSCARC_SCOPE_MAIN            !> Type of surrounding solver scope
INTEGER :: TYPE_SCOPE0      = NSCARC_SCOPE_MAIN            !> Type of surrounding solver scope, restore
INTEGER :: TYPE_METHOD      = NSCARC_METHOD_KRYLOV         !> Type of ScaRC method
INTEGER :: TYPE_METHOD_CORE = NSCARC_METHOD_KRYLOV         !> Core type of ScaRC method
INTEGER :: TYPE_TWOLEVEL    = NSCARC_TWOLEVEL_NONE         !> Type of two-level method
INTEGER :: TYPE_INTERPOL    = NSCARC_INTERPOL_CONSTANT     !> Type of interpolation method
INTEGER :: TYPE_PRECON      = NSCARC_PRECON_SSOR           !> Type of preconditioner for iterative solver
INTEGER :: TYPE_PRECON_CORE = NSCARC_PRECON_SSOR           !> Core type of preconditioner for iterative solver
INTEGER :: TYPE_KRYLOV      = NSCARC_KRYLOV_CG             !> Type of Krylov method (CG/BICG)
INTEGER :: TYPE_MULTIGRID   = NSCARC_MULTIGRID_GEOMETRIC   !> Type of multigrid method (GMG/AMG)
INTEGER :: TYPE_MKL         = NSCARC_UNDEFINED             !> Type of MKL method (PARDISO/CLUSTER_SPARSE_SOLVER)
INTEGER :: TYPE_ACCURACY    = NSCARC_ACCURACY_ABSOLUTE     !> Type of requested accuracy
INTEGER :: TYPE_SMOOTH      = NSCARC_SMOOTH_SSOR           !> Type of smoother for multigrid method
INTEGER :: TYPE_CYCLE       = NSCARC_CYCLE_V               !> Type of cycling for multigrid method
INTEGER :: TYPE_COARSENING  = NSCARC_UNDEFINED             !> Type of coarsening algorithm for AMG
INTEGER :: TYPE_COARSE      = NSCARC_COARSE_ITERATIVE      !> Type of coarse grid solver for multigrid method
INTEGER :: TYPE_CASE        = NSCARC_UNDEFINED             !> Type of predefined test case
INTEGER :: TYPE_DEBUG       = NSCARC_UNDEFINED             !> Type of debugging level
INTEGER :: TYPE_INITIAL     = NSCARC_UNDEFINED             !> Type of initial solution
INTEGER :: TYPE_EXCHANGE    = NSCARC_UNDEFINED             !> Type of data exchange
INTEGER :: TYPE_VECTOR         = NSCARC_UNDEFINED             !> Type of vector to point to
INTEGER :: TYPE_LAYER       = NSCARC_LAYER_ONE             !> Type of layers (for overlap)

INTEGER :: TYPE_MKL_LEVEL(NSCARC_LEVEL_MAX)=NSCARC_UNDEFINED

INTEGER :: NLEVEL, NLEVEL_MAX, NLEVEL_MIN                 !> Total, minimum and maximum number of multigrid levels
INTEGER :: N_REQ, N_EXCHANGE, TAG                               !> Variables for data exchange
INTEGER :: SNODE, RNODE                                   !> Process identifier for data exchange

INTEGER :: N_CELLS_GLOBAL(NSCARC_LEVEL_MAX)  = 0               !> number of global cells 
INTEGER :: NUNKH_GLOB(NSCARC_LEVEL_MAX) = 0               !> number of global cells 

INTEGER,  ALLOCATABLE, DIMENSION (:)  :: REQ              !> Request array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: COUNTS           !> Counter array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: DISPLS           !> Displacement array for data exchange

REAL(EB):: SP_GLOBAL                                      !> Global scalar product
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: SP_LOCAL         !> Local scalar procucts

!> ------------------------------------------------------------------------------------------------
!> Time measurement
!> ------------------------------------------------------------------------------------------------
TYPE SCARC_TIME_TYPE
REAL(EB) :: OVERALL              !> complete time
REAL(EB) :: SOLVER               !> time for solver (general version)
REAL(EB) :: CLUSTER              !> time for cluster solver
REAL(EB) :: PARDISO              !> time for pardiso solver
REAL(EB) :: KRYLOV               !> time for krylov solver
REAL(EB) :: MULTIGRID            !> time for multigrid solver
REAL(EB) :: MATVEC               !> time for matrix vector multiplication
REAL(EB) :: SCALPROD             !> time for scalar product 
REAL(EB) :: L2NORM               !> time for l2-norm
REAL(EB) :: PRECON               !> time for preconditioner
REAL(EB) :: SMOOTH               !> time for smoother
REAL(EB) :: COARSE               !> time for coarse grid solver
REAL(EB) :: EXCHANGE             !> time for data exchange
END TYPE SCARC_TIME_TYPE

TYPE (SCARC_TIME_TYPE), ALLOCATABLE, DIMENSION(:) :: TSETUP, TSUM, TSTEP

PRIVATE

!> --------------------------------------------------------------------------------------------
!> Scopes for different solution methods (CG/BICG/GMG/AMG/MKL)
!> --------------------------------------------------------------------------------------------
TYPE SCARC_SCOPE_TYPE
INTEGER :: X, F, Y, G, W, D, Z, AUX         !> local references of solver vectors
INTEGER :: ITE, NIT                         !> local references to iteration counters
INTEGER :: NTYPE                            !> type of scope (main/precon)
INTEGER :: ISTACK                           !> position in calling stack
REAL(EB) :: EPS, RES, RESIN, OMEGA          !> local references to solver parameters
CHARACTER(30) :: CSOLVER = 'null'
END TYPE SCARC_SCOPE_TYPE

!> --------------------------------------------------------------------------------------------
!> Face information related to wall cells and neighbors
!> --------------------------------------------------------------------------------------------
TYPE SCARC_FACE_TYPE
INTEGER :: NFC, NFW                                !> number of cells and wall cells along face
INTEGER :: NFX, NFY, NFZ                           !> local face dimensions
INTEGER :: NCPL                                    !> number of adjacent couplings
INTEGER :: NUM_NEIGHBORS = 0                       !> number of adjacent neighbors
INTEGER :: IWG_PTR                                 !> first (global) IW number to that face
INTEGER :: IOFFSET_WALL   = 0                      !> counter for wall cells over all faces
INTEGER , ALLOCATABLE, DIMENSION(:) :: NEIGHBORS   !> adjacent neighbors
REAL(EB), POINTER, DIMENSION(:) :: DH              !> adjacent grid sizes
END TYPE SCARC_FACE_TYPE

!> --------------------------------------------------------------------------------------------
!> Wall information related to neighbors and BC's
!> --------------------------------------------------------------------------------------------
TYPE SCARC_WALL_TYPE
INTEGER :: BTYPE                                  !> boundary type of wall cell (1:Dirichlet/2:Neumann/3:Internal)
INTEGER :: BOUNDARY_TYPE = 0                      !> boundary type of wall cell (1:Dirichlet/2:Neumann/3:Internal)
INTEGER :: IOR                                    !> orientation of wall cell
INTEGER :: NOM                                    !> neighbor at wall cell
INTEGER :: NCPL = 1                               !> number of couplings at wall cell (depending on resolution of neighbor)
INTEGER :: IWL                                    !> corresponding local wall cell number for neighbor NOM
INTEGER :: CGSC, UNKH                             !> variables of DISCRET-structure 
INTEGER :: ICW  = IS_UNDEFINED                    !> internal wall cell for IW
INTEGER :: ICO                                    !> overlapping cell for IW
INTEGER :: IXG, IYG, IZG                          !> x-, y- and z-indices of ghost cells
INTEGER :: IXW, IYW, IZW                          !> x-, y- and z-indices of (internal) wall cells
INTEGER :: IXN(2), IYN(2), IZN(2)                 !> x-, y- and z-indices of neighboring cells
INTEGER, ALLOCATABLE, DIMENSION(:) :: ICE, ICE2   !> extended cell for IW
INTEGER, ALLOCATABLE, DIMENSION(:) :: ICG, ICG2   !> ghost cell for IW
END TYPE SCARC_WALL_TYPE

!> --------------------------------------------------------------------------------------------
!> Administration of grids inkluding mappings between different informations
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MAPPING_TYPE
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICE_TO_IWG        !> mapping from ICE to IWG
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICE_TO_IWL        !> mapping from ICE to IWL
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICE_TO_ICG        !> mapping from ICE to ICG
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICG_TO_IWG        !> mapping from ICG to IWG
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICG_TO_ICE        !> mapping from ICG to ICE
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICG_TO_ICO        !> mapping from ICG to ICE
INTEGER, ALLOCATABLE, DIMENSION (:)  :: IWL_TO_IWG        !> mapping from IWL to IWG
INTEGER, ALLOCATABLE, DIMENSION (:)  :: IWL_TO_ICW        !> mapping from IWL to ICW
INTEGER, ALLOCATABLE, DIMENSION (:)  :: IWL_TO_ICO        !> mapping from IWL to ICO
INTEGER, ALLOCATABLE, DIMENSION (:,:):: IWL_TO_ICG        !> mapping from IWL to ICG
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICN_TO_ICE        !> mapping from ICN to ICE    ! ACHTUNG: RIESIG
END TYPE SCARC_MAPPING_TYPE

!> --------------------------------------------------------------------------------------------
!> Matrices including storage pointers and lengths 
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MATRIX_TYPE
INTEGER :: NAV, NAC, NAR, NAS, NAE                        !> number of elements for different parts of structure
INTEGER :: NSTENCIL                                       !> number of points in matrix stencil
REAL (EB), ALLOCATABLE, DIMENSION (:) :: VAL              !> system, prolongation, restriction and strength matrix
INTEGER,   ALLOCATABLE, DIMENSION (:) :: ROW              !> row pointer 
INTEGER,   ALLOCATABLE, DIMENSION (:) :: COL              !> column pointer
INTEGER,   ALLOCATABLE, DIMENSION (:) :: COL_GLOBAL       !> column pointer
INTEGER,   ALLOCATABLE, DIMENSION (:) :: POS              !> position
INTEGER,   ALLOCATABLE, DIMENSION (:) :: STENCIL          !> size information 
INTEGER,   ALLOCATABLE, DIMENSION (:) :: TAG              !> marking tag
END TYPE SCARC_MATRIX_TYPE

!> --------------------------------------------------------------------------------------------
!> Package of vectors used for different solvers
!> --------------------------------------------------------------------------------------------
TYPE SCARC_SOLVER_TYPE
REAL (EB), ALLOCATABLE, DIMENSION (:) :: X            !> solution vectors for CG, MG and MKL
REAL (EB), ALLOCATABLE, DIMENSION (:) :: F            !> right hand side vectors for CG, MG and MKL
REAL (EB), ALLOCATABLE, DIMENSION (:) :: D            !> defect vectors for CG and MG
REAL (EB), ALLOCATABLE, DIMENSION (:) :: Y            !> auxiliary vectors for CG and MG
REAL (EB), ALLOCATABLE, DIMENSION (:) :: G            !> auxiliary vectors for CG and MG
REAL (EB), ALLOCATABLE, DIMENSION (:) :: W            !> auxiliary vectors for CG and MG
REAL (EB), ALLOCATABLE, DIMENSION (:) :: Z            !> auxiliary vectors for CG and MG
END TYPE SCARC_SOLVER_TYPE

!> --------------------------------------------------------------------------------------------
!> Obstruction information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_OBST_TYPE
INTEGER :: I1, I2, J1, J2, K1, K2
END TYPE SCARC_OBST_TYPE

!> --------------------------------------------------------------------------------------------
!> Geometric multigrid information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_GMG_TYPE
INTEGER :: CYCLE_COUNT(2) = 0                          !> Counter for multigrid cycling
END TYPE SCARC_GMG_TYPE

!> --------------------------------------------------------------------------------------------
!> Algebraic multigrid information 
!> --------------------------------------------------------------------------------------------
TYPE SCARC_AMG_TYPE
REAL(EB):: MAX_ROWSUM=0.9_EB
INTEGER :: NCF, NCFE, NCFS                                !> number of internal and extended fine cells 
INTEGER :: NCC, NCCE, NCCS                                !> number of internal and extended coarse cells 
INTEGER :: NCE, NCCI                                      !> number of internal coarse cells
INTEGER :: NC, NCW, NCG                                   !> number of coarse wall and ghost cells 
INTEGER :: NP, NPE, NPS                                   !> number of elements for prolongation matrix
INTEGER :: NR, NRE, NRS                                   !> number of elements for restriction matrix
INTEGER :: ICG0                                           !> auxiliary number
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: MEASURE          !> measure of different cells
INTEGER,  ALLOCATABLE, DIMENSION (:,:):: CELL_MAP         !> cell mapping for coarsening
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: CELL_TYPE        !> cell types 
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: CELL_PTR         !> pointers for coarsening
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: GRAPH            !> graph vector 
END TYPE SCARC_AMG_TYPE

!> --------------------------------------------------------------------------------------------
!> Store scoping information of parent (calling) solver routing
!> --------------------------------------------------------------------------------------------
TYPE SCARC_PARENT_TYPE
INTEGER :: TYPE_METHOD           !> type of solution method
INTEGER :: TYPE_SCOPE            !> type of related scope
INTEGER :: TYPE_PRECON           !> type of preconditioning method
INTEGER :: TYPE_SMOOTH           !> type of smoothing
INTEGER :: TYPE_TWOLEVEL         !> type of twolevel method
INTEGER :: TYPE_INTERPOL         !> type of interpolation method
INTEGER :: TYPE_ACCURACY         !> type of accuracy requirements (relative/absolute)
INTEGER :: TYPE_CYCLE            !> type of multigrid cycle
END TYPE SCARC_PARENT_TYPE

!> --------------------------------------------------------------------------------------------
!> Preconditioners information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_FFT_TYPE
REAL (EB), ALLOCATABLE, DIMENSION (:, :)    :: BXS, BXF, BYS, BYF, BZS, BZF
REAL (EB), ALLOCATABLE, DIMENSION (:, :, :) :: WORK
END TYPE SCARC_FFT_TYPE

!> --------------------------------------------------------------------------------------------
!> MKL information
!> --------------------------------------------------------------------------------------------
#ifdef WITH_MKL
TYPE SCARC_MKL_TYPE
TYPE(MKL_PARDISO_HANDLE),               ALLOCATABLE :: PT_H(:), PT(:)
TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), ALLOCATABLE :: CT_H(:), CT(:)
INTEGER, ALLOCATABLE :: IPARM(:)                      
INTEGER :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR, MSGLVL
INTEGER :: PERM(1)
END TYPE SCARC_MKL_TYPE
#endif

!> --------------------------------------------------------------------------------------------
!> Data related to a single grid level
!> --------------------------------------------------------------------------------------------
TYPE SCARC_LEVEL_TYPE

INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: DISCRET             !> basic discretization information
INTEGER :: SUBDIVISION(3,-3:3)=0                                !> basic information related to single faces

INTEGER :: NOM                                                  !> number of adjacent neighbor
INTEGER :: IOR = 0                                              !> local orientations
INTEGER :: N_OBST                                               !> number of obstructions
INTEGER :: N_CELLS                                              !> number of cells
INTEGER :: N_WALL_CELLS                                         !> number of wall cells
INTEGER :: N_EXTERNAL_WALL_CELLS                                !> number of external cells
INTEGER :: N_INTERNAL_WALL_CELLS                                !> number of internal cells
INTEGER :: N_OFFSET_MESH                                        !> offset for global mesh numbering
INTEGER :: NUNKH_GLOB = 0                                       !> number of global cells 

INTEGER :: NX, NY, NZ                                           !> number of grid cells in x-, y- and z-direction
INTEGER :: NC, NCS                                              !> number of cells
INTEGER :: NCW, NA                                              !> number of coarse wall cells
INTEGER :: NW, NWL                                              !> number of global and local wall cells
INTEGER :: NC_GLOBAL = 0                                        !> number of global cells 
INTEGER :: NCG=0, NCG0=0                                        !> number of ghost cells
INTEGER :: NCE=0, NCE0=0                                        !> number of extended cells plus auxiliary variable
INTEGER :: NCO=0, NCO0=0                                        !> number of overlapping cells plus auxiliary variable
INTEGER :: NCPL=1, NCPL_MAX=-10                                 !> number of couplings
INTEGER :: NCPLS, NCPLR                                         !> number of couplings to send and read

REAL(EB) :: DX , DY , DZ                                        !> step sizes in x-, y- and z-direction
REAL(EB) :: DXI, DYI, DZI                                       !> inversed of step sizes in x-, y- and z-direction
REAL(EB) :: DXI2, DYI2, DZI2                                    !> squared and inversed step sizes in x-, y- and z-direction
REAL(EB) :: DH = 0.0_EB                                         !> local step sizes
REAL (EB), POINTER, DIMENSION (:)     :: XCOR, YCOR, ZCOR       !> coordinate vectors in x-, y- and z-direction
REAL (EB), POINTER, DIMENSION (:)     :: XMID, YMID, ZMID       !> midpoint vectors in x-, y- and z-direction
REAL (EB), ALLOCATABLE, DIMENSION (:) :: DXL, DYL, DZL          !> step size vectors in x-, y- and z-direction
  
INTEGER, ALLOCATABLE, DIMENSION(:)    :: NC_LOCAL               !> number of local cells
INTEGER, ALLOCATABLE, DIMENSION(:)    :: NC_OFFSET              !> offset in cell numbering 
INTEGER, ALLOCATABLE, DIMENSION(:)    :: NUNKH_LOC              !> number of local cells
INTEGER, ALLOCATABLE, DIMENSION(:)    :: UNKH_IND               !> offset in cell numbering 
INTEGER, ALLOCATABLE, DIMENSION(:)    :: INTERNAL_BDRY_CELL     !> index of internal boundary cells
INTEGER, ALLOCATABLE, DIMENSION(:)    :: EXT_PTR                !> index of external cells
INTEGER, ALLOCATABLE, DIMENSION(:,:,:):: CELL_INDEX             !> wall index
INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: WALL_INDEX             !> wall index

INTEGER :: IWL = 0                                              !> local wall cell numbers
INTEGER :: ICG = 0                                              !> ghost cell counter
INTEGER :: ICG0 = 0                                             !> auxiliary ghost cell counter
INTEGER :: ICG_PTR = 0                                          !> ghost cell pointer
INTEGER :: ICO_PTR = 0                                          !> overlapping cell pointer
INTEGER :: ICE_PTR = 0                                          !> extended cell pointer
INTEGER :: IWL_PTR = 0                                          !> local wall cell pointer
INTEGER :: IWG_PTR = 0                                          !> global wall cell pointer

TYPE (SCARC_WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: WALL       !> description of wall information
TYPE (SCARC_FACE_TYPE), ALLOCATABLE, DIMENSION(:) :: FACE       !> description of faces
TYPE (SCARC_OBST_TYPE), ALLOCATABLE, DIMENSION(:) :: OBST       !> description of obstructions

TYPE (SCARC_MAPPING_TYPE) :: MAP                                !> mappings between different grid description arrays
TYPE (SCARC_MATRIX_TYPE)  :: A, AS, P, R, S, ST                 !> system, system symmetric, prolongation, restriction, ...
                                                                !> strength, strength transpose
TYPE (SCARC_SOLVER_TYPE)  :: MAIN, PRECON, COARSE               !> different types of solvers
TYPE (SCARC_FFT_TYPE)     :: FFT                                !> description of preconditioner
TYPE (SCARC_GMG_TYPE)     :: GMG                                !> description of geometric multigrid
TYPE (SCARC_AMG_TYPE)     :: AMG                                !> description of algebraic multigrid
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE)     :: MKL                                !> MKL-type
#endif

END TYPE SCARC_LEVEL_TYPE

!> --------------------------------------------------------------------------------------------
!> Administration other mesh data needed for the coupling of adjacent neighbors
!> --------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE
REAL (EB), ALLOCATABLE, DIMENSION (:) :: SEND_REAL, RECV_REAL       !> main real send and receive buffers
INTEGER  , ALLOCATABLE, DIMENSION (:) :: SEND_INT, RECV_INT         !> main integer send and receive buffers
REAL(EB) :: SEND_REAL_BASIC(50), RECV_REAL_BASIC(50)                          !> initial real send and receive buffers
INTEGER  :: SEND_INT_BASIC(50)  , RECV_INT_BASIC(50)                          !> initial integer send and receive buffers
INTEGER  :: NICMAX_R=0, NICMAX_S=0, NIC_R=0, NIC_S=0
TYPE (SCARC_LEVEL_TYPE), ALLOCATABLE, DIMENSION(:) :: LEVEL
END TYPE OSCARC_TYPE


!> --------------------------------------------------------------------------------------------
!> Basic administration type for ScaRC-method
!> --------------------------------------------------------------------------------------------
TYPE SCARC_TYPE
INTEGER :: NUM_NEIGHBORS = 0                                 !> Number of adjacent neighbors of whole mesh
INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS              !> List of adjacent neighbors of whole mesh
INTEGER  :: IBAR, JBAR, KBAR                                 !> number of cells in x-, y- and z-direction on finest level
REAL(EB) :: XS, XF, YS, YF, ZS, ZF                           !> x-, y- and z-bounds of mesh (same for all levels)
TYPE (SCARC_SCOPE_TYPE),  ALLOCATABLE, DIMENSION(:) :: SCOPE     !> Scope type
TYPE (SCARC_LEVEL_TYPE) , ALLOCATABLE, DIMENSION(:) :: LEVEL     !> Grid level type
TYPE (OSCARC_TYPE)      , ALLOCATABLE, DIMENSION(:) :: OSCARC    !> ScaRC type on other mesh
END TYPE SCARC_TYPE


!> --------------------------------------------------------------------------------------------
!> Globally used types
!> --------------------------------------------------------------------------------------------
TYPE (SCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET ::  SCARC

CONTAINS


!> ------------------------------------------------------------------------------------------------
!> Initialize ScaRC structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP
REAL(EB):: TNOW

SCARC_ROUTINE = 'SCARC_SETUP'

CALL SCARC_SETUP_TIMING          !> allocate different arrays for time measurments
TNOW = CURRENT_TIME()

CALL SCARC_PARSE_INPUT           !> Parse input parameters for ScaRC

CALL SCARC_SETUP_DEBUGGING       !> setup debugging information if requested
CALL SCARC_SETUP_LEVEL           !> setup different grid levels
CALL SCARC_SETUP_TYPES           !> setup basic ScaRC-types for all used grid levels
CALL SCARC_SETUP_MESHES          !> setup mesh information
CALL SCARC_SETUP_DISCRET         !> setup discretization information
CALL SCARC_SETUP_INTERFACE       !> setup structures related to mesh interfaces
CALL SCARC_SETUP_GLOBAL_FINE     !> setup global variables on fine level
CALL SCARC_SETUP_WALLS           !> setup information along neighboring walls
CALL SCARC_SETUP_GLOBAL_COARSE   !> setup global variables for all coarser levels
CALL SCARC_SETUP_EXCHANGE        !> setup information for data exchange
CALL SCARC_SETUP_SYSTEM          !> setup linear system of equations
CALL SCARC_SETUP_COARSENING      !> setup coarsening on different grid levels if requested (AMG only)
CALL SCARC_SETUP_SOLVER          !> setup single solvers (allocate work and auxiliary vectors)

TSETUP(MYID+1)%OVERALL = TSETUP(MYID+1)%OVERALL + CURRENT_TIME() - TNOW

END SUBROUTINE SCARC_SETUP


!> ----------------------------------------------------------------------------------------------------
!> Determine types of input parameters
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PARSE_INPUT

SCARC_ROUTINE = 'SCARC_PARSE_INPUT'

ITERATE_PRESSURE = .TRUE.  ! Although there is no need to do pressure iterations to drive down velocity error
                           ! leave it .TRUE. to write out velocity error diagnostics.

!> 
!> ------------- set type of discretization
!> 
SELECT CASE (TRIM(SCARC_DISCRETIZATION))
   CASE ('STRUCTURED')
      TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED
   CASE ('UNSTRUCTURED')
      TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED
END SELECT
PRES_ON_WHOLE_DOMAIN = (TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED)

!> 
!> ------------ set type of global solver
!> 
SELECT CASE (TRIM(SCARC_METHOD))

   CASE ('KRYLOV')

      TYPE_METHOD      = NSCARC_METHOD_KRYLOV
      TYPE_METHOD_CORE = NSCARC_METHOD_KRYLOV

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
         CASE DEFAULT
            WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_TWOLEVEL)
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
      END SELECT

      !> set type of interpolation for two-level Krylov method
      SELECT CASE (TRIM(SCARC_KRYLOV_INTERPOL))
         CASE ('NONE')
            TYPE_INTERPOL = NSCARC_UNDEFINED
         CASE ('CONSTANT')
            TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
         CASE ('BILINEAR')
            TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
         CASE DEFAULT
            WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_KRYLOV_INTERPOL)
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
               CASE ('PARDISO')
                  TYPE_SMOOTH = NSCARC_SMOOTH_PARDISO
               CASE ('CLUSTER')
                  TYPE_SMOOTH = NSCARC_SMOOTH_CLUSTER
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
      TYPE_PRECON_CORE = TYPE_PRECON

   CASE ('MULTIGRID')

      TYPE_METHOD      = NSCARC_METHOD_MULTIGRID
      TYPE_METHOD_CORE = NSCARC_METHOD_MULTIGRID

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

      !> set type of smoother (JACOBI/SSOR)
      SELECT CASE (TRIM(SCARC_SMOOTH))                        !> use same parameters as for preconditioner
         CASE ('JACOBI')
            TYPE_PRECON = NSCARC_PRECON_JACOBI
            TYPE_SMOOTH = NSCARC_SMOOTH_JACOBI
         CASE ('SSOR')
            TYPE_PRECON = NSCARC_PRECON_SSOR
            TYPE_SMOOTH = NSCARC_SMOOTH_SSOR
         CASE ('FFT')
            TYPE_PRECON = NSCARC_PRECON_FFT
            TYPE_SMOOTH = NSCARC_SMOOTH_FFT
#ifdef WITH_MKL
         CASE ('PARDISO')
            TYPE_PRECON = NSCARC_PRECON_PARDISO
            TYPE_SMOOTH = NSCARC_SMOOTH_PARDISO
         CASE ('CLUSTER')
            TYPE_PRECON = NSCARC_PRECON_CLUSTER
            TYPE_SMOOTH = NSCARC_SMOOTH_CLUSTER
#endif
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

!> 
!> if a multigrid solver is used (either as main solver or as preconditioner)
!> set types for multigrid, coarse grid solver and cycling pattern
!> 
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_PRECON_MULTIGRID) THEN

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

   !> set type of coarse grid solver (CG/GE)
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

   !> set type of cycling pattern (F/V/W)
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

!> set type of coarse grid solver
SELECT CASE (TRIM(SCARC_COARSE))
   CASE ('ITERATIVE')
      TYPE_COARSE = NSCARC_COARSE_ITERATIVE
      TYPE_KRYLOV = NSCARC_KRYLOV_CG
   CASE ('DIRECT')
      TYPE_COARSE = NSCARC_COARSE_DIRECT
      TYPE_MKL    = NSCARC_MKL_COARSE
      IF (SCARC_COARSE_MTYPE == 'SYMMETRIC') THEN
         SCARC_MKL_MTYPE = 'SYMMETRIC'
      ELSE
         SCARC_MKL_MTYPE = 'NONSYMMETRIC'
      ENDIF
   CASE DEFAULT
      WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_COARSE)
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
END SELECT

!> 
!> set test case
!> 
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
      TYPE_CASE = NSCARC_UNDEFINED
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
      WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_ACCURACY)
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
END SELECT

!> 
!> set level of debugging (NONE/LESS/MEDIUM/MUCH)
!> 
SELECT CASE (TRIM(SCARC_DEBUG))
   CASE ('NONE')
      TYPE_DEBUG = NSCARC_UNDEFINED
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
   CASE ('MATLAB')
      TYPE_DEBUG = NSCARC_DEBUG_MATLAB
   CASE DEFAULT
      WRITE(SCARC_MESSAGE,'(3A)') TRIM(SCARC_ROUTINE),': Error with input parameter ',TRIM(SCARC_DEBUG)
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
END SELECT

!> 
!> set type of initial solution (not yet used, may be used to define own initial vector)
!> 
SELECT CASE (TRIM(SCARC_INITIAL))
   CASE ('null')
      TYPE_INITIAL = NSCARC_UNDEFINED
   CASE DEFAULT
      TYPE_INITIAL = NSCARC_UNDEFINED
END SELECT

!> Define some logical variables for notational convenience
BCG = TYPE_METHOD == NSCARC_METHOD_KRYLOV
BCGGMG  = BCG .AND. TYPE_PRECON == NSCARC_PRECON_MULTIGRID .AND. &
          TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC

BMG = TYPE_METHOD == NSCARC_METHOD_MULTIGRID
BGMG = BMG .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC
BAMG = BMG .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC

BTWOLEVEL   = BCG .AND. TYPE_PRECON /= NSCARC_PRECON_MULTIGRID .AND. &
              TYPE_TWOLEVEL > NSCARC_TWOLEVEL_NONE
BMULTILEVEL = BGMG .OR. BCGGMG .OR. BTWOLEVEL

BCGADD = BTWOLEVEL .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_ADD
BCGMUL = BTWOLEVEL .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_MUL

END SUBROUTINE SCARC_PARSE_INPUT


!> ------------------------------------------------------------------------------------------------
!> Setup debug file if requested
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DEBUGGING
INTEGER:: NM, LASTID

SCARC_ROUTINE = 'SCARC_SETUP_DEBUGGING'
IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) THEN
   LASTID = -99999
   DO NM=LOWER_MESH_INDEX, UPPER_MESH_INDEX
      IF (MYID == LASTID) CYCLE
      WRITE (SCARC_FN, '(A,A,i3.3)') TRIM(CHID),'.scarc',MYID+1
      WRITE(*,*) 'Opening ', SCARC_FN
      LU_SCARC = GET_FILE_NUMBER()
      OPEN (LU_SCARC, FILE=SCARC_FN, ACTION = 'readwrite')
      WRITE(*,*) 'MYID=',MYID,': LU_SCARC=',LU_SCARC
      LASTID = MYID
   ENDDO
ENDIF
END SUBROUTINE SCARC_SETUP_DEBUGGING


!> ------------------------------------------------------------------------------------------------
!> Determine number of grid levels  (1 for CG/BICG-method, NLEVEL for MG-method)
!> Note: NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LEVEL
#ifdef WITH_MKL
INTEGER :: NL
#endif

SCARC_ROUTINE = 'SCARC_SETUP_LEVEL'
SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !> ----------------------- Krylov method --------------------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV: SELECT CASE (TYPE_PRECON)

#ifdef WITH_MKL
         CASE (NSCARC_PRECON_PARDISO)
            TYPE_MKL_LEVEL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
            IF (BTWOLEVEL) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF

         CASE (NSCARC_PRECON_CLUSTER)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            TYPE_MKL_LEVEL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
#endif

         CASE (NSCARC_PRECON_MULTIGRID)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)

         CASE DEFAULT
            IF (BTWOLEVEL) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
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
               CASE (NSCARC_SMOOTH_PARDISO)
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     TYPE_MKL_LEVEL(NL) = NSCARC_MKL_LOCAL
                  ENDDO
               CASE (NSCARC_SMOOTH_CLUSTER)
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     TYPE_MKL_LEVEL(NL) = NSCARC_MKL_GLOBAL
                  ENDDO
            END SELECT SELECT_SMOOTHER
      
            IF (TYPE_MKL == NSCARC_MKL_COARSE) TYPE_MKL_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif
      
         !> first, only finest level is set, further levels are defined during coarsening process
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_AMG)
      
      END SELECT SELECT_MULTIGRID
      
   !> ----------------------- MKL method -----------------------------------------
   CASE (NSCARC_METHOD_MKL)

      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)

      SELECT_MKL: SELECT CASE (TYPE_MKL)
         CASE (NSCARC_MKL_LOCAL)
            TYPE_MKL_LEVEL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
         CASE (NSCARC_MKL_GLOBAL)
            TYPE_MKL_LEVEL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
      END SELECT SELECT_MKL

END SELECT SELECT_METHOD

END SUBROUTINE SCARC_SETUP_LEVEL


!> ------------------------------------------------------------------------------------------------
!> Setup single level in case of default Krylov method
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS(NTYPE)
INTEGER, INTENT(IN) :: NTYPE
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM

SCARC_ROUTINE = 'SCARC_GET_NUMBER_OF_LEVELS'
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
      IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
         NLEVEL_MAX  = NLEVEL_MIN + SCARC_MULTIGRID_LEVEL - 1
      ELSE
         NLEVEL_MAX  = NLEVEL
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

SCARC_ROUTINE = 'SCARC_GET_MAX_LEVEL'

IF (BMULTILEVEL .AND.  MOD(NC,2)/=0) THEN
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

!> Divide by 2 as often as possible or until user defined max-level is reached
NC0=NC
DO NL=1,NSCARC_LEVEL_MAX
   NC0=NC0/2
   IF (MOD(NC0,2)/=0) EXIT                !> NC no longer divisable by two
   IF (NL==SCARC_MULTIGRID_LEVEL) EXIT    !> max number of levels defined by user
   IF (NC0==1) EXIT                       !> NC is power of two, minimum has been reached
ENDDO

SCARC_GET_MAX_LEVEL=NL
RETURN
1000 FORMAT(A,'=',I3,' must be divisable by two for Geometric Multigrid or Two-Level Method!')
END FUNCTION SCARC_GET_MAX_LEVEL

!> ------------------------------------------------------------------------------------------------
!> Allocate basic ScaRC-structures for all needed levels
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TYPES
INTEGER :: NM

SCARC_ROUTINE = 'SCARC_SETUP_TYPES'

!> Basic information for all requested grid levels
ALLOCATE (SCARC(NMESHES), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERROR)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   !> Detailed information for single grid levels
   ALLOCATE (SCARC(NM)%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'LEVEL', IERROR)

   !> Needed information about other meshes
   ALLOCATE (SCARC(NM)%OSCARC(NMESHES), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSCARC', IERROR)

ENDDO MESHES_LOOP
END SUBROUTINE SCARC_SETUP_TYPES


!> ------------------------------------------------------------------------------------------------
!> Setup debug file if requested
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TIMING
SCARC_ROUTINE = 'SCARC_SETUP_TIMING'
ALLOCATE (TSETUP(LOWER_MESH_INDEX:UPPER_MESH_INDEX), STAT = IERROR)
CALL CHKMEMERR (SCARC_ROUTINE, 'TSETUP', IERROR)
ALLOCATE (TSTEP(LOWER_MESH_INDEX:UPPER_MESH_INDEX), STAT = IERROR)
CALL CHKMEMERR (SCARC_ROUTINE, 'TSTEP', IERROR)
ALLOCATE (TSUM(LOWER_MESH_INDEX:UPPER_MESH_INDEX), STAT = IERROR)
CALL CHKMEMERR (SCARC_ROUTINE, 'TSUM', IERROR)
END SUBROUTINE SCARC_SETUP_TIMING


!> ----------------------------------------------------------------------------------------------------
!> Setup geometry information for mesh NM
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MESHES
INTEGER  :: NL, NM, IX, IY, IZ, IO
INTEGER  :: NX, NY, NZ
TYPE (MESH_TYPE)       , POINTER :: M
TYPE (SCARC_TYPE)      , POINTER :: S
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

IERROR=0
SCARC_ROUTINE = 'SCARC_SETUP_MESHES'

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

      L%N_EXTERNAL_WALL_CELLS = M%N_EXTERNAL_WALL_CELLS
      L%N_INTERNAL_WALL_CELLS = M%N_INTERNAL_WALL_CELLS

      L%NW = L%N_EXTERNAL_WALL_CELLS + L%N_INTERNAL_WALL_CELLS 

      !> Number of local cells per mesh
      ALLOCATE(L%NC_LOCAL(NMESHES), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','NC_LOCAL',IERROR)
      L%NC_LOCAL = 0

      !> Offset of cell numbers related to global numbering
      ALLOCATE(L%NC_OFFSET(NMESHES), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','NC_OFFSET',IERROR)
      L%NC_OFFSET = 0

      !> Number of local cells per mesh for DISCRET
      ALLOCATE(L%NUNKH_LOC(NMESHES), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','NUNKH_LOC',IERROR)
      L%NUNKH_LOC = 0

      !> Offset of cell numbers related to global numbering for DISCRET
      ALLOCATE(L%UNKH_IND(NMESHES), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','UNKH_IND',IERROR)
      L%UNKH_IND = 0

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

      !> step widths in x-, y- and z-direction for level 'NL'
      L%DX = (S%XF-S%XS)/REAL(L%NX,EB)
      L%DY = (S%YF-S%YS)/REAL(L%NY,EB)
      L%DZ = (S%ZF-S%ZS)/REAL(L%NZ,EB)

      L%DXI = 1.0_EB/L%DX
      L%DYI = 1.0_EB/L%DY
      L%DZI = 1.0_EB/L%DZ

      L%DXI2 = L%DXI**2
      L%DYI2 = L%DYI**2
      L%DZI2 = L%DZI**2

      !> needed in case of GMG with multiple grid levels
      NX=NX/2
      IF (.NOT.TWO_D) NY=NY/2
      NZ=NZ/2

      !> Allocate vectors for coordinate information
      IF (NL == NLEVEL_MIN) THEN

         L%XCOR => M%X
         L%YCOR => M%Y
         L%ZCOR => M%Z

         L%XMID => M%XC
         L%YMID => M%YC
         L%ZMID => M%ZC

      ELSE

         ALLOCATE(L%XCOR(0:L%NX), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','XCOR',IERROR)
         ALLOCATE(L%YCOR(0:L%NY), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','YCOR',IERROR)
         ALLOCATE(L%ZCOR(0:L%NZ), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','ZCOR',IERROR)

         !> compute coordinates in x-, y- and z-direction
         DO IX = 0, L%NX
            L%XCOR(IX) = S%XS + IX*L%DX
         ENDDO
         DO IY = 0, L%NY
            L%YCOR(IY) = S%YS + IY*L%DY
         ENDDO
         DO IZ = 0, L%NZ
            L%ZCOR(IZ) = S%ZS + IZ*L%DZ
         ENDDO

         !> compute midpoints in x-, y- and z-direction
         ALLOCATE(L%XMID(0:L%NX+1), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','XMID',IERROR)
         ALLOCATE(L%YMID(0:L%NY+1), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','YMID',IERROR)
         ALLOCATE(L%ZMID(0:L%NZ+1), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_MESHES','ZMID',IERROR)

         L%XMID(0) = S%XS - 0.5_EB*L%DX
         DO IX = 1, L%NX
            L%XMID(IX) = 0.5_EB*(L%XCOR(IX-1) + L%XCOR(IX))
         ENDDO
         L%XMID(L%NX+1) = S%XF + 0.5_EB*L%DX

         L%YMID(0) = S%YS - 0.5_EB*L%DY
         DO IY = 1, L%NY
            L%YMID(IY) = 0.5_EB*(L%YCOR(IY-1) + L%YCOR(IY))
         ENDDO
         L%YMID(L%NY+1) = S%YF + 0.5_EB*L%DY

         L%ZMID(0) = S%ZS - 0.5_EB*L%DZ
         DO IZ = 1, L%NZ
            L%ZMID(IZ) = 0.5_EB*(L%ZCOR(IZ-1) + L%ZCOR(IZ))
         ENDDO
         L%ZMID(L%NZ+1) = S%ZF + 0.5_EB*L%DZ

      ENDIF

      !> Allocate vectors for step sizes in different directions
      ALLOCATE(L%DXL(0:L%NX), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','DXL',IERROR)
      L%DXL = 0.0_EB
      ALLOCATE(L%DYL(0:L%NY), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','DYL',IERROR)
      L%DYL = 0.0_EB
      ALLOCATE(L%DZL(0:L%NZ), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_MESHES','DZL',IERROR)
      L%DZL = 0.0_EB

      !> set step sizes between cell midpoints, use interior step sizes for ghost cells as initial values
      !> correct sizes for ghost cells are exchanged later
      DO IX = 1, L%NX-1
         L%DXL(IX) = L%XMID(IX+1) - L%XMID(IX)
      ENDDO
      L%DXL(0)     = L%DXL(1)
      L%DXL(L%NX) = L%DXL(L%NX-1)

      DO IY = 1, L%NY-1
         L%DYL(IY) = L%YMID(IY+1) - L%YMID(IY)
      ENDDO
      L%DYL(0)     = L%DYL(1)
      L%DYL(L%NY) = L%DYL(L%NY-1)

      DO IZ = 1, L%NZ-1
         L%DZL(IZ) = L%ZMID(IZ+1) - L%ZMID(IZ)
      ENDDO
      L%DZL(0)     = L%DZL(1)
      L%DZL(L%NZ) = L%DZL(L%NZ-1)

   ENDDO LEVEL_LEVEL_LOOP
ENDDO LEVEL_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MESHES


!> ------------------------------------------------------------------------------------------------
!> Setup communication structure for data exchange along mesh interfaces
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERFACE
INTEGER :: NM, NOM, NL
TYPE (MESH_TYPE)       ,  POINTER :: M
TYPE (SCARC_TYPE)      ,  POINTER :: S
TYPE (OMESH_TYPE)      ,  POINTER :: OM
TYPE (OSCARC_TYPE)     ,  POINTER :: OS
TYPE (SCARC_LEVEL_TYPE),  POINTER :: F, C, OF, OC
TYPE (SCARC_MATRIX_TYPE), POINTER :: OFA, OCA

SCARC_ROUTINE = 'SCARC_SETUP_INTERFACE'

!> Initialize communication counter for ScaRC, use same TAG for all communications
N_REQ  =  0
N_EXCHANGE =  0
TAG   = 99

!> Allocate basic WALL and FACE types on mesh NM for all requested grid levels
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   S => SCARC(NM)

   F => S%LEVEL(NLEVEL_MIN)                          !> fine mesh level

   ALLOCATE(F%WALL(F%NW), STAT=IERROR)
   CALL ChkMemErr('SCARC_SETUP_INTERFACE','WALL',IERROR)

   ALLOCATE(F%FACE(-3:3), STAT=IERROR)
   CALL ChkMemErr('SCARC_SETUP_INTERFACE','FACE',IERROR)

   IF (NLEVEL_MAX > NLEVEL_MIN .AND. .NOT.BAMG) THEN
      DO NL=NLEVEL_MIN+1,NLEVEL_MAX
         C => S%LEVEL(NL)                            !> coarser grid levels
         ALLOCATE(C%FACE(-3:3), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_INTERFACE','WALL',IERROR)
      ENDDO
   ENDIF

   !> Get communication lengths for other meshes
   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES
      IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
      OM => M%OMESH(NOM)                              !> other mesh structure
      OS => S%OSCARC(NOM)                             !> other scarc structure
      OS%NIC_R    = OM%NIC_R
      OS%NIC_S    = OM%NIC_S
      OS%NICMAX_R = OM%NIC_R
      OS%NICMAX_S = OM%NIC_S
      IF (OS%NICMAX_R==0 .AND. OS%NICMAX_S==0)  CYCLE OTHER_MESHES_LOOP
      N_EXCHANGE  = N_EXCHANGE+1
   ENDDO OTHER_MESHES_LOOP
ENDDO MESHES_LOOP

!> Initialize level structures on neighboring meshes
LEVEL_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   S => SCARC(NM)

   LEVEL_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      IF (NOM == NM) CYCLE LEVEL_OTHER_MESHES_LOOP

      OS => S%OSCARC(NOM)
      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0)  CYCLE LEVEL_OTHER_MESHES_LOOP

      ALLOCATE (S%OSCARC(NOM)%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_INTERFACE', 'OS%LEVEL', IERROR)

      OF  => S%OSCARC(NOM)%LEVEL(NLEVEL_MIN)            ! other fine scarc structure
      OFA => S%OSCARC(NOM)%LEVEL(NLEVEL_MIN)%A          ! matrix of other fine scarc structure

      OF%NX = MESHES(NOM)%IBAR
      OF%NY = MESHES(NOM)%JBAR
      OF%NZ = MESHES(NOM)%KBAR

      OF%N_EXTERNAL_WALL_CELLS = MESHES(NOM)%N_EXTERNAL_WALL_CELLS
      OF%N_INTERNAL_WALL_CELLS = MESHES(NOM)%N_INTERNAL_WALL_CELLS

      OF%NW  = OF%N_EXTERNAL_WALL_CELLS + OF%N_INTERNAL_WALL_CELLS
      OF%NC  = OF%NX*OF%NY*OF%NZ
      OF%NCG = 0

      OFA%NAV = 0
      OFA%NAC = 0
      OFA%NAR = 0

      IF (S%OSCARC(NOM)%NICMAX_S == 0 .AND. S%OSCARC(NOM)%NICMAX_R == 0) CYCLE LEVEL_OTHER_MESHES_LOOP

      !> In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN .AND. .NOT.BAMG) THEN

         DO NL=NLEVEL_MIN+1,NLEVEL_MAX

            OF  => S%OSCARC(NOM)%LEVEL(NL-1)            !> fine level of other meshes
            OC  => S%OSCARC(NOM)%LEVEL(NL)              !> coarse level of other meshes
            OCA => S%OSCARC(NOM)%LEVEL(NL)%A            !> matrix on coarse level of other meshes

            OC%NX = OF%NX/2                          
            IF (TWO_D) THEN                            
               OC%NY = 1
            ELSE
               OC%NY = OF%NY/2
            ENDIF
            OC%NZ = OF%NZ/2

            OC%N_EXTERNAL_WALL_CELLS = 2*OC%NX * 2*OC%NZ + &
                                        2*OC%NX * 2*OC%NY + &
                                        2*OC%NY * 2*OC%NZ 

            OC%NC  = OC%NX * OC%NY * OC%NZ
            OC%NW  = OC%N_EXTERNAL_WALL_CELLS 
            OC%NCG = 0

            OCA%NAV = 0
            OCA%NAC = 0
            OCA%NAR = 0

         ENDDO
      ENDIF

   ENDDO LEVEL_OTHER_MESHES_LOOP
ENDDO LEVEL_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_INTERFACE


!> ------------------------------------------------------------------------------------------------
!> Initialize arrays for data exchange
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE
INTEGER :: NL, NM, NOM, NLEN_SEND, NLEN_RECV
INTEGER :: NLMIN, NLMAX
INTEGER :: INBR
TYPE ( SCARC_TYPE)     , POINTER :: S
TYPE (OSCARC_TYPE)     , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: L, OL

SCARC_ROUTINE = 'SCARC_SETUP_EXCHANGE'

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
   DO INBR = 1, S%NUM_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)

      L  => SCARC(NM)%LEVEL(NLEVEL_MIN)
      OS => SCARC(NM)%OSCARC(NOM)
      OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)

      !> ---------------- real buffers
      NLEN_SEND = NSCARC_MAX_STENCIL*MAX(OL%NWL, OL%NCG)+10
      NLEN_RECV = NLEN_SEND

      CALL SCARC_ALLOCATE_REAL1(OS%SEND_REAL, 1, NLEN_SEND, 'OS%SEND_REAL', .TRUE.)
      CALL SCARC_ALLOCATE_REAL1(OS%RECV_REAL, 1, NLEN_RECV, 'OS%RECV_REAL', .TRUE.)

      !> ---------------- integer buffers
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

      CALL SCARC_ALLOCATE_INT1(OS%SEND_INT, 1, NLEN_SEND, 'OS%SEND_INT', .TRUE.)
      CALL SCARC_ALLOCATE_INT1(OS%RECV_INT, 1, NLEN_RECV, 'OS%RECV_INT', .TRUE.)

      !> ---------- neighboring wall (and face?) structures for common wall cells
      ALLOCATE (OL%WALL(OL%NCG), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_INTERFACE', 'OL%WALL', IERROR)

      !ALLOCATE (OL%FACE(-3:3), STAT=IERROR)
      !CALL CHKMEMERR ('SCARC_SETUP_INTERFACE', 'OL%FACE', IERROR)

      !> In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN .AND. .NOT.BAMG) THEN
         DO NL=NLEVEL_MIN+1,NLEVEL_MAX
            OL => S%OSCARC(NOM)%LEVEL(NL)                          !> pointer to coarser level
            ALLOCATE (OL%WALL(OL%NCG), STAT=IERROR)
            CALL CHKMEMERR ('SCARC_SETUP_INTERFACE', 'OL%WALL', IERROR)
            !ALLOCATE (OL%FACE(-3:3), STAT=IERROR)
            !CALL CHKMEMERR ('SCARC_SETUP_INTERFACE', 'OL%FACE', IERROR)
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
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_DISCRET, NL)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WIDTHINFO, NL)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WALLINFO, NL)
   ENDDO
ENDIF

!>
!> Correct boundary types for cells adjacent to obstructions on ghost cells
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NLEVEL_MIN)
   IF (.NOT.BAMG) THEN
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX
         CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
      ENDDO
   ENDIF
ENDDO

END SUBROUTINE SCARC_SETUP_EXCHANGE


!> ----------------------------------------------------------------------------------------------------
!> Setup neighborship structures and boundary conditions
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLS
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES
INTEGER :: NM, NL
INTEGER :: IREFINE, IFACE, III
INTEGER :: INBR_FACE, INBR_MESH
INTEGER :: IWG, IWL, IWC
INTEGER :: NOM_LAST , NOM
INTEGER :: NCPL_LAST, NCPL
INTEGER :: IOR_LAST, IOR0
INTEGER :: FACE_NEIGHBORS(-3:3, NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: MESH_NEIGHBORS(NSCARC_MAX_MESH_NEIGHBORS)
INTEGER :: NUM_FACE_NEIGHBORS(-3:3)
INTEGER :: NUM_MESH_NEIGHBORS
LOGICAL :: BKNOWN
TYPE (MESH_TYPE)        , POINTER :: M
TYPE (SCARC_TYPE)       , POINTER :: S
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L, F, C
TYPE (SCARC_LEVEL_TYPE) , POINTER :: OL, OF, OC

SCARC_ROUTINE = 'SCARC_SETUP_WALLS'

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   S => SCARC(NM)

   FACE_NEIGHBORS = -1
   MESH_NEIGHBORS = -1

   NUM_FACE_NEIGHBORS = 0
   NUM_MESH_NEIGHBORS = 0

   !> For all solvers: Determine array WALL and PRESSURE_BC_INDEX on finest level 
   L => S%LEVEL(NLEVEL_MIN)

   L%NCE = L%NCS
   L%NCO = L%NCS

   L%N_EXTERNAL_WALL_CELLS = M%N_EXTERNAL_WALL_CELLS
   L%N_INTERNAL_WALL_CELLS = M%N_INTERNAL_WALL_CELLS

   !> Determine basic data for single faces (orientation, dimensions)
   FACES_OF_MESH_INDEX1: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE FACES_OF_MESH_INDEX1

      !> information about face orientation and local dimensions
      SELECT CASE (ABS(IOR0))
         CASE (1)
            L%FACE(IOR0)%NFC  =  L%NX                    !> number of cells between opposite mesh faces
            L%FACE(IOR0)%NFX  =  1                       !> number of cells in x-direction
            L%FACE(IOR0)%NFY  =  L%NY                    !> number of cells in y-direction
            L%FACE(IOR0)%NFZ  =  L%NZ                    !> number of cells in z-direction
            L%FACE(IOR0)%NFW  =  L%NY*L%NZ               !> number of wall cells at that face
            L%FACE(IOR0)%DH   => L%DXL                   !> step size vector between opposite mesh faces
         CASE (2)
            L%FACE(IOR0)%NFC  =  L%NY                    !> see above
            L%FACE(IOR0)%NFX  =  L%NX
            L%FACE(IOR0)%NFY  =  1
            L%FACE(IOR0)%NFZ  =  L%NZ
            L%FACE(IOR0)%NFW  =  L%NX*L%NZ
            L%FACE(IOR0)%DH   => L%DYL
         CASE (3)
            L%FACE(IOR0)%NFC  =  L%NZ                    !> see above
            L%FACE(IOR0)%NFX  =  L%NX
            L%FACE(IOR0)%NFY  =  L%NY
            L%FACE(IOR0)%NFZ  =  1
            L%FACE(IOR0)%NFW  =  L%NX*L%NY
            L%FACE(IOR0)%DH   => L%DZL
      END SELECT

   ENDDO FACES_OF_MESH_INDEX1

   !> Store local IWG-number for each face
   IWG = 1
   FACE_ORDER_LOOP: DO IFACE = 1, 6
      IOR0 = FACE_ORDER_XYZ(IFACE)
      L%FACE(IOR0)%IWG_PTR = IWG
      IWG = IWG + L%FACE(IOR0)%NFW
   ENDDO FACE_ORDER_LOOP

   !> L loop over global IW's:
   !> store basic data and determine number of adajacent neighbors to each
   !> face with corresponding number of IW's
   IOR_LAST  =  0
   NOM_LAST  = -1
   NCPL_LAST = -1
   IWL = 0

   !> L process external wall cells
   EXTERNAL_WALL_CELLS_LOOP1: DO IWG = 1, L%N_EXTERNAL_WALL_CELLS

      !> Determine and store neighbors, orientation and number of couplings for a single wall cell
      NOM  =  M%EXTERNAL_WALL(IWG)%NOM
      IOR0 =  M%WALL(IWG)%ONE_D%IOR
      NCPL = (M%EXTERNAL_WALL(IWG)%IIO_MAX - M%EXTERNAL_WALL(IWG)%IIO_MIN + 1) * &
             (M%EXTERNAL_WALL(IWG)%JJO_MAX - M%EXTERNAL_WALL(IWG)%JJO_MIN + 1) * &
             (M%EXTERNAL_WALL(IWG)%KKO_MAX - M%EXTERNAL_WALL(IWG)%KKO_MIN + 1)

      L%WALL(IWG)%NOM  = NOM                            !> store number of neighbor in wall cell
      L%WALL(IWG)%IOR  = IOR0                           !> store orientation of that cell

      IWL = IWL + 1                                     !> count local wall cells for that face

      IF (NOM /= 0) THEN
         L%WALL(IWG)%NCPL = NCPL                        !> store number of couplings for that cell
         BKNOWN = .FALSE.
         DO INBR_FACE = 1, NUM_FACE_NEIGHBORS(IOR0)
            IF (FACE_NEIGHBORS(IOR0, INBR_FACE) == NOM) THEN
               BKNOWN = .TRUE.
               EXIT
            ENDIF
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

         IF (BKNOWN) THEN
            OL%NWL = OL%NWL + 1                         !> increase own counter for local wall cells
            OL%NCG = OL%NCG + NCPL                      !> increase counter for local ghost cells
         ELSE
            OL%NWL = 1                                  !> initialize own counter for local wall cells
            OL%NCG = NCPL                               !> initialize counter for local ghost cells
            L%NCPL_MAX  = MAX(L%NCPL_MAX, NCPL)         !> get max NCPL ever used on this mesh
         ENDIF
      ENDIF

      IF (NOM /= 0 .AND. .NOT.BKNOWN) THEN
         NUM_FACE_NEIGHBORS(IOR0) = NUM_FACE_NEIGHBORS(IOR0) + 1     !> increase neighbor counter for face
         NUM_MESH_NEIGHBORS       = NUM_MESH_NEIGHBORS       + 1     !> increase neighbor counter for mesh
         FACE_NEIGHBORS(IOR0, NUM_FACE_NEIGHBORS(IOR0)) = NOM        !> store number of neighbor for face
         MESH_NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0))       = NOM        !> store number of neighbor for mesh
      ENDIF

      IOR_LAST  = IOR0                                               !> save former values
      NCPL_LAST = NCPL
      NOM_LAST  = NOM

   ENDDO EXTERNAL_WALL_CELLS_LOOP1
 
   !> Then process internal wall cells
   INTERNAL_WALL_CELLS_LOOP1: DO IWG = L%N_EXTERNAL_WALL_CELLS+1, L%N_EXTERNAL_WALL_CELLS+L%N_INTERNAL_WALL_CELLS

      L%WALL(IWG)%IOR  = M%WALL(IWG)%ONE_D%IOR
      L%WALL(IWG)%NOM  = 0

      L%WALL(IWG)%BTYPE  = NEUMANN
      L%WALL(IWG)%BOUNDARY_TYPE  = M%WALL(IWG)%BOUNDARY_TYPE

      L%WALL(IWG)%IXG =  M%WALL(IWG)%ONE_D%II                        !> ghost cell indices
      L%WALL(IWG)%IYG =  M%WALL(IWG)%ONE_D%JJ
      L%WALL(IWG)%IZG =  M%WALL(IWG)%ONE_D%KK

      L%WALL(IWG)%IXW =  M%WALL(IWG)%ONE_D%IIG                       !> (internal) wall cell indices
      L%WALL(IWG)%IYW =  M%WALL(IWG)%ONE_D%JJG
      L%WALL(IWG)%IZW =  M%WALL(IWG)%ONE_D%KKG

   ENDDO INTERNAL_WALL_CELLS_LOOP1

   !> Allocate array which stores numbers of all neighboring meshes
   IF (NUM_MESH_NEIGHBORS /= 0) THEN
      ALLOCATE(S%NEIGHBORS(NUM_MESH_NEIGHBORS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_WALLS','MESH_NEIGHBORS',IERROR)
   ENDIF
   S%NUM_NEIGHBORS = NUM_MESH_NEIGHBORS

   !> Store information about adjacent neighbors on different faces
   !> Allocate corresponding index arrays in OSCARC-structures
   !> First allocate administrative mapping arrays for own mesh
   CALL SCARC_SETUP_MAPPINGS(NM, NLEVEL_MIN)

   INBR_MESH = 1
   NEIGHBORS_OF_FACE_LOOP: DO IOR0 = -3, 3
      IF (IOR0 == 0) CYCLE NEIGHBORS_OF_FACE_LOOP

      !> if there are neighbors at face IOR0 store information about them
      IF (NUM_FACE_NEIGHBORS(IOR0) /= 0) THEN

         L%FACE(IOR0)%NUM_NEIGHBORS = NUM_FACE_NEIGHBORS(IOR0)        !> store number of neighbors on face

         !> allocate array for storing the numbers of the single neighbors
         ALLOCATE(L%FACE(IOR0)%NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0)), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_WALLS','NOMS',IERROR)

         !> store every neighbor and allocate corresponding administration arrays
         DO INBR_FACE = 1, NUM_FACE_NEIGHBORS(IOR0)

            NOM = FACE_NEIGHBORS(IOR0, INBR_FACE)

            L%FACE(IOR0)%NEIGHBORS(INBR_FACE) = NOM                   !> store NOM as a neighbor of that face
            S%NEIGHBORS(INBR_MESH)            = NOM                   !> store NOM as a neighbor of that mesh
            INBR_MESH = INBR_MESH + 1

            !> allocate administrative arrays for neighboring meshes
            CALL SCARC_SETUP_OMAPPINGS(NM, NOM, NLEVEL_MIN)

         ENDDO
      ENDIF
   ENDDO NEIGHBORS_OF_FACE_LOOP

   !> Second loop over external wall cells:
   !> Store detailed coordinate and cell data and get type of boundary condition
   L%ICE_PTR = L%NCS
   L%ICO_PTR = L%NCS

   IOR_LAST = 0
   WALL_CELLS_LOOP2: DO IWG = 1, L%N_EXTERNAL_WALL_CELLS

      !> Determine neighbors and orientation again
      NOM  = L%WALL(IWG)%NOM
      IOR0 = L%WALL(IWG)%IOR
      NCPL = L%WALL(IWG)%NCPL

      !> Determine boundary type for IW
      IF (M%EXTERNAL_WALL(IWG)%NOM /= 0) THEN
         L%WALL(IWG)%BTYPE = INTERNAL
      ELSE IF (M%WALL(IWG)%PRESSURE_BC_INDEX == DIRICHLET) THEN
         L%WALL(IWG)%BTYPE = DIRICHLET
      ELSE
         L%WALL(IWG)%BTYPE = NEUMANN
      ENDIF

      L%WALL(IWG)%BOUNDARY_TYPE  = M%WALL(IWG)%BOUNDARY_TYPE

      L%WALL(IWG)%IXG =  M%WALL(IWG)%ONE_D%II                                 !> ghost cell indices
      L%WALL(IWG)%IYG =  M%WALL(IWG)%ONE_D%JJ
      L%WALL(IWG)%IZG =  M%WALL(IWG)%ONE_D%KK

      L%WALL(IWG)%IXW =  M%WALL(IWG)%ONE_D%IIG                                !> (internal) wall cell indices
      L%WALL(IWG)%IYW =  M%WALL(IWG)%ONE_D%JJG
      L%WALL(IWG)%IZW =  M%WALL(IWG)%ONE_D%KKG

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

!> 
!> Set DISCRET information on finest level and if requested also on coarser levels
!> 
CALL SCARC_SETUP_DISCRET_GLOBAL(NLEVEL_MIN)
IF (.NOT.BAMG) THEN
   DO NL = NLEVEL_MIN+1, NLEVEL_MAX
      CALL SCARC_SETUP_DISCRET_LEVEL(NL)
      CALL SCARC_SETUP_DISCRET_GLOBAL(NL)
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

         F => S%LEVEL(NL-1)                   !> fine level
         C => S%LEVEL(NL)                     !> coarse level

         IREFINE=IREFINE*2

         CALL SCARC_CHECK_DIVISIBILITY(F%NCE-F%NCS, 'F%NCE')
         CALL SCARC_CHECK_DIVISIBILITY(F%NCO-F%NCS, 'F%NCO')

         C%NCE = C%NCS + (F%NCE-F%NCS)/2
         C%NCO = C%NCS + (F%NCO-F%NCS)/2

         C%ICE_PTR = C%NCS
         C%ICO_PTR = C%NCS

         C%N_EXTERNAL_WALL_CELLS = SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
         C%N_INTERNAL_WALL_CELLS = SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)

         C%NW = C%N_EXTERNAL_WALL_CELLS + C%N_INTERNAL_WALL_CELLS
         ALLOCATE(C%WALL(C%NW), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_INTERFACE','WALL',IERROR)

         CALL SCARC_SETUP_MAPPINGS(NM, NL)

         !> compute FACE and WALL information for all faces of coarser level
         IWC = 1
         IWG = 1
         DO IFACE = 1, 6

            IOR0 = FACE_ORDER_XYZ(IFACE)

            !> compute mesh dimensions of coarser mesh level
            CALL SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)

            !> for every neighbor do:
            IF (F%FACE(IOR0)%NUM_NEIGHBORS /= 0) THEN
               DO INBR_FACE = 1, F%FACE(IOR0)%NUM_NEIGHBORS

                  NOM = F%FACE(IOR0)%NEIGHBORS(INBR_FACE)

                  OC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  OF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL-1)

                  OC%NCPL = OF%NCPL

                  CALL SCARC_CHECK_DIVISIBILITY(OF%NWL, 'OF%NWL')
                  CALL SCARC_CHECK_DIVISIBILITY(OF%NCG, 'OF%NCG')

                  IF (.NOT.TWO_D) THEN
                     OC%NWL = OF%NWL/4
                     OC%NCG = OF%NCG/4
                  ELSE
                     OC%NWL = OF%NWL/2
                     OC%NCG = OF%NCG/2
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

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALLINFO , NLEVEL_MIN,  'WALL')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACEINFO , NLEVEL_MIN,  'FACE')

END SUBROUTINE SCARC_SETUP_WALLS


!> -----------------------------------------------------------------------------------------
!> --- Setup CELL_INDEX array on coarser grid levels in case of MG-method
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CELL_INDEX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: I, J, K, IZERO, NO
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_OBST_TYPE) , POINTER :: OB

SCARC_ROUTINE = 'SCARC_SETUP_CELL_INDEX'

L => SCARC(NM)%LEVEL(NL)

ALLOCATE(L%CELL_INDEX(0:L%NX+1,0:L%NY+1,0:L%NZ+1),STAT=IZERO)
CALL ChkMemErr('SCARC_SETUP_CELL_INDEX','CELL_INDEX',IZERO)  

L%CELL_INDEX = 0
L%N_CELLS = 0

!>
!> Preset it for all grid cells
!>
DO K=0,L%NZ+1
   DO J=0,L%NY+1
      DO I=0,1
         IF (L%CELL_INDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL_INDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
      DO I=L%NX,L%NX+1
         IF (L%CELL_INDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL_INDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO K=0,L%NZ+1
   DO I=0,L%NX+1
      DO J=0,1
         IF (L%CELL_INDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL_INDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
      DO J=L%NY,L%NY+1
         IF (L%CELL_INDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL_INDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO J=0,L%NY+1
   DO I=0,L%NX+1
      DO K=0,1
         IF (L%CELL_INDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL_INDEX(I,J,K) = L%N_CELLS
         ENDIF
      ENDDO
      DO K=L%NZ,L%NZ+1
         IF (L%CELL_INDEX(I,J,K)==0) THEN
            L%N_CELLS = L%N_CELLS + 1
            L%CELL_INDEX(I,J,K) = L%N_CELLS
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
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELLS = L%N_CELLS + 1
               L%CELL_INDEX(I,J,K) = L%N_CELLS
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_CELL_INDEX


!> -----------------------------------------------------------------------------------------
!> --- Setup WALL_INDEX array on coarser grid levels in case of MG-method
!> --- corresponding to cell index information
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_INDEX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: I, J, K, ICG, IZERO, IW, IOR0
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

SCARC_ROUTINE = 'SCARC_SETUP_WALL_INDEX'
L => SCARC(NM)%LEVEL(NL)

ALLOCATE(L%WALL_INDEX(L%N_CELLS, -3:3), STAT=IZERO) 
CALL ChkMemErr('SCARC_SETUP_WALL_INDEX','WALL_INDEX',IZERO)  
L%WALL_INDEX = 0

DO IW = 1, L%NW

   I = L%WALL(IW)%IXW
   J = L%WALL(IW)%IYW
   K = L%WALL(IW)%IZW

   IOR0 = L%WALL(IW)%IOR
   ICG  = L%CELL_INDEX(I,J,K)
  
   L%WALL_INDEX(ICG,-IOR0) = IW

ENDDO

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
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: C, F
TYPE (SCARC_OBST_TYPE), POINTER :: OB

SCARC_ROUTINE = 'SCARC_SETUP_INTERNAL_WALL_COORDS'
M   => MESHES(NM)

F  => SCARC(NM)%LEVEL(NL-1)
C  => SCARC(NM)%LEVEL(NL)

IWC = C%N_EXTERNAL_WALL_CELLS + 1

DO IO = 1, C%N_OBST

   OB => C%OBST(IO)

   !> Analyze IOR = 1
   I = OB%I1
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = C%DISCRET(I, J, K, UNKH)
         !IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH)  WRITE(LU_SCARC,*) 'A: Processing ', IC, I,J,K, 1, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         C%WALL(IWC)%IXW = I+1
         C%WALL(IWC)%IYW = J
         C%WALL(IWC)%IZW = K
         C%WALL(IWC)%IXG = I
         C%WALL(IWC)%IYG = J
         C%WALL(IWC)%IZG = K
         C%WALL(IWC)%IOR = 1
         C%WALL(IWC)%BTYPE = NEUMANN
         C%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = C%DISCRET(I+1, J, K, UNKH)
         !IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'B: Processing ', IC, I,J,K,-1, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         C%WALL(IWC)%IXW = I
         C%WALL(IWC)%IYW = J
         C%WALL(IWC)%IZW = K
         C%WALL(IWC)%IXG = I+1
         C%WALL(IWC)%IYG = J
         C%WALL(IWC)%IZG = K
         C%WALL(IWC)%IOR =-1
         C%WALL(IWC)%BTYPE = NEUMANN
         C%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 2
   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = C%DISCRET(I, J, K, UNKH)
         !IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'C: Processing ', IC, I,J,K,2, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         C%WALL(IWC)%IXW = I
         C%WALL(IWC)%IYW = J+1
         C%WALL(IWC)%IZW = K
         C%WALL(IWC)%IXG = I
         C%WALL(IWC)%IYG = J
         C%WALL(IWC)%IZG = K
         C%WALL(IWC)%IOR = 2
         C%WALL(IWC)%BTYPE = NEUMANN
         C%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -2
   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = C%DISCRET(I, J+1, K, UNKH)
         !IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'D: Processing ', IC, I,J,K,-2, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         C%WALL(IWC)%IXW = I
         C%WALL(IWC)%IYW = J
         C%WALL(IWC)%IZW = K
         C%WALL(IWC)%IXG = I
         C%WALL(IWC)%IYG = J+1
         C%WALL(IWC)%IZG = K
         C%WALL(IWC)%IOR =-2
         C%WALL(IWC)%BTYPE = NEUMANN
         C%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 3
   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = C%DISCRET(I, J, K, UNKH)
         !IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'E: Processing ', IC, I,J,K,3, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         C%WALL(IWC)%IXW = I
         C%WALL(IWC)%IYW = J
         C%WALL(IWC)%IZW = K+1
         C%WALL(IWC)%IXG = I
         C%WALL(IWC)%IYG = J
         C%WALL(IWC)%IZG = K
         C%WALL(IWC)%IOR = 3
         C%WALL(IWC)%BTYPE = NEUMANN
         C%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -3
   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = C%DISCRET(I, J, K+1, UNKH)
         !IF (TYPE_DEBUG >= NSCARC_DEBUG_MUCH) WRITE(LU_SCARC,*) 'F: Processing ', IC, I,J,K,-3, IWC
         IF (IC == IS_UNDEFINED) CYCLE
         C%WALL(IWC)%IXW = I
         C%WALL(IWC)%IYW = J
         C%WALL(IWC)%IZW = K
         C%WALL(IWC)%IXG = I
         C%WALL(IWC)%IYG = J
         C%WALL(IWC)%IZG = K+1
         C%WALL(IWC)%IOR =-3
         C%WALL(IWC)%BTYPE = NEUMANN
         C%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO
         
ENDDO

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
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

SCARC_ROUTINE = 'SCARC_SETUP_INTERNAL_NEUMANNS'
M => MESHES(NM)
L => SCARC(NM)%LEVEL(NL)

!IF (NL == NLEVEL_MIN) THEN
!   CELL_INDEX => M%CELL_INDEX
!   WALL_INDEX => M%WALL_INDEX
!ELSE
!   CELL_INDEX => L%CELL_INDEX
!   WALL_INDEX => L%WALL_INDEX
!ENDIF

DO IWG = 1, L%N_EXTERNAL_WALL_CELLS

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
  ! IF (L%DISCRET(IX, IY, IZ, CGSC) == IS_SOLID) L%WALL(IWG)%BTYPE=NEUMANN

   IF (L%WALL(IWG)%BOUNDARY_TYPE == SOLID_BOUNDARY) L%WALL(IWG)%BTYPE=NEUMANN

ENDDO

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
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: F, C

SCARC_ROUTINE = 'SCARC_SETUP_EXTERNAL_WALL_CELLS'

M   => MESHES(NM)
C => SCARC(NM)%LEVEL(NL)
F => SCARC(NM)%LEVEL(NL-1)

IF (NL == NLEVEL_MIN+1) THEN
   CELL_INDEX => M%CELL_INDEX
   WALL_INDEX => M%WALL_INDEX
ELSE
   CELL_INDEX => F%CELL_INDEX
   WALL_INDEX => F%WALL_INDEX
ENDIF

IWC = 0
ICF = 0
IWF = 0

IF (TWO_D) THEN
   IYC = 1
   IYF = 1

   !> IOR = 1
   IOR0 = 1
   DO IZC = 1, C%NZ 
      IZF = 2*IZC - 1
      ICF(1) = CELL_INDEX(1  , IYF  , IZF  )
      ICF(2) = CELL_INDEX(1  , IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
   ENDDO

   !> IOR = -1
   IOR0 = -1
   DO IZC = 1, C%NZ 
      IZF = 2*IZC - 1
      ICF(1) = CELL_INDEX(F%NX, IYF  , IZF  )
      ICF(2) = CELL_INDEX(F%NX, IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
   ENDDO

   !> IOR = 2
   IOR0 = 2
   DO IZC = 1, C%NZ 
      IZF = 2*IZC - 1
      DO IXC = 1, C%NX 
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
   DO IZC = 1, C%NZ 
      IZF = 2*IZC - 1
      DO IXC = 1, C%NX 
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
   DO IXC = 1, C%NX 
      IXF = 2*IXC - 1
      ICF(1) = CELL_INDEX(IXF    , IYF  , 1)
      ICF(2) = CELL_INDEX(IXF+1  , IYF  , 1)
      IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
   ENDDO

   !> IOR = -3
   IOR0 = -3
   DO IXC = 1, C%NX 
      IXF = 2*IXC - 1
      ICF(1) = CELL_INDEX(IXF  , IYF  , F%NZ)
      ICF(2) = CELL_INDEX(IXF+1, IYF  , F%NZ)
      IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 2, NM, NL-1)) IWC = IWC + 1
   ENDDO

ELSE

   !> IOR = 1
   IOR0 = 1
   DO IZC = 1, C%NZ 
      IZF = 2*IZC - 1
      DO IYC = 1, C%NY 
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
   DO IZC = 1, C%NZ 
      IZF = 2*IZC - 1
      DO IYC = 1, C%NY 
         IYF = 2*IYC - 1
         ICF(1) = CELL_INDEX(F%NX, IYF  , IZF  )
         ICF(2) = CELL_INDEX(F%NX, IYF+1, IZF  )
         ICF(3) = CELL_INDEX(F%NX, IYF  , IZF+1)
         ICF(4) = CELL_INDEX(F%NX, IYF+1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = 2
   IOR0 = 2
   DO IZC = 1, C%NZ 
      IZF = 2*IZC - 1
      DO IXC = 1, C%NX 
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
   DO IZC = 1, C%NZ 
      IZF = 2*IZC - 1
      DO IXC = 1, C%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF    , F%NY, IZF  )
         ICF(2) = CELL_INDEX(IXF+1  , F%NY, IZF  )
         ICF(3) = CELL_INDEX(IXF    , F%NY, IZF+1)
         ICF(4) = CELL_INDEX(IXF+1  , F%NY, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = 3
   IOR0 = 3
   DO IYC = 1, C%NY 
      IYF = 2*IYC - 1
      DO IXC = 1, C%NX 
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
   DO IYC = 1, C%NY 
      IYF = 2*IYC - 1
      DO IXC = 1, C%NX 
         IXF = 2*IXC - 1
         ICF(1) = CELL_INDEX(IXF  , IYF  , F%NZ)
         ICF(2) = CELL_INDEX(IXF+1, IYF  , F%NZ)
         ICF(3) = CELL_INDEX(IXF  , IYF+1, F%NZ)
         ICF(4) = CELL_INDEX(IXF+1, IYF+1, F%NZ)
         IF (IS_EXTERNAL_WALLCELL(WALL_INDEX, IOR0, ICF, 4, NM, NL-1)) IWC = IWC + 1
      ENDDO
   ENDDO

ENDIF

SCARC_COUNT_EXTERNAL_WALL_CELLS = IWC
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
TYPE (SCARC_LEVEL_TYPE), POINTER :: F, C
TYPE (SCARC_OBST_TYPE),  POINTER :: OB

F => SCARC(NM)%LEVEL(NL-1)
C => SCARC(NM)%LEVEL(NL)

C%N_OBST = F%N_OBST                   !> Number of obstructions is the same on all levels

ALLOCATE(C%OBST(C%N_OBST), STAT=IERROR)
CALL ChkMemErr('SCARC_COUNT_INTERNAL_WALL_CELLS','OBST',IERROR)

NW_INT = 0
IWC = C%N_EXTERNAL_WALL_CELLS + 1

DO IO = 1, F%N_OBST

   OB => C%OBST(IO)                    !> obstruction pointer on coarse level

   OB%I1 = (F%OBST(IO)%I1+1)/2
   OB%I2 =  F%OBST(IO)%I2/2

   IF (TWO_D) THEN
      OB%J1 = 0
      OB%J2 = 1
   ELSE
      OB%J1 = (F%OBST(IO)%J1+1)/2
      OB%J2 =  F%OBST(IO)%J2/2
   ENDIF

   OB%K1 = (F%OBST(IO)%K1+1)/2
   OB%K2 =  F%OBST(IO)%K2/2
   
   !> Analyze IOR = 1
   I = OB%I1
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = C%DISCRET(I, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = C%DISCRET(I+1, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 2
   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = C%DISCRET(I, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -2
   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = C%DISCRET(I, J+1, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 3
   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = C%DISCRET(I, J, K, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -3
   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = C%DISCRET(I, J, K+1, UNKH)
         IF (IC == IS_UNDEFINED) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO
         
ENDDO

SCARC_COUNT_INTERNAL_WALL_CELLS = NW_INT

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

IS_EXTERNAL_WALLCELL = .FALSE.

L => SCARC(NM)%LEVEL(NL)

DO I = 1, NCNT
   IWF(I) = WALL_INDEX(ICF(I), -IOR0)
ENDDO

BSUM = 0.0_EB
IWF_LAST = 0

DO I = 1, NCNT
   IF (IWF(I)>0) THEN
      BSUM = BSUM + REAL(L%WALL(IWF(I))%BOUNDARY_TYPE,EB)
      IWF_LAST = IWF(I)
   ENDIF
ENDDO

IF (IWF_LAST == 0) RETURN

IF (ABS(BSUM/REAL(NCNT,EB) - REAL(L%WALL(IWF_LAST)%BOUNDARY_TYPE,EB)) < 1E-12) THEN
   IS_EXTERNAL_WALLCELL = .TRUE.
   RETURN
ELSE
   WRITE(*,*) 'Wrong sum of BOUNDARY_TYPE for IOR =  1'
   STOP
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
TYPE (SCARC_MAPPING_TYPE), POINTER :: MAP, OMAP

SCARC_ROUTINE = 'SCARC_SETUP_WALLCELL_NEIGHBOR'

L    => SCARC(NM)%LEVEL(NL)
MAP  => L%MAP

OL   => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
OMAP => OL%MAP

!> store basic information about neighbor and orientation
OL%IOR  = IOR0
OL%NOM  = NOM

ICO  = L%ICO_PTR
ICE  = L%ICE_PTR

ICG  = OL%ICG_PTR
IWL  = OL%IWL_PTR

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
ALLOCATE(L%WALL(IWG)%ICE(OL%NCPL), STAT=IERROR)
CALL ChkMemErr('SCARC_SETUP_WALLS','ICE',IERROR)
L%WALL(IWG)%ICE = NSCARC_UNDEFINED

ALLOCATE(L%WALL(IWG)%ICG(OL%NCPL), STAT=IERROR)
CALL ChkMemErr('SCARC_SETUP_WALLS','ICG',IERROR)
L%WALL(IWG)%ICG = NSCARC_UNDEFINED

IWL = IWL + 1
ICO = ICO + 1

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

         L%WALL(IWG)%ICE(ICPL) = ICE                      !> number of extended grid cell
         L%WALL(IWG)%ICG(ICPL) = ICG                      !> number of ghost grid cell
   
         MAP%ICE_TO_IWG(ICE)  = IWG                       !> map extended cell to global wall cell
         MAP%ICE_TO_IWL(ICE)  = IWL                       !> map extended cell to local wall cell
         MAP%ICE_TO_ICG(ICE)  = ICG                       !> map extended cell to ghost cell

         OMAP%ICG_TO_IWG(ICG) = IWG                       !> map ghost cell to global wall cell
         OMAP%ICG_TO_ICO(ICG) = ICO                       !> map ghost cell to extended grid cell
         OMAP%ICG_TO_ICE(ICG) = L%WALL(IWG)%ICE(ICPL)     !> map ghost cell to extended grid cell

      ENDDO
   ENDDO
ENDDO

IF (BAMG) THEN
   DO IZ = NZ1, NZ2
      DO IY = NY1, NY2
         DO IX = NX1, NX2
            OMAP%IWL_TO_ICG(IWL, ICPL) = ICG              !> map extended cell to ghost cell
         ENDDO
      ENDDO
   ENDDO
ENDIF

L%WALL(IWG)%ICO = ICO                                     !> number of overlapping cell
L%WALL(IWG)%IWL = IWL                                     !> number of local wall cell

L%ICO_PTR  = ICO                                          !> store overlapping cell counter
L%ICE_PTR  = ICE                                          !> store extended cell counter

OL%IWL_PTR = IWL                                          !> store local wall cell pointer
OL%ICG_PTR = ICG                                          !> store ghost cell counter

OMAP%IWL_TO_IWG(IWL) = IWG                                !> map local wall cell to global wall cell
OMAP%IWL_TO_ICO(IWL) = ICO                                !> map local wall cell to internal grid cell
OMAP%IWL_TO_ICW(IWL) = L%WALL(IWG)%ICW                    !> map local wall cell to internal grid cell

END SUBROUTINE SCARC_SETUP_WALLCELL_NEIGHBOR


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MAPPINGS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: NCS, NCO, NCE, NW
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_MAPPING_TYPE), POINTER :: MAP

SCARC_ROUTINE = 'SCARC_SETUP_MAPPINGS'

L   => SCARC(NM)%LEVEL(NL)
MAP => SCARC(NM)%LEVEL(NL)%MAP

NCS = L%NCS
NCE = L%NCE
NCO = L%NCO
NW  = L%NW

!> Allocate mappings for ICE
ALLOCATE(MAP%ICE_TO_IWG(NCS+1:NCE), STAT=IERROR)
CALL ChkMemErr(SCARC_ROUTINE,'ICE_TO_IWG',IERROR)
MAP%ICE_TO_IWG=0

ALLOCATE(MAP%ICE_TO_IWL(NCS+1:NCE), STAT=IERROR)
CALL ChkMemErr(SCARC_ROUTINE,'ICE_TO_IWL',IERROR)
MAP%ICE_TO_IWL=0

ALLOCATE(MAP%ICE_TO_ICG(NCS+1:NCE), STAT=IERROR)
CALL ChkMemErr(SCARC_ROUTINE,'ICE_TO_ICG',IERROR)
MAP%ICE_TO_ICG=0

!> Allocate remaining mappings (outdated)
ALLOCATE(L%INTERNAL_BDRY_CELL(NCS), STAT=IERROR)
CALL ChkMemErr(SCARC_ROUTINE,'INTERNAL_BDRY_CELL',IERROR)
L%INTERNAL_BDRY_CELL=0

END SUBROUTINE SCARC_SETUP_MAPPINGS


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OMAPPINGS(NM, NOM, NL)
INTEGER, INTENT(IN) :: NM, NOM, NL
INTEGER:: NWL, NCG, NCPL
TYPE (SCARC_LEVEL_TYPE), POINTER :: OL
TYPE (SCARC_MAPPING_TYPE), POINTER :: OMAP

SCARC_ROUTINE = 'SCARC_SETUP_OMAPPINGS'

!> Get corresponding dimensions
OL   => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
OMAP => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%MAP

NWL  = OL%NWL
NCG  = OL%NCG
NCPL = OL%NCPL

!> initialize some counters
OL%IWL = 0
OL%ICG = 0

!> Allocate mapping array for local to global wall cells
ALLOCATE(OMAP%IWL_TO_IWG(NWL), STAT=IERROR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','IWL_TO_IWG',IERROR)
OMAP%IWL_TO_IWG = 0

ALLOCATE(OMAP%IWL_TO_ICW(NWL), STAT=IERROR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','IWL_TO_ICW',IERROR)
OMAP%IWL_TO_ICW = 0

ALLOCATE(OMAP%IWL_TO_ICO(NWL), STAT=IERROR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','IWL_TO_ICO',IERROR)
OMAP%IWL_TO_ICO = 0

IF (BAMG) THEN
   ALLOCATE(OMAP%IWL_TO_ICG(NWL, NCPL), STAT=IERROR)
   CALL ChkMemErr('SCARC_SETUP_MAPPINGS','L%IWL_TO_ICG',IERROR)
   OMAP%IWL_TO_ICG = 0
ENDIF

IF (PRES_ON_WHOLE_DOMAIN.AND.OL%NCS>0) THEN
   ALLOCATE(OMAP%ICN_TO_ICE(1:OL%NCS), STAT=IERROR)             !> ACHTUNG: kann verbessert werden??? riesig!
   CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICN_TO_ICE',IERROR)
   OMAP%ICN_TO_ICE = 0
ENDIF


!> Allocate mapping arrays for ghost cells
ALLOCATE(OMAP%ICG_TO_IWG(NCG), STAT=IERROR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_IWG',IERROR)
OMAP%ICG_TO_IWG = 0

ALLOCATE(OMAP%ICG_TO_ICO(NCG), STAT=IERROR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_ICO',IERROR)
OMAP%ICG_TO_ICO = 0

ALLOCATE(OMAP%ICG_TO_ICE(NCG), STAT=IERROR)
CALL ChkMemErr('SCARC_SETUP_OMAPPINGS','ICG_TO_ICE',IERROR)
OMAP%ICG_TO_ICE = 0

END SUBROUTINE SCARC_SETUP_OMAPPINGS


!> ----------------------------------------------------------------------------------------------------
!> Check divisibility by 2 of a given number of elements (in one grid direction)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_DIVISIBILITY(NN, CDIR)
INTEGER, INTENT(IN) :: NN
CHARACTER(*) , INTENT(IN) :: CDIR
IF (MOD(NN,2) /= 0) THEN
   WRITE(SCARC_MESSAGE,'(4A,I8,A)') TRIM(SCARC_ROUTINE),': Parameter ',CDIR,' = ',NN,' not divisable by 2'
   CALL SHUTDOWN(SCARC_MESSAGE); RETURN
ENDIF
END SUBROUTINE SCARC_CHECK_DIVISIBILITY


!> ----------------------------------------------------------------------------------------------------
!> Set wall cell information on coarse level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)
INTEGER, INTENT(IN)    :: IOR0, NM, NL
INTEGER, INTENT(INOUT) :: IWG
INTEGER:: INBR_FACE
TYPE (SCARC_LEVEL_TYPE), POINTER :: C, F            !> LEVEL types for coarse and fine grid
TYPE (SCARC_FACE_TYPE) , POINTER :: FC, FF          !> FACE types on coarse and fine levels

SCARC_ROUTINE = 'SCARC_SETUP_FACE_DIMENSIONS'

!> reference coarse and fine LEVEL type
C => SCARC(NM)%LEVEL(NL)
F => SCARC(NM)%LEVEL(NL-1)

!> reference coarse and fine FACE type
FC => SCARC(NM)%LEVEL(NL)%FACE(IOR0)
FF => SCARC(NM)%LEVEL(NL-1)%FACE(IOR0)

!> initialize FACE type for coarser mesh
FC%IWG_PTR = IWG

FC%NUM_NEIGHBORS = FF%NUM_NEIGHBORS
IF (FC%NUM_NEIGHBORS /= 0) THEN
   ALLOCATE(FC%NEIGHBORS(FC%NUM_NEIGHBORS), STAT=IERROR)
   CALL ChkMemErr('SCARC_FACE_DIMENSIONS','FC%NEIGHBORS',IERROR)
ENDIF
DO INBR_FACE= 1, FC%NUM_NEIGHBORS
   FC%NEIGHBORS(INBR_FACE) = FF%NEIGHBORS(INBR_FACE)
ENDDO

SELECT CASE (ABS(IOR0))

   CASE (1)

      FC%DH => C%DXL
      FC%NFX = 1                                              !> no extension in x-direction
      FC%NFC = C%NX                                           !> number of cells between opposite faces
      IF (.NOT.TWO_D) THEN                                    !> only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z')              !> number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2

   CASE (2)

      FC%DH => C%DYL
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X')              !> number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2
      FC%NFY = 1                                              !> no extension in y-direction
      FC%NFC = C%NY                                           !> number of cells between opposite faces
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z')              !> number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2

   CASE (3)

      FC%DH => C%DZL
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X')             !> number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2
      IF (.NOT.TWO_D) THEN                                   !> only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF
      FC%NFZ = 1                                             !> no extension in y-direction
      FC%NFC = C%NZ                                          !> number of cells between opposite faces
END SELECT

FC%NFW = FC%NFX * FC%NFY * FC%NFZ                            !> get number of wall cells for that face
IWG = IWG + FC%NFW                                           !> increase global wall cell counter

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
TYPE (SCARC_WALL_TYPE) , POINTER, DIMENSION(:) :: WC, WF        !> coarse and fine WALL types
TYPE (SCARC_FACE_TYPE) , POINTER               :: FF, FC        !> coarse and fine FACE types
TYPE (SCARC_LEVEL_TYPE), POINTER               :: C, F, OC      !> coarse and fine LEVEL types

SCARC_ROUTINE = 'SCARC_SETUP_FACE'

!> reference coarse and fine LEVEL type
C => SCARC(NM)%LEVEL(NL)
F => SCARC(NM)%LEVEL(NL-1)

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
         NX1 = C%NX+1
         NX2 = C%NX+1
      ENDIF
      NY1 = 1
      NY2 = C%NY
      NZ1 = 1
      NZ2 = C%NZ
   CASE (2)
      NX1 = 1                                                      !> set dimensions for wall cell counting
      NX2 = C%NX
      IF (IOR0 > 0) THEN
         NY1 = 0
         NY2 = 0
      ELSE
         NY1 = C%NY+1
         NY2 = C%NY+1
      ENDIF
      NZ1 = 1
      NZ2 = C%NZ
   CASE (3)
      NX1 = 1                                                      !> set dimensions for wall cell counting
      NX2 = C%NX
      NY1 = 1
      NY2 = C%NY
      IF (IOR0 > 0) THEN
         NZ1 = 0
         NZ2 = 0
      ELSE
         NZ1 =C%NZ+1
         NZ2 =C%NZ+1
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
               WC(IWC)%ICW = (IZ-1)*C%NX*C%NY + (IY-1)*C%NX + IX + 1
            CASE (-1)
               WC(IWC)%ICW = (IZ-1)*C%NX*C%NY + (IY-1)*C%NX + IX - 1
            CASE (2)
               WC(IWC)%ICW = (IZ-1)*C%NX*C%NY +  IY   *C%NX + IX
            CASE (-2)
               WC(IWC)%ICW = (IZ-1)*C%NX*C%NY + (IY-2)*C%NX + IX
            CASE (3)
               WC(IWC)%ICW =  IZ   *C%NX*C%NY + (IY-1)*C%NX + IX
            CASE (-3)
               WC(IWC)%ICW = (IZ-2)*C%NX*C%NY + (IY-1)*C%NX + IX
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
                  IWF(1) = FF%IWG_PTR + 2*(IZ-1)*F%NX + 2*(IX - 1)
               CASE ( 3)
                  IWF(1) = FF%IWG_PTR + 2*(IX-1)
            END SELECT
            IWF(2) = IWF(1)+1

            !> set fine cell neighbors (they must be the same for all fine IW's)
            NOMF(1) = WF(IWF(1))%NOM
            NOMF(2) = WF(IWF(2))%NOM
            IF (NOMF(1) /= NOMF(2)) THEN
               WRITE(SCARC_MESSAGE,'(2A,2I4)') TRIM(SCARC_ROUTINE),': Inconsistent neighbors ',NOMF(1),NOMF(2)
               CALL SHUTDOWN(SCARC_MESSAGE); RETURN
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

            !> in case of an internal boundary set neighboring WALL cells
            IF (NOMF(1) > 0) THEN

               OC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)

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

               !>!
               !> Allocate and specify ICN and ICE arrays for OC
               !>!
               NCPL = (IZ2-IZ1+1)*(IX2-IX1+1)
               OC%NCPL = NCPL

               CALL SCARC_SETUP_WALLCELL_NEIGHBOR(IX1, IX2, 1, 1, IZ1, IZ2, IWC, IOR0, NM, NOMF(1), NL)

            ENDIF


         !> ------------------------------------------------------------
         !> 3D-version
         !> ------------------------------------------------------------
         ELSE

            !> determine fine IW's, which must be merged to one coarse IW
            SELECT CASE (ABS(IOR0))
               CASE (1)
                  IWF(1) = FF%IWG_PTR + (2*IZ-2)*F%NY + 2*IY - 2
                  IWF(3) = FF%IWG_PTR + (2*IZ-1)*F%NY + 2*IY - 2
               CASE (2)
                  IWF(1) = FF%IWG_PTR + (2*IZ-2)*F%NX + 2*IX - 2
                  IWF(3) = FF%IWG_PTR + (2*IZ-1)*F%NX + 2*IX - 2
               CASE (3)
                  IWF(1) = FF%IWG_PTR + (2*IY-2)*F%NX + 2*IX - 2
                  IWF(3) = FF%IWG_PTR + (2*IY-1)*F%NX + 2*IX - 2
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

            !> in case of an internal boundary set WALL(10:15,IWC)
            IF (NOMF(1) > 0) THEN

               OC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)

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

               !> Allocate and specify ICN and ICE arrays for OC
               NCPL = (IZ2-IZ1+1)*(IY2-IY1+1)*(IX2-IX1+1)
               OC%NCPL = NCPL

               CALL SCARC_SETUP_WALLCELL_NEIGHBOR(IX1, IX2, IY1, IY2, IZ1, IZ2, IWC, IOR0, NM, NOMF(1), NL)

            ENDIF
         ENDIF
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
TYPE (SCARC_LEVEL_TYPE), POINTER :: L, OL

L  => SCARC(NM)%LEVEL(NL)
OS => SCARC(NM)%OSCARC(NOM)
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

OS%NIC_R = 0
FOUND = .FALSE.

SEARCH_LOOP: DO IW=1, OL%NCG

   !> neighborship structure already known from finest level
   IF (L%WALL(IW)%NOM/=NOM) CYCLE SEARCH_LOOP
   OS%NIC_R = OS%NIC_R + 1
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
SUBROUTINE SCARC_SETUP_GLOBAL_FINE
INTEGER :: NM, NM2, NL, NP
INTEGER :: IREFINE, IOFF
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

SCARC_ROUTINE = 'SCARC_SETUP_GLOBAL_FINE'
!>
!> Allocate and preset counter and displacement vector for global data exchanges
!>
ALLOCATE(COUNTS(0:N_MPI_PROCESSES-1), STAT=IERROR)
CALL CHKMEMERR (SCARC_ROUTINE, 'COUNTS', IERROR)

ALLOCATE(DISPLS(0:N_MPI_PROCESSES-1), STAT=IERROR)
CALL CHKMEMERR (SCARC_ROUTINE, 'COUNTS', IERROR)

COUNTS = 0
DO NP = 0, N_MPI_PROCESSES-1
   DO NM = 1, NMESHES
      IF (PROCESS(NM)==NP) COUNTS(NP) = COUNTS(NP) + 1
   ENDDO
ENDDO

DISPLS(0) = 0
DO NP = 1, N_MPI_PROCESSES-1
   DISPLS(NP) = COUNTS(NP-1) + DISPLS(NP-1)
ENDDO


!>
!> Allocate arrays for global scalar products
!>
ALLOCATE(SP_LOCAL(NMESHES), STAT=IERROR)
CALL CHKMEMERR (SCARC_ROUTINE, 'SP_LOCAL', IERROR)
SP_LOCAL = 0.0_EB


!>
!> Compute global number of cells for all levels
!>
IREFINE = 0
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MIN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      L => SCARC(NM)%LEVEL(NL)

      DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L%NC_LOCAL(NM2) = SCARC(NM2)%LEVEL(NL)%NCS
      ENDDO

      IF (PROCESS(NM) /= MYID) THEN
         WRITE(SCARC_MESSAGE,'(3A)') TRIM('SCARC_SETUP_GLOBAL_FINE'),': Wrong process assignment'
         CALL SHUTDOWN(SCARC_MESSAGE); RETURN
      ENDIF

   ENDDO

   !> Determine global number of cells for all levels
   IF (N_MPI_PROCESSES > 1) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,L%NC_LOCAL,COUNTS,DISPLS,&
                          MPI_INTEGER,MPI_COMM_WORLD,IERROR)
   ENDIF

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      L => SCARC(NM)%LEVEL(NL)
      L%NC_GLOBAL = SUM(L%NC_LOCAL(1:NMESHES))
      N_CELLS_GLOBAL(NL) = L%NC_GLOBAL


      IF (NMESHES > 1) THEN
         IOFF = 0
         DO NM2 = 2, NMESHES
            IOFF = IOFF + L%NC_LOCAL(NM2-1)
            L%NC_OFFSET(NM2) = L%NC_OFFSET(NM2) + IOFF
         ENDDO
      ENDIF

   ENDDO

ENDDO LEVEL_LOOP

END SUBROUTINE SCARC_SETUP_GLOBAL_FINE


!> ----------------------------------------------------------------------------------------------------
!> Allocate several global structures for data exchange
!> --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_SETUP_GLOBAL_COARSE
INTEGER :: NM, NM2, NL
INTEGER :: IREFINE, IOFF
!INTEGER, DIMENSION(:), ALLOCATABLE :: NC_AUX
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

!>
!> Compute global number of cells for all higher levels levels
!>
IREFINE = 0
LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      L => SCARC(NM)%LEVEL(NL)

      DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L%NC_LOCAL(NM2) = SCARC(NM2)%LEVEL(NL)%NCS
      ENDDO

      IF (PROCESS(NM) /= MYID) THEN
         WRITE(SCARC_MESSAGE,'(3A)') TRIM('SCARC_SETUP_GLOBAL_FINE'),': Wrong process assignment'
         CALL SHUTDOWN(SCARC_MESSAGE); RETURN
      ENDIF

   ENDDO

   !> Determine global number of cells for all levels
   IF (N_MPI_PROCESSES > 1) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,L%NC_LOCAL,COUNTS,DISPLS,&
                          MPI_INTEGER,MPI_COMM_WORLD,IERROR)
   ENDIF

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      L => SCARC(NM)%LEVEL(NL)
      L%NC_GLOBAL = SUM(L%NC_LOCAL(1:NMESHES))
      N_CELLS_GLOBAL(NL) = L%NC_GLOBAL


      IF (NMESHES > 1) THEN
         IOFF = 0
         DO NM2 = 2, NMESHES
            IOFF = IOFF + L%NC_LOCAL(NM2-1)
            L%NC_OFFSET(NM2) = L%NC_OFFSET(NM2) + IOFF
         ENDDO
      ENDIF

   ENDDO

ENDDO LEVEL_LOOP

END SUBROUTINE SCARC_SETUP_GLOBAL_COARSE


!> ----------------------------------------------------------------------------------------------------
!> Setup system of equation:
!> Define matrix stencils and initialize matrices and boundary conditions on all needed levels
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM
INTEGER :: NM, NL

SCARC_ROUTINE = 'SCARC_SETUP_SYSTEM'

!> ---------------------------------------------------------------------------------------------
!> Setup sizes for system matrices 
!> ---------------------------------------------------------------------------------------------
SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   CASE (NSCARC_METHOD_KRYLOV)
      CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NLEVEL_MIN)
      IF (BTWOLEVEL) CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NLEVEL_MAX)
      IF (BCGGMG) THEN
         DO NL=NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NL)
         ENDDO
      ENDIF

   CASE (NSCARC_METHOD_MULTIGRID)
      SELECT CASE (TYPE_MULTIGRID)
         CASE (NSCARC_MULTIGRID_GEOMETRIC)
            DO NL=NLEVEL_MIN, NLEVEL_MAX
               CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NL)
            ENDDO
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)
            CALL SCARC_SETUP_SIZES (NSCARC_SIZE_MATRIX , NLEVEL_MIN)
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
            CASE (NSCARC_PRECON_MULTIGRID)

               SELECT_PRECON_MG: SELECT CASE (TYPE_MULTIGRID)

                  !> Geometric multigrid:
                  !>    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
                  CASE (NSCARC_MULTIGRID_GEOMETRIC)

                     DO NL = NLEVEL_MIN, NLEVEL_MAX
                        CALL SCARC_SETUP_MATRIX  (NM, NL)
                        CALL SCARC_SETUP_BOUNDARY(NM, NL)
#if defined(WITH_MKL)
                        IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_LOCAL)  CALL SCARC_SETUP_MATRIX_MKL(NM, NL)
                        IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) CALL SCARC_SETUP_MATRIX_MKL(NM, NL)
#endif
                     ENDDO

#if defined(WITH_MKL)
                     IF (TYPE_MKL == NSCARC_MKL_COARSE) CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MAX)
#endif

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
            CASE (NSCARC_PRECON_PARDISO)

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
            CASE (NSCARC_PRECON_CLUSTER)

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
                  IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_LOCAL   .OR. &
                      TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
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
      CASE (NSCARC_METHOD_MKL)

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

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NLEVEL_MIN, 'MATRIX000')
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
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE),  POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A
!#ifdef WITH_MKL
!TYPE (SCARC_MKL_TYPE), POINTER :: LMKL
!#endif

SCARC_ROUTINE = 'SCARC_SETUP_MATRIX'

M => MESHES(NM)
L => SCARC(NM)%LEVEL(NL)
A => SCARC(NM)%LEVEL(NL)%A

IF (TWO_D) THEN
   NSTENCIL = 5
ELSE
   NSTENCIL = 7
ENDIF
NAV = L%NCS*NSTENCIL + 10
NAC = L%NCS*NSTENCIL + 10
NAR = L%NCS+10
CALL SCARC_ALLOCATE_MATRIX(A, NAV, NAC, NAR, NSTENCIL, 'A')

#ifdef WITH_MKL
IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
     TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
   CALL SCARC_ALLOCATE_INT1 (A%COL_GLOBAL, 1, A%NAV+10, 'A%COL_GLOBAL', .TRUE.)
ENDIF
#endif

IF (NL == NLEVEL_MIN) THEN
   CELL_INDEX => M%CELL_INDEX
   WALL_INDEX => M%WALL_INDEX
   !NCELLS     =  N_CELLS(NM)             ! ACHTUNG: Kann weg ?
ELSE
   CELL_INDEX => L%CELL_INDEX
   WALL_INDEX => L%WALL_INDEX
   !NCELLS     =  L%N_CELLS             ! ACHTUNG: Kann weg ?
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

         IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. L%DISCRET(IX, IY, IZ, CGSC)/=IS_GASPHASE) CYCLE

         IWTYPE = 0
         IWNUM  = 0

         IC  = L%DISCRET(IX, IY, IZ, UNKH)
         ICI = CELL_INDEX(IX, IY, IZ)

         CALL SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)

         !> IOR = 1
         IOR0 = 1
         IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
         IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
         IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
         !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            ICNUM  = L%DISCRET(IX-1, IY, IZ, UNKH)
            ICTYPE = L%DISCRET(IX-1, IY, IZ, CGSC)
            CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
         ENDIF

         !> IOR = -1
         IOR0 = -1
         IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
         IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
         IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
         !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            ICNUM  = L%DISCRET(IX+1, IY, IZ, UNKH)
            ICTYPE = L%DISCRET(IX+1, IY, IZ, CGSC)
            CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
         ENDIF

         !> IOR = 3
         IOR0 = 3
         IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
         IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
         IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
         !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            ICNUM  = L%DISCRET(IX, IY, IZ-1, UNKH)
            ICTYPE = L%DISCRET(IX, IY, IZ-1, CGSC)
            CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
         ENDIF

         !> IOR = -3
         IOR0 = -3
         IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
         IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
         IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
         !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            ICNUM  = L%DISCRET(IX, IY, IZ+1, UNKH)
            ICTYPE = L%DISCRET(IX, IY, IZ+1, CGSC)
            CALL SCARC_SETUP_MATRIX_SUBDIAG  (IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
         ENDIF

      ENDDO
   ENDDO

ELSE   ! 3D

   IP  = 1
   DO IZ = 1, L%NZ
      DO IY = 1, L%NY
         DO IX = 1, L%NX
 
            IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. L%DISCRET(IX, IY, IZ, CGSC)/=IS_GASPHASE) CYCLE
 
            IWTYPE = 0
            IWNUM  = 0
 
            IC  = L%DISCRET(IX, IY, IZ, UNKH)
            ICI = CELL_INDEX(IX, IY, IZ)
 
            CALL SCARC_SETUP_MATRIX_MAINDIAG (IC, IP, IX, IY, IZ, NM, NL)
 
            !> set subdiagonal entries depending on type of neighboring cell
            !> IOR = 1
            IOR0 = 1
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)          ! get wall index of related cell
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE    ! if there is a boundary, get corresponding type
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN    ! do, if neighbor is internal or on adjacent mesh
               ICNUM  = L%DISCRET(IX-1, IY, IZ, UNKH)
               ICTYPE = L%DISCRET(IX-1, IY, IZ, CGSC)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = -1
            IOR0 = -1
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%DISCRET(IX+1, IY, IZ, UNKH)
               ICTYPE = L%DISCRET(IX+1, IY, IZ, CGSC)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IX, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = 2
            IOR0 = 2
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%DISCRET(IX, IY-1, IZ, UNKH)
               ICTYPE = L%DISCRET(IX, IY-1, IZ, CGSC)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IY, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = -2
            IOR0 = -2
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%DISCRET(IX, IY+1, IZ, UNKH)
               ICTYPE = L%DISCRET(IX, IY+1, IZ, CGSC)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IY, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = 3
            IOR0 = 3
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%DISCRET(IX, IY, IZ-1, UNKH)
               ICTYPE = L%DISCRET(IX, IY, IZ-1, CGSC)
               CALL SCARC_SETUP_MATRIX_SUBDIAG(IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
            !> IOR = -3
            IOR0 = -3
            IF (ICI   /= 0) IWNUM  = WALL_INDEX(ICI, -IOR0)
            IF (IWNUM /= 0) IWTYPE = L%WALL(IWNUM)%BOUNDARY_TYPE
            IF (PRES_ON_WHOLE_DOMAIN .OR. IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
            !IF (IWNUM == 0 .OR. IWTYPE==INTERPOLATED_BOUNDARY) THEN
               ICNUM  = L%DISCRET(IX, IY, IZ+1, UNKH)
               ICTYPE = L%DISCRET(IX, IY, IZ+1, CGSC)
               CALL SCARC_SETUP_MATRIX_SUBDIAG  (IOR0, IC, IP, IZ, ICTYPE, ICNUM, NM, NL)
            ENDIF
 
         ENDDO
      ENDDO
   ENDDO

ENDIF
A%ROW(L%NCS+1) = IP
A%NAV          = IP -1                         !> set correct number of matrix entries
A%NAS          = IP -1                                   

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
TYPE (SCARC_LEVEL_TYPE),  POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A

L => SCARC(NM)%LEVEL(NL)
A => SCARC(NM)%LEVEL(NL)%A

A%VAL(IP) = A%VAL(IP) - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D) A%VAL(IP) = A%VAL(IP) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
A%VAL(IP) = A%VAL(IP) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

A%ROW(IC) = IP
A%COL(IP) = IC
      
#ifdef WITH_MKL
IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
     TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
   A%COL_GLOBAL(IP) = A%COL(IP) + L%NC_OFFSET(NM)
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
TYPE (MESH_TYPE)        , POINTER :: M
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A

M => MESHES(NM)
L => SCARC(NM)%LEVEL(NL)
A => SCARC(NM)%LEVEL(NL)%A

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

   IF (PRES_ON_WHOLE_DOMAIN .OR. NBR_TYPE == IS_GASPHASE) THEN
      A%VAL(IP) = A%VAL(IP) + DSCAL
      A%COL(IP) = NBR_NUM
#ifdef WITH_MKL
      IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
          (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
           TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
         A%COL_GLOBAL(IP)= A%COL(IP) + L%NC_OFFSET(NM)
      ENDIF
#endif
      IP = IP + 1
   ELSE
      WRITE(*,*) 'MISSING, NOTHING TO DO?'
   ENDIF
         
!> if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
ELSE IF (L%FACE(IOR0)%NUM_NEIGHBORS /= 0) THEN
   
   IF (HAS_NEIGHBOR(L%DISCRET, L%WALL, IC, IW, L%FACE(IOR0)%IWG_PTR, L%FACE(IOR0)%NFW)) THEN
      A%VAL(IP) = A%VAL(IP) + DSCAL
      A%COL(IP) = L%WALL(IW)%ICE(1)
      !A%COL(IP) = L%WALL(IW)%ICO
         
#ifdef WITH_MKL
      IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
          (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
           TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
         IX = L%WALL(IW)%IXG
         IY = L%WALL(IW)%IYG
         IZ = L%WALL(IW)%IZG
         A%COL_GLOBAL(IP) = L%DISCRET(IX, IY, IZ, UNKH) + L%NC_OFFSET(L%WALL(IW)%NOM)
      ENDIF
#endif
         
      IP = IP + 1
   ENDIF
   
ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG


!> ------------------------------------------------------------------------------------------------
!> Determine if cell IC has a neighbor and, if yes, return corresponding IW-value
!> ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION HAS_NEIGHBOR(DISCRET0, WALL0, IC, IW, IW0, IL0)
TYPE (SCARC_WALL_TYPE), DIMENSION(:), INTENT(IN):: WALL0
INTEGER, DIMENSION(0:,0:,0:,:), INTENT(IN):: DISCRET0
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

  IC0 = DISCRET0(IXW0, IYW0, IZW0, UNKH) 
  IF (IC == IC0 .AND. (WALL0(IW)%NOM /= 0.AND.DISCRET0(IXG0, IYG0, IZG0, CGSC)/=IS_SOLID)) THEN
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
REAL(EB) :: VAL, VALS, DIFF
LOGICAL  :: BSYM
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, AS
INTEGER, DIMENSION(:), ALLOCATABLE :: ICOL_AUX, JC_AUX

SCARC_ROUTINE = 'SCARC_SETUP_MATRIX_MKL'

L  => SCARC(NM)%LEVEL(NL)
A  => SCARC(NM)%LEVEL(NL)%A
AS => SCARC(NM)%LEVEL(NL)%AS

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
                  IF (ABS(VAL - VALS) < 1E-10) THEN
                     BSYM=.TRUE.
                     CYCLE COLUMN_LOOP
                  ENDIF
               ENDIF
            ENDDO
            IF (.NOT.BSYM) THEN
               WRITE(SCARC_MESSAGE,'(3A)') SCARC_ROUTINE,': System matrix not symmetric'
               CALL SHUTDOWN(SCARC_MESSAGE); RETURN
            ENDIF
         ENDIF
      ENDDO COLUMN_LOOP
   
   ENDDO

   !> Compute number of entries in symmetric matrix
   A%NAS = 0
   DO IC = 1, L%NCS
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_LOCAL) THEN
            JC = A%COL(ICOL)
            IF (JC >= IC .AND. JC <= L%NCS) A%NAS = A%NAS+1   !really neccesary ???
         ELSE IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
            JC = A%COL_GLOBAL(ICOL)
            IF (JC >= IC + L%NC_OFFSET(NM)) A%NAS = A%NAS+1
         ELSE
            WRITE(SCARC_MESSAGE,'(3A)') SCARC_ROUTINE,'ERROR IN SETUP_MATRIX_MKL'
            WRITE(*,*) 
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
CALL SCARC_ALLOCATE_MATRIX(AS, A%NAV, A%NAC, A%NAR, A%NSTENCIL, 'AS')

IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
     TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
   ALLOCATE (ICOL_AUX(A%NSTENCIL), STAT=IERROR)
   CALL CHKMEMERR (SCARC_ROUTINE, 'ICOL_AUX', IERROR)
   ALLOCATE (JC_AUX(A%NSTENCIL), STAT=IERROR)
   CALL CHKMEMERR (SCARC_ROUTINE, 'JC_AUX', IERROR)
ENDIF


!> ------------------------------------------------------------------------------------------------
!> extract symmetric matrix part from usual system matrix
!> ------------------------------------------------------------------------------------------------
IAS = 1
DO IC = 1, L%NCS
   AS%ROW(IC) = IAS

   !> blockwise use of MKL solver
   IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_LOCAL) THEN    !really neccesary ??

      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         JC = A%COL(ICOL)

            IF (JC >= IC .AND. JC <= L%NCS) THEN
               AS%COL(IAS) = A%COL(ICOL) 
               AS%VAL(IAS)     = A%VAL(ICOL)
               IAS = IAS + 1
            ENDIF
      ENDDO

   !> global use of MKL solver
   ELSE IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN   

      !> store indices of all diagonal and upper-diagonal entries
      ICOL_AUX = 0
      JC_AUX   = 99999999
      ISYM = 1
      JC0 = A%COL_GLOBAL(A%ROW(IC))
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         !JC = A%COL(ICOL)
         JC = A%COL_GLOBAL(ICOL)
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
               AS%VAL(IAS) = A%VAL(ICOL)
               !AS%COL(IAS) = A%COL(ICOL)
               AS%COL(IAS) = A%COL_GLOBAL(ICOL)
               JC_AUX(ISYM) = 99999999            ! mark entry as already used
               IAS  = IAS  + 1
            ENDIF
         ENDDO
         JSYM = JSYM + 1
      ENDDO SORT_LOOP
   ENDIF
ENDDO

AS%ROW(L%NCS+1) = IAS
AS%COL = AS%COL 

IF ((TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
    (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
     TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
   DEALLOCATE (ICOL_AUX, STAT=IERROR)
   DEALLOCATE (JC_AUX, STAT=IERROR)
   DEALLOCATE (A%COL_GLOBAL, STAT=IERROR)
ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_MKL
#endif


!> ----------------------------------------------------------------------------------------------------
!> Extract overlapping matrix data
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OVERLAPS (NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM, LAST_NOM
INTEGER :: IW, IWG, IWC, IG, IC, JC, JCA, ICW, ICN, ICE, ICE0, ICE2, ICG, ICG2
INTEGER :: IX, IY, IZ
INTEGER :: ICOL, ICCE, ICPL, NCOL
INTEGER :: ICC, ICOLG, ICOLE, IOR0, NCPL, NCE0
TYPE (SCARC_TYPE)       , POINTER :: S
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L, OL, F, C, OF, OC
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, OA
TYPE (SCARC_AMG_TYPE), POINTER :: AMG

SCARC_ROUTINE = 'SCARC_SETUP_OVERLAPS'

SELECT_TYPE: SELECT CASE (NTYPE)

   !> --------------- set stencil sizes on overlapping parts  -----------------------------
   CASE (NSCARC_MATRIX_STENCIL)

      STENCIL_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         S  => SCARC(NM)
         L  => S%LEVEL(NL)
         A  => S%LEVEL(NL)%A

         STENCIL_EXTERNAL_CELL_LOOP: DO ICE = L%NCS+1, L%NCE

           IW = ABS(L%MAP%ICE_TO_IWG(ICE))

           IF (IW > 0) THEN
              ICG = L%WALL(IW)%ICG(1)
           ELSE
              IW  = ABS(IW)
              ICG = L%WALL(ABS(IW))%ICG2(1)
           ENDIF

           NOM = L%WALL(IW)%NOM
           OA => S%OSCARC(NOM)%LEVEL(NL)%A

           A%ROW(ICE+1) = A%ROW(ICE) + OA%STENCIL(ICG)

         ENDDO STENCIL_EXTERNAL_CELL_LOOP

      ENDDO STENCIL_MESHES_LOOP

   !> -------------------- after exchange of system matrix ------------------------
   CASE (NSCARC_MATRIX_SYSTEM)

      SYSTEM_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         S => SCARC(NM)
         L => S%LEVEL(NL)
         AMG => S%LEVEL(NL)%AMG

         NCE0 = L%NCE+1
         DO ICE = L%NCS+1, L%NCE

            IC = ICE
            IWG = L%MAP%ICE_TO_IWG(ICE)
            IX  = L%WALL(IWG)%IXG
            IY  = L%WALL(IWG)%IYG
            IZ  = L%WALL(IWG)%IZG
            ICN = L%DISCRET(IX, IY, IZ, UNKH)

            NOM = L%WALL(IWG)%NOM
            OL => S%OSCARC(NOM)%LEVEL(NL)

            ICPL = 0
            DO ICOL = A%ROW(ICE),A%ROW(ICE+1)-1
               JC  = A%COL(ICOL)
               JCA = ABS(A%COL(ICOL))
WRITE(*,*) 'ACHTUNG HIER WEGEN ICN SCHAUEN!'
               IF (JC < 0 .AND. JCA <= OL%NC) THEN
                  !IF (OL%MAP%ICN_TO_ICE(JCA)==0) THEN
                  IF (JC==-1) THEN
                     OL%MAP%ICN_TO_ICE(JCA) = NCE0
                     AMG%CELL_MAP(2, ICE) = NCE0
                     NCE0 = NCE0 + 1
                  ELSE
                     AMG%CELL_MAP(2, ICE) = OL%MAP%ICN_TO_ICE(JCA)
                  ENDIF
                  AMG%CELL_MAP(3, ICE) = JC
                  A%COL(ICOL) = AMG%CELL_MAP(2,ICE)
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

         S  => SCARC(NM)
         F => S%LEVEL(NL)
         C => S%LEVEL(NL+1)

         IG   = 0
         IWC  = 0
         ICCE = F%AMG%NCC
         PROLONGATION_WALL_LOOP1: DO IW = 1, F%NW

            NOM  = F%WALL(IW)%NOM
            IF (NOM == 0) CYCLE PROLONGATION_WALL_LOOP1

            OF => S%OSCARC(NOM)%LEVEL(NL)
            OC => S%OSCARC(NOM)%LEVEL(NL+1)

            IOR0 = F%WALL(IW)%IOR
            ICW  = F%WALL(IW)%ICW
            ICC  = F%AMG%CELL_TYPE(ICW)

            ICE   = F%WALL(IW)%ICE(1)
            ICG   = F%WALL(IW)%ICG(1)
            NCOL  = OF%P%ROW(ICG +1)-OF%P%ROW(ICG)         !> wozu wird das nochmal gebraucht ??

            IF (ICC >= NSCARC_CELL_TYPE_COARSE) THEN
               IWC = IWC + 1
               !C%INTERNAL_BDRY_CELL(ICC) = NSCARC_LAYER_ONE
               C%WALL_INDEX(ICC, IOR0)   = IWC
               C%WALL(IWC)%NOM = NOM
               C%WALL(IWC)%ICW = ICC
               C%WALL(IWC)%IOR = IOR0

            ENDIF

            IF (F%AMG%CELL_TYPE(ICE) >= NSCARC_CELL_TYPE_COARSE) THEN
               ICCE = ICCE + 1
               F%AMG%CELL_TYPE(ICE) = ICCE
               !C%EXT_PTR(ICCE) = OF%AMG%CELL_TYPE(ICG)
               ICOL = OF%P%ROW(ICG)
               IG  = IG + 1
               !OF%ICG0 = OF%ICG0 + 1
               OF%ICG0 = IG
               OC%MAP%ICN_TO_ICE(OF%P%COL(ICOL))   = ICCE
               OC%MAP%ICG_TO_ICE(OF%ICG0) = ICCE
            ENDIF

            IF (NL/=NLEVEL_MIN .OR. TYPE_LAYER/=NSCARC_LAYER_TWO) CYCLE PROLONGATION_WALL_LOOP1

            ICE2=0
            ICG2=0
            WRITE(*,*) 'ACHTUNG: ICE2 UND ICG2 NICHT RICHTIG INITIALISIERT!'
            IF (F%AMG%CELL_TYPE(ICE2) >= NSCARC_CELL_TYPE_COARSE) THEN
               F%AMG%CELL_TYPE(ICE2) = ICCE
               C%EXT_PTR(ICCE) = OF%AMG%CELL_TYPE(ICG2)
               ICOL = OF%P%ROW(ICG2)
               IG  = IG + 1
               OF%ICG0 = OF%ICG0 + 1
               OC%MAP%ICN_TO_ICE(OF%P%COL(ICOL)) = ICCE
               ICCE = ICCE + 1
            ENDIF

         ENDDO PROLONGATION_WALL_LOOP1

         !> Replace negative cell numbers of neighboring cells for internal wall cells
         !> Note that this only holds true for fine grid cells
         !> (weights of coarse cells don't reach in neighboring mesh)
         IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) THEN

         PROLONGATION_WALL_LOOP2: DO IW = 1, F%NW

            NOM  = F%WALL(IW)%NOM
            IF (NOM == 0) CYCLE PROLONGATION_WALL_LOOP2

            OF => S%OSCARC(NOM)%LEVEL(NL)

            ICW  = F%WALL(IW)%ICW
            ICC  = F%AMG%CELL_TYPE(ICW)

            IF (ICC < NSCARC_CELL_TYPE_COARSE) THEN

               DO ICOL = F%P%ROW(ICW), F%P%ROW(ICW+1)-1
                  JC   = F%P%COL(ICOL)

                  !> Additionally identify coarse cells from second layer
                  !> adjacent to internal boundary
WRITE(*,*) 'ACHTUNG HIER WEGEN ICN SCHAUEN'
                  IF (JC < 0) THEN
                     IC = OF%MAP%ICN_TO_ICE(ABS(JC))
                     F%P%COL(ICOL) = F%AMG%CELL_TYPE(IC)

                  ELSE IF (JC <= F%AMG%NCC) THEN
                     IF (C%INTERNAL_BDRY_CELL(JC) > NSCARC_LAYER_ONE) THEN
                        IWC = IWC + 1
                        C%INTERNAL_BDRY_CELL(JC) = NSCARC_LAYER_TWO
                        C%WALL_INDEX(ICC, IOR0)  = -IWC
                        C%WALL(IWC)%NOM = NOM
                        C%WALL(IWC)%ICW = JC
                        C%WALL(IWC)%IOR = F%WALL(IW)%IOR
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO PROLONGATION_WALL_LOOP2

         PROLONGATION_ECELL_LOOP: DO ICE0 = F%NCS+1, F%NCE

            IW  = L%MAP%ICE_TO_IWG(ICE0)
            NOM = F%WALL(ABS(IW))%NOM
            OF => S%OSCARC(NOM)%LEVEL(NL)
            OC => S%OSCARC(NOM)%LEVEL(NL+1)

            !> Replace positive cell numbers of neighboring cells in first layer

            !WRITE(*,*) 'ACHTUNG: HIER WALL%ICN aendern auf OFI!'
            !WRITE(*,*) 'ACHTUNG: MUSS NOCH UEBERARBEITET WERDEN !>!'

            IF (IW > 0) THEN
               ICG = F%WALL(IW)%ICG(1)
               ICE = F%WALL(IW)%ICE(1)
               ICOLE = F%P%ROW(ICE)

               DO ICOLG = OF%P%ROW(ICG), OF%P%ROW(ICG+1)-1
                  F%P%COL(ICOLE) = OF%P%COL(ICOLG)
                  F%P%COL(ICOLE) = OF%MAP%ICN_TO_ICE(OF%P%COL(ICOLG))
                  F%P%VAL(ICOLE) = OF%P%VAL(ICOLG)
                  ICOLE = ICOLE + 1

!>                 IF (JC > 0 .AND. JC <= OF%AMG%NCC) THEN
!>                    IC = OC%MAP%ICN_TO_ICE(JC)
!>                    !OF%P%COL(ICOLG) = IC
!>                    F%P%COL(ICOLE)  = IC
!>                 ELSE
!>                    !OF%P%COL(ICOLG) = F%AMG%CELL_TYPE(ABS(JC))
!>                    F%P%COL(ICOLE) = F%AMG%CELL_TYPE(ABS(JC))
!>                 ENDIF
!>                 F%P%VAL(ICOLE) = OF%P%VAL(ICOLG)
!>                 ICOLE = ICOLE + 1
               ENDDO
               F%P%ROW(ICE+1) = ICOLE

            !> Replace positive cell numbers of neighboring cells in second layer (if requested)
            ELSE
               IW = ABS(IW)

               WRITE(*,*) 'D: HIER WALL%ICN aendern auf OFI!'
               WRITE(*,*) 'D: MUSS NOCH UEBERARBEITET WERDEN !>'

               IF (NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) THEN
                  ICOLE = F%P%ROW(ICE2)
                  DO ICOLG = OF%P%ROW(ICG2), OF%P%ROW(ICG2+1)-1
                     JC   = OF%P%COL(ICOLG)
                     IF (JC > 0 .AND. JC <= OF%AMG%NCC) THEN
                        IC = OC%MAP%ICN_TO_ICE(JC)
                        IF (IC > 0) THEN
                           OF%P%COL(ICOLG) = IC
                           F%P%COL(ICOLE) = IC
                        ELSE
                           OF%P%COL(ICOLG) = -OF%P%COL(ICOLG)
                           F%P%COL(ICOLE)  =  OF%P%COL(ICOLG)
                        ENDIF
                     ELSE
                        OF%P%COL(ICOLG) = F%AMG%CELL_TYPE(ABS(JC))
                        F%P%COL(ICOLE) = F%AMG%CELL_TYPE(ABS(JC))
                     ENDIF
                     F%P%VAL(ICOLE) = OF%P%VAL(ICOLG)
                     ICOLE = ICOLE + 1
                  ENDDO
                  F%P%ROW(ICE2+1) = ICOLE
               ENDIF
            ENDIF

         ENDDO PROLONGATION_ECELL_LOOP

         C%NW = IWC
         ENDIF     !GMG-ENDIF

      ENDDO PROLONGATION_MESHES_LOOP

   !> -------------------- after exchange of GMG-matrix -------------------------------
   CASE (NSCARC_MATRIX_GMG)

      GMG_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         S => SCARC(NM)
         F => S%LEVEL(NL)
         C => S%LEVEL(NL+1)

         IWC  = 0
         ICCE = F%AMG%NCC
         LAST_NOM = 0
         GMG_WALL_LOOP1: DO IW = 1, F%NW

            NOM  = F%WALL(IW)%NOM
            IF (NOM == 0) CYCLE GMG_WALL_LOOP1
            IF (NOM /= LAST_NOM) IG = 0
            LAST_NOM = NOM

            OF => S%OSCARC(NOM)%LEVEL(NL)
            OC => S%OSCARC(NOM)%LEVEL(NL+1)

            IOR0 = F%WALL(IW)%IOR
            ICW  = F%WALL(IW)%ICW
            ICC  = F%AMG%CELL_TYPE(ICW)

            WRITE(*,*) 'ACHTUNG: HIER WALL%ICN aendern auf OFI!'
            WRITE(*,*) 'ACHTUNG: MUSS NOCH UEBERARBEITET WERDEN !>'
            WRITE(*,*) 'ACHTUNG: ICG NICHT RICHTIG INITIALISIERT!'
            ICG=0

            !ICE   = F%WALL(IW)%ICE(1)
            !ICG   = F%WALL(IW)%ICG(1)
            NCOL  = OF%P%ROW(ICG +1)-OF%P%ROW(ICG)

            IF (ICC >= NSCARC_CELL_TYPE_COARSE) THEN

               IWC  = IWC + 1
               ICCE = ICCE + 1
               NCPL = 1

               C%INTERNAL_BDRY_CELL(ICC) = NSCARC_LAYER_ONE
               C%WALL_INDEX(ICC, IOR0)   = IWC

               WRITE(*,*) 'ACHTUNG: HIER WALL%ICN aendern auf OFI!'
               WRITE(*,*) 'ACHTUNG: MUSS NOCH UEBERARBEITET WERDEN !>'

               C%WALL(IWC)%NOM  = NOM
               !C%WALL(IWC)%ICW  = ICC
               C%WALL(IWC)%IOR  = IOR0
               !C%WALL(IWC)%NCPL = NCPL

               !IF (F%AMG%CELL_TYPE(ICE) < NSCARC_CELL_TYPE_COARSE) THEN
               !   WRITE(SCARC_MESSAGE,'(2A)') TRIM(SCARC_ROUTINE),': Wrong cominbation of coarse cells in OVERLAP_GMG, stop!'
               !   CALL SHUTDOWN(SCARC_MESSAGE); RETURN
               !NDIF

               WRITE(*,*) 'ACHTUNG: HIER WALL%ICN aendern auf OFI!'
               WRITE(*,*) 'ACHTUNG: MUSS NOCH UEBERARBEITET WERDEN !>'
               WRITE(*,*) 'ACHTUNG: ICE NICHT RICHTIG INIALISIERT'
               ICE=0

               !ALLOCATE(C%WALL(IWC)%ICG(NCPL), STAT=IERROR)
               !CALL ChkMemErr('SCARC_SETUP_OVERLAP_GMG','ICG',IERROR)

               !ALLOCATE(C%WALL(IWC)%ICE(NCPL), STAT=IERROR)
               !CALL ChkMemErr('SCARC_SETUP_SETUP_OVERLAP_GMG','ICE',IERROR)

               !ALLOCATE(C%WALL(IWC)%ICN(NCPL), STAT=IERROR)
               !CALL ChkMemErr('SCARC_SETUP_SETUP_OVERLAP_GMG','ICN',IERROR)

               F%AMG%CELL_TYPE(ICE) = ICCE
               C%EXT_PTR(ICCE) = OF%AMG%CELL_TYPE(ICG)
               ICOL = OF%P%ROW(ICG)

               IG  = IG + 1
               OC%MAP%ICN_TO_ICE(OF%P%COL(ICOL)) = ICCE
               OC%MAP%ICG_TO_ICE(IG)             = ICCE

               !C%WALL(IWC)%ICE(1) = ICCE
               !C%WALL(IWC)%ICG(1) = IG
               !C%WALL(IWC)%ICN(1) = OF%AMG%CELL_TYPE(ICG)

            ENDIF

         ENDDO GMG_WALL_LOOP1
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
TYPE (SCARC_LEVEL_TYPE),  POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A
TYPE (MESH_TYPE), POINTER :: M

SCARC_ROUTINE = 'SCARC_SETUP_BOUNDARY'

M => MESHES(NM)
L => SCARC(NM)%LEVEL(NL)
A => SCARC(NM)%LEVEL(NL)%A

IF (PRES_ON_WHOLE_DOMAIN) THEN
   NW = L%N_EXTERNAL_WALL_CELLS
ELSE
   NW = L%N_EXTERNAL_WALL_CELLS + L%N_INTERNAL_WALL_CELLS
ENDIF

IF (TWO_D) THEN

   WALL_CELLS_LOOP2D: DO IW = 1, NW

      IOR0 = L%WALL(IW)%IOR
      IF (ABS(IOR0) == 2) CYCLE          !> 2D: cycle boundaries in y-direction

      I    = L%WALL(IW)%IXW
      J    = 1
      K    = L%WALL(IW)%IZW

      NOM  = L%WALL(IW)%NOM

      ITYPE = L%DISCRET(I,J,K,CGSC)

      IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.ITYPE/=IS_GASPHASE) CYCLE

      L%WALL(IW)%ICW = L%DISCRET(I, J, K, UNKH)
      IC = L%DISCRET(I, J, K, UNKH)

      SELECT CASE (ABS(IOR0))
         CASE (1)
            DBC= L%DXI2
         CASE (3)
            DBC= L%DZI2
      END SELECT

      IP = A%ROW(IC)
      A_OLD = A%VAL(IP)
      SELECT CASE (L%WALL(IW)%BTYPE)
         CASE (DIRICHLET)                      !> set Dirichlet BC's along open boundary cells
            A%VAL(IP) = A%VAL(IP) - DBC
         !CASE (INTERNAL)                      !> do nothing along internal boundaries (only debugging)
         CASE (NEUMANN)                        !> set Neumann BC's at all other nodes
            A%VAL(IP) = A%VAL(IP) + DBC
      END SELECT

   ENDDO WALL_CELLS_LOOP2D

ELSE

   WALL_CELLS_LOOP3D: DO IW = 1, NW

      IOR0 = L%WALL(IW)%IOR

      I    = L%WALL(IW)%IXW
      J    = L%WALL(IW)%IYW
      K    = L%WALL(IW)%IZW

      NOM  = L%WALL(IW)%NOM
      ITYPE = L%DISCRET(I,J,K,CGSC)

      IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.ITYPE/=IS_GASPHASE) CYCLE

      L%WALL(IW)%ICW = L%DISCRET(I, J, K, UNKH)
      IC = L%DISCRET(I, J, K, UNKH)

      SELECT CASE (ABS(IOR0))
         CASE (1)
            DBC= L%DXI2           ! Achtung: Wirklich richtig oder Mittelwert?
         CASE (2)
            DBC= L%DYI2
         CASE (3)
            DBC= L%DZI2
      END SELECT
      IP = A%ROW(IC)
      SELECT CASE (L%WALL(IW)%BTYPE)
         CASE (DIRICHLET)                      !> set Dirichlet BC's at open and null boundary cells
            A%VAL(IP) = A%VAL(IP) - DBC
         !CASE (INTERNAL)                      !> do nothing along internal boundaries (only debugging)
         CASE (NEUMANN)                        !> set Neumann BC's at all other cells
            A%VAL(IP) = A%VAL(IP) + DBC
      END SELECT

   ENDDO WALL_CELLS_LOOP3D

ENDIF

END SUBROUTINE SCARC_SETUP_BOUNDARY

!> ----------------------------------------------------------------------------------------------------
!> Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SOLVER
INTEGER :: NM, NL

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   SELECT_METHOD: SELECT CASE (TRIM(SCARC_METHOD))

      !> ------------------ Krylov method -------------------------------------
      CASE ('KRYLOV')

         !> allocate necessary working vectors for finest level
         CALL SCARC_SETUP_KRYLOV(NSCARC_SCOPE_MAIN, NM, NLEVEL_MIN)

         !> allocate additional workspace for different preconditioners
         SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

            !> Preconditioning by blockwise FFT
            CASE (NSCARC_PRECON_FFT)                 
               CALL SCARC_SETUP_FFT(NM, NLEVEL_MIN)

            !> Preconditioning by blockwise Multigrid
            CASE (NSCARC_PRECON_MULTIGRID)          
               DO NL = NLEVEL_MIN, NLEVEL_MAX
                  CALL SCARC_SETUP_MULTIGRID(NSCARC_SCOPE_PRECON, NM, NL)
               ENDDO
               IF (TYPE_SMOOTH == NSCARC_SMOOTH_FFT) THEN    
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     CALL SCARC_SETUP_FFT(NM, NL)
                  ENDDO
               ELSE IF (TYPE_SMOOTH == NSCARC_SMOOTH_PARDISO) THEN    
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     IF (TYPE_SMOOTH == NSCARC_SMOOTH_PARDISO) CALL SCARC_SETUP_PARDISO(NM, NL)
                  ENDDO
               ENDIF
               CALL SCARC_SETUP_COARSESOLVER(NSCARC_SCOPE_PRECON, NM, NLEVEL_MAX)

#ifdef WITH_MKL
            !> Preconditioning by blockwise Pardiso
            CASE (NSCARC_PRECON_PARDISO)
               CALL SCARC_SETUP_PARDISO(NM, NLEVEL_MIN)

            !> Preconditioning by global LU-decomposition (Cluster sparse solver)
            CASE (NSCARC_PRECON_CLUSTER)         
               CALL SCARC_SETUP_CLUSTER(NM, NLEVEL_MIN)
#endif
   
         END SELECT SELECT_KRYLOV_PRECON

         !> if Twolevel-method is chosen, allocate intermediate structures for interpolation
         !> and workspace for coarse grid solver
         IF (BTWOLEVEL) THEN
            DO NL = NLEVEL_MIN+1, NLEVEL_MAX
               CALL SCARC_SETUP_INTERPOLATION(NM, NL)
            ENDDO
            CALL SCARC_SETUP_COARSESOLVER(NSCARC_SCOPE_MAIN, NM, NLEVEL_MAX)
         ENDIF
         

      !> ------------------ Multigrid method -------------------------------------
      CASE ('MULTIGRID')

         DO NL = NLEVEL_MIN, NLEVEL_MAX
            CALL SCARC_SETUP_MULTIGRID(NSCARC_SCOPE_MAIN, NM, NL)
            IF (TYPE_SMOOTH == NSCARC_SMOOTH_FFT) CALL SCARC_SETUP_FFT(NM, NL)
#ifdef WITH_MKL
            IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_LOCAL)  CALL SCARC_SETUP_PARDISO(NM, NL)
            IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) CALL SCARC_SETUP_CLUSTER(NM, NL)
#endif
         ENDDO
         CALL SCARC_SETUP_COARSESOLVER(NSCARC_SCOPE_MAIN, NM, NLEVEL_MAX)


#ifdef WITH_MKL
      !> ------------------ MKL method -------------------------------------
      CASE ('MKL')

         CALL SCARC_SETUP_MKL(NM, NLEVEL_MIN)

         SELECT_MKL: SELECT CASE (TYPE_MKL)
            CASE (NSCARC_MKL_GLOBAL)
               CALL SCARC_SETUP_CLUSTER(NM, NL)
            CASE (NSCARC_MKL_LOCAL)
               CALL SCARC_SETUP_PARDISO(NM, NL)
         END SELECT SELECT_MKL
#endif

   END SELECT SELECT_METHOD
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_SOLVER


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_KRYLOV(NSCOPE, NM, NL)
INTEGER, INTENT(IN) :: NSCOPE, NM, NL
INTEGER :: NC
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_SOLVER_TYPE) , POINTER :: S

SCARC_ROUTINE = 'SCARC_SETUP_KRYLOV'

L => SCARC(NM)%LEVEL(NL)
NC = L%NCE

SELECT CASE(NSCOPE)
   CASE (NSCARC_SCOPE_MAIN)              !> Used as main solver 
      S => L%MAIN
   CASE (NSCARC_SCOPE_PRECON)            !> Used as preconditioner solver 
      S => L%PRECON
END SELECT

CALL SCARC_ALLOCATE_REAL1(S%X, 1, NC, 'X', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%F, 1, NC, 'F', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%D, 1, NC, 'D', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%G, 1, NC, 'G', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%W, 1, NC, 'W', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%Y, 1, NC, 'Y', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%Z, 1, NC, 'Z', .TRUE.)

END SUBROUTINE SCARC_SETUP_KRYLOV

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Geometric Multigrid method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MULTIGRID(NSCOPE, NM, NL)
INTEGER, INTENT(IN) :: NSCOPE, NM, NL
INTEGER :: NC
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_SOLVER_TYPE), POINTER :: S

SCARC_ROUTINE = 'SCARC_SETUP_MULTIGRID'

L => SCARC(NM)%LEVEL(NL)
NC = L%NCE

SELECT CASE(NSCOPE)
   CASE (NSCARC_SCOPE_MAIN)              !> Used as main solver 
      S => L%MAIN
   CASE (NSCARC_SCOPE_PRECON)            !> Used as preconditioner solver 
      S => L%PRECON
END SELECT

CALL SCARC_ALLOCATE_REAL1(S%X, 1, NC, 'X', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%F, 1, NC, 'F', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%D, 1, NC, 'D', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%Z, 1, NC, 'Z', .TRUE.)

END SUBROUTINE SCARC_SETUP_MULTIGRID

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for MKL-methods
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSESOLVER(NSCOPE, NM, NL)
INTEGER, INTENT(IN) :: NSCOPE, NM, NL

SCARC_ROUTINE = 'SCARC_SETUP_COARSESOLVER'

SELECT_COARSE: SELECT CASE (TYPE_COARSE)

   CASE (NSCARC_COARSE_ITERATIVE)
      CALL SCARC_SETUP_KRYLOV(NSCOPE, NM, NL)

#ifdef WITH_MKL
   CASE (NSCARC_COARSE_DIRECT)
      IF (N_MPI_PROCESSES > 1) THEN
         CALL SCARC_SETUP_CLUSTER(NM, NL)
      ELSE
         CALL SCARC_SETUP_PARDISO(NM, NL)
      ENDIF
#endif

   CASE DEFAULT
      WRITE(SCARC_MESSAGE,'(3A)') SCARC_ROUTINE, 'Error with input parameter ',TYPE_COARSE
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN

END SELECT SELECT_COARSE

END SUBROUTINE SCARC_SETUP_COARSESOLVER


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for MKL-methods
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MKL(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: NC
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_SOLVER_TYPE) , POINTER :: S

SCARC_ROUTINE = 'SCARC_SETUP_MKL'

L => SCARC(NM)%LEVEL(NL)
S => SCARC(NM)%LEVEL(NL)%MAIN
NC = L%NCE

CALL SCARC_ALLOCATE_REAL1(S%X, 1, NC, 'X', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%F, 1, NC, 'F', .TRUE.)

END SUBROUTINE SCARC_SETUP_MKL


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for additive or multiplicative coarse grid 
!> (corresponding to Schwarz domain decomposition method)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERPOLATION(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: NC
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_SOLVER_TYPE) , POINTER :: S

SCARC_ROUTINE = 'SCARC_SETUP_INTERPOLATION'

L => SCARC(NM)%LEVEL(NL)
S => SCARC(NM)%LEVEL(NL)%MAIN
NC = L%NCE

CALL SCARC_ALLOCATE_REAL1(S%X, 1, NC, 'X', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%F, 1, NC, 'F', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%W, 1, NC, 'W', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%G, 1, NC, 'G', .TRUE.)
CALL SCARC_ALLOCATE_REAL1(S%Z, 1, NC, 'Z', .TRUE.)

END SUBROUTINE SCARC_SETUP_INTERPOLATION

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for blockwise FFT methods
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFT(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: NC
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_FFT_TYPE), POINTER :: F

SCARC_ROUTINE = 'SCARC_SETUP_FFT'

L => SCARC(NM)%LEVEL(NL)
F => SCARC(NM)%LEVEL(NL)%FFT
NC = L%NCE

IF (L%NY == 1) THEN
   CALL SCARC_ALLOCATE_REAL3(F%WORK, 1, L%NX+1, 1, 1,      1, L%NZ+1, 'FFT', .TRUE.)
ELSE
   CALL SCARC_ALLOCATE_REAL3(F%WORK, 1, L%NX+1, 1, L%NY+1, 1, L%NZ+1, 'FFT', .TRUE.)
ENDIF

IF (L%NZ>1) THEN
   IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(F%BXS, 1, L%NY+1, 1, L%NZ+1, 'BXS', .TRUE.)
   IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BXS, 1,      1, 1, L%NZ+1, 'BXS', .TRUE.)
ELSE
   CALL SCARC_ALLOCATE_REAL2(F%BXS, 1, L%NY+1, 1, 1, 'BXS', .TRUE.)
ENDIF

IF (L%NZ>1) THEN
   IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(F%BXF, 1, L%NY+1, 1, L%NZ+1, 'BXF', .TRUE.)
   IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BXF, 1,      1, 1, L%NZ+1, 'BXF', .TRUE.)
ELSE
   CALL SCARC_ALLOCATE_REAL2(F%BXF, 1, L%NY+1, 1, 1,'BXF', .TRUE.)
ENDIF

IF (L%NZ > 1) THEN
   CALL SCARC_ALLOCATE_REAL2(F%BYS, 1, L%NX+1,1, L%NZ+1, 'BYS', .TRUE.)
ELSE
   CALL SCARC_ALLOCATE_REAL2(F%BYS, 1, L%NX+1,1,      1, 'BYS', .TRUE.)
ENDIF

IF (L%NZ > 1) THEN
   CALL SCARC_ALLOCATE_REAL2(F%BYF, 1, L%NX+1,1, L%NZ+1, 'BYF', .TRUE.)
ELSE
   CALL SCARC_ALLOCATE_REAL2(F%BYF, 1, L%NX+1,1,      1, 'BYF', .TRUE.)
ENDIF

IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(F%BZS, 1, L%NX+1, 1, L%NY+1, 'BZS', .TRUE.)
IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BZS, 1, L%NX+1, 1,      1, 'BZS', .TRUE.)

IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(F%BZF, 1, L%NX+1, 1, L%NY+1, 'BZF', .TRUE.)
IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(F%BZF, 1, L%NX+1, 1,      1, 'BZF', .TRUE.)

END SUBROUTINE SCARC_SETUP_FFT


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CLUSTER(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I !, IC, IP
REAL (EB) :: TNOW, DDUMMY(1)=0.0_EB
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_MKL_TYPE)    , POINTER :: M
TYPE (SCARC_MATRIX_TYPE) , POINTER :: AS

SCARC_ROUTINE = 'SCARC_SETUP_CLUSTER'

TNOW = CURRENT_TIME()

L  => SCARC(NM)%LEVEL(NL)
AS => SCARC(NM)%LEVEL(NL)%AS
M  => L%MKL

!> Allocate workspace for parameters needed in MKL-routine
IF (.NOT.ALLOCATED(M%IPARM)) THEN
   ALLOCATE(M%IPARM(64), STAT=IERROR) 
   CALL CHKMEMERR ('SCARC', 'IPARM', IERROR)
   M%IPARM = 0
ENDIF

!> Allocate workspace for pointers needed in MKL-routine
IF (.NOT.ALLOCATED(M%CT)) THEN
   ALLOCATE(M%CT(64), STAT=IERROR) 
   CALL CHKMEMERR ('SCARC', 'CT', IERROR)
   DO I=1,64
      M%CT(I)%DUMMY = 0
   ENDDO
ENDIF

!> Define corresponding parameters
M%NRHS   =  1         ! one right hand side
M%MAXFCT =  1         ! one matrix
M%MNUM   =  1         ! number of matrix to be factorized
M%ERROR  =  0         ! initialize error flag
M%MSGLVL =  0         ! print statistical information

IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
   M%MTYPE  = -2         ! Matrix type real and symmetric indefinite
ELSE
   M%MTYPE  = 11         ! Matrix type real and non-symmetric
ENDIF

! Define control parameter vector iparm:
!M%IPARM(1) = 1   ! no solver default
!M%IPARM(2) = 3   ! Parallel fill-in reordering from METIS
!M%IPARM(4) = 0   ! no iterative-direct algorithm
!M%IPARM(5) = 0   ! no user fill-in reducing permutation
!M%IPARM(6) = 2   ! =0 solution on the first n components of x
!M%IPARM(8) = 2   ! numbers of iterative refinement steps
!M%IPARM(10) = 13 ! perturb the pivot elements with 1E-13
!M%IPARM(11) = 1  ! use nonsymmetric permutation and scaling MPS  !>!>! was 1
!M%IPARM(13) = 1  ! maximum weighted matching algorithm is switched-off
!               !(default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
!M%IPARM(14) = 0  ! Output: number of perturbed pivots
!M%IPARM(18) = 0 !-1 ! Output: number of nonzeros in the factor LU
!M%IPARM(19) = 0 !-1 ! Output: Mflops for LU factorization
!M%IPARM(20) = 0  ! Output: Numbers of CG Iterations
!M%IPARM(21) = 1  ! 1x1 diagonal pivoting for symmetric indefinite matrices.
!M%IPARM(24) = 0
!M%IPARM(27) = 1 ! Check matrix
!M%IPARM(37) = 0 ! CSR-format
!M%IPARM(40) = 2 ! Matrix, solution and rhs provided in distributed assembled matrix input format.

M%IPARM(1)  =  1      ! supply own parameters
M%IPARM(2)  =  3      ! supply own parameters
M%IPARM(4)  =  0      ! supply own parameters
M%IPARM(5)  =  0      ! supply own parameters
M%IPARM(6)  =  0      ! write solution to x
M%IPARM(8)  =  2      ! automatic setting of iterative refinement steps
M%IPARM(10) = 13      ! pivoting perturbation
M%IPARM(11) =  1      ! pivoting perturbation
M%IPARM(13) =  1      ! pivoting perturbation
M%IPARM(14) =  0      ! pivoting perturbation
M%IPARM(18) =  0      ! pivoting perturbation
M%IPARM(19) =  0      ! pivoting perturbation
M%IPARM(20) =  0      ! pivoting perturbation
M%IPARM(21) =  1      ! Bunch-Kaufman pivoting which is default in case of IPARM(0)=0
M%IPARM(24) =  0      ! Bunch-Kaufman pivoting which is default in case of IPARM(0)=0
M%IPARM(27) =  1      ! use matrix checker
M%IPARM(40) = 2       ! provide matrix in distributed format
M%IPARM(41) = L%NC_OFFSET(NM) + 1                  ! first global cell number for mesh NM
M%IPARM(42) = L%NC_OFFSET(NM) + L%NC_LOCAL(NM)    ! last global cell number for mesh NM
!M%IPARM(39) = 2                                    ! provide matrix in distributed format
!M%IPARM(40) = L%NC_OFFSET(NM)+1                   ! first global cell number for mesh NM
!M%IPARM(41) = L%NC_OFFSET(NM)+L%NC_LOCAL(NM)     ! last global cell number for mesh NM


! perform only reordering and symbolic factorization
M%PHASE = 11
CALL CLUSTER_SPARSE_SOLVER(M%CT, M%MAXFCT, M%MNUM, M%MTYPE, M%PHASE, L%NC_GLOBAL, &
                           AS%VAL, AS%ROW, AS%COL, M%PERM, M%NRHS, M%IPARM, &
                           M%MSGLVL, DDUMMY, DDUMMY, MPI_COMM_WORLD, M%ERROR)
               
! perform only factorization
M%PHASE = 22 
CALL CLUSTER_SPARSE_SOLVER(M%CT, M%MAXFCT, M%MNUM, M%MTYPE, M%PHASE, L%NC_GLOBAL, &
                           AS%VAL, AS%ROW, AS%COL, M%PERM, M%NRHS, M%IPARM, &
                           M%MSGLVL, DDUMMY, DDUMMY, MPI_COMM_WORLD, M%ERROR)
         
TSETUP(MYID+1)%CLUSTER=TSETUP(MYID+1)%CLUSTER+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SETUP_CLUSTER

!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PARDISO(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, IDUMMY(1)=0 !, IP, IC
TYPE (SCARC_LEVEL_TYPE)  , POINTER :: L
TYPE (SCARC_MKL_TYPE)    , POINTER :: M
TYPE (SCARC_MATRIX_TYPE) , POINTER :: AS
REAL (EB) :: TNOW, DDUMMY(1)=0.0_EB

SCARC_ROUTINE = 'SCARC_SETUP_PARDISO'

TNOW = CURRENT_TIME()

L  => SCARC(NM)%LEVEL(NL)
AS => SCARC(NM)%LEVEL(NL)%AS
M  => L%MKL

!> Allocate workspace for parameters needed in MKL-routine
IF (.NOT.ALLOCATED(M%IPARM)) THEN
   ALLOCATE(M%IPARM(64), STAT=IERROR) 
   CALL CHKMEMERR ('SCARC', 'IPARM', IERROR)
   M%IPARM = 0
ENDIF

!> Allocate workspace for pointers needed in MKL-routine
IF (.NOT.ALLOCATED(M%PT)) THEN
   ALLOCATE(M%PT(64), STAT=IERROR) 
   CALL CHKMEMERR ('SCARC', 'PT', IERROR)
   DO I=1,64
      M%PT(I)%DUMMY = 0
   ENDDO
ENDIF

!> Define corresponding parameters
M%NRHS   = 1
M%MAXFCT = 1
M%MNUM   = 1

!M%IPARM(1) = 1   ! no solver default
!M%IPARM(2) = 2   ! fill-in reordering from METIS
!M%IPARM(4) = 0   ! no iterative-direct algorithm
!M%IPARM(5) = 0   ! no user fill-in reducing permutation
!M%IPARM(6) = 2   ! =0 solution on the first n components of x
!M%IPARM(8) = 2   ! numbers of iterative refinement steps
!M%IPARM(10) = 13 ! perturb the pivot elements with 1E-13
!M%IPARM(11) = 1  ! use nonsymmetric permutation and scaling MPS  !>!>! was 1
!M%IPARM(13) = 1  ! maximum weighted matching algorithm is switched-off
!                  !(default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
!M%IPARM(14) = 0  ! Output: number of perturbed pivots
!M%IPARM(18) = 0  !-1 ! Output: number of nonzeros in the factor LU
!M%IPARM(19) = 0  !-1 ! Output: Mflops for LU factorization
!M%IPARM(20) = 0  ! Output: Numbers of CG Iterations
!M%IPARM(21) = 1  ! 1x1 diagonal pivoting for symmetric indefinite matrices.
!M%IPARM(24) = 0
!M%IPARM(27) = 1  ! Check matrix
!M%IPARM(37) = 0  ! Matrix, solution and rhs provided in distributed assembled matrix input format.  ???
!M%IPARM(40) = 2  ! Matrix, solution and rhs provided in distributed assembled matrix input format.   ???

M%IPARM(1)  =  1      ! no solver default
M%IPARM(2)  =  2      ! nested dissection algorithm
M%IPARM(4)  =  0      ! factorization computed as required by phase
M%IPARM(5)  =  0      ! user permutation ignored
M%IPARM(6)  =  0      ! write solution on x
M%IPARM(8)  =  2      ! numbers of iterative refinement steps
M%IPARM(10) = 13      ! perturb the pivot elements with 1E-13
M%IPARM(11) =  0      ! disable scaling (default for SPD)
M%IPARM(13) =  0      ! disable matching
M%IPARM(18) = -1      ! Output: number of nonzeros in the factor LU
M%IPARM(19) = -1      ! Output: number of floating points operations
M%IPARM(20) =  1      ! Output: Numbers of CG Iterations
M%IPARM(27) =  1      ! use matrix checker
M%IPARM(37) =  0      ! matrix storage in CSR-format

M%ERROR  =  0       ! initialize error flag
M%MSGLVL =  0       ! print statistical information
M%MTYPE  = -2       ! Matrix type real non-symmetric

! perform only reordering and symbolic factorization
M%PHASE = 11
CALL PARDISO_D(M%PT, M%MAXFCT, M%MNUM, M%MTYPE, M%PHASE, L%NCS, AS%VAL, AS%ROW, AS%COL, &
               IDUMMY, M%NRHS, M%IPARM, M%MSGLVL, DDUMMY, DDUMMY, M%ERROR)
               
! perform only Factorization
M%PHASE = 22 
CALL PARDISO_D(M%PT, M%MAXFCT, M%MNUM, M%MTYPE, M%PHASE, L%NCS, AS%VAL, AS%ROW, AS%COL, &
               IDUMMY, M%NRHS, M%IPARM, M%MSGLVL, DDUMMY, DDUMMY, M%ERROR)

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
TYPE (SCARC_LEVEL_TYPE),  POINTER :: L, OL, F, C, OF, OC
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, OA
TYPE (SCARC_AMG_TYPE),    POINTER :: AMG

SCARC_ROUTINE = 'SCARC_SETUP_SIZES'

SELECT CASE (NTYPE)

   !> --------------------------------------------------------------------------------------------------
   !> Define sizes for system matrix A (including extended regions related to overlaps)
   !> --------------------------------------------------------------------------------------------------
   CASE (NSCARC_SIZE_MATRIX)

      LEVEL_SYSTEM_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L   => SCARC(NM)%LEVEL(NL)
         A   => SCARC(NM)%LEVEL(NL)%A
         AMG => SCARC(NM)%LEVEL(NL)%AMG

         IF (TWO_D) THEN
            A%NSTENCIL = 5
         ELSE
            A%NSTENCIL = 7
         ENDIF
         A%NAV = L%NCS * A%NSTENCIL
         A%NAC = A%NAV
         A%NAR = L%NCS
         A%NAS = L%NCS * A%NSTENCIL

         !> Determine sizes of overlapped parts for later communication with corresponding neighbors
         DO IW = 1, L%NW
            NOM = L%WALL(IW)%NOM
            IF (NOM /= 0) THEN
               OA => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%A
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
            A => SCARC(NM)%LEVEL(NL)%A
            A%NAE = A%NAV
   
            ALLOCATE (AMG%CELL_MAP(3,L%NCS+1:2*L%NCE), STAT=IERROR)
            CALL CHKMEMERR ('SCARC', 'AMG%CELL_MAP', IERROR)
            AMG%CELL_MAP = 0
   
            LEVEL_OTHER_MESHES_LOOP2: DO NOM = 1, NMESHES
   
               IF (NOM == NM) CYCLE LEVEL_OTHER_MESHES_LOOP2
               OS => SCARC(NM)%OSCARC(NOM)
               IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE LEVEL_OTHER_MESHES_LOOP2
   
               OL => OS%LEVEL(NL)
               A%NAE  = A%NAE  + OL%A%NAV
   
               ALLOCATE (OL%A%STENCIL(1:OL%NCG), STAT=IERROR)
               CALL CHKMEMERR ('SCARC', 'STENCIL', IERROR)
               OL%A%STENCIL = 0
   
            ENDDO LEVEL_OTHER_MESHES_LOOP2
   
         ENDDO LEVEL_SYSTEM_MESHES_LOOP2
      ENDIF

   !> -------------------------------------------------------------------------------------------
   !> Define sizes for transfer matrices P and R (including extended regions related to overlaps)
   !> -------------------------------------------------------------------------------------------
   CASE (NSCARC_SIZE_TRANSFER)

      TRANSFER_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         F => SCARC(NM)%LEVEL(NL)
         C => SCARC(NM)%LEVEL(NL+1)

         !> Determine dimensions of restriction and prolongation matrices in own mesh
         F%AMG%NCF  = 0
         F%AMG%NCC  = 0
         F%AMG%NP   = 0
         F%AMG%NR   = 0
         F%AMG%NCCI = 0

         DO IC = 1, F%NC
            IF (F%AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
               F%AMG%NCC = F%AMG%NCC + 1
               F%AMG%NP  = F%AMG%NP  + 1
               F%AMG%NR  = F%AMG%NP
    !>         F%AMG%CELL_TYPE(IC)  = F%AMG%NCC
            ELSE
               F%AMG%NCF = F%AMG%NCF + 1
               F%AMG%NP  = F%AMG%NP  + A%ROW(IC+1)-A%ROW(IC) - 1
               F%AMG%NR  = F%AMG%NP
            ENDIF
         ENDDO

    !>   F%AMG%AMG%NCCE = F%NCC
    !>   F%AMG%NCFE = F%NCF
    !>   DO IC = F%NC+1, F%NCE
    !>      IF (F%AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
    !>         F%AMG%NCCE = F%AMG%NCCE + 1
    !>         F%AMG%CELL_TYPE(IC)  = F%AMG%NCCE
    !>      ELSE
    !>         F%AMG%NCFE = F%AMG%NCFE + 1
    !>      ENDIF
    !>   ENDDO

         F%AMG%NCCI = F%AMG%NCC
         F%AMG%NCCE = F%AMG%NCC
         !F%NCE0 = F%NC+1      !> really needed? OUTDATED ?

         !> Determine number of coarse and fine cells and check correctness of computation
         IF (F%AMG%NCC + F%AMG%NCF /= F%NC) THEN
            WRITE(SCARC_MESSAGE,'(2A,I8,A,I8,A,I4)') TRIM(SCARC_ROUTINE),&
                 ': N_CELLS_COARSE + N_CELLS_FINE = ', F%AMG%NCC + F%AMG%NCF, &
                 '  differs from number of cells = ', F%NC, ' on level ', NL
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
         ENDIF


         !> define variables for overlapping parts
         TRANSFER_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

            IF (NOM == NM) CYCLE TRANSFER_OTHER_MESHES_LOOP
            OS => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE TRANSFER_OTHER_MESHES_LOOP

            OF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            OF%AMG%NP   = 0
            OF%AMG%NR   = 0
            OF%AMG%NPS  = 0
            OF%AMG%NRS  = 0
            OF%AMG%NCC  = 0
            OF%AMG%NCF  = 0
            OF%AMG%NCCS = 0
            OF%AMG%NCFS = 0
            OF%AMG%ICG0 = 0

         ENDDO TRANSFER_OTHER_MESHES_LOOP

         !> Determine sizes of overlapped parts for later communication with corresponding neighbors
         ICCE = F%AMG%NCC
         IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) F%AMG%NCW = 0

         DO IW = 1, F%NW
            NOM = F%WALL(IW)%NOM
            IF (NOM /= 0) THEN

               OF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

               !WRITE(*,*) 'ACHTUNG: HIER WALL%ICN aendern auf OFI!'
               !WRITE(*,*) 'ACHTUNG: MUSS NOCH UEBERARBEITET WERDEN !>'

               IC  = F%WALL(IW)%ICW
               !ICE = F%WALL(IW)%ICE(1)

               IF (TYPE_COARSENING<NSCARC_COARSENING_GMG.AND.F%AMG%CELL_TYPE(IC)>0) F%NCW = F%NCW + 1

               !IF (F%AMG%CELL_TYPE(ICE) >= NSCARC_CELL_TYPE_COARSE) THEN
               IF (F%AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
                  OF%AMG%NPS  = OF%AMG%NPS  + 1
                  OF%AMG%ICG0 = OF%AMG%ICG0  + 1
                  OF%AMG%NRS  = OF%AMG%NPS
                  OF%AMG%NCCS = OF%AMG%NCCS + 1
               ELSE
                  OF%AMG%NPS  = OF%AMG%NPS  + A%ROW(IC+1)-A%ROW(IC) - 1
                  OF%AMG%NRS  = OF%AMG%NPS
                  OF%AMG%NCFS = OF%AMG%NCFS + 1
               ENDIF

               OF%AMG%NCC=OF%AMG%NCCS
               OF%AMG%NCF=OF%AMG%NCFS

            ENDIF
         ENDDO

         !>!
         !> Determine new numbering for coarse cells in interior of mesh
         !>!
         IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) THEN
            ICP   = 0
            DO IC = 1, F%NCE
               IF (F%AMG%CELL_TYPE(IC) == NSCARC_CELL_TYPE_COARSE) THEN
                  ICP = ICP + 1
                  F%AMG%CELL_TYPE(IC) = ICP
               ENDIF
            ENDDO
         ENDIF

      ENDDO TRANSFER_MESHES_LOOP

      !> Exchange sizes for AMG matrices and vectors
      IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_TRANSFER_SIZE, NL)

      !> Determine extended sizes for extended prolongation and restriction matrix
      TRANSFER_MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         F => SCARC(NM)%LEVEL(NL)

         F%AMG%NPE  = F%AMG%NP
         F%AMG%NRE  = F%AMG%NR
         F%AMG%NCCE = F%AMG%NCC
         F%AMG%NCFE = F%AMG%NCF

         TRANSFER_OTHER_MESHES_LOOP2: DO NOM = 1, NMESHES

            IF (NOM == NM) CYCLE TRANSFER_OTHER_MESHES_LOOP2

            OS => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE TRANSFER_OTHER_MESHES_LOOP2

            OF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            OC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

            OC%AMG%NC  = OF%AMG%NCCS
            OF%AMG%NCC = OF%AMG%NCCS
            OC%AMG%NCG = OF%AMG%ICG0

            ALLOCATE (OC%MAP%ICG_TO_ICE(1:OC%NC), STAT=IERROR)
            CALL CHKMEMERR ('SCARC', 'OC%MAP%ICG_TO_ICE', IERROR)
            OC%MAP%ICG_TO_ICE = 0

WRITE(*,*) 'ACHTUNG HIER WEGEN ICN SCHAUEN'
            ALLOCATE (OC%MAP%ICN_TO_ICE(1:OF%AMG%NCCI), STAT=IERROR)
            CALL CHKMEMERR ('SCARC', 'OC%MAP%ICN_TO_ICE', IERROR)
            OC%MAP%ICN_TO_ICE = 0

            DO IG = 1, OC%NC
               OC%MAP%ICG_TO_ICE(IG) = F%AMG%NCCE + IG
            ENDDO

            F%AMG%NPE  = F%AMG%NPE  + OF%AMG%NP
            F%AMG%NRE  = F%AMG%NRE  + OF%AMG%NR
            F%AMG%NCCE = F%AMG%NCCE + OF%AMG%NCC
            F%AMG%NCFE = F%AMG%NCFE + OF%AMG%NCF
         ENDDO TRANSFER_OTHER_MESHES_LOOP2
      ENDDO TRANSFER_MESHES_LOOP2

      !> Determine extended sizes for extended prolongation and restriction matrix
      TRANSFER_MESHES_LOOP3: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         F => SCARC(NM)%LEVEL(NL)
         C => SCARC(NM)%LEVEL(NL+1)

         C%NC  = F%AMG%NCC
         C%NCE = F%AMG%NCCE
         C%NW  = F%AMG%NCW

         ALLOCATE (C%WALL(1:C%NW), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'C%WALL', IERROR)

         ALLOCATE (C%AMG%CELL_TYPE(1:C%NCE), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'C%AMG%CELL_TYPE', IERROR)
         C%AMG%CELL_TYPE = NSCARC_UNDEFINED

         ALLOCATE (C%AMG%MEASURE(1:C%NCE), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'C%AMG%MEASURE', IERROR)
         C%AMG%MEASURE = NSCARC_MEASURE_NONE

         ALLOCATE(C%INTERNAL_BDRY_CELL(C%NC), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_WALLS','INTERNAL_BDRY_CELL',IERROR)
         C%INTERNAL_BDRY_CELL = NSCARC_UNDEFINED

         ALLOCATE(C%WALL_INDEX(C%NC, -3:3), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_WALLS','WALL_INDEX',IERROR)
         C%WALL_INDEX = 0

         IF (C%NCE > C%NC) THEN
            ALLOCATE(C%MAP%ICE_TO_IWG(C%NC+1:C%NCE), STAT=IERROR)
            CALL ChkMemErr('SCARC_SETUP_WALLS','ICE_TO_IWG',IERROR)

            ALLOCATE(C%EXT_PTR(C%NC+1:C%NCE), STAT=IERROR)
            CALL ChkMemErr('SCARC_SETUP_WALLS','EXT_PTR',IERROR)
         ENDIF

      ENDDO TRANSFER_MESHES_LOOP3

END SELECT

END SUBROUTINE SCARC_SETUP_SIZES



!> ------------------------------------------------------------------------------------------------
!> Initialize global 3D-solver methods (cg/mg)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSENING
INTEGER :: NL, NM, NOM
TYPE (SCARC_LEVEL_TYPE), POINTER :: L, F, C
TYPE (SCARC_LEVEL_TYPE), POINTER :: OL
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, S, ST, P, R, OP, OR
TYPE (SCARC_AMG_TYPE), POINTER :: AMG
TYPE (OSCARC_TYPE), POINTER :: OS

IF (.NOT.BAMG) RETURN
SCARC_ROUTINE = 'SCARC_SETUP_COARSENING'

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

      L  => SCARC(NM)%LEVEL(NL)
      A  => SCARC(NM)%LEVEL(NL)%A
      S  => SCARC(NM)%LEVEL(NL)%S
      ST => SCARC(NM)%LEVEL(NL)%ST

      ALLOCATE (AMG%MEASURE(1:L%NCE), STAT=IERROR)
      CALL CHKMEMERR (SCARC_ROUTINE, 'AMG%MEASURE', IERROR)
      AMG%MEASURE = NSCARC_MEASURE_NONE

      ALLOCATE (AMG%CELL_TYPE(1:L%NCE), STAT=IERROR)
      CALL CHKMEMERR (SCARC_ROUTINE, 'AMG%CELL_TYPE', IERROR)
      AMG%CELL_TYPE = NSCARC_UNDEFINED

      ALLOCATE (AMG%CELL_PTR(1:L%NCE), STAT=IERROR)
      CALL CHKMEMERR (SCARC_ROUTINE, 'AMG%CELL_PTR', IERROR)
      AMG%CELL_PTR = -1

      ALLOCATE (AMG%GRAPH(1:L%NCE), STAT=IERROR)
      CALL CHKMEMERR (SCARC_ROUTINE, 'AMG%GRAPH', IERROR)
      AMG%GRAPH = NSCARC_UNDEFINED

      CALL SCARC_COPY_MATRIX(A, S, 'S')
      CALL SCARC_ALLOCATE_MATRIX(ST, A%NAC, A%NAR, A%NAC, A%NSTENCIL, 'ST')

      OTHER_MESH_INDEX1: DO NOM = 1, NMESHES

         IF (NOM == NM) CYCLE OTHER_MESH_INDEX1

         OS  => SCARC(NM)%OSCARC(NOM)
         IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE OTHER_MESH_INDEX1

         OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

         ALLOCATE (OL%AMG%MEASURE(1:OL%NCG), STAT=IERROR)
         CALL CHKMEMERR (SCARC_ROUTINE, 'OL%AMG%MEASURE', IERROR)
         OL%AMG%MEASURE = NSCARC_MEASURE_NONE

         ALLOCATE (OL%AMG%CELL_TYPE(1:OL%NCG), STAT=IERROR)
         CALL CHKMEMERR (SCARC_ROUTINE, 'OL%AMG%CELL_TYPE', IERROR)
         OL%AMG%CELL_TYPE = NSCARC_UNDEFINED

         ALLOCATE (OL%AMG%GRAPH(1:OL%NCG), STAT=IERROR)
         CALL CHKMEMERR (SCARC_ROUTINE, 'OL%AMG%GRAPH', IERROR)
         OL%AMG%GRAPH = NSCARC_UNDEFINED

      ENDDO OTHER_MESH_INDEX1

   ENDDO

   !ENDIF

   CALL SCARC_SETUP_STRENGTH_MATRIX  (NSCARC_COARSENING_RS3, NL)

   !> Then set measures and CELL_TYPEs on internal cells due to chosen coarsening strategy
   SELECT CASE (TYPE_COARSENING)
      CASE (NSCARC_COARSENING_GMG)
         CALL SCARC_SETUP_CELL_TYPES (NSCARC_COARSENING_GMG, NL)
      CASE (NSCARC_COARSENING_RS3)
         CALL SCARC_SETUP_COLORING  (NSCARC_COARSENING_RS3, NL)
      CASE (NSCARC_COARSENING_FALGOUT)
         CALL SCARC_SETUP_COLORING  (NSCARC_COARSENING_RS3, NL)
         CALL SCARC_SETUP_COLORING  (NSCARC_COARSENING_FALGOUT, NL)
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
      A => SCARC(NM)%LEVEL(NL)%A
      P => SCARC(NM)%LEVEL(NL)%P
      R => SCARC(NM)%LEVEL(NL)%R

      !> allocate prolongation and restriction matrix 
      CALL SCARC_ALLOCATE_MATRIX(P, AMG%NPE+10, AMG%NPE+10, AMG%NCE +10, -1, 'P')
      CALL SCARC_ALLOCATE_MATRIX(R, AMG%NRE+10, AMG%NRE+10, AMG%NCCE+10, -1, 'P')

      !> allocate auxiliary tag arrays to mark relevant positions in A and P
      ALLOCATE (A%TAG(L%NCE), STAT=IERROR)
      CALL CHKMEMERR (SCARC_ROUTINE, 'A%TAG', IERROR)
      A%TAG = 0

      ALLOCATE (P%TAG(AMG%NCCE), STAT=IERROR)
      CALL CHKMEMERR (SCARC_ROUTINE, 'P_TAG', IERROR)
      P%TAG = 0

      OTHER_MESH_INDEX2: DO NOM = 1, NMESHES

         IF (NOM == NM) CYCLE OTHER_MESH_INDEX2
         OS => SCARC(NM)%OSCARC(NOM)

         IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE OTHER_MESH_INDEX2
         OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
         OP => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%P
         OR => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%R

         CALL SCARC_ALLOCATE_MATRIX(OP, OL%AMG%NP+10, OL%AMG%NP+10, OL%NC+10, -1, 'OP')
         CALL SCARC_ALLOCATE_MATRIX(OR, OL%AMG%NP+10, OL%AMG%NP+10, OL%NC+10, -1, 'OR')   ! wirklich OL%NP statt OL%NR?

      ENDDO OTHER_MESH_INDEX2
   ENDDO

   !>
   !> determine prolongation and restriction matrix
   !> set corresponding overlap information between neighboring meshes for coarser grid level
   !>
   IF (TYPE_COARSENING >= NSCARC_COARSENING_GMG) THEN
      CALL SCARC_SETUP_RESTRICTION(NL)
      CALL SCARC_SETUP_PROLONGATION(NL)
      CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_GMG, NL)
   ELSE
      CALL SCARC_SETUP_PROLONGATION(NL)
      CALL SCARC_SETUP_OVERLAPS (NSCARC_MATRIX_PROLONGATION, NL)
      CALL SCARC_MATRIX_TRANSPOSE(F%P%VAL, F%P%ROW, F%P%COL, &
                                  F%R%VAL, F%R%ROW, F%R%COL, F%NCE, F%AMG%NCCE )
   ENDIF

   !> -----------------------------------------------------------------------------------------
   !> Allocate coarse grid matrix including pointer arrays
   !> Note: number of cells on coarse level corresponds to number of c-points on fine level
   !> Compute coarse grid matrix by multiplication with restriction and prolongation matrix:
   !>  A_coarse := R * A_fine * P
   !> -----------------------------------------------------------------------------------------
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      F => SCARC(NM)%LEVEL(NL)               !> Pointer to fine level
      C => SCARC(NM)%LEVEL(NL+1)             !> Pointer to coarse level

      C%AMG%NC  = F%AMG%NCC
      C%AMG%NCE = F%AMG%NCCE
      C%NW      = F%NCW

      !> Allocate WALL-information for coarser grids (own and neighboring grids)
      ALLOCATE(C%MAP%ICG_TO_ICE(C%NC+1:C%NCE), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_COARSENING','ICG_TO_ICE',IERROR)
      C%MAP%ICG_TO_ICE = 0

      IF (TYPE_COARSENING >= NSCARC_COARSENING_GMG) THEN
         ALLOCATE(C%XCOR(0:C%NX), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_CELL_TYPE','C%XCOR',IERROR)

         ALLOCATE(C%YCOR(0:C%NY), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_CELL_TYPE','C%YCOR',IERROR)

         ALLOCATE(C%ZCOR(0:C%NZ), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_CELL_TYPE','C%ZCOR',IERROR)

         CALL SCARC_SETUP_COORDINATES_AMG(NM, NL)
      ENDIF

      !DO IC = F%NC+1,F%NCE
      !>   ICC = F%AMG%CELL_TYPE(IC)
      !ENDDO


   ENDDO
   IF (NMESHES>1) CALL SCARC_SETUP_WALLS_AMG (NL)          !ACHTUNG: HIER NOCHMAL CHECKEN ----

   CALL SCARC_SETUP_SUBDIVISION_AMG(NL)                       !HIER NOCHMAL CHECKEN ----

   CALL SCARC_TRANSFER_MATRIX (NL)
   CALL SCARC_SETUP_SYSTEM_AMG (NL)

   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SUBDIAG, NL+1)

   !IF (NMESHES > 1) THEN
   !>   IF (TYPE_COARSENING == NSCARC_COARSENING_GMG)  THEN
   !>      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SYSTEM, NL+1)
   !>   ENDIF
   !ENDIF

!>! ACHTUNG : Deallocate auxiliary arrays ----!>!
   IF (C%NC <= 4) THEN
      NLEVEL_MAX = NL + 1
      EXIT LEVEL_LOOP
   ENDIF

ENDDO LEVEL_LOOP

!DO NL = NLEVEL_MIN, NLEVEL_MAX
!>  DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!>     F => SCARC(NM)%LEVEL(NL)
!>     DEALLOCATE(F%AMG%MEASURE)
!>     DEALLOCATE(F%AMG%CELL_TYPE)
!>  ENDDO
!ENDDO

END SUBROUTINE SCARC_SETUP_COARSENING



!> ----------------------------------------------------------------------------------------------------
!> Store subdivision information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COORDINATES_AMG(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER IXC, IYC, IZC
TYPE (SCARC_LEVEL_TYPE), POINTER :: C, F

SCARC_ROUTINE = 'SCARC_SETUP_COORDINATES_AMG'

F => SCARC(NM)%LEVEL(NL)
C => SCARC(NM)%LEVEL(NL+1)

IF (TWO_D) THEN
   DO IXC = 0, C%NX-1
      C%XCOR(IXC)  = F%XCOR(2*IXC)
   ENDDO
   C%XCOR(IXC)  = F%XCOR(F%NX)
   DO IZC = 0, C%NZ-1
      C%ZCOR(IZC)  = F%ZCOR(2*IZC)
   ENDDO
   C%ZCOR(IZC)  = F%ZCOR(F%NZ)
ELSE
   DO IZC = 1, C%NZ
      DO IYC = 1, C%NY
         DO IXC = 1, C%NX
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_COORDINATES_AMG


!> ----------------------------------------------------------------------------------------------------
!> Store subdivision information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIVISION_AMG(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IOR_LAST, NOM, INBR
INTEGER :: NEIGHBORS(20,-3:3)
TYPE (SCARC_LEVEL_TYPE), POINTER :: C

SCARC_ROUTINE = 'SCARC_SETUP_SUBDIVISION_AMG'

IOR_LAST    = 0
NEIGHBORS   = 0

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   C => SCARC(NM)%LEVEL(NL+1)
   C%SUBDIVISION = 0

   WALL_CELLS_LOOP: DO IW = 1, C%NW

      IOR0 = C%WALL(IW)%IOR

      IF (IOR_LAST /= IOR0) C%SUBDIVISION(1,IOR0) = IW
      C%SUBDIVISION(2,IOR0) = C%SUBDIVISION(2,IOR0) + 1

      NOM= C%WALL(IW)%NOM

      IF (NOM /= 0) THEN
         NEIGHBOR_LOOP: DO INBR = 1, 20
            IF (NOM == NEIGHBORS(INBR, IOR0)) THEN
               EXIT NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) /= 0) THEN
               CYCLE NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) == 0) THEN
               NEIGHBORS(INBR, IOR0) = NOM
               C%SUBDIVISION(3,IOR0) = C%SUBDIVISION(3,IOR0) + 1
               EXIT NEIGHBOR_LOOP
            ELSE
               WRITE(SCARC_MESSAGE,'(2A)') TRIM(SCARC_ROUTINE),': More than 20 neighbors at one face not allowed yet!'
               CALL SHUTDOWN(SCARC_MESSAGE); RETURN
            ENDIF
         ENDDO NEIGHBOR_LOOP
      ENDIF

      IOR_LAST = IOR0

   ENDDO WALL_CELLS_LOOP
ENDDO MESHES_LOOP1

END SUBROUTINE SCARC_SETUP_SUBDIVISION_AMG



!> ------------------------------------------------------------------------------------------------
!> Determine if a cell IC is strongly coupled to another cell JC
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLS_AMG(NL)
INTEGER , INTENT(IN) :: NL
INTEGER :: NM, NOM, ICOL, IW, ICW, NCPL, JC
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: C, OC

SCARC_ROUTINE = 'SCARC_SETUP_WALLS_AMG'

!> -------------------------------------------------------------------------
!> Loop over all boundary cells IW of fine grid
!> Get corresponding adjacent and ghost cell
!> -------------------------------------------------------------------------
IF (TYPE_COARSENING < NSCARC_COARSENING_GMG) THEN
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   C => SCARC(NM)%LEVEL(NL+1)

   DO IW = 1, C%NW

      NOM = C%WALL(IW)%NOM
      OC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

      ICW = C%WALL(IW)%ICW

      IF (TYPE_COARSENING == NSCARC_COARSENING_GMG) NCPL = 1

      WRITE(*,*) 'L: HIER WALL%ICN �ndern auf OFI!'
      WRITE(*,*) 'L: MUSS NOCH UEBERARBEITET WERDEN !>'

!>     C%WALL(IW)%NCPL = NCPL

!>     ALLOCATE(C%WALL(IW)%ICN(NCPL), STAT=IERROR)
!>     CALL ChkMemErr('SCARC_SETUP_WALLS_AMG','ICN',IERROR)

!>     ALLOCATE(C%WALL(IW)%ICE(NCPL), STAT=IERROR)
!>     CALL ChkMemErr('SCARC_SETUP_WALLS_AMG','ICE',IERROR)

!>     ALLOCATE(C%WALL(IW)%ICG(NCPL), STAT=IERROR)
!>     CALL ChkMemErr('SCARC_SETUP_WALLS_AMG','ICG',IERROR)

      NCPL = 0
      DO ICOL = C%A%ROW(ICW)+1, C%A%ROW(ICW+1)-1
         JC = C%A%COL(ICOL)
         IF (JC > C%NC) THEN
            NCPL = NCPL + 1
!>           C%WALL(IW)%ICE(NCPL) = JC
!>           C%WALL(IW)%ICN(NCPL) = C%EXT_PTR(JC)
         ENDIF
      ENDDO

    ENDDO

ENDDO MESHES_LOOP

ENDIF

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MESHINFO, NL+1)

MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
      OS => SCARC(NM)%OSCARC(NOM)

      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE OTHER_MESHES_LOOP
      OC  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL+1)

      ALLOCATE(OC%WALL(1:OC%NW), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_WALLS_AMG','OC%WALL',IERROR)

   ENDDO OTHER_MESHES_LOOP
ENDDO MESHES_LOOP2

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WALLINFO, NL+1)

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
TYPE (SCARC_MATRIX_TYPE), POINTER :: A

SCARC_ROUTINE = 'SCARC_SETUP_STRENGTH_MATRIX'

!> Only dummy (NTYPE really used ?)
IC = NTYPE

STRENGTH_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)
   A => SCARC(NM)%LEVEL(NL)%A

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
      L%S%COL(IDIAG) = -1

      IF ((ROW_SUM > MAX_ROW_SUM) .AND. (MAX_ROW_SUM < 1.0_EB)) THEN
         !> set all dependencies to be weak
         DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
            JC = A%COL(ICOL)
            IF (JC <= L%NC) THEN
               L%S%COL(ICOL) = -1
            ENDIF
         ENDDO
      ELSE
         !> set dependencies to be weak related to threshold
         DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
            JC = A%COL(ICOL)
            IF (JC <= L%NC) THEN
               IF (A%VAL(ICOL) <= THRESHOLD * ROW_SCALE) L%S%COL(ICOL) = -1
            !ELSE
            !>   L%S%COL(ICOL) = -1
            ENDIF
          ENDDO
      ENDIF

   ENDDO STRENGTH_CELL_LOOP

   !> Compress strength matrix
   IS = 1
   STRENGTH_CELL_LOOP2: DO IC = 1, L%NC
      L%S%ROW(IC) = IS
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         !IF (L%S%COL(ICOL) > -1.AND.L%S%COL(ICOL)<=L%NC) THEN
         IF (L%S%COL(ICOL) > -1) THEN
            L%S%COL(IS) = L%S%COL(ICOL)
            IS = IS + 1
         ENDIF
      ENDDO
   ENDDO STRENGTH_CELL_LOOP2
   L%S%ROW(L%NC+1) = IS

   DO IC = 1, L%NCE+1
      L%ST%ROW(IC) = 0
   ENDDO

   IS = L%S%ROW(L%NC+1)-1
   DO ICOL = 1, IS
      L%ST%ROW(L%S%COL(ICOL)+1) = L%ST%ROW(L%S%COL(ICOL)+1) + 1
   ENDDO
   L%ST%ROW(1) = 1

   DO IC = 1, L%NCE
      L%ST%ROW(IC+1)= L%ST%ROW(IC+1) + L%ST%ROW(IC)
   ENDDO
   DO IC = 1, L%NC
      DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1
         IPOS = L%S%COL(ICOL)
         L%ST%COL(L%ST%ROW(IPOS)) = IC
         L%ST%ROW(IPOS) = L%ST%ROW(IPOS) + 1
      ENDDO
   ENDDO
   DO IC = L%NCE+1, 2, -1
      L%ST%ROW(IC) = L%ST%ROW(IC-1)
   ENDDO
   L%ST%ROW(1) = 1

ENDDO STRENGTH_MESHES_LOOP

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
LOGICAL :: BEMPTY=.FALSE., BNONEMPTY=.FALSE., BNEIGHBOR, BREAK = .TRUE., BCYCLE = .TRUE.
REAL(EB) :: RAND_NUM, MEASURE, NEW_MEASURE
REAL(EB) :: MEASURE_MAX, EPS
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_AMG_TYPE), POINTER :: AMG

SCARC_ROUTINE = 'SCARC_SETUP_COLORING'

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

         L => SCARC(NM)%LEVEL(NL)
         AMG => SCARC(NM)%LEVEL(NL)%AMG

         REMAINING_CELLS = 0

         !> Currently the measures are computed as row sums of ST (number of influences for IC)
         RS3_MEASURE_LOOP0: DO IC = 1, L%NC
            AMG%MEASURE(IC) = L%ST%ROW(IC+1)-L%ST%ROW(IC)
         ENDDO RS3_MEASURE_LOOP0

         !> Subdivide in fine and coarse grid cells
         RS3_MEASURE_LOOP1: DO IC = 1, L%NC

            IF (L%S%ROW(IC+1)-L%S%ROW(IC) == 0) THEN
               AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE0
               AMG%MEASURE(IC) = 0.0_EB
             !> IF (AGGRESSIVE2) AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE0
            ELSE
               AMG%CELL_TYPE(IC) = NSCARC_UNDEFINED
               REMAINING_CELLS   = REMAINING_CELLS + 1
            ENDIF
         ENDDO RS3_MEASURE_LOOP1

         RS3_MEASURE_LOOP2: DO IC = 1, L%NC
            MEASURE = AMG%MEASURE(IC)
            IF (AMG%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_FINE0 .AND. AMG%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_COARSE0) THEN
               IF (AMG%MEASURE(IC) > 0.0_EB) THEN
                  WRITE(*,*) 'ACHTUNG: KORREKT?'
                  !AMG%MEASURE(IC) = AMG%MEASURE(IC)
               ELSE
                  IF (AMG%MEASURE(IC) < 0.0_EB) WRITE(*,*) 'SCARC_SETUP_MEASURE: Negative measure !>'
                  AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
                  DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1
                     JC = L%S%COL(ICOL)
                     IF (AMG%CELL_TYPE(JC) /= NSCARC_CELL_TYPE_COARSE0 .AND. AMG%CELL_TYPE(JC) /= NSCARC_CELL_TYPE_FINE0) THEN
                        IF (JC < IC) THEN
                           NEW_MEASURE = AMG%MEASURE(JC)
                           IF (NEW_MEASURE > 0.0_EB) AMG%MEASURE(JC) = 0.0_EB
                           NEW_MEASURE = AMG%MEASURE(JC)+1
                           AMG%MEASURE(JC) = NEW_MEASURE
                        ELSE
                           NEW_MEASURE = AMG%MEASURE(JC)+1
                        ENDIF
                     ENDIF
                  ENDDO
                  REMAINING_CELLS = REMAINING_CELLS - 1
               ENDIF
            ENDIF
         ENDDO RS3_MEASURE_LOOP2

         RS3_CYCLE_LOOP: DO WHILE (REMAINING_CELLS > 0)


            !> get maximum (remaining) measure for all cells
            MEASURE_MAX = MAXVAL(AMG%MEASURE(1:L%NC))
            IF (MEASURE_MAX <= EPS) EXIT RS3_CYCLE_LOOP

            RS3_CELL_LOOP: DO IC = 1, L%NC

             !> Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, AMG%MEASURE(IC))) THEN

                  AMG%MEASURE(IC)  = NSCARC_MEASURE_NONE
                  AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
                  REMAINING_CELLS = REMAINING_CELLS - 1

                  !> Determine set of fine cells
                  DO ICOL = L%ST%ROW(IC), L%ST%ROW(IC+1)-1

                   !> IF JC hasn't been marked yet, set it to be a fine cell which is no longer measured
                     JC = L%ST%COL(ICOL)

                     IF (JC > L%NC) CYCLE
                     IF (AMG%CELL_TYPE(JC) == NSCARC_UNDEFINED) THEN

                        AMG%MEASURE(JC)  = NSCARC_MEASURE_NONE
                        AMG%CELL_TYPE(JC) = NSCARC_CELL_TYPE_FINE
                        REMAINING_CELLS = REMAINING_CELLS - 1

                        !>  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                        DO JCOL = L%S%ROW(JC), L%S%ROW(JC+1)-1
                           KC = L%S%COL(JCOL)
                           IF (AMG%CELL_TYPE(KC)==NSCARC_UNDEFINED) THEN
                              AMG%MEASURE(KC) = AMG%MEASURE(KC) + 1.0_EB
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
                  EXIT RS3_CELL_LOOP
               ENDIF
            ENDDO RS3_CELL_LOOP

            DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1
               JC = L%S%COL(ICOL)
               IF (JC > L%NC) CYCLE
               IF (AMG%CELL_TYPE(JC) == NSCARC_UNDEFINED) THEN
                  MEASURE = AMG%MEASURE(JC) - 1
                  AMG%MEASURE(JC) = MEASURE
                  IF (MEASURE > 0.0_EB) THEN
                     AMG%MEASURE(JC) = MEASURE
                  ELSE
                     AMG%CELL_TYPE(JC) = NSCARC_CELL_TYPE_FINE
                     REMAINING_CELLS = REMAINING_CELLS - 1
                     DO JCOL = L%S%ROW(JC), L%S%ROW(JC+1)-1
                        KC = L%S%COL(JCOL)
                        IF (AMG%CELL_TYPE(KC)==NSCARC_UNDEFINED) THEN
                           AMG%MEASURE(KC) = AMG%MEASURE(KC) + 1.0_EB
                           MEASURE_MAX = MAX(MEASURE_MAX, AMG%MEASURE(KC))
                        ENDIF
                     ENDDO
                  ENDIF

               ENDIF
            ENDDO

         ENDDO RS3_CYCLE_LOOP
         L%NCW = 0

         DO IC = 1, L%NC
            IF (AMG%CELL_TYPE(IC) == NSCARC_CELL_TYPE_COARSE0) AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
         ENDDO

         !> exchange information with neighboring meshes
         IF (NMESHES > 1) THEN

            DO IC = 1, L%NCE
               AMG%GRAPH(IC) = -1
            ENDDO

            IC0 = 1
            DO IC = 1, L%NC
               IF (ICT2 /= IC) ICT = -1
               IF (AMG%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE) THEN

                  DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1
                     JC = L%S%COL(ICOL)
                     IF (JC <= L%NC) THEN
                        IF (AMG%CELL_TYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                           AMG%GRAPH(JC) = IC
                        ENDIF
                     ENDIF
                  ENDDO

                  !> Hier fehlt noch die Abfrage nach Nachbarmesh !>!

                  DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1
                     JC = L%S%COL(ICOL)
                     IF (JC <= L%NC) THEN
                     IF (AMG%CELL_TYPE(JC) == NSCARC_CELL_TYPE_FINE) THEN

                        !> DIESER PART WIRD ERST FUER ANDERE VERFEINERUNGEN AKTIV
                        !> ACHTUNG: DANN NOCHMAL UEBERPRUEFEN!
                        BEMPTY = .TRUE.
                        DO JCOL = L%S%ROW(JC), L%S%ROW(JC+1)-1
                           KC = L%S%COL(JCOL)
                           IF (AMG%GRAPH(KC) == IC) THEN
                              BEMPTY = .FALSE.
                              EXIT
                           ENDIF
                           IF (BEMPTY) THEN
                              IF (BNONEMPTY) THEN
                                 AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
                                 IF (ICT > -1) THEN
                                    AMG%CELL_TYPE(ICT) = NSCARC_CELL_TYPE_FINE
                                    ICT = -1
                                 ENDIF
                                 !> Hier fehlt noch Nachbaranteil
                                 BNONEMPTY = .FALSE.
                                 EXIT
                              ELSE
                                 ICT  = JC
                                 ICT2 = IC
                                 AMG%CELL_TYPE(JC) = NSCARC_CELL_TYPE_COARSE
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
      ENDDO RS3_MESH_INDEX

      !> Exchange CELL_TYPE-data along internal boundaries
      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)

      !> Third pass: identify fine cells along boundary and get their coarse neighbors
      RS3_MESH_INDEX2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         DO IC = 1, L%NC
            AMG%GRAPH(IC) = -1
         ENDDO
      ENDDO RS3_MESH_INDEX2


!> ----------------------------------------------------------------------------------------
!> Falgout coarsening
!> ----------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_FALGOUT)

      FCELL = ZCELL
      FALGOUT_MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)

         !> Reset the measures as row sums of ST (number of influences for IC)
         !> plus random number to make them unique
         FALGOUT_CELL_LOOP1: DO IC = 1, L%NC
            RAND_NUM = 0.01_EB
            DO IRAND = 1, 5
               CALL RANDOM_NUMBER(RAND_NUM)
               RAND_NUM = RAND_NUM + RAND_NUM/10**(IRAND-1)
               AMG%MEASURE(IC) = AMG%MEASURE(IC) + RAND_NUM
            ENDDO
            AMG%MEASURE(IC) = L%S%ROW(IC+1)-L%S%ROW(IC) + RAND_NUM
         ENDDO FALGOUT_CELL_LOOP1

      ENDDO FALGOUT_MESHES_LOOP1

      !> Initial exchange of measure array
      IF (NMESHES > 1) THEN
         TYPE_VECTOR = NSCARC_VECTOR_MEASURE
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR , NL)
      ENDIF

      FALGOUT_INIT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         !> reset CELL_TYPE for cells with neighbors in other meshes
         FALGOUT_INTERNAL_CELL_LOOP1: DO IC = 1, L%NC
            BNEIGHBOR = .FALSE.
            DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1
               JC = L%S%COL(ICOL)
               IF (JC > L%NC) THEN
                  BNEIGHBOR = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (BNEIGHBOR .OR. AMG%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE) AMG%CELL_TYPE(IC)=NSCARC_UNDEFINED
         ENDDO FALGOUT_INTERNAL_CELL_LOOP1

         !> initialize GRAPH and reset CELL_TYPE on ghost cells
         FALGOUT_EXTENDED_CELL_LOOP: DO IC = L%NC+1, L%NCE
            AMG%GRAPH(IC) = IC
            AMG%CELL_TYPE(IC) = NSCARC_UNDEFINED
         ENDDO FALGOUT_EXTENDED_CELL_LOOP

         !> reset CELL_TYPE on internal wall cells
         DO IW = 1, L%NW
            IF (L%WALL(IW)%NOM /= 0) THEN
               IC = L%WALL(IW)%ICW
               AMG%CELL_TYPE(IC) = NSCARC_UNDEFINED
            ENDIF
         ENDDO

         !> reset CELL_TYPE on internal fine cells
         ICG = 1
         FALGOUT_INTERNAL_CELL_LOOP2: DO IC = 1, L%NC
            IF (AMG%CELL_TYPE(IC)<NSCARC_UNDEFINED) THEN
               AMG%CELL_TYPE(IC)=NSCARC_UNDEFINED
            ENDIF
            IF (AMG%CELL_TYPE(IC) == NSCARC_CELL_TYPE_ZPNT) THEN
               IF (AMG%MEASURE(IC) >= 1.0_EB .OR. (L%S%ROW(IC+1)-L%S%ROW(IC)) > 0) THEN
                  AMG%CELL_TYPE(IC) = NSCARC_UNDEFINED
                  AMG%GRAPH(ICG) = IC
                  ICG = ICG + 1
               ELSE
                  AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
               ENDIF
            ELSE IF (AMG%CELL_TYPE(IC) == NSCARC_CELL_TYPE_SFPNT) THEN
               AMG%MEASURE(IC) = NSCARC_MEASURE_NONE
            ELSE
               AMG%GRAPH(ICG) = IC
               ICG = ICG + 1
            ENDIF
         ENDDO FALGOUT_INTERNAL_CELL_LOOP2

         IGRAPH  = ICG-1
         IGRAPHE = L%NCE-L%NC

      ENDDO FALGOUT_INIT_LOOP

      FALGOUT_EXTERNAL_LOOP: DO IC = L%NC+1, L%NCE
         AMG%MEASURE(IC) = NSCARC_MEASURE_NONE
      ENDDO FALGOUT_EXTERNAL_LOOP

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

                  IC = AMG%GRAPH(ICG)

                  !> if cell isn't marked as coarse cell yet and has measure
                  !> less than 1, mark it as fine cell
                  !> take care that all dependencies have been taken into account
                  IF (AMG%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_COARSE .AND. &
                      AMG%MEASURE(IC)  <  NSCARC_MEASURE_ONE) THEN
                     AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
                     DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1
                        JC = L%S%COL(ICOL)
                        IF (JC < 0) CYCLE
                        AMG%CELL_TYPE(IC) = NSCARC_UNDEFINED
                     ENDDO
                  ENDIF

                  !> if cell is already marked as fine or coarse, set its measure to zero
                  !> and extract it from the graph (put it at the end of the graph
                  !> array and decrease number of relevant graph entries)
                  IF (AMG%CELL_TYPE(IC) /= NSCARC_UNDEFINED) THEN
                     AMG%MEASURE(IC) = NSCARC_MEASURE_NONE
                     AMG%GRAPH(ICG) = AMG%GRAPH(IGRAPH)
                     AMG%GRAPH(IGRAPH) = IC
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
               IC = AMG%GRAPH(ICG)

               IF (IC < 0) CYCLE
               IF (AMG%CELL_TYPE(IC) < NSCARC_UNDEFINED) THEN
                  AMG%GRAPH(ICG) = AMG%GRAPH(L%NC+IGRAPHE)
                  AMG%GRAPH(L%NC+IGRAPHE) = IC
                  IGRAPHE = IGRAPHE - 1
               ENDIF
               ICG = ICG + 1
            ENDDO FALGOUT_GRAPH_LOOP2

            AMG%MEASURE(L%NC+1:L%NCE) = NSCARC_MEASURE_NONE

            FALGOUT_GRAPH_LOOP3: DO ICG = 1, IGRAPH

               IC = AMG%GRAPH(ICG)

               !> Coarse cells don't interpolate from influencing neighbors
               IF (AMG%CELL_TYPE(IC) > NSCARC_UNDEFINED) THEN

                  AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE

                  DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1
                     JC = L%S%COL(ICOL)

                     !> remove edge from S and decrement measures of unmarked neighbors
                     IF (JC > -1) THEN
                        L%S%COL(ICOL) = - L%S%COL(ICOL) - 1
                        IF (AMG%CELL_TYPE(JC) == NSCARC_UNDEFINED) AMG%MEASURE(JC) = AMG%MEASURE(JC) - 1
                     ENDIF
                  ENDDO

               ELSE

                  !> dependencies which are already marked
                  DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1

                     JC = L%S%COL(ICOL)
                     IF (JC < 0) JC = -JC - 1

                     !> remove edge from S and temporarily reset CELL_TYPE
                     IF (AMG%CELL_TYPE(JC) > NSCARC_UNDEFINED) THEN
                        IF (L%S%COL(ICOL) > -1) L%S%COL(ICOL) = -L%S%COL(ICOL) - 1
                        AMG%CELL_TYPE(JC) = NSCARC_CELL_TYPE_CPNT
                     ELSE IF (AMG%CELL_TYPE(JC) == NSCARC_CELL_TYPE_SFPNT) THEN    !> necessary ??
                        IF (L%S%COL(ICOL) > -1) L%S%COL(ICOL) = -L%S%COL(ICOL) - 1
                     ENDIF
                  ENDDO

                  !> dependencies which aren't marked yet
                  DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1

                     JC = L%S%COL(ICOL)
                     IF (JC > -1 .AND. JC<=L%NC) THEN

                        BREAK = .TRUE.

                        !> check if there are common C-points
                        DO JCOL = L%S%ROW(JC), L%S%ROW(JC+1)-1
                           KC = L%S%COL(JCOL)
                           IF (KC < 0) KC = -KC - 1
                           IF (KC <= L%NC) THEN
                              !> remove edge from S and update measure
                              IF (AMG%CELL_TYPE(KC) == NSCARC_CELL_TYPE_CPNT) THEN
                                 L%S%COL(ICOL) = - L%S%COL(ICOL) - 1
                                 AMG%MEASURE(JC) = AMG%MEASURE(JC) - 1
                                 BREAK = .FALSE.
                                 EXIT
                              ENDIF
                           ENDIF
                        ENDDO

                        IF (BREAK) THEN
                           DO JCOL = L%S%ROW(JC), L%S%ROW(JC+1)-1
                              KC = L%S%COL(JCOL)
                              IF (KC < 0) KC = -KC - 1
                              IF (KC > L%NC) THEN
                                 !> remove edge from S and update measure
                                 IF (AMG%CELL_TYPE(KC) == NSCARC_CELL_TYPE_CPNT) THEN
                                    L%S%COL(ICOL) = - L%S%COL(ICOL) - 1
                                    AMG%MEASURE(JC) = AMG%MEASURE(JC) - 1
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
               DO ICOL = L%S%ROW(IC), L%S%ROW(IC+1)-1
                  JC = L%S%COL(ICOL)
                  IF (JC < 1) JC = -JC - 1
                  IF (AMG%CELL_TYPE(JC) == NSCARC_CELL_TYPE_CPNT) AMG%CELL_TYPE(JC) = NSCARC_CELL_TYPE_COARSE
               ENDDO

            ENDDO FALGOUT_GRAPH_LOOP3

         ENDDO FALGOUT_GRAPH_LOOP

       !> reset S-matrix
         DO ICOL = 1, L%S%ROW(L%NC+1)
            IF (L%S%COL(ICOL) < 0) L%S%COL(ICOL) = -L%S%COL(ICOL)-1
         ENDDO

      ENDDO FALGOUT_COLORING_LOOP

END SELECT


END SUBROUTINE SCARC_SETUP_COLORING


!> ----------------------------------------------------------------------------------------
!> Select an independent set from the graph
!> ----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRAPHSET(IGRAPH, NM, NL)
INTEGER, INTENT(IN):: IGRAPH, NM, NL
INTEGER :: IG, IC, ICOL, JC
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A
TYPE (SCARC_AMG_TYPE), POINTER :: AMG

SCARC_ROUTINE = 'SCARC_SETUP_GRAPHSET'

L => SCARC(NM)%LEVEL(NL)
A => SCARC(NM)%LEVEL(NL)%A
AMG => SCARC(NM)%LEVEL(NL)%AMG

!> First mark every cell from the 'internal' graphset with measure bigger
!> than one as a coarse cell
DO IG = 1, IGRAPH
   IC = AMG%GRAPH(IG)
   IF (AMG%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
   ENDIF
ENDDO

!> Do the same with cells on the overlapping areas
DO IG = L%NC+1, L%NCE
   IC = AMG%GRAPH(IG)
   IF (IC < 0) CYCLE
   IF (AMG%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
   ENDIF
ENDDO

!> remove nodes from the initial independent set depending on their measure
!> For each cell consider every connection in stencil and set cell with
!> 1G/smaller measure to Zero-Type
DO IG = 1, IGRAPH
   IC = AMG%GRAPH(IG)

   IF (IC < 0) CYCLE

   IF (AMG%MEASURE(IC) > NSCARC_MEASURE_ONE) THEN
      DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
         JC = A%COL(ICOL)
         IF (JC < 0) CYCLE
         IF (AMG%MEASURE(JC) > NSCARC_MEASURE_ONE) THEN
            IF (AMG%MEASURE(IC) > AMG%MEASURE(JC)) THEN
               AMG%CELL_TYPE(JC) = NSCARC_UNDEFINED
            ELSE IF (AMG%MEASURE(JC) > AMG%MEASURE(IC)) THEN
               AMG%CELL_TYPE(IC) = NSCARC_UNDEFINED
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
         CELL_TYPE(PC(INX,INZ)) = NSCARC_CELL_TYPE_SFINE
      ENDDO
   ENDDO
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
LOGICAL  :: BEVENX, BEVENZ, BTWOLEVEL_X, BTHREE_X, BTWOLEVEL_Z, BTHREE_Z, BREMOVE=.FALSE.
TYPE (SCARC_LEVEL_TYPE), POINTER :: L, F, C
TYPE (SCARC_MATRIX_TYPE), POINTER :: A
TYPE (SCARC_AMG_TYPE), POINTER :: AMG

SCARC_ROUTINE = 'SCARC_SETUP_CELL_TYPES'

EPS = 1.0E-12
MEASURE_MAX = 0.0_EB


!> Define CELL_TYPEs for corresponding coarsening strategy
SELECT CASE (NTYPE)

   !> ---------------------------------------------------------------------------------------------
   !> standard Rue-Stueben coarsening
   !>  - identify first point IC with maximum measure and set it to be a C-point
   !>  - identify surrounding strongly connected points JC and set them to be F-points
   !>  - increase measures of strongly connected neighbours KC of JC
   !>  - decrease measures of strongly connected neighbours of IC
   !> ---------------------------------------------------------------------------------------------

   CASE (NSCARC_COARSENING_RS3)

      MEASURE_TYPE = 1
      RS3_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L   => SCARC(NM)%LEVEL(NL)
         A   => SCARC(NM)%LEVEL(NL)%A
         AMG => SCARC(NM)%LEVEL(NL)%AMG

         REMAINING_CELLS = L%NC

         RS3_CYCLE_LOOP: DO WHILE (REMAINING_CELLS > 0)

            !> get maximum (remaining) measure for all cells
            MEASURE_MAX = MAXVAL(AMG%MEASURE(1:L%NC))
            IF (MEASURE_MAX <= EPS) EXIT RS3_CYCLE_LOOP

            RS3_CELL_LOOP: DO IC = 1, L%NC

               !> Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, AMG%MEASURE(IC))) THEN

                  AMG%MEASURE(IC)  = NSCARC_MEASURE_NONE
                  AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE
                  REMAINING_CELLS = REMAINING_CELLS - 1

                  !> Determine set of fine cells
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1

                     IF (STRONGLY_COUPLED(A, IC, ICOL)) THEN

                        !> IF JC hasn't been marked yet, set it to be a fine cell which is no longer measured
                        JC = A%COL(ICOL)
                        IF (AMG%CELL_TYPE(JC) == NSCARC_UNDEFINED) THEN

                           AMG%MEASURE(JC)  = NSCARC_MEASURE_NONE
                           AMG%CELL_TYPE(JC) = NSCARC_CELL_TYPE_FINE
                           REMAINING_CELLS = REMAINING_CELLS - 1

                           !>  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                           DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
                              IF (STRONGLY_COUPLED(A, JC, JCOL)) THEN
                                 KC = A%COL(JCOL)
                                 IF (KC /= 0) THEN
                                    IF (KC /= IC .AND. (AMG%CELL_TYPE(KC)==NSCARC_UNDEFINED)) THEN
                                       AMG%MEASURE(KC) = AMG%MEASURE(KC) + 1.0_EB
                                       MEASURE_MAX = MAX(MEASURE_MAX, AMG%MEASURE(KC))
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDDO

                        ENDIF
                     ENDIF
                  ENDDO
                  EXIT RS3_CELL_LOOP
               ENDIF
            ENDDO RS3_CELL_LOOP
         ENDDO RS3_CYCLE_LOOP
         L%NCW = 0
      ENDDO RS3_LOOP

   !> ---------------------------------------------------------------------------------------------
   !> set CELL_TYPEs for Falgout coarsening
   !> Care: ZPOINT is equal to FINE in this case !>!
   !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_FALGOUT)

      FALGOUT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)

         !> copy column pointers for matrix
         DO IA = 1, L%NA
            L%S%COL(IA) = A%COL(IA)
         ENDDO

         !> this part may be changed for other coarsening strategies 
         ICOUNT = 1
         FALGOUT_CELL_LOOP: DO IC = 1, L%NC
            IF (AMG%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE) THEN
               IF (AMG%MEASURE(IC) >= 1.0_EB .AND. A%ROW(IC+1)-A%ROW(IC)+1 > 0) THEN
                  AMG%CELL_TYPE(IC)  = NSCARC_UNDEFINED
                  AMG%GRAPH(ICOUNT) = IC
                  ICOUNT = ICOUNT + 1
                  IF (ICOUNT > L%NC) THEN
                     WRITE(SCARC_MESSAGE,'(2A,I8)') TRIM(SCARC_ROUTINE),': Graph size must be increased'
                     CALL SHUTDOWN(SCARC_MESSAGE); RETURN
                  ENDIF
               ELSE
                  AMG%CELL_TYPE(IC)  = NSCARC_CELL_TYPE_FINE
               ENDIF
            ELSEIF (AMG%CELL_TYPE(IC) == NSCARC_CELL_TYPE_FINE0) THEN
               AMG%MEASURE(IC) = NSCARC_MEASURE_NONE
            ELSE
               AMG%GRAPH(ICOUNT) = IC
               ICOUNT = ICOUNT + 1
            ENDIF
         ENDDO FALGOUT_CELL_LOOP

         IGRAPH = ICOUNT
         ILOOP  = 0
         FALGOUT_COARSENING_LOOP: DO

            !> Care: Is graph_size always L%NC ?
            IF (ILOOP > 0) THEN

               DO IG = 1, IGRAPH

                  IC = AMG%GRAPH(IG)
                  JC = IG

                  !> make IC a fine cell and look for all dependencies
                  IF ((AMG%CELL_TYPE(IC) /= NSCARC_CELL_TYPE_COARSE) .AND. (AMG%MEASURE(IC) < 1)) THEN
                     AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_FINE
                     DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
                        IF (L%S%COL(IROW) > -1) AMG%CELL_TYPE(IC) = NSCARC_UNDEFINED
                     ENDDO
                  ENDIF

                  !> remove cell from graph
                  IF (AMG%CELL_TYPE(IC) /= NSCARC_UNDEFINED) THEN
                     AMG%MEASURE(JC) = NSCARC_MEASURE_NONE
                     IGRAPH = IGRAPH - 1
                     AMG%GRAPH(JC) = AMG%GRAPH(IGRAPH)
                     AMG%GRAPH(IGRAPH) = IC
                     JC = JC - 1               !> correct ???
                  ENDIF

               ENDDO

            ENDIF
            ILOOP = ILOOP + 1

            FALGOUT_GRAPH_LOOP: DO IG = 1, IGRAPH

               IC = AMG%GRAPH(IG)

               !> C-points are not influenced by neighbors (no interpolation)
               IF (AMG%CELL_TYPE(IC) > 0) THEN

                  AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE            !> define it as C-point

                  DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
                     JC = L%S%COL(IROW)
                     IF (JC > -1) THEN
                        L%S%COL(IROW) = -L%S%COL(IROW) - 1  !> remove edge
                        IF (AMG%CELL_TYPE(JC) == NSCARC_UNDEFINED) THEN
                           AMG%MEASURE(JC) = AMG%MEASURE(JC) - 1.0_EB   !> decrement measure of unmarked neigh.
                        ENDIF
                     ENDIF
                  ENDDO

               ELSE

                  !> dependencies which have already been marked
                  DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
                     JC = L%S%COL(IROW)
                     IF (JC < 0) JC = -JC-1
                     IF (AMG%CELL_TYPE(JC) > 0) THEN
                        IF (L%S%COL(IROW) > -1) THEN
                           L%S%COL(IROW) = -L%S%COL(IROW) -1   !> remove edge
                           AMG%CELL_TYPE(JC)    = NSCARC_CELL_TYPE_COMMON   !> temporarily set CELL_TYPE to common
                        ENDIF
                     ELSEIF (AMG%CELL_TYPE(JC) == NSCARC_CELL_TYPE_FINE0) THEN
                        IF (L%S%COL(IROW) > -1) THEN
                           L%S%COL(IROW) = -L%S%COL(IROW) -1   !> remove edge
                        ENDIF
                     ENDIF
                  ENDDO

                  !> dependencies which haven't been marked yet
                  DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
                     IF (L%S%COL(IROW) > -1) THEN
                        BREMOVE = .TRUE.
                        JC = L%S%COL(IROW)
                        DO JROW = A%ROW(JC)+1, A%ROW(JC+1)-1
                           KC = L%S%COL(JROW)
                           IF (KC < 0) KC = -KC-1                        !> check for all dependencies !>
                           IF (AMG%CELL_TYPE(KC) == NSCARC_CELL_TYPE_COMMON) THEN
                              L%S%COL(IROW) = -L%S%COL(IROW)-1
                              AMG%MEASURE(JC) = AMG%MEASURE(JC) - 1.0_EB
                              BREMOVE = .FALSE.
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO

               ENDIF

               !> reset CELL_TYPES
               DO IROW = A%ROW(IC)+1, A%ROW(IC+1)-1
                  JC = L%S%COL(IROW)
                  IF (JC < 0) JC = -JC-1
                  IF (AMG%CELL_TYPE(JC) == NSCARC_CELL_TYPE_COMMON) AMG%CELL_TYPE(JC)=NSCARC_CELL_TYPE_COARSE
               ENDDO

            ENDDO FALGOUT_GRAPH_LOOP
         ENDDO FALGOUT_COARSENING_LOOP

      ENDDO FALGOUT_LOOP

   !> ---------------------------------------------------------------------------------------------
   !> set CELL_TYPEs for GMG-like interpolation
   !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_GMG)

      IF (TWO_D) THEN

         GMG_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

            F => SCARC(NM)%LEVEL(NL)
            C => SCARC(NM)%LEVEL(NL+1)

            !> Analyze grid sizes
            INCRX = 2
            C%NX = F%NX/2
            IF (MOD(F%NX,2) == 0) THEN
               BEVENX=.TRUE.
            ELSE
               BEVENX=.FALSE.
            ENDIF

            INCRZ  = 2
            C%NZ = F%NZ/2
            IF (MOD(F%NZ,2) == 0) THEN
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

            C%NY = 1
            C%NC = C%NX * C%NZ

            SELECT CASE(ICASE)
               CASE (0)
                  DO IZ = 1, F%NZ ,2
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, 1,F%NX  ,1,1,IZ,2,1,2)
                  ENDDO
               CASE (1)
                  DO IZ = 1, F%NZ-3,2
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, 1,F%NX  ,1,1,IZ,2,1,2)
                  ENDDO
                  CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, 1,F%NX  ,1,1,F%NZ-2,2,1,3)
               CASE (2)
                  DO IZ = 1, F%NZ,2
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, 1       ,F%NX-3,1,1,IZ,2,1,2)
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, F%NX-2,F%NX-2,1,1,IZ,3,1,2)
                  ENDDO
               CASE (3)
                  DO IZ = 1, F%NZ-3,2
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, 1       ,F%NX-3,1,1,IZ,2,1,2)
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, F%NX-2,F%NX-2,1,1,IZ,3,1,2)
                  ENDDO
                  CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, 1       ,F%NX-3,1,1,F%NZ-2,2,1,3)
                  CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, F%NX-2,F%NX-2,1,1,F%NZ-2,3,1,3)
            END SELECT

            IF (NMESHES > 1) THEN
               F%NCW = 0
               DO IC = 1, F%NC
                  IF (F%AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     DO IOR0 = -3, 3
                        IW = F%WALL_INDEX(IC, IOR0)
                        IF (F%WALL(IW)%NOM /= 0) F%NCW = F%NCW + 1
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF

         ENDDO GMG_LOOP
      ELSE
         WRITE(*,*) 'GMG_COARSENING FOR 3D not yet implemented'
      ENDIF

   !> ---------------------------------------------------------------------------------------------
   !> set CELL_TYPEs for GMG-like interpolation
   !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_GMG3)

      IF (TWO_D) THEN

         GMG3_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

            F => SCARC(NM)%LEVEL(NL)
            C => SCARC(NM)%LEVEL(NL+1)

            !> Analyze grid sizes
            IF (MOD(F%NX,3) == 0) THEN
               INCRX    = 3
               C%NX   = F%NX/3
               BTWOLEVEL_X   = .FALSE.
               BTHREE_X = .TRUE.
            ELSE IF (MOD(F%NX,2) == 0) THEN
               INCRX    = 2
               C%NX   = F%NX/2
               BTWOLEVEL_X   = .TRUE.
               BTHREE_X = .FALSE.
            ELSE
               INCRX    = 2
               C%NX   = F%NX/2 + 1
               BTWOLEVEL_X   = .FALSE.
               BTHREE_X = .FALSE.
            ENDIF

            IF (MOD(F%NZ,3) == 0) THEN
               INCRZ    = 3
               C%NZ   = F%NZ/3
               BTWOLEVEL_Z   = .FALSE.
               BTHREE_Z = .TRUE.
            ELSE IF (MOD(F%NZ,2) == 0) THEN
               INCRZ    = 2
               C%NZ   = F%NZ/2
               BTWOLEVEL_Z   = .TRUE.
               BTHREE_Z = .FALSE.
            ELSE
               INCRZ    = 2
               C%NZ   = F%NZ/2 + 1
               BTWOLEVEL_Z   = .FALSE.
               BTHREE_Z = .FALSE.
            ENDIF

            IF (.NOT.(BTWOLEVEL_X.OR.BTHREE_X).OR..NOT.(BTWOLEVEL_Z.OR.BTHREE_Z)) THEN
               WRITE(SCARC_MESSAGE,'(2A,I8)') TRIM(SCARC_ROUTINE),': Cell numbers not divisable by 2 or 3!'
               CALL SHUTDOWN(SCARC_MESSAGE); RETURN
            ENDIF

            IF (BTHREE_X) THEN
               ICASE = 0
            ELSE IF (BTWOLEVEL_X) THEN
               ICASE = 1
            ELSE
               ICASE = 2
            ENDIF

            C%NY = 1
            C%NC = C%NX * C%NZ

            SELECT CASE(ICASE)
               CASE (0)
                  DO IZ = 1, F%NZ, INCRZ
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, 1,F%NX  ,1,1,IZ,3,1,INCRZ)
                  ENDDO
               CASE (1)
                  DO IZ = 1, F%NZ, INCRZ
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, 1,F%NX  ,1,1,IZ,2,1,INCRZ)
                  ENDDO
               CASE (3)
                  DO IZ = 1, F%NZ, INCRZ
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, 1       ,F%NX-3,1,1,IZ,2,1,INCRZ)
                     CALL  SCARC_PATCH_CELL_TYPES(F%AMG%CELL_TYPE,F%NX,F%NZ, F%NX-2,F%NX-2,1,1,IZ,3,1,INCRZ)
                  ENDDO
            END SELECT

            IF (NMESHES > 1) THEN
               F%NCW = 0
               DO IC = 1, F%NC
                  IF (F%AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     DO IOR0 = -3, 3
                        IW = F%WALL_INDEX(IC, IOR0)
                        IF (F%WALL(IW)%NOM /= 0) F%NCW = F%NCW + 1
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF

         ENDDO GMG3_LOOP

      ELSE
         WRITE(*,*) 'GMG3_COARSENING FOR 3D not yet implemented'
      ENDIF

   !> ---------------------------------------------------------------------------------------------
   !> first set CELL_TYPEs for cells in internal boundary layers (adjacent to neighboring meshes)
   !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_BDRY)

      BDRY_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)

         !> First define coarse cells on boundary layer(adjacent to neighboring meshes)
         BDRY_CYCLE_LOOP: DO

            !>!
            !> get maximum (remaining) measure for all cells
            !>!
            MEASURE_MAX = 0.0_EB
            BDRY_LOOP1: DO IW = 1, L%NW
               IF (L%WALL(IW)%NOM == 0) CYCLE BDRY_LOOP1
               IC = L%WALL(IW)%ICW
               MEASURE_MAX = MAX(MEASURE_MAX, AMG%MEASURE(IC))
            ENDDO BDRY_LOOP1
            IF (MEASURE_MAX <= EPS) EXIT BDRY_CYCLE_LOOP

            BDRY_LOOP2: DO IW = 1, L%NW

               IF (L%WALL(IW)%NOM == 0) CYCLE BDRY_LOOP2
               IC = L%WALL(IW)%ICW

               !> Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, AMG%MEASURE(IC))) THEN

                  AMG%MEASURE(IC)  = 0.0_EB
                  AMG%CELL_TYPE(IC) = NSCARC_CELL_TYPE_COARSE

                  !> Determine set of fine cells
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1

                   !> JC is set to be a fine cell which is no longer measured
                     JC = A%COL(ICOL)

                     AMG%MEASURE(JC)  = 0.0_EB
                     AMG%CELL_TYPE(JC) = NSCARC_CELL_TYPE_SFINE

                     !>  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                     DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
                        KC = A%COL(JCOL)
                        IF (KC > 0 .AND. KC /= IC .AND. (AMG%CELL_TYPE(KC)==NSCARC_UNDEFINED)) THEN
                           AMG%MEASURE(KC) = AMG%MEASURE(KC) + 1.0_EB
                           MEASURE_MAX = MAX(MEASURE_MAX, AMG%MEASURE(KC))
                        ENDIF
                     ENDDO

                  ENDDO

                  EXIT BDRY_LOOP2
               ENDIF

            ENDDO BDRY_LOOP2
         ENDDO BDRY_CYCLE_LOOP
      ENDDO BDRY_LOOP



END SELECT

!> -----------------------------------------------------------------------------------------------------
!> Exchange CELL_TYPEs along internal boundaries
!> -----------------------------------------------------------------------------------------------------
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_TYPE, NL)

END SUBROUTINE SCARC_SETUP_CELL_TYPES


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
INTEGER  :: IW, IW0, IW2, IDIAG, JDIAG, JCO, IC2, ICP, ICP2, ICOL0, ICOL_FIRST, ICOL_LAST
REAL(EB) :: SUM_COUPLED, SUM_CPOINTS, SCAL, SUM_COARSE, SUM_DIAG
REAL(EB) :: DATA_SUM, DATA_DIAG, DATA_INTERPOL
REAL(EB) :: VALUES(20), WEIGHTS(20)
INTEGER  :: NEIGHBOR(20), NWEIGHTS, NWEIGHTS2
INTEGER  :: COARSE_CELL(20), COARSE_INDEX(20)
LOGICAL  :: BFIRST
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, P, R
TYPE (SCARC_AMG_TYPE), POINTER :: AMG

SCARC_ROUTINE = 'SCARC_SETUP_PROLONGATION'

SELECT_INTERPOLATION: SELECT CASE (TYPE_INTERPOL)

   !> ------------------------------------------------------------------------------------------
   !> Standard interpolation
   !> ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_STANDARD)

      MESHES_LOOP_STANDARD: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(*,*) 'Standard interpolatiion not yet implemented'
      ENDDO MESHES_LOOP_STANDARD
      stop

   !> ------------------------------------------------------------------------------------------
   !> GMG-like interpolation
   !> ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_GMG, NSCARC_INTERPOL_GMG3)
      IF (TWO_D) THEN
         SCAL = 4.0_EB
      ELSE
         SCAL = 4.0_EB
      ENDIF
      !> temporarily
      SCAL = 2.0_EB

      MESHES_LOOP_GMG: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L   => SCARC(NM)%LEVEL(NL)
         P   => SCARC(NM)%LEVEL(NL)%P
         R   => SCARC(NM)%LEVEL(NL)%R
         AMG => SCARC(NM)%LEVEL(NL)%AMG

         ICP2 = 1
         DO IC = 1, L%NCE
            BFIRST = .TRUE.
            DO ICP = 1, AMG%NCCE

               ROW_LOOP: DO IC2 = R%ROW(ICP),R%ROW(ICP+1)-1
                  IF (R%COL(IC2) == IC) THEN
                     P%VAL(ICP2) = SCAL*R%VAL(IC2)
                     IF (BFIRST) THEN
                        P%ROW(IC) = ICP2
                     ENDIF
                     P%COL(ICP2) = ICP
                     ICP2 = ICP2 + 1
                     BFIRST = .FALSE.
                     EXIT ROW_LOOP
                  ENDIF
               ENDDO ROW_LOOP

            ENDDO

         ENDDO
         P%ROW(L%NC+1)=ICP2

      ENDDO MESHES_LOOP_GMG

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)

   !> ------------------------------------------------------------------------------------------
   !> Classical interpolation
   !> ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_CLASSICAL)

      MESHES_LOOP_CLASSICAL: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%LEVEL(NL)%A
         P => SCARC(NM)%LEVEL(NL)%P
         R => SCARC(NM)%LEVEL(NL)%R

         IP = 1
         INTERNAL_CELL_LOOP_CLASSICAL: DO IC = 1, L%NC

            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0

            !>!
            !> If IC is a coarse cell, its value is taken
            !>!
            IF (AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               NEIGHBOR(1)= AMG%CELL_TYPE(IC)
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
                  IF (AMG%CELL_TYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
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

                  SELECT CASE (AMG%CELL_TYPE(JC))

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
                  NEIGHBOR(IW) = AMG%CELL_TYPE(NEIGHBOR(IW))
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

         !> Build Interpolation
         ICOL0 = 1
         INTERNAL_CELL_LOOP2_CLASSICAL2: DO IC = 1, L%NC

            !> For a C-point IC, use identity and set mapping
            IF (AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               P%ROW(IC) = ICOL0
               P%COL(ICOL0) = AMG%CELL_TYPE(IC)
               P%VAL(ICOL0) = 1.0_EB

               ICOL0 = ICOL0 + 1

            !> For a F-point IC, build interpolation
            ELSE

               P%ROW(IC)= ICOL0               !> diagonal part of P
               ICOL_FIRST = ICOL0

               DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                  JC = A%COL(ICOL)

                  !> If JC is a C-point, initialize interpolation to zero and set column number in P%COL
                  IF (AMG%CELL_TYPE(JC) >= NSCARC_CELL_TYPE_COARSE) THEN
                     AMG%CELL_PTR(JC) = ICOL0
                     P%COL(ICOL0) = AMG%CELL_TYPE(JC)
                     P%VAL(ICOL0) = 0.0_EB
                     ICOL0 = ICOL0 + 1
                  !> If JC is a F-point, set it to be a strong F-point with relevant connections
                  ELSE IF (AMG%CELL_TYPE(JC) /= NSCARC_CELL_TYPE_WFINE) THEN
                     AMG%CELL_PTR(JC) = NSCARC_CELL_TYPE_SFINE
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
                  IF (AMG%CELL_PTR(JC) >= ICOL_FIRST) THEN
                     JCOL = AMG%CELL_PTR(JC)
                     P%VAL(JCOL) = P%VAL(JCOL) + A%VAL(ICOL)

                  !> If JC is a strongly coupled F-point to IC, the value a_(IC, JC) must be distributed
                  !> to C-points which are strongly coupled to IC (no distribution to the diagonal part)
                  ELSE IF (AMG%CELL_PTR(JC) >= NSCARC_CELL_TYPE_SFINE) THEN
                     DATA_SUM  = 0.0_EB
                     SIGN0 = 1
                     IF (A%VAL(IDIAG) > 0) SIGN0=-1

                     !> search JC in the row of A and get sum of the couplings to C-points strongly coupled to IC
                     DO JCOL = A%ROW(JC), A%ROW(JC+1)-1
                        KC = A%COL(JCOL)
                        IF (AMG%CELL_PTR(KC) > ICOL_FIRST .AND. SIGN0*A%VAL(JCOL) >0) THEN
                           DATA_SUM = DATA_SUM + A%VAL(JCOL)
                        ENDIF
                     ENDDO

                     IF (DATA_SUM .NE. 0) THEN

                        DATA_INTERPOL = A%VAL(ICOL)/DATA_SUM

                        !> loop over row of A for JC and spread data
                        DO JCOL = A%ROW(JC), A%ROW(JC+1)-1
                           KC = A%COL(JCOL)
                           IF (AMG%CELL_PTR(KC) >= ICOL_FIRST .AND. SIGN0*A%VAL(JCOL) > 0) THEN
                              KCOL = AMG%CELL_PTR(KC)
                              P%VAL(KCOL) = P%VAL(KCOL) + DATA_INTERPOL * A%VAL(JCOL)
                           ENDIF
                        ENDDO


                     ENDIF

                  !> If JC is a weakly coupled F-point to IC, the value a_(IC, JC) must be
                  !> distributed to the diagonal
                  ELSE IF (AMG%CELL_TYPE(JC) .NE. NSCARC_CELL_TYPE_WFINE) THEN
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

         IP = 1
         INTERNAL_CELL_LOOP: DO IC = 1, L%NC

            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0

            !>!
            !> If IC is a coarse cell, its value is taken
            !>!
            IF (AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

               NEIGHBOR(1)= AMG%CELL_TYPE(IC)
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
               SELECT_FPOINT_TYPE: SELECT CASE(AMG%CELL_TYPE(IC))

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
                        IF (AMG%CELL_TYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + A%VAL(ICOL)
                     ENDDO

                     SCAL = - SUM_COUPLED/SUM_CPOINTS

                     !> for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                     IW = 1
                     DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                        JC = A%COL(ICOL)
                        IF (AMG%CELL_TYPE(JC) > 0 ) THEN
                           NEIGHBOR(IW) = AMG%CELL_TYPE(JC)
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
                           IF (AMG%CELL_TYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + A%VAL(JCOL)
                        ENDDO

                        IF (SUM_CPOINTS == 0.0_EB) THEN
                           WRITE(SCARC_MESSAGE,'(2A,I8,A,I8,A,I8,A,I8)') TRIM(SCARC_ROUTINE),&
                              ': Error in SCARC_INTERPOLATION on mesh ',NM,' for cell ',IC,' with neighbor ', JC, &
                              '  on level ', NL
                           CALL SHUTDOWN(SCARC_MESSAGE); RETURN
                        ENDIF
                        SCAL =  A%VAL(ICOL)/A%VAL(IDIAG) * SUM_COUPLED/SUM_CPOINTS

                        !> Get diagonal matrix a_jj for point JC
                        JDIAG = A%ROW(JC)

                        !> Compute interpolation weights for all strong coarse cell couplings KC of JC
                        !> note that a coarse cell KC may be considered several times for different JC's
                        DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
                           KC = A%COL(JCOL)
                           IF (AMG%CELL_TYPE(KC) > 0) THEN
                             NEIGHBOR(IW) = AMG%CELL_TYPE(KC)
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
                 P%COL(ICOL) = AMG%CELL_TYPE(JC)
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
                    JCO = L%DISCRET(L%WALL(IW)%IXG,L%WALL(IW)%IYG,L%WALL(IW)%IZG,UNKH)  ! ACHTUNG: SCHAUEN!
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

         IP = 1
         INTERNAL_CELL_BDRY_LOOP: DO IC = 1, L%NC

            DIRECT_BDRY_IF: IF (L%INTERNAL_BDRY_CELL(IC)/=0) THEN

               VALUES   = 0.0_EB
               WEIGHTS  = 0.0_EB
               NEIGHBOR = 0

               !>!
               !> If IC is a coarse cell, its value is taken
               !>!
               IF (AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

                  NEIGHBOR(1)= AMG%CELL_TYPE(IC)
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
                     IF (L%INTERNAL_BDRY_CELL(JC)/=0.AND.AMG%CELL_TYPE(JC)>0) SUM_CPOINTS = SUM_CPOINTS + A%VAL(ICOL)
                  ENDDO

                  SCAL = - SUM_COUPLED/SUM_CPOINTS

                  !> for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                  IW = 1
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     JC = A%COL(ICOL)
                     IF (JC > L%NC) CYCLE
                     IF (L%INTERNAL_BDRY_CELL(JC)/=0.AND.AMG%CELL_TYPE(JC) > 0 ) THEN
                        NEIGHBOR(IW) = AMG%CELL_TYPE(JC)
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
               IF (AMG%CELL_TYPE(IC) >= NSCARC_CELL_TYPE_COARSE) THEN

                  NEIGHBOR(1)= AMG%CELL_TYPE(IC)
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
                  SELECT_FPOINT_BDRY_TYPE: SELECT CASE(AMG%CELL_TYPE(IC))

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
                           IF (AMG%CELL_TYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + A%VAL(ICOL)
                        ENDDO

                        SCAL = - SUM_COUPLED/SUM_CPOINTS

                        !> for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                        IW = 1
                        DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                           JC = A%COL(ICOL)
                           IF (AMG%CELL_TYPE(JC) > 0 ) THEN
                              NEIGHBOR(IW) = AMG%CELL_TYPE(JC)
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
                              IF (AMG%CELL_TYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + A%VAL(JCOL)
                           ENDDO

                           IF (SUM_CPOINTS == 0.0_EB) THEN
                              WRITE(SCARC_MESSAGE,'(2A,I8,A,I8,A,I8,A,I8)') TRIM(SCARC_ROUTINE),&
                                 ': Error in SCARC_INTERPOLATION on mesh ',NM,' for cell ',IC,' with neighbor ', JC, &
                                 '  on level ', NL
                              CALL SHUTDOWN(SCARC_MESSAGE); RETURN
                           ENDIF
                           SCAL =  A%VAL(ICOL)/A%VAL(IDIAG) * SUM_COUPLED/SUM_CPOINTS

                           !> Get diagonal matrix a_jj for point JC
                           JDIAG = A%ROW(JC)

                           !> Compute interpolation weights for all strong coarse cell couplings KC of JC
                           !> note that a coarse cell KC may be considered several times for different JC's
                           DO JCOL = A%ROW(JC)+1, A%ROW(JC+1)-1
                              KC = A%COL(JCOL)
                              IF (AMG%CELL_TYPE(KC) > 0) THEN
                                NEIGHBOR(IW) = AMG%CELL_TYPE(KC)
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
                 P%COL(ICOL) = AMG%CELL_TYPE(JC)
              ENDIF
            ENDDO
         ENDDO INTERNAL_CELL_BDRY_LOOP2


      ENDDO MESHES_LOOP_DIRECT_BDRY

      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_PROLONGATION, NL)


   !> ---------------------------------------------------------------------------------------------
   !> Multipass interpolation
   !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_MULTIPASS)

      MESHES_LOOP_MULTIPASS: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(*,*) 'Multipass interpolation not yet implemented'
      ENDDO MESHES_LOOP_MULTIPASS
      stop

END SELECT SELECT_INTERPOLATION


END SUBROUTINE SCARC_SETUP_PROLONGATION

!> ------------------------------------------------------------------------------------------------
!> Numbering of single patches
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PATCH_WEIGHTS(R, ROW, COL, CELL_TYPE, WEIGHTS, NX, NX1, NX2, NY1, NY2, IZ, &
                               INCRX, INCRY, INCRZ, IP, ICC)
INTEGER, INTENT(IN)    :: NX, NX1, NX2, NY1, NY2, IZ, INCRX, INCRY, INCRZ
INTEGER, INTENT(INOUT) :: IP, ICC
INTEGER,  DIMENSION(:)  , INTENT(OUT) :: ROW, COL, CELL_TYPE
REAL(EB), DIMENSION(:,:), INTENT(IN)  :: WEIGHTS
REAL(EB), DIMENSION(:)  , INTENT(OUT) :: R
INTEGER :: IX, INX, INZ, PX(3), PY(3), PZ(3), PC(3,3), ISAVE

ISAVE = NY1       ! Dummy call, only temporarily (only 2D version currently)
ISAVE = NY2       ! Dummy call, only temporarily (only 2D version currently)
ISAVE = INCRY     ! Dummy call, only temporarily (only 2D version currently)

PY=1
DO INZ = 1, INCRZ
   PZ(INZ) = IZ + INZ - 1
ENDDO

DO IX = NX1, NX2, INCRX

   DO INX = 1, INCRX
      PX(INX) = IX + INX - 1
   ENDDO

   ROW(ICC) = IP

   DO INZ = 1, INCRZ
      DO INX = 1, INCRX
         PC(INX,INZ) = (PZ(INZ)-1)*NX + PX(INX)
         IF (CELL_TYPE(PC(INX,INZ))== NSCARC_CELL_TYPE_COARSE) CELL_TYPE(PC(INX,INZ))=ICC
         COL(IP)   = PC(INX,INZ)
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

SCARC_ROUTINE = 'SCARC_SETUP_MATRIX_TRANSPOSE'

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
TYPE (SCARC_LEVEL_TYPE), POINTER :: F, C
TYPE (SCARC_MATRIX_TYPE), POINTER :: A

SCARC_ROUTINE = 'SCARC_SETUP_SYSTEM_AMG'

MESHES_LOOP_SYSTEM_AMG: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   F => SCARC(NM)%LEVEL(NL)
   C => SCARC(NM)%LEVEL(NL+1)

   IROW = 1

!> 
!> Still experimental code:
!> First look at all exterior entries
!> 
!>  IF (C%NW > 0) THEN
!>
!>     LOOP_NBR_CELLS: DO IC0 = 1, L%NCC
!>
!>        IF (C%INTERNAL_BDRY(IC) < 0) CYCLE LOOP_NBR_CELLS
!>
!>        IROW_INIT = IROW
!>
!>      !> loop over all entries in row IC0 of R
!>        LOOP_NBR_R_ENTRIES: DO ICOL1 = F%R%ROW(IC0), F%R%ROW(IC0+1)-1
!>
!>           IC1 = F%R%COL(ICOL1)
!>           IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '--------- ICOL1 =', ICOL1,'   : IC1=',IC1
!>
!>         !> loop over all entries in row IC1 of A
!>           LOOP_NBR_A_ENTRIES: DO ICOL2 = A%ROW(IC1), A%ROW(IC1+1)-1
!>
!>              IC2 = A%COL(ICOL2)
!>              IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '   ------ ICOL2 =', ICOL2,'   : IC2=',IC2
!>
!>              IF (F%A%TAG(IC2) /= IC0) THEN
!>
!>               !> mark IC2 as already considered
!>                 F%A%TAG(IC2) = IC0
!>                 IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '          A%TAG(',IC2,') =', F%A%TAG(IC2)
!>
!>               !> loop over all entries in row IC2 of P
!>                 LOOP_NBR_P_ENTRIES: DO ICOL3 = F%P%ROW(IC2), F%P%ROW(IC2+1)-1
!>
!>                    IC3 = F%P%COL(ICOL3)
!>                    IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '    !> ---- ICOL3 =', ICOL3,'   : IC3=',IC3
!>                    !IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
!>                  !>   '          P_TAG(',IC3,') =', F%P%TAG(IC3),': IROW_INIT=',IROW_INIT,' : IROW=',IROW
!>
!>                  !> verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
!>                    IF (F%P%TAG(IC3) < IROW_INIT) THEN
!>                       F%P%TAG(IC3) = IROW
!>                       IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,*) '          P_TAG(',IC3,') =', F%P%TAG(IC3)
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


   !> Determine size of matrix RAP
   !> loop over interior c-cells
   LOOP1_OWN_CELLS: DO IC0 = 1, F%AMG%NCCE

      IROW_INIT = IROW

      IF (BSQUARE) THEN
         F%P%TAG(IC0) = IROW
         IROW = IROW + 1
      ENDIF

      !> loop over all entries in row IC0 of R
      LOOP1_R_ENTRIES: DO ICOL1 = F%R%ROW(IC0), F%R%ROW(IC0+1)-1
         IC1 = F%R%COL(ICOL1)

         !> loop over all entries in row IC1 of A
         LOOP1_A_ENTRIES: DO ICOL2 = F%A%ROW(IC1), F%A%ROW(IC1+1)-1
            IC2 = A%COL(ICOL2)

            IF (F%A%TAG(IC2) /= IC0) THEN
               F%A%TAG(IC2) = IC0                           !> mark IC2 as already considered

               !> loop over all entries in row IC2 of P
               LOOP1_P_ENTRIES: DO ICOL3 = F%P%ROW(IC2), F%P%ROW(IC2+1)-1

                  IC3 = F%P%COL(ICOL3)

                  !> verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
                  IF (IC3 >0) THEN   !> ONLY TEMPORARILY
                  IF (F%P%TAG(IC3) < IROW_INIT) THEN
                     F%P%TAG(IC3) = IROW
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

   !> Determine size of matrix RAP  == C%A
   C%A%NAV = IROW
   C%A%NAE = C%A%NAV + 100  !> ONLY TEMPORARILY

   CALL SCARC_ALLOCATE_MATRIX(C%A, C%A%NAE+1, C%A%NAE+1, C%NCE+1, -1, 'C%A')
   C%A%ROW(F%AMG%NCC+1) = IROW

   F%A%TAG = 0
   F%P%TAG = 0

   IROW = 1

   !> loop over interior c-cells
   LOOP2_C_CELLS: DO IC0 = 1, F%AMG%NCC
      IROW_INIT = IROW
      C%A%ROW(IC0) = IROW_INIT

      IF (BSQUARE) THEN
         F%P%TAG(IC0) = IROW
         C%A%COL(IROW) = IC0
         C%A%VAL(IROW) = 0.0_EB
         IROW = IROW + 1
      ENDIF

      !> loop over all entries in row IC0 of R
      LOOP2_R_ENTRIES: DO ICOL1 = F%R%ROW(IC0), F%R%ROW(IC0+1)-1

         IC1 = F%R%COL(ICOL1)
         R_VALUE = F%R%VAL(ICOL1)

         !> loop over all entries in row IC1 of A
         LOOP2_A_ENTRIES: DO ICOL2 = A%ROW(IC1), A%ROW(IC1+1)-1

            IC2 = A%COL(ICOL2)
            RA_VALUE = R_VALUE * A%VAL(ICOL2)

            !> Hasn't cell IC2 been considered before? (new values for RAP can only be done for unmarked cells)
            IF (F%A%TAG(IC2) /= IC0) THEN

               !> mark IC2 as already considered
               F%A%TAG(IC2) = IC0

               !> loop over all entries in row IC2 of P
               LOOP2_P_ENTRIES1: DO ICOL3 = F%P%ROW(IC2), F%P%ROW(IC2+1)-1

                  IC3 = F%P%COL(ICOL3)
                  RAP_VALUE = RA_VALUE * F%P%VAL(ICOL3)

                  !> verify that P_TAG for entry RAP_(IC0, IC3) hasn't been considered yet; if not, mark it
                  IF (IC3 >0) THEN     !> ONLY TEMPORARILY

                  IF (F%P%TAG(IC3) < IROW_INIT) THEN
                     F%P%TAG(IC3)  = IROW
                     C%A%COL(IROW) = F%P%COL(ICOL3)
                     C%A%VAL(IROW)     = RAP_VALUE
                     IROW = IROW + 1
                  ELSE
                     C%A%VAL(F%P%TAG(IC3)) = C%A%VAL(F%P%TAG(IC3)) + RAP_VALUE
                  ENDIF

                  ENDIF  !> ONLY TEMPORARILY
               ENDDO LOOP2_P_ENTRIES1

            !> or has it been already considered
            ELSE

               !> loop over all entries in row IC2 of P
               LOOP2_P_ENTRIES2: DO ICOL3 = F%P%ROW(IC2), F%P%ROW(IC2+1)-1
                  IC3 = F%P%COL(ICOL3)
                  RAP_VALUE = RA_VALUE * F%P%VAL(ICOL3)
                  C%A%VAL(F%P%TAG(IC3)) = C%A%VAL(F%P%TAG(IC3)) + RAP_VALUE
               ENDDO LOOP2_P_ENTRIES2
            ENDIF

         ENDDO LOOP2_A_ENTRIES
      ENDDO LOOP2_R_ENTRIES

      !> Store counters
      IROW_SAVE = IROW

   ENDDO LOOP2_C_CELLS
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
LOGICAL  :: BFIRST, BEVENX, BEVENZ, BTWOLEVEL_X, BTWOLEVEL_Z, BTHREE_X, BTHREE_Z
REAL(EB) :: WEIGHTS2X2(2,2), WEIGHTS2X3(2,3), WEIGHTS3x2(3,2), WEIGHTS3X3(3,3)
TYPE (OSCARC_TYPE), POINTER :: OS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_LEVEL_TYPE), POINTER :: OL
TYPE (SCARC_MATRIX_TYPE), POINTER :: P, R, OP, OR
TYPE (SCARC_AMG_TYPE), POINTER :: AMG

SCARC_ROUTINE = 'SCARC_SETUP_RESTRICTION'

SELECT CASE (TYPE_INTERPOL)

 !> ------------------------------------------------------------------------------------------
 !> GMG-like interpolation
 !> ------------------------------------------------------------------------------------------
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

      MESHES_LOOP_GMG: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         R => SCARC(NM)%LEVEL(NL)%R
         AMG => SCARC(NM)%LEVEL(NL)%AMG
         IP = 1
         ICC = 1

         !> Analyze grid sizes
         INCRX = 2
         IF (MOD(L%NX,2) == 0) THEN
            BEVENX=.TRUE.
         ELSE
            BEVENX=.FALSE.
         ENDIF

         INCRZ  = 2
         IF (MOD(L%NZ,2) == 0) THEN
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
               DO IZ = 1, L%NZ ,2
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS2x2, L%NX,&
                                           1      , L%NX  , 1, 1, IZ , 2, 1, 2, IP, ICC)
               ENDDO
            CASE (1)
               DO IZ = 1, L%NZ-3,2
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS2x2, L%NX,&
                                           1      , L%NX  , 1, 1, IZ , 2, 1, 2, IP, ICC)
               ENDDO
               CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS2x3, L%NX,&
                                           1      , L%NX  , 1, 1, L%NZ-2, 2, 1, 3, IP, ICC)
            CASE (2)
               DO IZ = 1, L%NZ,2
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS2x2, L%NX,&
                                           1      , L%NX-3, 1, 1, IZ, 2, 1, 2, IP, ICC)
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS3x2, L%NX,&
                                           L%NX-2, L%NX-2, 1, 1, IZ, 3, 1, 2, IP, ICC)
               ENDDO
            CASE (3)
               DO IZ = 1, L%NZ-3,2
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS2x2, L%NX,&
                                           1      , L%NX-3, 1, 1, IZ, 2, 1, 2, IP, ICC)
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS3x2, L%NX,&
                                           L%NX-2, L%NX-2, 1, 1, IZ, 3, 1, 2, IP, ICC)
               ENDDO
               CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS2x3, L%NX,&
                                        1      , L%NX-3, 1, 1, L%NZ-2, 2, 1, 3, IP, ICC)
               CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS3x3, L%NX,&
                                        L%NX-2, L%NX-2, 1, 1, L%NZ-2, 3, 1, 3, IP, ICC)
         END SELECT
         R%ROW(AMG%NCC+1) = IP

      ENDDO MESHES_LOOP_GMG

 !> ------------------------------------------------------------------------------------------
 !> GMG-like interpolation
 !> ------------------------------------------------------------------------------------------
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

      MESHES_LOOP_GMG3: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         AMG => SCARC(NM)%LEVEL(NL)%AMG
         IP = 1
         ICC = 1

         !> Analyze grid sizes
         IF (MOD(L%NX,3) == 0) THEN
            INCRX   = 3
            BTWOLEVEL_X   = .FALSE.
            BTHREE_X = .TRUE.
         ELSE IF (MOD(L%NX,2) == 0) THEN
            INCRX    = 2
            BTWOLEVEL_X   = .TRUE.
            BTHREE_X = .FALSE.
         ELSE
            INCRX  = 2
            BTWOLEVEL_X   = .FALSE.
            BTHREE_X = .FALSE.
         ENDIF

         IF (MOD(L%NZ,3) == 0) THEN
            INCRZ   = 3
            BTWOLEVEL_Z   = .FALSE.
            BTHREE_Z = .TRUE.
         ELSE IF (MOD(L%NZ,2) == 0) THEN
            INCRZ    = 2
            BTWOLEVEL_Z   = .TRUE.
            BTHREE_Z = .FALSE.
         ELSE
            INCRZ  = 2
            BTWOLEVEL_Z   = .FALSE.
            BTHREE_Z = .FALSE.
         ENDIF

         IF (.NOT.(BTWOLEVEL_X.OR.BTHREE_X).OR..NOT.(BTWOLEVEL_Z.OR.BTHREE_Z)) THEN
            WRITE(SCARC_MESSAGE,'(2A)') TRIM(SCARC_ROUTINE),': Cell numbers not divisable by 2 or 3!'
            CALL SHUTDOWN(SCARC_MESSAGE); RETURN
         ENDIF

         IF (BTHREE_X) THEN
            ICASE = 0
         ELSE IF (BTWOLEVEL_X) THEN
            ICASE = 1
         ELSE
            ICASE = 2
         ENDIF

         SELECT CASE(ICASE)
            CASE (0)
               DO IZ = 1, L%NZ ,INCRZ
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS3x3, L%NX,&
                                           1      , L%NX  , 1, 1, IZ , 3, 1, INCRZ, IP, ICC)
               ENDDO
            CASE (1)
               DO IZ = 1, L%NZ,INCRZ
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS2x3, L%NX,&
                                           1      , L%NX  , 1, 1, IZ , 2, 1, INCRZ, IP, ICC)
               ENDDO
            CASE (2)
               DO IZ = 1, L%NZ,INCRZ
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS2x2, L%NX,&
                                           1      , L%NX-3, 1, 1, IZ, 3, 1, INCRZ, IP, ICC)
                  CALL SCARC_PATCH_WEIGHTS(R%VAL, R%ROW, R%COL, AMG%CELL_TYPE, WEIGHTS3x2, L%NX,&
                                           L%NX-2, L%NX-2, 1, 1, IZ, 2, 1, INCRZ, IP, ICC)
               ENDDO
         END SELECT
         R%ROW(AMG%NCC+1) = IP

      ENDDO MESHES_LOOP_GMG3

 !> ------------------------------------------------------------------------------------------
 !> DEFAULT
 !> ------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_CLASSICAL2)

      MESHES_LOOP_CLASSICAL2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         P => SCARC(NM)%LEVEL(NL)%P
         R => SCARC(NM)%LEVEL(NL)%R
         AMG => SCARC(NM)%LEVEL(NL)%AMG
         CALL SCARC_MATRIX_TRANSPOSE(P%VAL, P%ROW, P%COL, R%VAL, R%ROW, R%COL, L%NCE, AMG%NCCE)
      ENDDO MESHES_LOOP_CLASSICAL2


 !> ------------------------------------------------------------------------------------------
 !> DEFAULT
 !> ------------------------------------------------------------------------------------------
   CASE DEFAULT

      MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         R => SCARC(NM)%LEVEL(NL)%R
         P => SCARC(NM)%LEVEL(NL)%P

         IC2 = 1
         DO ICP = 1, AMG%NCCE
            BFIRST = .TRUE.
            DO IC = 1, L%NCE

               ROW_LOOP: DO ICP2 = P%ROW(IC),P%ROW(IC+1)-1
                  IF (P%COL(ICP2) == ICP) THEN
                     R%VAL(IC2) = P%VAL(ICP2)
                     IF (BFIRST) THEN
                        R%ROW(ICP) = IC2
                     ENDIF
                     R%COL(IC2) = IC
                     IC2 = IC2 + 1
                     BFIRST = .FALSE.
                     EXIT ROW_LOOP
                  ENDIF
               ENDDO ROW_LOOP

            ENDDO

         ENDDO
         R%ROW(AMG%NCCE+1)=IC2

         OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

            IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
            OS => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE OTHER_MESHES_LOOP
            OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            OP => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%P
            OR => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%R

            IC2 = 1
            ICC = 0
            OTHER_GHOSTCELL_LOOP: DO IG0 = 1, OL%NCG

               ICP = OL%AMG%CELL_TYPE(IG0)
               IF (ICP < NSCARC_CELL_TYPE_COARSE) CYCLE OTHER_GHOSTCELL_LOOP
               ICC = ICC + 1

               BFIRST = .TRUE.
               DO IG = 1, OL%NCG

                  OTHER_ROW_LOOP: DO ICP2 = OP%ROW(IG),OP%ROW(IG+1)-1
                     IF (OP%COL(ICP2) == ICC) THEN
                        OR%VAL(IC2) = OP%VAL(ICP2)
                        IF (BFIRST) THEN
                           OR%ROW(ICC) = IC2
                        ENDIF
                        OR%COL(IC2) = IG
                        IC2 = IC2 + 1
                        BFIRST = .FALSE.
                        EXIT OTHER_ROW_LOOP
                     ENDIF
                  ENDDO OTHER_ROW_LOOP

               ENDDO

            ENDDO OTHER_GHOSTCELL_LOOP
            OR%ROW(ICC+1)=IC2

         ENDDO OTHER_MESHES_LOOP

      ENDDO MESHES_LOOP

END SELECT

END SUBROUTINE SCARC_SETUP_RESTRICTION


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
INTEGER :: NM, IC, ICO2, ICP1, ICP2, IP, ICOL, IDIAG, IOR0, IW
REAL (EB), ALLOCATABLE, DIMENSION(:) :: VAL
!REAL(EB):: MAUX1(50,50)=0.0_EB!, MAUX2(30,30)=0.0_EB
REAL(EB) :: AUX1, AUX2, PW !, MATRIX(16,16)
TYPE (SCARC_LEVEL_TYPE), POINTER :: F, C

WRITE(*,*) 'ACHTUNG, MUSS ALTERNATIV ZUR FOLGEROUTINE GENUTZT WERDEN'
SCARC_ROUTINE = 'SCARC_TRANSFER_MATRIX'

SELECT CASE (TYPE_COARSENING)

 !> ---------------------------------------------------------------------------------
 !> GMG-like coarsening
 !> ---------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_GMG, NSCARC_INTERPOL_GMG3)

      MESHES_GMG: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         F => SCARC(NM)%LEVEL(NL)             !> pointer to fine level
         C => SCARC(NM)%LEVEL(NL+1)           !> pointer to coarse level

         IF (NMESHES==1) THEN
            ALLOCATE (VAL(F%NC), STAT=IERROR)
         ELSE
            ALLOCATE (VAL(F%NCE), STAT=IERROR)
         ENDIF
         CALL CHKMEMERR ('SCARC_TRANSFER_MATRIX', 'VAL', IERROR)
         VAL = 0.0_EB

         IP = 1
         ICP1_GMG_LOOP: DO ICP1 = 1, F%AMG%NCC

            ICP2_GMG_LOOP1: DO ICP2 = 1, F%AMG%NCC

               AUX2 = 0.0_EB

               DO IC = 1, F%NC

                  IDIAG = F%A%ROW(IC)
                  !PW    = 0.5_EB* P_WEIGHT(IC, ICP2, NM, NL)
                  PW    = P_WEIGHT(IC, ICP2, NM, NL)
                  AUX1  = F%A%VAL(IDIAG) * PW

                  COUPLINGS_GMG_LOOP: DO ICOL = F%A%ROW(IC)+1, F%A%ROW(IC+1)-1
                     ICO2  = F%A%COL(ICOL)
                     !PW = 0.5_EB* P_WEIGHT(ICO2, ICP2, NM, NL)
                     PW =  P_WEIGHT(ICO2, ICP2, NM, NL)
                     AUX1 = AUX1 + F%A%VAL(ICOL) * PW
                  ENDDO COUPLINGS_GMG_LOOP
!MATRIX(IC, ICP2 ) = AUX1

                  AUX2 = AUX2 + R_WEIGHT(IC, ICP1, NM, NL) * AUX1

               ENDDO
               VAL(ICP2) = AUX2

            ENDDO ICP2_GMG_LOOP1

            !> analyze new matrix line and store it corresponding to compact storage technique:
            !> (diagonal entry first)
            C%A%VAL(IP)   = VAL(ICP1)
            C%A%ROW(ICP1) = IP
            C%A%COL(IP)   = ICP1

            IP  = IP + 1
            ICP2_GMG_LOOP3: DO ICP2 = 1, F%AMG%NCC
               IF (ICP2 /= ICP1 .AND. ABS(VAL(ICP2)) >= 1.0E-12_EB) THEN
                  C%A%VAL(IP) = VAL(ICP2)
                  C%A%COL(IP) = ICP2
                  IP  = IP + 1
               ENDIF
            ENDDO ICP2_GMG_LOOP3
            DO IOR0 = -3, 3
               IW = C%WALL_INDEX(ICP1, IOR0)
               IF (IW > 0) THEN
                  C%A%VAL(IP) = 0.0_EB
                  C%A%COL(IP) = -IW
                  IP  = IP + 1
               ENDIF
            ENDDO

            VAL = 0.0_EB
         ENDDO ICP1_GMG_LOOP

         C%A%ROW(C%NC+1) = IP
         C%A%NAV = IP - 1

         DEALLOCATE(VAL)

      ENDDO MESHES_GMG

 !> ---------------------------------------------------------------------------------
 !> Default case
 !> ---------------------------------------------------------------------------------
   CASE DEFAULT

      MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         F => SCARC(NM)%LEVEL(NL)             !> pointer to fine level
         C => SCARC(NM)%LEVEL(NL+1)           !> pointer to coarse level

         IF (NMESHES==1) THEN
            ALLOCATE (VAL(F%NC), STAT=IERROR)
         ELSE
            ALLOCATE (VAL(F%NCE), STAT=IERROR)
         ENDIF
         CALL CHKMEMERR ('SCARC_TRANSFER_MATRIX', 'VAL', IERROR)
         VAL = 0.0_EB

         IP = 1
         ICP1_LOOP: DO ICP1 = 1, F%AMG%NCC


            !> ---------------------------------------------------------------------------------
            !> First: loop over the internal coarse cells
            !> ---------------------------------------------------------------------------------
            ICP2_LOOP1: DO ICP2 = 1, F%AMG%NCCE

               AUX2 = 0.0_EB

               DO IC = 1, F%NC

                  IDIAG = F%A%ROW(IC)
                  PW    = P_WEIGHT(IC, ICP2, NM, NL)
                  AUX1  = F%A%VAL(IDIAG) * PW

                  COUPLINGS_LOOP: DO ICOL = F%A%ROW(IC)+1, F%A%ROW(IC+1)-1
                     ICO2 = F%A%COL(ICOL)
                     PW = P_WEIGHT(ICO2, ICP2, NM, NL)
                     AUX1 = AUX1 + F%A%VAL(ICOL) * PW
                  ENDDO COUPLINGS_LOOP

                  AUX2 = AUX2 + P_WEIGHT(IC, ICP1, NM, NL) * AUX1

               ENDDO
               VAL(ICP2) = AUX2

            ENDDO ICP2_LOOP1

            !> analyze new matrix line and store it corresponding to compact storage technique:
            !> (diagonal entry first)
            C%A%VAL(IP)       = VAL(ICP1)
            C%A%ROW(ICP1) = IP
            C%A%COL(IP)   = ICP1

            IP  = IP + 1
            ICP2_LOOP2: DO ICP2 = 1, F%AMG%NCCE
            !DO ICP2 = 1, F%NCC
               IF (ICP2 /= ICP1 .AND. ABS(VAL(ICP2)) >= 1.0E-12_EB) THEN
                  C%A%VAL(IP) = VAL(ICP2)
                  C%A%COL(IP) = ICP2
                  IP  = IP + 1
               ENDIF
            ENDDO ICP2_LOOP2

            VAL = 0.0_EB
         ENDDO ICP1_LOOP

         C%A%ROW(C%NC+1) = IP
         C%A%NAV = IP - 1

         DEALLOCATE(VAL)

      ENDDO MESHES_LOOP

END SELECT

IF (NMESHES > 1) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_MATRIX_SUBDIAG, NL)

END SUBROUTINE SCARC_TRANSFER_MATRIX


!> ------------------------------------------------------------------------------------------------
!> Determine if corresponding coarse point contributes a non-zero interpolation weight or not
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION P_WEIGHT(IC, ICP, NM, NL)
INTEGER, INTENT(IN):: IC, ICP, NM, NL
INTEGER :: ICOL
REAL(EB) :: VAL
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: P

L => SCARC(NM)%LEVEL(NL)
P => SCARC(NM)%LEVEL(NL)%P

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
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: R

L => SCARC(NM)%LEVEL(NL)
R => SCARC(NM)%LEVEL(NL)%R

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

SCARC_ROUTINE = 'SCARC_SOLVER'
TNOW = CURRENT_TIME()

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !> ---------------- Krylov method (CG/BICG) --------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV: SELECT CASE (TYPE_KRYLOV)
         CASE (NSCARC_KRYLOV_CG)
            CALL SCARC_METHOD_CG  (NSCARC_SCOPE_MAIN)
         CASE (NSCARC_KRYLOV_BICG)
            CALL SCARC_METHOD_BICG(NSCARC_SCOPE_MAIN)
      END SELECT SELECT_KRYLOV

   !> ---------------- Multigrid method ---------------------------------------
   CASE (NSCARC_METHOD_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID(NSCARC_SCOPE_MAIN)

   !> ---------------- MKL method ---------------------------------------------
#ifdef WITH_MKL
   CASE (NSCARC_METHOD_MKL)

      SELECT_MKL: SELECT CASE (TYPE_MKL)
         CASE (NSCARC_MKL_GLOBAL) 
            CALL SCARC_METHOD_CLUSTER(NSCARC_SCOPE_MAIN)
         CASE (NSCARC_MKL_LOCAL) 
            CALL SCARC_METHOD_PARDISO(NSCARC_SCOPE_MAIN)
      END SELECT SELECT_MKL
#endif

END SELECT SELECT_METHOD

TSTEP(MYID+1)%SOLVER=MAX(TSTEP(MYID+1)%SOLVER,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%SOLVER =TSUM(MYID+1)%SOLVER+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SOLVER


!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for compact storage technique
!> ------------------------------------------------------------------------------------------------
SUBROUTINE POINT_TO_VECTOR(NVECTOR, NM, NL, VEC)
REAL(EB), INTENT(OUT), POINTER, DIMENSION(:) :: VEC
TYPE (SCARC_SOLVER_TYPE), POINTER :: SOLVER=>NULL()
INTEGER, INTENT(IN):: NVECTOR, NM, NL

SELECT CASE (TYPE_SCOPE)
   CASE (NSCARC_SCOPE_MAIN) 
      SOLVER => SCARC(NM)%LEVEL(NL)%MAIN
   CASE (NSCARC_SCOPE_PRECON) 
      SOLVER => SCARC(NM)%LEVEL(NL)%PRECON
END SELECT

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_X)
      VEC => SOLVER%X
   CASE (NSCARC_VECTOR_F)
      VEC => SOLVER%F
   CASE (NSCARC_VECTOR_D)
      VEC => SOLVER%D
   CASE (NSCARC_VECTOR_G)
      VEC => SOLVER%G
   CASE (NSCARC_VECTOR_W)
      VEC => SOLVER%W
   CASE (NSCARC_VECTOR_Y)
      VEC => SOLVER%Y
   CASE (NSCARC_VECTOR_Z)
      VEC => SOLVER%Z
END SELECT

END SUBROUTINE POINT_TO_VECTOR



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
!> Set pointer to chosen vector for compact storage technique
!> ------------------------------------------------------------------------------------------------
SUBROUTINE POINT_TO_SOLVER(SOLVER, NM, NL, NTYPE)
TYPE (SCARC_SOLVER_TYPE), POINTER, INTENT(OUT):: SOLVER
INTEGER, INTENT(IN):: NM, NL, NTYPE

SELECT CASE (NTYPE)
   CASE (NSCARC_SOLVER_MAIN)
      SOLVER => SCARC(NM)%LEVEL(NL)%MAIN
   CASE (NSCARC_SOLVER_PRECON)
      SOLVER => SCARC(NM)%LEVEL(NL)%PRECON
   CASE (NSCARC_SOLVER_COARSE)
      SOLVER => SCARC(NM)%LEVEL(NL)%COARSE
END SELECT

END SUBROUTINE POINT_TO_SOLVER


!> ------------------------------------------------------------------------------------------------
!> Compute global matrix-vector product (including data exchange along internal boundaries)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATVEC_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), POINTER, DIMENSION(:)     :: V1, V2
REAL(EB) :: VCO, TNOW
INTEGER , POINTER :: NC, NCE
INTEGER :: NM, IC, JC, ICOL
TYPE (SCARC_MATRIX_TYPE), POINTER:: A

SCARC_ROUTINE = 'SCARC_MATVEC_PRODUCT'
TNOW = CURRENT_TIME()

!> 
!> Exchange internal boundary values of vector1 such that the ghost values contain the corresponding
!> overlapped values of adjacent neighbor
!>
TYPE_VECTOR = NVECTOR1
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

!> 
!> Perform global matrix-vector product:
!> Note: - matrix already contains subdiagonal values from neighbor along internal boundaries
!>       - if vector1 contains neighboring values, then correct values of global matvec are achieved
!> 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL POINT_TO_VECTOR (NVECTOR1, NM, NL, V1)
   CALL POINT_TO_VECTOR (NVECTOR2, NM, NL, V2)

   NC => SCARC(NM)%LEVEL(NL)%NCS                                      !> number of cells
   NCE => SCARC(NM)%LEVEL(NL)%NCE                                     !> number of cells

   A => SCARC(NM)%LEVEL(NL)%A                                         !> system matrix

   DO IC = 1, NC

      ICOL = A%ROW(IC)                                                !> diagonal entry
      JC   = A%COL(ICOL)

      V2(IC) = A%VAL(ICOL)* V1(JC)
      VCO = V2(IC)

      DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1                            !> subdiagonal entry
         JC = A%COL(ICOL)
         V2(IC) =  V2(IC) + A%VAL(ICOL)* V1(JC)
         VCO = V2(IC)
      ENDDO
   ENDDO
ENDDO

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

SP_LOCAL = 0.0_EB

!> 
!> Compute local scalar products on single meshes und process group
!> 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

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

ENDDO

!> 
!> Compute global scalar product as sum of local scalar products
!> 
NL0  = NL
SP_GLOBAL   = 0.0_EB

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,SP_LOCAL,COUNTS,DISPLS,&
                       MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
ENDIF
SP_GLOBAL = SUM(SP_LOCAL(1:NMESHES))

SCARC_SCALAR_PRODUCT = SP_GLOBAL

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

SP_GLOBAL = SCARC_SCALAR_PRODUCT(NVECTOR1, NVECTOR1, NL)
SP_GLOBAL = SQRT (SP_GLOBAL/REAL(N_CELLS_GLOBAL(NL), EB))

SCARC_L2NORM = SP_GLOBAL

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
SUBROUTINE SCARC_COPY_INT(X, Y, SCAL, NLEN)
INTEGER,  INTENT(IN) , DIMENSION(:) :: X
INTEGER,  INTENT(OUT), DIMENSION(:) :: Y
INTEGER,  INTENT(IN) :: NLEN
REAL(EB) :: SCAL
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
REAL(EB), DIMENSION(:)    , POINTER ::  V1, V2
INTEGER  :: NM

#if defined(WITH_MKL)
INTEGER , POINTER :: NC
EXTERNAL :: DCOPY, DSCAL
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

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

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL POINT_TO_VECTOR(NVECTOR, NM, NL, VC)
   VC =  0.0_EB
ENDDO

END SUBROUTINE SCARC_VECTOR_CLEAR


!> ------------------------------------------------------------------------------------------------
!> Perform preconditioning
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_BLOCK_SOLVER (NSOL, NRHS, NSCOPE, NL)
USE POIS, ONLY: H2CZSS, H3CZSS
INTEGER, INTENT(IN):: NSOL, NRHS, NSCOPE, NL
REAL(EB), DIMENSION(:)    , POINTER ::  V1, V2
INTEGER , POINTER:: NC
INTEGER  :: NPRECON, NM, I, J, K, IC, ICOL!, IROW, JCOL
REAL(EB) :: AUX, OMEGA=1.5_EB
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE),  POINTER :: L
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, AS
TYPE (SCARC_FFT_TYPE), POINTER :: FFT
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: MKL
#endif

REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

!> Choose preconditioner depending on given scope
SELECT CASE (NSCOPE)
   CASE (NSCARC_SCOPE_MAIN)
      NPRECON = TYPE_PRECON
   CASE (NSCARC_SCOPE_PRECON)
      NPRECON = TYPE_SMOOTH
   CASE (NSCARC_SCOPE_COARSE)
     NPRECON = NSCARC_PRECON_SSOR
END SELECT

SELECT CASE (NPRECON)

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise Jacobi
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_JACOBI)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL POINT_TO_VECTOR(NRHS, NM, NL, V2)

         NC => SCARC(NM)%LEVEL(NL)%NCS
         A  => SCARC(NM)%LEVEL(NL)%A

         DO IC = 1, NC
            V2(IC) = V2(IC) / A%VAL(A%ROW(IC))
         ENDDO

      ENDDO

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise SSOR
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_SSOR)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL POINT_TO_VECTOR(NRHS, NM, NL, V2)

         NC => SCARC(NM)%LEVEL(NL)%NCS
         A  => SCARC(NM)%LEVEL(NL)%A

         !> forward SOR step
         FORWARD_CELL_LOOP: DO IC = 1, NC
            AUX = 0.0_EB
            LOWER_DIAG_LOOP: DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
               IF (A%COL(ICOL) >= IC) CYCLE LOWER_DIAG_LOOP
               IF (A%COL(ICOL) <= NC) AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))
            ENDDO LOWER_DIAG_LOOP
            V2(IC) = (V2(IC) - AUX * OMEGA) / A%VAL(A%ROW(IC))
         ENDDO FORWARD_CELL_LOOP

         !> backward SOR step
         BACKWARD_CELL_LOOP: DO IC = NC-1, 1, -1
            AUX = 0.0_EB
            UPPER_DIAG_LOOP: DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
               IF (A%COL(ICOL) <= IC) CYCLE UPPER_DIAG_LOOP
               IF (A%COL(ICOL) <= NC) AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))
            ENDDO UPPER_DIAG_LOOP
            V2(IC) = V2(IC) - AUX * OMEGA / A%VAL(A%ROW(IC))
         ENDDO BACKWARD_CELL_LOOP

      ENDDO

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise Geometric Multigrid 
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID (NSCARC_SCOPE_PRECON)


   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise FFT based on Crayfishpak
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_FFT)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL POINT_TO_VECTOR(NSOL, NM, NL, V1)
         CALL POINT_TO_VECTOR(NRHS, NM, NL, V2)

         M   => MESHES(NM)
         L   => SCARC(NM)%LEVEL(NL)
         FFT => SCARC(NM)%LEVEL(NL)%FFT

         DO K = 1, L%NZ
            DO J = 1, L%NY
               DO I = 1, L%NX
                  IC = (K-1) * L%NX * L%NY + (J-1) * L%NX + I
                  FFT%WORK(I, J, K) = V1(IC)
               ENDDO
            ENDDO
         ENDDO
         IF (TWO_D) THEN
            FFT%BXS=0.0_EB
            FFT%BXF=0.0_EB
            FFT%BZS=0.0_EB
            FFT%BZF=0.0_EB
            CALL H2CZSS (FFT%BXS, FFT%BXF, FFT%BZS, FFT%BZF,&
                         L%NX+1, FFT%WORK, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
         ELSE
            FFT%BXS=0.0_EB
            FFT%BXF=0.0_EB
            FFT%BYS=0.0_EB
            FFT%BYF=0.0_EB
            FFT%BZS=0.0_EB
            FFT%BZF=0.0_EB
            CALL H3CZSS (FFT%BXS, FFT%BXF, FFT%BYS, FFT%BYF, FFT%BZS, FFT%BZF, &
                         L%NX+1, L%NY+1, FFT%WORK, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
         ENDIF
         DO K = 1, L%NZ
            DO J = 1, L%NY
               DO I = 1, L%NX
                  IC = (K-1) * L%NX * L%NY + (J-1) * L%NX + I
                  V2(IC) = FFT%WORK(I, J, K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

#ifdef WITH_MKL
   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by Cluster Sparse Solver from MKL
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_CLUSTER)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      
         L   => SCARC(NM)%LEVEL(NL)
         MKL => SCARC(NM)%LEVEL(NL)%MKL
      
         MKL%MSGLVL =  0         ! do not print statistical information
         MKL%PHASE  = 33         ! only solving
      
         CALL POINT_TO_VECTOR (NSOL, NM, NL, V1)
         CALL POINT_TO_VECTOR (NRHS, NM, NL, V2)
         
         NC => SCARC(NM)%LEVEL(NL)%NCS
         AS => SCARC(NM)%LEVEL(NL)%AS

         CALL CLUSTER_SPARSE_SOLVER(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC_GLOBAL, &
                                    AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                    MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)

         IF (MKL%ERROR /= 0) WRITE(LU_SCARC,*) 'The following ERROR was detected: ', MKL%ERROR
      
      ENDDO

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by Pardiso Solver from MKL
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_PRECON_PARDISO)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      
         L   => SCARC(NM)%LEVEL(NL)
         MKL => SCARC(NM)%LEVEL(NL)%MKL
      
         MKL%MSGLVL =  0         ! do not print statistical information
         MKL%PHASE  = 33         ! only solving
      
         CALL POINT_TO_VECTOR (NSOL, NM, NL, V1)
         CALL POINT_TO_VECTOR (NRHS, NM, NL, V2)
         
         NC => SCARC(NM)%LEVEL(NL)%NCS
         AS => SCARC(NM)%LEVEL(NL)%AS

         CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, NC, AS%VAL, AS%ROW, AS%COL, &
                        MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1, V2, MKL%ERROR)
      
         IF (MKL%ERROR /= 0) WRITE(LU_SCARC,*) 'The following ERROR was detected: ', MKL%ERROR
      
      ENDDO

#endif

END SELECT

TSTEP(MYID+1)%PRECON=MAX(TSTEP(MYID+1)%PRECON,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%PRECON =TSUM(MYID+1)%PRECON+CURRENT_TIME()-TNOW

END SUBROUTINE SCARC_BLOCK_SOLVER

!> ------------------------------------------------------------------------------------------------
!> Save settings of calling ENV-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SAVE_PARENT(PARENT)
TYPE (SCARC_PARENT_TYPE), INTENT(OUT):: PARENT
PARENT%TYPE_METHOD   = TYPE_METHOD
PARENT%TYPE_SCOPE    = TYPE_SCOPE
PARENT%TYPE_PRECON   = TYPE_PRECON
PARENT%TYPE_SMOOTH   = TYPE_SMOOTH
PARENT%TYPE_INTERPOL = TYPE_INTERPOL
PARENT%TYPE_TWOLEVEL = TYPE_TWOLEVEL
PARENT%TYPE_CYCLE    = TYPE_CYCLE
PARENT%TYPE_ACCURACY = TYPE_ACCURACY
END SUBROUTINE SCARC_SAVE_PARENT


!> ------------------------------------------------------------------------------------------------
!> Reset settings of calling CURRENT-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_PARENT(PARENT)
TYPE (SCARC_PARENT_TYPE), INTENT(IN):: PARENT
TYPE_METHOD   = PARENT%TYPE_METHOD
TYPE_SCOPE    = PARENT%TYPE_SCOPE
TYPE_PRECON   = PARENT%TYPE_PRECON
TYPE_TWOLEVEL = PARENT%TYPE_TWOLEVEL
TYPE_INTERPOL = PARENT%TYPE_INTERPOL
TYPE_SMOOTH   = PARENT%TYPE_SMOOTH
TYPE_CYCLE    = PARENT%TYPE_CYCLE
TYPE_ACCURACY = PARENT%TYPE_ACCURACY
END SUBROUTINE SCARC_RESTORE_PARENT


!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
#ifdef WITH_MKL
SUBROUTINE SCARC_METHOD_CLUSTER(NSCOPE)
INTEGER, INTENT(IN) :: NSCOPE
INTEGER ::  NM, NL! , IP, IX, IY, IZ, IC
REAL (EB) :: TNOW
REAL(EB), POINTER, DIMENSION(:) :: XS, FS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MKL_TYPE)   , POINTER :: MKL
TYPE (SCARC_MATRIX_TYPE), POINTER :: AS
TYPE (SCARC_PARENT_TYPE) :: PARENT
TYPE (SCARC_SCOPE_TYPE) :: INTEL

SCARC_ROUTINE = 'SCARC_METHOD_CLUSTER'
TNOW = CURRENT_TIME()

CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(INTEL, NSCOPE, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L     => SCARC(NM)%LEVEL(NL)
   AS    => SCARC(NM)%LEVEL(NL)%AS                      !> point to symmetric system matrix
   MKL   => SCARC(NM)%LEVEL(NL)%MKL                     !> point to MKL structure

   CALL POINT_TO_VECTOR (INTEL%X, NM, NL, XS)
   CALL POINT_TO_VECTOR (INTEL%F, NM, NL, FS)

   MKL%MSGLVL =  0            ! do not print statistical information any more
   MKL%PHASE  = 33            ! only solving
   CALL CLUSTER_SPARSE_SOLVER(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC_GLOBAL, &
                              AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                              MKL%MSGLVL, FS, XS, MPI_COMM_WORLD, MKL%ERROR)

   IF (MKL%ERROR /= 0) WRITE(LU_SCARC,*) 'The following ERROR was detected: ', MKL%ERROR

ENDDO MESHES_LOOP

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS  (NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS (NLEVEL_MIN)
ENDIF

CALL SCARC_RESTORE_PARENT(PARENT)

TSTEP(MYID+1)%CLUSTER=MAX(TSTEP(MYID+1)%CLUSTER,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%CLUSTER =TSUM(MYID+1)%CLUSTER+CURRENT_TIME()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_CLUSTER
#endif


!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
#ifdef WITH_MKL
SUBROUTINE SCARC_METHOD_PARDISO(NSCOPE)
INTEGER, INTENT(IN) :: NSCOPE
INTEGER ::  NM, NL !, IP, IC
REAL (EB) :: TNOW
REAL(EB), POINTER, DIMENSION(:) :: XS, FS
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L
TYPE (SCARC_MKL_TYPE)   , POINTER :: MKL
TYPE (SCARC_MATRIX_TYPE), POINTER :: AS
TYPE (SCARC_PARENT_TYPE) :: PARENT
TYPE (SCARC_SCOPE_TYPE) :: INTEL

SCARC_ROUTINE = 'SCARC_METHOD_PARDISO'
TNOW = CURRENT_TIME()

CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(INTEL, NSCOPE, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L   => SCARC(NM)%LEVEL(NL)
   AS  => SCARC(NM)%LEVEL(NL)%AS
   MKL => SCARC(NM)%LEVEL(NL)%MKL

   MKL%MSGLVL =  0         ! do not print statistical information
   MKL%PHASE  = 33         ! only solving

   CALL POINT_TO_VECTOR (INTEL%X, NM, NL, XS)
   CALL POINT_TO_VECTOR (INTEL%F, NM, NL, FS)

   CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NCS, &
                  AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                  MKL%MSGLVL, FS, XS, MKL%ERROR)

   IF (MKL%ERROR /= 0) WRITE(LU_SCARC,*) 'The following ERROR was detected: ', MKL%ERROR

ENDDO MESHES_LOOP

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RESTORE_PARENT(PARENT)

TSTEP(MYID+1)%PARDISO=MAX(TSTEP(MYID+1)%PARDISO,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%PARDISO =TSUM(MYID+1)%PARDISO+CURRENT_TIME()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_PARDISO
#endif


!> ------------------------------------------------------------------------------------------------
!> Perform global CG-method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CG(NSCOPE)
INTEGER, INTENT(IN) :: NSCOPE
INTEGER   :: NL
INTEGER   :: ISTATE, ITE
REAL (EB) :: SIGMA0, SIGMA1, ALPHA0, GAMMA0
REAL (EB) :: TNOW
TYPE (SCARC_SCOPE_TYPE) :: CG
TYPE (SCARC_PARENT_TYPE) :: PARENT

SCARC_ROUTINE = 'SCARC_METHOD_CG'
TNOW = CURRENT_TIME()

!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Save SETTINGS (in case that subsequent solvers with different settings are called)
!>   - Define parameters for current scope (note: NL denotes the finest level)
!>   - Initialize solution, right hand side vector
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(CG, NSCOPE, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

!SCARC(1)%LEVEL(1)%MAIN%F = 0.0_EB
!SCARC(1)%LEVEL(1)%MAIN%F(25) = 1000.0_EB
!SCARC(1)%LEVEL(1)%MAIN%F(33) = 1000.0_EB

!> ------------------------------------------------------------------------------------------------
!> Compute initial residual and perform initial preconditioning
!> ------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (CG%X, CG%W, NL)                         !>  W := A*X
CALL SCARC_VECTOR_SUM     (CG%F, CG%W, -1.0_EB, 1.0_EB, NL)        !>  W := W - F

CG%RES = SCARC_L2NORM (CG%W, NL)                                   !>  RESIN := ||W||
CG%RESIN = CG%RES

call scarc_debug_level (CG%X, 'X INIT', NL)
call scarc_debug_level (CG%F, 'F INIT', NL)
call scarc_debug_level (CG%W, 'W INIT', NL)

ISTATE = SCARC_CONVERGENCE_STATE (CG, 0, NL)                       !>  RES < TOL already ??

IF (ISTATE /= NSCARC_STATE_CONV0) THEN                             !>  if no convergence yet, start precon
   CALL SCARC_PRECONDITIONER(CG, NSCOPE, NL)
   SIGMA0 = SCARC_SCALAR_PRODUCT(CG%W, CG%G, NL)                   !>  SIGMA0 := (W,G)
   CALL SCARC_VECTOR_COPY (CG%G, CG%D, -1.0_EB, NL)                !>  D := -G
ELSE
   CG%NIT=0                                                        !>  if already convergence, don't iterate
ENDIF

call scarc_debug_level (CG%G, 'G INIT ', NL)
call scarc_debug_level (CG%D, 'D INIT ', NL)

!> ------------------------------------------------------------------------------------------------
!> Perform conjugate gradient looping
!> ------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, CG%NIT

   ITE_CG = ITE

   CALL SCARC_MATVEC_PRODUCT (CG%D, CG%Y, NL)                      !>  Y := A*D

   ALPHA0 = SCARC_SCALAR_PRODUCT (CG%D, CG%Y, NL)                  !>  ALPHA0 := (D,Y)
   ALPHA0 = SIGMA0/ALPHA0

   CALL SCARC_VECTOR_SUM (CG%D, CG%X, ALPHA0, 1.0_EB, NL)          !>  X := ALPHA0*D + X
   CALL SCARC_VECTOR_SUM (CG%Y, CG%W, ALPHA0, 1.0_EB, NL)          !>  W := ALPHA0*Y + W

call scarc_debug_level (CG%D, 'D ITE FINE', NL)
call scarc_debug_level (CG%X, 'X ITE FINE', NL)
call scarc_debug_level (CG%Y, 'Y ITE FINE', NL)
call scarc_debug_level (CG%W, 'W ITE FINE', NL)

   CG%RES = SCARC_L2NORM (CG%W, NL)                                !>  RES := ||W||
   ISTATE = SCARC_CONVERGENCE_STATE (CG, ITE, NL)                  !>  RES < TOL ??
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP

   CALL SCARC_PRECONDITIONER(CG, NSCOPE, NL)

   SIGMA1 = SCARC_SCALAR_PRODUCT (CG%W, CG%G, NL)                  !>  SIGMA1 := (W,G)
   GAMMA0 = SIGMA1/SIGMA0
   SIGMA0 = SIGMA1

   CALL SCARC_VECTOR_SUM (CG%G, CG%D, -1.0_EB, GAMMA0, NL)         !>  D := -G + GAMMA0*D

CALL SCARC_DEBUG_LEVEL (CG%G, 'G ITE FINE2', NL)
CALL SCARC_DEBUG_LEVEL (CG%D, 'D ITE FINE2', NL)

ENDDO CG_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of CG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(CG, ITE, ISTATE)

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS  (NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS (NLEVEL_MIN)
ENDIF

call scarc_debug_level (CG%F, 'F FINAL', NL)
call scarc_debug_level (CG%X, 'X FINAL', NL)

CALL SCARC_RESTORE_PARENT(PARENT)

TSTEP(MYID+1)%KRYLOV=MAX(TSTEP(MYID+1)%KRYLOV,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%KRYLOV =TSUM(MYID+1)%KRYLOV+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_METHOD_CG


!> ------------------------------------------------------------------------------------------------
!> Perform global BICGstab-method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_BICG(NSCOPE)
INTEGER, INTENT(IN) :: NSCOPE
INTEGER   :: NL = NSCARC_UNDEFINED
INTEGER   :: ISTATE, ITE
REAL (EB) :: ALPHA0, ALPHA1, ALPHA2, RHO0, RHO1, DTHETA, DBETA
REAL (EB) :: TNOW
TYPE (SCARC_SCOPE_TYPE)  :: BICG
TYPE (SCARC_PARENT_TYPE) :: PARENT

SCARC_ROUTINE = 'SCARC_METHOD_BICG'
TNOW = CURRENT_TIME()

!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Save SETTINGS (in case that subsequent solvers with different settings are called)
!>   - Define parameters for current scope (note: NL denotes the finest level)
!>   - Initialize solution, right hand side vector
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(BICG, NSCOPE, NL)
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
!CALL SCARC_BLOCK_SOLVER (BICG%W, BICG%W, NSCOPE, NL)                        !>  W := PRECON(W)
CALL SCARC_MATVEC_PRODUCT  (BICG%X, BICG%W, NL)                                 !>  W := A*X
CALL SCARC_VECTOR_SUM      (BICG%F, BICG%W, 1.0_EB, -1.0_EB, NL)                !>  W := F - W
CALL SCARC_BLOCK_SOLVER (BICG%W, BICG%W, NSCOPE, NL)                         !>  W := PRECON(W)

BICG%RESIN = SCARC_L2NORM (BICG%W, NL)                                          !>  RESIN := ||W||
CALL SCARC_CONVERGENCE_INFO (BICG%RESIN, 0, NL, BICG%CSOLVER)

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
   CALL SCARC_BLOCK_SOLVER (BICG%Y, BICG%Y, NSCOPE, NL)                      !> Z := PRECON(Z)

   DTHETA = SCARC_SCALAR_PRODUCT (BICG%G, BICG%Y, NL)                           !> DTHETA := (G,Y)
   DTHETA = RHO1/DTHETA

   CALL SCARC_VECTOR_SUM      (BICG%Y, BICG%W, -DTHETA, 1.0_EB, NL)             !> W := -DTHETA*Y + W
   CALL SCARC_MATVEC_PRODUCT  (BICG%W, BICG%D, NL)                              !> D := A*W
   CALL SCARC_BLOCK_SOLVER (BICG%D, BICG%D, NSCOPE, NL)                      !> D := PRECON(D)

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
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS(NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RESTORE_PARENT(PARENT)

TSTEP(MYID+1)%KRYLOV=MAX(TSTEP(MYID+1)%KRYLOV,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%KRYLOV =TSUM(MYID+1)%KRYLOV+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_METHOD_BICG

!> -----------------------------------------------------------------------------------------------
!> Preconditioning method
!> -----------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRECONDITIONER(SC, NSCOPE, NL)
TYPE (SCARC_SCOPE_TYPE), INTENT(IN):: SC
INTEGER, INTENT(IN) :: NSCOPE, NL

SCARC_ROUTINE = 'SCARC_PRECONDITIONER'
SELECT CASE (TYPE_TWOLEVEL)

   !> classical one-level preconditioning
   CASE (NSCARC_TWOLEVEL_NONE)

call scarc_debug_level (SC%W, 'NONE: W BEFORE', NL)
      CALL SCARC_VECTOR_COPY (SC%W, SC%G, 1.0_EB, NL)           !>  G := W 
      CALL SCARC_BLOCK_SOLVER (SC%W, SC%G, NSCOPE, NL)          !>  G := PRECON(W)

call scarc_debug_level (SC%G, 'NONE: G AFTER BLOCK', NL) 

   !> additive two-level preconditioning 
   CASE (NSCARC_TWOLEVEL_ADD)

      CALL SCARC_RESTRICTION (SC%W, SC%F, NL, NL+1)             !>  F_coarse := rest(W_fine)
      CALL SCARC_COARSE_SOLVER                                  !>  X_coarse := A_coarse^{-1}(F_coarse)
      CALL SCARC_PROLONGATION (SC%X, SC%Z, NL+1, NL)            !>  Z := prol(X_coarse)
      CALL SCARC_VECTOR_COPY (SC%W, SC%G, 1.0_EB, NL)           !>  G := W
      CALL SCARC_BLOCK_SOLVER (SC%W, SC%G, NSCOPE, NL)          !>  G := PRECON(W)

call scarc_debug_level (SC%W, 'ADD: W BEFORE', NL)
call scarc_debug_level (SC%F, 'ADD: F INIT COARSE', NL+1)
call scarc_debug_level (SC%X, 'ADD: X AFTER COARSE', NL+1)
call scarc_debug_level (SC%Z, 'ADD: Z AFTER PROL', NL)
call scarc_debug_level (SC%G, 'ADD: G AFTER BLOCK', NL)

      CALL SCARC_VECTOR_SUM (SC%Z, SC%G, 1.0_EB, 1.0_EB, NL)    !>  G := Z + G 

call scarc_debug_level (SC%G, 'ADD: G SUM', NL)

   !> multiplicative two-level preconditioning (coarse first, fine second)
   CASE (NSCARC_TWOLEVEL_MUL)

      CALL SCARC_RESTRICTION (SC%W, SC%F, NL, NL+1)             !>  F_coarse := rest(W_fine)
call scarc_debug_level (SC%W, 'MUL: W BEFORE COARSE', NL)
call scarc_debug_level (SC%F, 'MUL: F INIT COARSE', NL+1)

      CALL SCARC_COARSE_SOLVER                                  !>  X_coarse := A_coarse^{-1}(F_coarse)
call scarc_debug_level (SC%X, 'MUL: X AFTER COARSE', NL+1)
      CALL SCARC_PROLONGATION (SC%X, SC%Y, NL+1, NL)            !>  G := prol(X_coarse)
call scarc_debug_level (SC%Y, 'MUL: Y AFTER PROL', NL)
      CALL SCARC_MATVEC_PRODUCT (SC%Y, SC%Z, NL)                !>  Z := A_fine*G

call scarc_debug_level (SC%Z, 'MUL: Z NEW MATVEC', NL)

      CALL SCARC_VECTOR_SUM (SC%W, SC%Z, 1.0_EB, -1.0_EB, NL)   !>  Z := W - Z
call scarc_debug_level (SC%Z, 'MUL: Z NEW DEFECT', NL)
      CALL SCARC_VECTOR_COPY (SC%Z, SC%G, 1.0_EB, NL)           !>  G := Z
      CALL SCARC_BLOCK_SOLVER (SC%Z, SC%G, NSCOPE, NL)          !>  G := PRECON(Z)
call scarc_debug_level (SC%G, 'MUL: G AFTER BLOCK', NL)
      CALL SCARC_VECTOR_SUM (SC%Y, SC%G, 1.0_EB, 1.0_EB, NL)   !>  Z := W - Z
call scarc_debug_level (SC%G, 'MUL: G AFTER ADD TO Y', NL)


   !> multiplicative two-level preconditioning (fine first, coarse second)
   CASE (NSCARC_TWOLEVEL_MUL2)

      CALL SCARC_VECTOR_COPY (SC%W, SC%G, 1.0_EB, NL)           !>  G := W
      CALL SCARC_BLOCK_SOLVER (SC%W, SC%G, NSCOPE, NL)          !>  G := PRECON(W)
      CALL SCARC_MATVEC_PRODUCT (SC%G, SC%Z, NL)                !>  Z := A_fine*G

call scarc_debug_level (SC%W, 'ADD: W BEFORE', NL)
call scarc_debug_level (SC%G, 'ADD: G AFTER BLOCK', NL)
call scarc_debug_level (SC%Z, 'MUL: Z NEW MATVEC', NL)

      CALL SCARC_VECTOR_SUM (SC%W, SC%Z, 1.0_EB, -1.0_EB, NL)   !>  Z := W - Z

call scarc_debug_level (SC%Z, 'MUL: Z NEW DEFECT', NL)
      CALL SCARC_RESTRICTION (SC%Z, SC%F, NL, NL+1)             !>  F_coarse := rest(W_fine)
      CALL SCARC_COARSE_SOLVER                                  !>  X_coarse := A_coarse^{-1}(F_coarse)
      CALL SCARC_PROLONGATION (SC%X, SC%Z, NL+1, NL)            !>  G := prol(X_coarse)
      CALL SCARC_VECTOR_SUM (SC%Z, SC%G, 1.0_EB, 1.0_EB, NL)   !>  Z := W - Z

call scarc_debug_level (SC%F, 'MUL: F INIT COARSE', NL+1)
call scarc_debug_level (SC%X, 'MUL: X AFTER COARSE', NL+1)
call scarc_debug_level (SC%G, 'MUL: G AFTER BLOCK', NL)

END SELECT 

END SUBROUTINE SCARC_PRECONDITIONER


!> ------------------------------------------------------------------------------------------------
!> Call requested coarse grid solver (iterative/direct)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_COARSE_SOLVER

SELECT CASE (TYPE_COARSE)
   CASE (NSCARC_COARSE_ITERATIVE)
      CALL SCARC_METHOD_CG (NSCARC_SCOPE_COARSE)
   CASE (NSCARC_COARSE_DIRECT)
#if defined(WITH_MKL)
      IF (N_MPI_PROCESSES > 1) THEN
         CALL SCARC_METHOD_CLUSTER (NSCARC_SCOPE_COARSE)             !> call CLUSTER_SPARSE_SOLVER
      ELSE
         CALL SCARC_METHOD_PARDISO (NSCARC_SCOPE_COARSE)             !> call PARDISO
      ENDIF
#else
      WRITE(*,*) 'SCARC_METHOD_DIRECT not working yet '
#endif
END SELECT

END SUBROUTINE SCARC_COARSE_SOLVER


!> ------------------------------------------------------------------------------------------------
!> Perform geometric multigrid method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MULTIGRID(NSCOPE)
INTEGER, INTENT(IN) :: NSCOPE
INTEGER   :: NL = NSCARC_UNDEFINED
INTEGER   :: ISTATE, ICYCLE, ITE
REAL (EB) :: TNOW, TNOW_COARSE
TYPE (SCARC_SCOPE_TYPE)  :: MG
TYPE (SCARC_PARENT_TYPE) :: PARENT

SCARC_ROUTINE = 'SCARC_METHOD_MULTIGRID'
TNOW = CURRENT_TIME()

!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Save SETTINGS (in case that subsequent solvers with different settings are called)
!>   - Define parameters for current scope (note: NL denotes the finest level)
!>   - Initialize solution, right hand side vector
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SAVE_PARENT(PARENT)
CALL SCARC_SETUP_SCOPE(MG, NSCOPE, NL)
CALL SCARC_SETUP_WORKSPACE(NL)

!SCARC(1)%LEVEL(1)%MAIN%F = 0.0_EB
!SCARC(1)%LEVEL(1)%MAIN%F(25) = 1000.0_EB
!SCARC(1)%LEVEL(1)%MAIN%F(33) = 1000.0_EB

!> ------------------------------------------------------------------------------------------------
!> Compute initial defect:  RESIN := || F - A*X ||
!>   - Initialize cycle counts for MG-iteration
!>   - Perform initial matrix-vector product on finest level
!>   - calculate norm of initial residual on finest level
!> ------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (MG%X, MG%D, NL)                                       !>  D := A*X
CALL SCARC_VECTOR_SUM     (MG%F, MG%D, 1.0_EB, -1.0_EB, NL)                      !>  D := F - D

call scarc_debug_level (MG%X, 'X BEFORE', NL)
call scarc_debug_level (MG%F, 'F BEFORE', NL)
call scarc_debug_level (MG%D, 'D BEFORE', NL)

ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_SETUP, NL)
MG%RESIN = SCARC_L2NORM (MG%D, NL)                                               !>  RESIN := ||D||

CALL SCARC_CONVERGENCE_INFO(MG%RESIN, 0, NL, MG%CSOLVER)

!> ------------------------------------------------------------------------------------------------
!> Perform multigrid-looping (start each iteration on finest level)
!> ------------------------------------------------------------------------------------------------
MULTIGRID_LOOP: DO ITE = 1, MG%NIT

   ITE_MG = ITE

   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLE_EXIT)

      !> presmoothing  (smoothing/restriction till coarsest level is reached)
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)
         CALL SCARC_SMOOTHER (MG, NSCARC_CYCLE_PRESMOOTH, NSCOPE, NL)            !> D_fine   := smooth(defect)
call scarc_debug_level (MG%D, 'D AFTER SMOOTHING', NL)
call scarc_debug_level (MG%F, 'F AFTER SMOOTHING', NL+1)
         CALL SCARC_RESTRICTION (MG%D, MG%F, NL, NL+1)                           !> F_coarse := rest(D_fine)
call scarc_debug_level (MG%D, 'D AFTER RESTRICTION', NL)
call scarc_debug_level (MG%F, 'F AFTER RESTRICTION', NL+1)
         CALL SCARC_VECTOR_CLEAR (MG%X, NL+1)                                    !> X_coarse := 0.0
         NL = NL + 1                                                             !> set coarser level
      ENDDO PRESMOOTHING_LOOP

      !> coarse grid solver
      TNOW_COARSE = CURRENT_TIME()
      CALL SCARC_COARSE_SOLVER                                                   !> X_coarse := exact_sol(.)
call scarc_debug_level (MG%X, 'X AFTER COARSE', NLEVEL_MAX)
      TSTEP(MYID+1)%COARSE=MAX(TSTEP(MYID+1)%COARSE,CURRENT_TIME()-TNOW_COARSE)
      TSUM(MYID+1)%COARSE =TSUM(MYID+1)%COARSE+CURRENT_TIME()-TNOW_COARSE

      !> postsmoothing (smoothing/restriction till finest level is reached again)
      POSTSMOOTHING_LOOP: DO WHILE (NL > NLEVEL_MIN)
         NL=NL-1
         CALL SCARC_PROLONGATION (MG%X, MG%D, NL+1, NL)                          !> D_fine := prol(X_coarse)
call scarc_debug_level (MG%D, 'D AFTER PROLONGATION', NL)
         CALL SCARC_VECTOR_SUM (MG%D, MG%X, 1.0_EB, 1.0_EB, NL)                  !> X_fine := D_fine + X_fine
call scarc_debug_level (MG%X, 'X AFTER SUMMING', NL)
         CALL SCARC_SMOOTHER (MG, NSCARC_CYCLE_POSTSMOOTH, NSCOPE, NL)           !> D_fine := smooth(defect)
call scarc_debug_level (MG%D, 'D AFTER PROLONGATION', NL)
         ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_PROCEED, NL)                  !> perform requested cycle
         IF (ICYCLE /= NSCARC_CYCLE_POSTSMOOTH) CYCLE CYCLE_LOOP
      ENDDO POSTSMOOTHING_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) THEN
      WRITE(SCARC_MESSAGE,'(2A,I8,A)') TRIM(SCARC_ROUTINE),': Wrong level ', NL,' for multigrid method'
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   ENDIF

 !> ---------------------------------------------------------------------------------------------
 !> Compute norm of new residual on finest level and  leave loop correspondingly
 !> ---------------------------------------------------------------------------------------------
   CALL SCARC_MATVEC_PRODUCT (MG%X, MG%D, NL)                                    !> D := A*X
   CALL SCARC_VECTOR_SUM     (MG%F, MG%D, 1.0_EB, -1.0_EB, NL)                   !> D := F - D

   MG%RES = SCARC_L2NORM (MG%D, NL)                                              !> RES := ||D||
   ISTATE = SCARC_CONVERGENCE_STATE(MG, ITE, NL)                                 !> convergence ?

call scarc_debug_level (MG%D, 'D NEW DEFECT ', NL)

   IF (TYPE_DEBUG>=NSCARC_DEBUG_INFO1.AND.MYID==0) &
      WRITE(LU_OUTPUT,'(a,i3,a,e14.5,a,e14.5)') '       SCARC_MG-Iteration  =',ITE,': Residuum=',SCARC_RESIDUAL
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
      CALL SCARC_UPDATE_PRESSURE_MAINCELLS(NLEVEL_MIN)
      CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
   CASE (NSCARC_SCOPE_PRECON)
      CALL SCARC_UPDATE_PRECONDITIONER(NLEVEL_MIN)
END SELECT

call scarc_debug_level (MG%X, 'X FINAL ', NL)

CALL SCARC_RESTORE_PARENT(PARENT)

TSTEP(MYID+1)%MULTIGRID=MAX(TSTEP(MYID+1)%MULTIGRID,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%MULTIGRID =TSUM(MYID+1)%MULTIGRID+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_METHOD_MULTIGRID


!> ------------------------------------------------------------------------------------------------
!> Control multigrid cycling (F/V/W)
!> Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CYCLE_CONTROL(NSCOPE, NL)
INTEGER, INTENT(IN) :: NSCOPE, NL
INTEGER :: NM, NL0, ICYCLE

SELECT CASE (NSCOPE)

 !> ---------------------------------------------------------------------------------------------
 !> initialize cycle counts at beginning of multigrid method
 !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_SETUP)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         SCARC(NM)%LEVEL(NL)%GMG%CYCLE_COUNT(2)=1
         DO NL0 = NLEVEL_MIN+1, NLEVEL_MAX - 1
            IF (TYPE_CYCLE==NSCARC_CYCLE_F) THEN
               SCARC(NM)%LEVEL(NL0)%GMG%CYCLE_COUNT(2)=2
            ELSE
               SCARC(NM)%LEVEL(NL0)%GMG%CYCLE_COUNT(2)=TYPE_CYCLE
            ENDIF
         ENDDO
      ENDDO

      ICYCLE = NSCARC_CYCLE_NEXT

 !> ---------------------------------------------------------------------------------------------
 !> reset cycle counts at beginning of each new multigrid iteration
 !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_RESET)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         DO NL0 = NLEVEL_MIN, NLEVEL_MAX
            SCARC(NM)%LEVEL(NL0)%GMG%CYCLE_COUNT(1)=SCARC(NM)%LEVEL(NL0)%GMG%CYCLE_COUNT(2)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLE_NEXT

 !> ---------------------------------------------------------------------------------------------
 !> determine where to proceed with cycling
 !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_PROCEED)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         SCARC(NM)%LEVEL(NL)%GMG%CYCLE_COUNT(1)=SCARC(NM)%LEVEL(NL)%GMG%CYCLE_COUNT(1)-1

         IF (SCARC(NM)%LEVEL(NL)%GMG%CYCLE_COUNT(1)==0) THEN
            IF (TYPE_CYCLE==NSCARC_CYCLE_F) THEN
               SCARC(NM)%LEVEL(NL)%GMG%CYCLE_COUNT(1)=1
            ELSE
               SCARC(NM)%LEVEL(NL)%GMG%CYCLE_COUNT(1)=SCARC(NM)%LEVEL(NL)%GMG%CYCLE_COUNT(2)
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
SUBROUTINE SCARC_SMOOTHER(MG, NTYPE, NSCOPE, NL)
INTEGER , INTENT(IN) :: NTYPE, NL, NSCOPE
TYPE (SCARC_SCOPE_TYPE), INTENT(INOUT) :: MG
INTEGER :: ITE, ISTATE=0
REAL(EB):: TNOW
LOGICAL :: BMATVEC, BL2NORM
TYPE (SCARC_SCOPE_TYPE) :: SM

SCARC_ROUTINE = 'SCARC_SMOOTHER'
TNOW = CURRENT_TIME()

!> ------------------------------------------------------------------------------------------------
!> Initialization
!> ------------------------------------------------------------------------------------------------
SM%X = MG%X
SM%F = MG%F
SM%D = MG%D
SM%Z = MG%Z

SM%NIT   = SCARC_SMOOTH_ITERATIONS
SM%EPS   = SCARC_SMOOTH_ACCURACY
SM%OMEGA = SCARC_SMOOTH_OMEGA

SELECT CASE (NTYPE)
   CASE (NSCARC_CYCLE_PRESMOOTH)
      SM%CSOLVER = 'SCARC_PRESMOOTHER'
   CASE (NSCARC_CYCLE_POSTSMOOTH)
      SM%CSOLVER = 'SCARC_POSTSMOOTHER'
END SELECT

BL2NORM  = .TRUE.
IF (NTYPE == NSCARC_CYCLE_PRESMOOTH.AND.NL==1) THEN
   BMATVEC = .FALSE.
ELSE
   BMATVEC = .TRUE.
ENDIF
BMATVEC = .TRUE.
BL2NORM = .FALSE.

IF (TYPE_SMOOTH == NSCARC_SMOOTH_CLUSTER) THEN
  SM%OMEGA=1.0_EB
  SM%NIT = 1
ENDIF
IF (TYPE_SMOOTH == NSCARC_SMOOTH_PARDISO) SM%OMEGA = 1.0_EB

!> ------------------------------------------------------------------------------------------------
!> Calculate initial defect on level NL (only if BMATVEC = .TRUE.)
!> Because initial vector is set to zero, this defect corresponds to F
!> ------------------------------------------------------------------------------------------------
IF (BMATVEC) THEN
   CALL SCARC_MATVEC_PRODUCT (SM%X, SM%D, NL)                                  !>  D := A*X
   CALL SCARC_VECTOR_SUM     (SM%F, SM%D, 1.0_EB, -1.0_EB, NL)                 !>  D := F - D
ENDIF

CALL SCARC_DEBUG_LEVEL (SM%X, 'SMOOTH: X INIT', NL)
CALL SCARC_DEBUG_LEVEL (SM%F, 'SMOOTH: F INIT', NL)
CALL SCARC_DEBUG_LEVEL (SM%D, 'SMOOTH: D INIT', NL)


IF (BL2NORM.AND.BMATVEC) THEN
   SM%RESIN = SCARC_L2NORM (SM%D, NL)                                          !>  RESIN := ||D||
ELSE
   SM%RESIN = SCARC_RESIDUAL
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Smoothing loop
!> ------------------------------------------------------------------------------------------------
SMOOTH_LOOP: DO ITE=1, SM%NIT

   ITE_SM = ITE

   IF (TYPE_SMOOTH == NSCARC_SMOOTH_PARDISO .OR. TYPE_SMOOTH == NSCARC_SMOOTH_CLUSTER) THEN
      CALL SCARC_VECTOR_COPY(SM%D, SM%Z, 1.0_EB, NL)
      CALL SCARC_BLOCK_SOLVER (SM%Z, SM%D, NSCOPE, NL)                      !>  D := PRECON (D)
CALL SCARC_DEBUG_LEVEL (SM%D, 'SMOOTH: D ITE CLUSTER', NL)
   ELSE
      CALL SCARC_BLOCK_SOLVER (SM%D, SM%D, NSCOPE, NL)                      !>  D := PRECON (D)
CALL SCARC_DEBUG_LEVEL (SM%D, 'SMOOTH: D ITE ', NL)
   ENDIF

   CALL SCARC_VECTOR_SUM      (SM%D, SM%X, SM%OMEGA, 1.0_EB, NL)               !>  X := OMEGA*D + X
   CALL SCARC_MATVEC_PRODUCT  (SM%X, SM%D, NL)                                 !>  D := A*X

   CALL SCARC_VECTOR_SUM      (SM%F, SM%D, 1.0_EB, -1.0_EB, NL)                !>  D := F - D

CALL SCARC_DEBUG_LEVEL (SM%X, 'SMOOTH: X ITE ', NL)
CALL SCARC_DEBUG_LEVEL (SM%D, 'SMOOTH: D ITE ', NL)

   IF (BL2NORM.OR.ITE==SM%NIT) THEN
      SM%RES = SCARC_L2NORM (SM%D, NL)                                         !>  RES := ||D||
      ISTATE = SCARC_CONVERGENCE_STATE(SM, ITE, NL)
      IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP                     !>  RES < TOL ?
   ENDIF

ENDDO SMOOTH_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> ------------------------------------------------------------------------------------------------
IF (BL2NORM) CALL SCARC_CONVERGENCE_RATE(SM, ITE, ISTATE)

TSTEP(MYID+1)%SMOOTH=MAX(TSTEP(MYID+1)%SMOOTH,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%SMOOTH =TSUM(MYID+1)%SMOOTH+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SMOOTHER


!> ----------------------------------------------------------------------------------------------------
!> Save rhs for using ScaRC as postprocessing after an usual FFT-run
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SAVE_RHS
INTEGER :: NM, IC, I, J, K
REAL(EB), POINTER, DIMENSION(:) :: F
TYPE (MESH_TYPE), POINTER ::  M
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   M => MESHES(NM)
   L => SCARC(NM)%LEVEL(NLEVEL_MIN)

   SELECT CASE (TYPE_METHOD)
      CASE (NSCARC_METHOD_KRYLOV)
         F => L%MAIN%F
      CASE (NSCARC_METHOD_MULTIGRID)
         F => L%PRECON%F
   END SELECT

   DO K = 1, M%KBAR
      DO J = 1, M%JBAR
         DO I = 1, M%IBAR
            IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%DISCRET(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
            IC = L%DISCRET(I,J,K,UNKH)
            F(IC) = M%PRHS (I, J, K)
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SAVE_RHS


!> ------------------------------------------------------------------------------------------------
!> Setup environement in every solver CALL (i.e. set pointers to used vectors)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SCOPE(SCOPE, NSCOPE, NL)
INTEGER, INTENT(IN)  :: NSCOPE
INTEGER, INTENT(OUT) :: NL
TYPE (SCARC_SCOPE_TYPE), INTENT(OUT):: SCOPE

TYPE_SCOPE  = NSCOPE

!> ---------------------------------------------------------------------------------------------------
!> change current solution method if called from preconditioner or coarse grid solver
!> if call from main solver, everything is already set
!> ---------------------------------------------------------------------------------------------------
SELECT CASE (NSCOPE)
   CASE (NSCARC_SCOPE_COARSE)
      SELECT CASE (TYPE_COARSE)
         CASE( NSCARC_COARSE_ITERATIVE)
            TYPE_METHOD   = NSCARC_METHOD_KRYLOV
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE
            TYPE_PRECON   = NSCARC_PRECON_SSOR
         CASE( NSCARC_COARSE_DIRECT)
            TYPE_METHOD   = NSCARC_METHOD_MKL
      END SELECT
   CASE (NSCARC_SCOPE_PRECON)
      TYPE_METHOD = NSCARC_METHOD_MULTIGRID
END SELECT

SCOPE%RESIN = 1.0_EB
SCOPE%RES   = 0.0_EB

SELECT CASE (TYPE_METHOD)

   !> ---------------------------------------------------------------------------------------------------
   !> Krylov method
   !> ---------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT CASE (NSCOPE)
         CASE (NSCARC_SCOPE_MAIN)

            NL = NLEVEL_MIN

            SCOPE%CSOLVER = 'SCARC_GLOBAL_KRYLOV'
            SCOPE%EPS     =  SCARC_KRYLOV_ACCURACY
            SCOPE%NIT     =  SCARC_KRYLOV_ITERATIONS
            SCOPE%NTYPE   =  NSCARC_SCOPE_MAIN
            SCOPE%X = NSCARC_VECTOR_X
            SCOPE%F = NSCARC_VECTOR_F
            SCOPE%D = NSCARC_VECTOR_D
            SCOPE%G = NSCARC_VECTOR_G
            SCOPE%Y = NSCARC_VECTOR_Y
            SCOPE%W = NSCARC_VECTOR_W
            SCOPE%Z = NSCARC_VECTOR_Z

         CASE (NSCARC_SCOPE_COARSE)

            NL = NLEVEL_MAX

            SCOPE%EPS      =  SCARC_COARSE_ACCURACY
            SCOPE%NIT      =  SCARC_COARSE_ITERATIONS

            IF (TYPE_PRECON_CORE == NSCARC_PRECON_MULTIGRID) THEN   ! used in preconditioning MG
               SCOPE%CSOLVER = 'SCARC_COARSE_KRYLOV_PRECON'
               SCOPE%NTYPE   =  NSCARC_SCOPE_PRECON
               SCOPE%X = NSCARC_VECTOR_X
               SCOPE%F = NSCARC_VECTOR_F
               SCOPE%D = NSCARC_VECTOR_D
               SCOPE%G = NSCARC_VECTOR_G
               SCOPE%Y = NSCARC_VECTOR_Y
               SCOPE%W = NSCARC_VECTOR_W
               SCOPE%Z = NSCARC_VECTOR_Z
            ELSE                                                    ! used in twolevel Krylov
               SCOPE%CSOLVER = 'SCARC_COARSE_KRYLOV_MAIN'
               SCOPE%NTYPE   =  NSCARC_SCOPE_MAIN
               SCOPE%X = NSCARC_VECTOR_X
               SCOPE%F = NSCARC_VECTOR_F
               SCOPE%D = NSCARC_VECTOR_D
               SCOPE%G = NSCARC_VECTOR_G
               SCOPE%Y = NSCARC_VECTOR_Y
               SCOPE%W = NSCARC_VECTOR_W
               SCOPE%Z = NSCARC_VECTOR_Z

            ENDIF

      END SELECT

   !> ---------------------------------------------------------------------------------------------------
   !> Multigrid method
   !> ---------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_MULTIGRID)

      NL = NLEVEL_MIN

      SCOPE%EPS   = SCARC_MULTIGRID_ACCURACY
      SCOPE%NIT   = SCARC_MULTIGRID_ITERATIONS
      SCOPE%OMEGA = SCARC_SMOOTH_OMEGA

      !> select scope (multigrid as main solver or preconditioner)
      SELECT CASE (NSCOPE)
         CASE (NSCARC_SCOPE_MAIN)
            SCOPE%CSOLVER = 'SCARC_GLOBAL_MULTIGRID'
            SCOPE%NTYPE   =  NSCARC_SCOPE_MAIN
            SCOPE%X = NSCARC_VECTOR_X
            SCOPE%F = NSCARC_VECTOR_F
            SCOPE%D = NSCARC_VECTOR_D
            SCOPE%Z = NSCARC_VECTOR_Z
         CASE (NSCARC_SCOPE_PRECON)
            SCOPE%CSOLVER = 'SCARC_PRECON_MULTIGRID'
            SCOPE%NTYPE   =  NSCARC_SCOPE_PRECON
            SCOPE%X = NSCARC_VECTOR_X
            SCOPE%F = NSCARC_VECTOR_F
            SCOPE%D = NSCARC_VECTOR_D
            SCOPE%Z = NSCARC_VECTOR_Z
      END SELECT

   !> ---------------------------------------------------------------------------------------------------
   !> MKL method
   !> ---------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_MKL)

      NL = NLEVEL_MAX

      SCOPE%CSOLVER = 'SCARC_COARSE_MKL'
      SCOPE%NTYPE   =  NSCARC_SCOPE_MAIN
      SCOPE%X = NSCARC_VECTOR_X
      SCOPE%F = NSCARC_VECTOR_F

END SELECT

END SUBROUTINE SCARC_SETUP_SCOPE


!> ----------------------------------------------------------------------------------------------------
!> Set initial solution corresponding to boundary data in BXS, BXF, ...
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WORKSPACE(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, I, J, K, IC
REAL(EB), POINTER, DIMENSION(:) :: X, F
TYPE (MESH_TYPE), POINTER ::  M
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
TYPE (SCARC_WALL_TYPE), POINTER :: WC

SCARC_ROUTINE = 'SCARC_SETUP_WORKSPACE'

SELECT CASE (TYPE_SCOPE)

   !> --------------- IF used as main solver use values from pres-routine as initialization 
   CASE (NSCARC_SCOPE_MAIN) 

      !> 2D: Initialize solution and right hand side vector corresponding to boundary conditions
      IF (TWO_D) THEN
   
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
            M => MESHES(NM)
            L => SCARC(NM)%LEVEL(NL)
   
            X => L%MAIN%X
            F => L%MAIN%F
   
            !> get right hand side (PRHS from pres.f90) and initial vector (H or HS from last time step)
            DO K = 1, M%KBAR
               DO I = 1, M%IBAR
                  IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%DISCRET(I,1,K,CGSC) /= IS_GASPHASE) CYCLE
                  IC = L%DISCRET(I,1,K,UNKH)
                  !IF (.NOT.(SCARC_FFT.AND.SCARC_PRESSURE_ITERATIONS=1)) F(IC) = M%PRHS (I, J, K)
                  F(IC) = M%PRHS (I, 1, K)
                  IF (PREDICTOR) THEN
                     X(IC) = M%H (I, 1, K)
                  ELSE
                     X(IC) = M%HS(I, 1, K)
                  ENDIF
               ENDDO
            ENDDO
   
            !> set boundary conditions at exterior boundaries (corresponding to pois.f90)
            LEVEL_WALL_CELLS_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
   
               WC => L%WALL(IW)

               I    = WC%IXW
               J    = WC%IYW
               K    = WC%IZW
   
               IF (J /= 1) THEN
                  WRITE(SCARC_MESSAGE,'(2A,I8)') TRIM(SCARC_ROUTINE),': Wrong index for J=',J
                  CALL SHUTDOWN(SCARC_MESSAGE)
               ENDIF
   
               IOR0 = WC%IOR
   
               IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. L%DISCRET(I,1,K,CGSC) /= IS_GASPHASE) CYCLE
               IC = L%DISCRET(I,1,K,UNKH)
   
               SELECT CASE (IOR0)
                  CASE (1)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DXI2 * M%BXS(1,K)         !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) + L%DXI * M%BXS(1,K)                   !> Neumann
                     ENDIF
                  CASE (-1)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DXI2 *M%BXF(1,K)          !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) - L%DXI *M%BXF(1,K)                    !> Neumann
                     ENDIF
                  CASE (3)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DZI2 * M%BZS(I,1)         !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) + L%DZI * M%BZS(I,1)                   !> Neumann
                     ENDIF
                  CASE (-3)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DZI2 * M%BZF(I,1)         !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) - L%DZI  * M%BZF(I,1)                  !> Neumann
                     ENDIF
               END SELECT
   
            ENDDO LEVEL_WALL_CELLS_LOOP2D
         ENDDO

      !> 3D: Initialize solution and right hand side vector corresponding to boundary conditions
      ELSE
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
            M => MESHES(NM)
            L => SCARC(NM)%LEVEL(NL)
   
            X => L%MAIN%X
            F => L%MAIN%F
   
            !> get right hand side (PRHS from pres.f90) and initial vector (H or HS from last time step)
            DO K = 1, M%KBAR
               DO J = 1, M%JBAR
                  DO I = 1, M%IBAR
                     IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%DISCRET(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
                     IC = L%DISCRET(I,J,K,UNKH)
                     !IF (.NOT.(SCARC_FFT.AND.SCARC_PRESSURE_ITERATIONS=1)) F(IC) = M%PRHS (I, J, K)
                     F(IC) = M%PRHS (I, J, K)
                     IF (PREDICTOR) THEN
                        X(IC) = M%H (I, J, K)
                     ELSE
                        X(IC) = M%HS(I, J, K)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
   
            LEVEL_WALL_CELLS_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
   
               WC => L%WALL(IW)

               I    = WC%IXW
               J    = WC%IYW
               K    = WC%IZW

               IOR0 = WC%IOR
   
               IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. L%DISCRET(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
               IC = L%DISCRET(I,J,K,UNKH)
   
               SELECT CASE (IOR0)
                  CASE (1)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DXI2 * M%BXS(J,K)         !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) + L%DXI * M%BXS(J,K)                   !> Neumann
                     ENDIF
                  CASE (-1)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DXI2 *M%BXF(J,K)          !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) - L%DXI *M%BXF(J,K)                    !> Neumann
                     ENDIF
                  CASE (2)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DYI2 * M%BYS(I,K)         !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) + L%DYI * M%BYS(I,K)                   !> Neumann
                     ENDIF
                  CASE (-2)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DYI2 *M%BYF(I,K)          !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) - L%DYI *M%BYF(I,K)                    !> Neumann
                     ENDIF
                  CASE (3)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DZI2 * M%BZS(I,J)         !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) + L%DZI * M%BZS(I,J)                   !> Neumann
                     ENDIF
                  CASE (-3)
                     IF (WC%BTYPE==DIRICHLET) THEN
                        F(IC) = F(IC) - 2.0_EB * L%DZI2 * M%BZF(I,J)         !> Dirichlet
                     ELSE IF (WC%BTYPE==NEUMANN.AND.WC%NOM==0) THEN
                        F(IC) = F(IC) - L%DZI * M%BZF(I,J)                   !> Neumann
                     ENDIF
               END SELECT
   
            ENDDO LEVEL_WALL_CELLS_LOOP3D
   
         ENDDO
   
      ENDIF

   !> --------------- If MG is used as preconditioner, RHS for MG is transferred in NSCARC_VECTOR_G_MAIN
   CASE (NSCARC_SCOPE_PRECON)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         SCARC(NM)%LEVEL(NL)%PRECON%F = SCARC(NM)%LEVEL(NL)%MAIN%G
      ENDDO

   !> --------------- If used as coarse grid solver start with zero initialization 
   CASE (NSCARC_SCOPE_COARSE)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         IF (BCGGMG) THEN                         !> used within GMG as preconditioner
            SCARC(NM)%LEVEL(NL)%PRECON%X = 0.0_EB
         ELSE                                     !> used within Twolevel-Krylov
            SCARC(NM)%LEVEL(NL)%MAIN%X = 0.0_EB
         ENDIF
      ENDDO

END SELECT


!> In case of a Krylov method clear overlapping parts of auxiliary vectors
IF (BCG) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM)%LEVEL(NL)
      L%MAIN%D(L%NCS+1:L%NCE) = 0.0_EB
      L%MAIN%G(L%NCS+1:L%NCE) = 0.0_EB
      L%MAIN%Y(L%NCS+1:L%NCE) = 0.0_EB
      L%MAIN%W(L%NCS+1:L%NCE) = 0.0_EB
      L%MAIN%Z(L%NCS+1:L%NCE) = 0.0_EB
   ENDDO
   IF (TYPE_KRYLOV == NSCARC_KRYLOV_BICG) THEN
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         L%MAIN%Z(L%NCS+1:L%NCE) = 0.0_EB
      ENDDO
   ENDIF
ENDIF
!> In case of a twolevel-Krylov also clear vectors on max level
IF (BTWOLEVEL) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM)%LEVEL(NL)
      L%MAIN%D(L%NCS+1:L%NCE) = 0.0_EB
      L%MAIN%G(L%NCS+1:L%NCE) = 0.0_EB
      L%MAIN%Y(L%NCS+1:L%NCE) = 0.0_EB
      L%MAIN%W(L%NCS+1:L%NCE) = 0.0_EB
      L%MAIN%Z(L%NCS+1:L%NCE) = 0.0_EB
   ENDDO
ENDIF

!> In case of a multigrid method (as MA solver or preconditioner) also 
!> clear overlapping parts of corresponding auxiliary vectors
IF (BGMG) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM)%LEVEL(NL)
      L%MAIN%D(L%NCS+1:L%NCE) = 0.0_EB
      IF (TYPE_COARSE == NSCARC_COARSE_ITERATIVE) THEN
         L%MAIN%G = 0.0_EB
         L%MAIN%W = 0.0_EB
      ENDIF
   ENDDO
ENDIF
IF (BCGGMG) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM)%LEVEL(NL)
      L%PRECON%D(L%NCS+1:L%NCE) = 0.0_EB
      IF (TYPE_COARSE == NSCARC_COARSE_ITERATIVE) THEN
         L => SCARC(NM)%LEVEL(NLEVEL_MAX)
         L%PRECON%G = 0.0_EB
         L%PRECON%W = 0.0_EB
      ENDIF
   ENDDO
ENDIF


END SUBROUTINE SCARC_SETUP_WORKSPACE


!> ------------------------------------------------------------------------------------------------
!> Print out residual information for loop ITE
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_INFO(RES, ITE, NL, CSOLVER)
INTEGER, INTENT(IN) :: ITE, NL
REAL(EB), INTENT(IN) :: RES
CHARACTER(*), INTENT(IN) :: CSOLVER
INTEGER:: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (TYPE_DEBUG>NSCARC_DEBUG_INFO1.AND.MYID==0) write(LU_OUTPUT,1000) TRIM(CSOLVER), NL, ITE,  RES
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

IF (TYPE_DEBUG>=NSCARC_DEBUG_INFO2.AND.MYID==0) WRITE(LU_OUTPUT,1000) TRIM(SCOPE%CSOLVER), NL, ITE, SCOPE%RES
IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) WRITE(LU_SCARC,1000) TRIM(SCOPE%CSOLVER), NL, ITE, SCOPE%RES

SELECT CASE (TYPE_ACCURACY)
   CASE (NSCARC_ACCURACY_RELATIVE)
      IF (SCOPE%RES <= SCOPE%RESIN*SCOPE%EPS)  ISTATE = NSCARC_STATE_CONV
      IF (SCOPE%RES <= 1.0E-15)                ISTATE = NSCARC_STATE_CONV
   CASE (NSCARC_ACCURACY_ABSOLUTE)
      !IF (SCOPE%RES <= SCOPE%EPS .AND. SCOPE%RES <= SCOPE%RESIN*SCARC_ACCURACY_RELATIVE) THEN
      IF (SCOPE%RES <= SCOPE%EPS .AND. SCOPE%RES <= SCOPE%RESIN) THEN
         IF (ITE == 0) THEN
            ISTATE = NSCARC_STATE_CONV0
         ELSE
            ISTATE = NSCARC_STATE_CONV
         ENDIF
      ENDIF
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
   IF (ISTATE == NSCARC_STATE_CONV0) THEN
     ITERATIONS = 0
   ELSE IF (ISTATE == NSCARC_STATE_CONV) THEN
     ITERATIONS = ITE
   ELSE
     ITERATIONS = ITE-1
   ENDIF
   IF (SCOPE%RESIN >= 1.0E-70_EB) THEN
      IF (ITERATIONS == 0) THEN
         !CAPPA = (SCOPE%RES/SCOPE%RESIN)
         CAPPA = 0.0_EB
      ELSE
         IF (ISTATE == NSCARC_STATE_CONV0) THEN
            CAPPA = 0.0E0
         ELSE
            CAPPA = (SCOPE%RES/SCOPE%RESIN) ** (1.0_EB/ITERATIONS)
         ENDIF
      ENDIF
   ELSE
      CAPPA = 0.0_EB
   ENDIF
ENDIF

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   SCARC_CAPPA      = CAPPA
   SCARC_RESIDUAL   = RESIDUAL
   SCARC_ITERATIONS = ITERATIONS
ENDIF

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(LU_SCARC,2000) SCOPE%CSOLVER, ITERATIONS, CAPPA
IF (MYID==0.AND.TYPE_DEBUG>=NSCARC_DEBUG_INFO1) THEN
   IF (TRIM(SCOPE%CSOLVER) /=  'SCARC_COARSE_CG') WRITE(LU_OUTPUT,2000) SCOPE%CSOLVER, ITERATIONS, CAPPA
ENDIF

2000 FORMAT (/,7X,A25,': iterations: ',i6,':  convergence rate =',e14.6,/)
END SUBROUTINE SCARC_CONVERGENCE_RATE


!> ------------------------------------------------------------------------------------------------
!> Perform restriction from finer to coarser grid in multigrid method
!>    - 'FI' corresponds to finer   grid
!>    - 'CO' corresponds to coarser grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTRICTION (NVECTORF, NVECTORC, NLF, NLC)
INTEGER, INTENT(IN) :: NVECTORF, NVECTORC, NLF, NLC
INTEGER :: NM, ICOL, IC
REAL(EB), POINTER, DIMENSION(:)     :: FCC, DCF, R
INTEGER , POINTER, DIMENSION(:)     :: ROW, COL
INTEGER , POINTER :: NXC, NYC, NZC, NCC
INTEGER  :: NXF, NYF, NZF
INTEGER  :: IXF, IYF, IZF, ICF(8)=0, ICFB(-2:2,-2:2)=0
INTEGER  :: IXC, IYC, IZC, ICC
REAL(EB) :: AUX, AUX0, SCAL, W1, W3, W4, W9, W12, W16
TYPE (SCARC_LEVEL_TYPE), POINTER :: C, F

SCARC_ROUTINE = 'SCARC_RESTRICTION'

SCAL = 0.015625_EB
W16  = 16.0_EB
W12  = 12.0_EB
W9   =  9.0_EB
W4   =  4.0_EB
W3   =  3.0_EB
W1   =  1.0_EB

IF (BAMG.AND.TYPE_COARSENING < NSCARC_COARSENING_GMG) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_VECTOR, NLF)

call scarc_debug_level (NVECTORF, 'NVECTORF', NLF)
call scarc_debug_level (NVECTORC, 'NVECTORC', NLC)
!>
!> ------------------ Twolevel-CG or Geometric multigrid (as main solver or preconditioner) --------------
!>
IF (BMULTILEVEL) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      C => SCARC(NM)%LEVEL(NLC)                 !> pointer to coarse grid
      F => SCARC(NM)%LEVEL(NLF)                 !> pointer to fine grid

      NXC => C%NX
      NYC => C%NY
      NZC => C%NZ

      NXF = 2*NXC
      NYF = 2*NYC
      NZF = 2*NZC

      CALL POINT_TO_VECTOR(NVECTORF, NM, NLF, DCF)
      CALL POINT_TO_VECTOR(NVECTORC, NM, NLC, FCC)

      IF (TWO_D) THEN

         SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)
  
            CASE (NSCARC_INTERPOL_CONSTANT)

               DO IZC = 1, NZC
                  DO IXC = 1, NXC

                     IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. C%DISCRET(IXC,1,IZC,CGSC)/=IS_GASPHASE) CYCLE

                     IXF = 2*IXC
                     IZF = 2*IZC

                     ICC = C%DISCRET(IXC, 1, IZC, UNKH) 

                     ICF(1) = F%DISCRET(IXF-1, 1, IZF-1, UNKH) 
                     ICF(2) = F%DISCRET(IXF-1, 1, IZF  , UNKH) 
                     ICF(3) = F%DISCRET(IXF  , 1, IZF-1, UNKH) 
                     ICF(4) = F%DISCRET(IXF  , 1, IZF  , UNKH) 

                     FCC(ICC) = 0.25_EB * (  DCF(ICF(1)) &
                                           + DCF(ICF(2)) &
                                           + DCF(ICF(3)) &
                                           + DCF(ICF(4)) )
                  ENDDO
               ENDDO

            CASE (NSCARC_INTERPOL_BILINEAR)
    
               FCC=0.0_EB
      
               DO IZC = 1, NZC
                  DO IXC = 1, NXC
                     
                     IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. C%DISCRET(IXC,1,IZC,CGSC)/=IS_GASPHASE) CYCLE
         
                     IXF = 2*IXC
                     IZF = 2*IZC
         
                     ICC = C%DISCRET(IXC, 1, IZC, UNKH) 
         
                     ICFB(-2,-2) = F%DISCRET(IXF-2, 1, IZF-2, UNKH) 
                     ICFB(-1,-2) = F%DISCRET(IXF-1, 1, IZF-2, UNKH) 
                     ICFB( 1,-2) = F%DISCRET(IXF  , 1, IZF-2, UNKH) 
                     ICFB( 2,-2) = F%DISCRET(IXF+1, 1, IZF-2, UNKH) 
      
                     ICFB(-2,-1) = F%DISCRET(IXF-2, 1, IZF-1, UNKH) 
                     ICFB(-1,-1) = F%DISCRET(IXF-1, 1, IZF-1, UNKH) 
                     ICFB( 1,-1) = F%DISCRET(IXF  , 1, IZF-1, UNKH) 
                     ICFB( 2,-1) = F%DISCRET(IXF+1, 1, IZF-1, UNKH) 
      
                     ICFB(-2, 1) = F%DISCRET(IXF-2, 1, IZF, UNKH) 
                     ICFB(-1, 1) = F%DISCRET(IXF-1, 1, IZF, UNKH) 
                     ICFB( 1, 1) = F%DISCRET(IXF  , 1, IZF, UNKH) 
                     ICFB( 2, 1) = F%DISCRET(IXF+1, 1, IZF, UNKH) 
      
                     ICFB(-2, 2) = F%DISCRET(IXF-2, 1, IZF+1, UNKH) 
                     ICFB(-1, 2) = F%DISCRET(IXF-1, 1, IZF+1, UNKH) 
                     ICFB( 1, 2) = F%DISCRET(IXF  , 1, IZF+1, UNKH) 
                     ICFB( 2, 2) = F%DISCRET(IXF+1, 1, IZF+1, UNKH) 
         
      !IF (TYPE_DEBUG > NSCARC_DEBUG_EXTREME) THEN
      !   WRITE(LU_SCARC,'(A,2i6)') 'RESTRICTION_NEW: ----------- ', IXC, IZC
      !   WRITE(LU_SCARC,'(A,4i6)') ' 2:', ICFB(-2, 2), ICFB(-1, 2), ICFB(1, 2), ICFB(2, 2)
      !   WRITE(LU_SCARC,'(A,4i6)') ' 1:', ICFB(-2, 1), ICFB(-1, 1), ICFB(1, 1), ICFB(2, 1)
      !   WRITE(LU_SCARC,'(A,4i6)') '-1:', ICFB(-2,-1), ICFB(-1,-1), ICFB(1,-1), ICFB(2,-1)
      !   WRITE(LU_SCARC,'(A,4i6)') '-2:', ICFB(-2,-2), ICFB(-1,-2), ICFB(1,-2), ICFB(2,-2)
      !ENDIF
                     IF (IXC==1.AND.IZC==1) THEN
                        FCC(ICC) = SCAL*( &
                                     W4 *DCF(ICFB(-1, 2)) + W3 *DCF(ICFB(1, 2)) + W1*DCF(ICFB(2, 2)) + &
                                     W12*DCF(ICFB(-1, 1)) + W9 *DCF(ICFB(1, 1)) + W3*DCF(ICFB(2, 1)) + &
                                     W16*DCF(ICFB(-1,-1)) + W12*DCF(ICFB(1,-1)) + W4*DCF(ICFB(2,-1)) )
                     ELSE IF (IXC==NXC.AND.IZC==  1) THEN
                        FCC(ICC) = SCAL*( &
                                     W1 *DCF(ICFB(-2, 2)) + W3 *DCF(ICFB(-1, 2)) + W4 *DCF(ICFB(1, 2)) + &
                                     W3 *DCF(ICFB(-2, 1)) + W9 *DCF(ICFB(-1, 1)) + W12*DCF(ICFB(1, 1)) + &
                                     W4 *DCF(ICFB(-2,-1)) + W12*DCF(ICFB(-1,-1)) + W16*DCF(ICFB(1,-1)) )
                     ELSE IF (IXC==  1.AND.IZC==NZC) THEN
                        FCC(ICC) = SCAL*( &
                                     W16*DCF(ICFB(-1, 1)) + W12*DCF(ICFB(1, 1)) + W4*DCF(ICFB(2, 1)) + &
                                     W12*DCF(ICFB(-1,-1)) + W9 *DCF(ICFB(1,-1)) + W3*DCF(ICFB(2,-1)) + &
                                     W4 *DCF(ICFB(-1,-2)) + W3 *DCF(ICFB(1,-2)) + W1*DCF(ICFB(2,-2)) )
                     ELSE IF (IXC==NXC.AND.IZC==NZC) THEN
                        FCC(ICC) = SCAL*( &
                                     W4 *DCF(ICFB(-2, 1)) + W12*DCF(ICFB(-1, 1)) + W16*DCF(ICFB(1, 1)) + &
                                     W3 *DCF(ICFB(-2,-1)) + W9 *DCF(ICFB(-1,-1)) + W12*DCF(ICFB(1,-1)) + &
                                     W1 *DCF(ICFB(-2,-2)) + W3 *DCF(ICFB(-1,-2)) + W4 *DCF(ICFB(1,-2)) )
                     ELSE IF (IZC==  1) THEN
                        FCC(ICC) = SCAL*( &
                                     W1*DCF(ICFB(-2, 2)) + W3 *DCF(ICFB(-1, 2)) + W3 *DCF(ICFB(1, 2)) + W1*DCF(ICFB(2, 2)) + &
                                     W3*DCF(ICFB(-2, 1)) + W9 *DCF(ICFB(-1, 1)) + W9 *DCF(ICFB(1, 1)) + W3*DCF(ICFB(2, 1)) + &
                                     W4*DCF(ICFB(-2,-1)) + W12*DCF(ICFB(-1,-1)) + W12*DCF(ICFB(1,-1)) + W4*DCF(ICFB(2,-1)) )
                     ELSE IF (IZC==NZC) THEN
                        FCC(ICC) = SCAL*( &
                                     W4*DCF(ICFB(-2, 1)) + W12*DCF(ICFB(-1, 1)) + W12*DCF(ICFB(1, 1)) + W4*DCF(ICFB(2, 1)) + &
                                     W3*DCF(ICFB(-2,-1)) + W9 *DCF(ICFB(-1,-1)) + W9 *DCF(ICFB(1,-1)) + W3*DCF(ICFB(2,-1)) + &
                                     W1*DCF(ICFB(-2,-2)) + W3 *DCF(ICFB(-1,-2)) + W3 *DCF(ICFB(1,-2)) + W1*DCF(ICFB(2,-2)) )
                     ELSE IF (IXC==  1) THEN
                        FCC(ICC) = SCAL*( &
                                     W4 *DCF(ICFB(-1, 2)) + W3*DCF(ICFB(1, 2)) + W1*DCF(ICFB(2, 2)) +&
                                     W12*DCF(ICFB(-1, 1)) + W9*DCF(ICFB(1, 1)) + W3*DCF(ICFB(2, 1)) +&
                                     W12*DCF(ICFB(-1,-1)) + W9*DCF(ICFB(1,-1)) + W3*DCF(ICFB(2,-1)) +&
                                     W4 *DCF(ICFB(-1,-2)) + W3*DCF(ICFB(1,-2)) + W1*DCF(ICFB(2,-2)) )
                     ELSE IF (IXC==NXC) THEN
                        FCC(ICC) = SCAL*( &
                                     W1*DCF(ICFB(-2, 2)) + W3*DCF(ICFB(-1, 2)) + W4 *DCF(ICFB(1, 2)) + &
                                     W3*DCF(ICFB(-2, 1)) + W9*DCF(ICFB(-1, 1)) + W12*DCF(ICFB(1, 1)) +&
                                     W3*DCF(ICFB(-2,-1)) + W9*DCF(ICFB(-1,-1)) + W12*DCF(ICFB(1,-1)) +&
                                     W1*DCF(ICFB(-2,-2)) + W3*DCF(ICFB(-1,-2)) + W4 *DCF(ICFB(1,-2)) )
                     ELSE 
                        FCC(ICC) = SCAL*( &
                                     W1*DCF(ICFB(-2,-2)) + W3*DCF(ICFB(-1,-2)) + W3*DCF(ICFB(1,-2)) + W1*DCF(ICFB(2,-2)) +&
                                     W3*DCF(ICFB(-2,-1)) + W9*DCF(ICFB(-1,-1)) + W9*DCF(ICFB(1,-1)) + W3*DCF(ICFB(2,-1)) +&
                                     W3*DCF(ICFB(-2, 1)) + W9*DCF(ICFB(-1, 1)) + W9*DCF(ICFB(1, 1)) + W3*DCF(ICFB(2, 1)) +&
                                     W1*DCF(ICFB(-2, 2)) + W3*DCF(ICFB(-1, 2)) + W3*DCF(ICFB(1, 2)) + W1*DCF(ICFB(2, 2)) )
                     ENDIF
                  ENDDO
               ENDDO

         END SELECT SELECT_INTERPOL
      
      ELSE

         DO IZC = 1, NZC
            DO IYC = 1, NYC
               DO IXC = 1, NXC

                  IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. C%DISCRET(IXC,IYC,IZC,CGSC)/=IS_GASPHASE) CYCLE

                  IXF = 2*IXC
                  IYF = 2*IYC
                  IZF = 2*IZC

                  ICC = C%DISCRET(IXC, IYC, IZC, UNKH) 

                  ICF(1) = F%DISCRET(IXF-1, IYF-1, IZF-1, UNKH) 
                  ICF(2) = F%DISCRET(IXF-1, IYF-1, IZF  , UNKH) 
                  ICF(3) = F%DISCRET(IXF-1, IYF  , IZF-1, UNKH) 
                  ICF(4) = F%DISCRET(IXF-1, IYF  , IZF  , UNKH) 
                  ICF(5) = F%DISCRET(IXF  , IYF-1, IZF-1, UNKH) 
                  ICF(6) = F%DISCRET(IXF  , IYF-1, IZF  , UNKH) 
                  ICF(7) = F%DISCRET(IXF  , IYF  , IZF-1, UNKH) 
                  ICF(8) = F%DISCRET(IXF  , IYF  , IZF  , UNKH) 

                  FCC(ICC) = 0.125_EB * (  DCF(ICF(1)) &
                                         + DCF(ICF(2)) &
                                         + DCF(ICF(3)) &
                                         + DCF(ICF(4)) &
                                         + DCF(ICF(5)) &
                                         + DCF(ICF(6)) &
                                         + DCF(ICF(7)) &
                                         + DCF(ICF(8)) )

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

      C => SCARC(NM)%LEVEL(NLC)           !> pointer to coarse grid
      F => SCARC(NM)%LEVEL(NLF)           !> pointer to fine grid

      CALL POINT_TO_VECTOR(NVECTORF, NM, NLF, DCF)
      CALL POINT_TO_VECTOR(NVECTORC, NM, NLC, FCC)

      NCC => C%NCS

      R   => F%R%VAL
      ROW => F%R%ROW
      COL => F%R%COL

      DO ICC = 1, NCC
         AUX = 0.0_EB
         DO ICOL = ROW(ICC), ROW(ICC+1)-1
            IC = COL(ICOL)
            AUX0 = AUX
            AUX = AUX + DCF(IC) * R(ICOL)
         ENDDO
         FCC(ICC) = AUX
      ENDDO
   ENDDO

ENDIF

call scarc_debug_level (NVECTORC, 'NVECTORC AFTER', NLC)
END SUBROUTINE SCARC_RESTRICTION


!> ------------------------------------------------------------------------------------------------
!> Perform prolongation from coarser to finer grid 
!>    - 'CO' corresponds to coarser grid
!>    - 'FI' corresponds to finer   grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PROLONGATION (NVECTORC, NVECTORF, NLC, NLF)
INTEGER, INTENT(IN) :: NVECTORC, NVECTORF, NLC, NLF
INTEGER :: NM, ICOL, IC, I
REAL(EB), POINTER, DIMENSION(:)     :: XCC, DCF, P
INTEGER , POINTER, DIMENSION(:)     :: ROW, COL
INTEGER , POINTER :: NXC, NYC, NZC, NCF, NCEF
INTEGER  :: NXF, NYF, NZF
INTEGER  :: IXF, IYF, IZF, ICF(8)=0, ICFB(-1:1,-1:1)=0
INTEGER  :: IXC, IYC, IZC, ICC
REAL(EB) :: AUX, SCAL, W1, W3, W4, W9, W12
TYPE (SCARC_LEVEL_TYPE), POINTER :: C, F

SCARC_ROUTINE = 'SCARC_PROLONGATION'

SCAL = 0.0625_EB
W12  = 12.0_EB
W9   =  9.0_EB
W4   =  4.0_EB
W3   =  3.0_EB
W1   =  1.0_EB

!>
!> ------------------ Twolevel CG or Geometric Multigrid -------------------------------------------
!>
IF (BMULTILEVEL) THEN

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         C => SCARC(NM)%LEVEL(NLC)
         F => SCARC(NM)%LEVEL(NLF)

         NXC => C%NX
         NYC => C%NY
         NZC => C%NZ

         NXF = 2*NXC
         NYF = 2*NYC
         NZF = 2*NZC

         CALL POINT_TO_VECTOR(NVECTORC, NM, NLC, XCC)
         CALL POINT_TO_VECTOR(NVECTORF, NM, NLF, DCF)

         IF (TWO_D) THEN

            SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)

               CASE (NSCARC_INTERPOL_CONSTANT)

                  DO IZC = 1, NZC
                     DO IXC = 1, NXC
      
                        IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.C%DISCRET(IXC,1,IZC,CGSC)/=IS_GASPHASE) CYCLE
      
                        IXF = 2*IXC
                        IYF = 1
                        IZF = 2*IZC
      
                        ICC = C%DISCRET(IXC, 1, IZC, UNKH) 
      
                        ICF(1) = F%DISCRET(IXF-1, 1, IZF-1, UNKH) 
                        ICF(2) = F%DISCRET(IXF-1, 1, IZF  , UNKH) 
                        ICF(3) = F%DISCRET(IXF  , 1, IZF-1, UNKH) 
                        ICF(4) = F%DISCRET(IXF  , 1, IZF  , UNKH) 
      
                        DO I = 1, 4
                           DCF(ICF(I)) = XCC(ICC)
                        ENDDO
                     ENDDO
                  ENDDO

               CASE (NSCARC_INTERPOL_BILINEAR)

                  DO IZC = 1, NZC
                     DO IXC = 1, NXC
         
                        IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.C%DISCRET(IXC,1,IZC,CGSC)/=IS_GASPHASE) CYCLE
         
                        IXF = 2*IXC
                        IZF = 2*IZC
         
                        ICC = C%DISCRET(IXC, 1, IZC, UNKH) 
         
                        ICFB(-1,-1) = F%DISCRET(IXF-1, 1, IZF-1, UNKH) 
                        ICFB(-1, 1) = F%DISCRET(IXF-1, 1, IZF  , UNKH) 
                        ICFB( 1,-1) = F%DISCRET(IXF  , 1, IZF-1, UNKH) 
                        ICFB( 1, 1) = F%DISCRET(IXF  , 1, IZF  , UNKH) 
         
         !IF (TYPE_DEBUG > NSCARC_DEBUG_EXTREME) THEN
         !   WRITE(LU_SCARC,'(A,2i6,A,4i6)') 'PROLONGATION_NEW:', IXC, IZC, ':', ICFB(-1,-1), ICFB(-1,1), ICFB(1,-1), ICFB(1,1)
         !ENDIF
                        IF (IXC==1.AND.IZC==1) THEN
                           DCF(ICFB(-1,-1)) = XCC(ICC)
                           DCF(ICFB(-1, 1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC+NXC))
                           DCF(ICFB( 1,-1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC+1))
                           DCF(ICFB( 1, 1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC+1)+W3*XCC(ICC+NXC)+W1*XCC(ICC+NXC+1))
                        ELSE IF (IXC==1 .AND. IZC==NZC) THEN
                           DCF(ICFB(-1,-1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC-NXC))
                           DCF(ICFB(-1, 1)) = XCC(ICC)
                           DCF(ICFB( 1,-1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC+1)+W3*XCC(ICC-NXC)+W1*XCC(ICC-NXC+1))
                           DCF(ICFB( 1, 1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC+1))
                        ELSE IF (IXC==NXC .AND. IZC==1) THEN
                           DCF(ICFB(-1,-1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC-1))
                           DCF(ICFB(-1, 1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC-1)+W3*XCC(ICC+NXC)+W1*XCC(ICC+NXC-1))
                           DCF(ICFB( 1,-1)) = XCC(ICC)
                           DCF(ICFB( 1, 1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC+NXC))
                        ELSE IF (IXC==NXC .AND. IZC==NZC) THEN
                           DCF(ICFB(-1,-1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC-1)+W3*XCC(ICC-NXC)+W1*XCC(ICC-NXC-1))
                           DCF(ICFB(-1, 1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC-1))
                           DCF(ICFB( 1,-1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC-NXC))
                           DCF(ICFB( 1, 1)) = XCC(ICC)
                        ELSE IF (IZC==1) THEN
                           DCF(ICFB(-1,-1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC-1))
                           DCF(ICFB(-1, 1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC-1)+W3*XCC(ICC+NXC)+W1*XCC(ICC+NXC-1))
                           DCF(ICFB( 1,-1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC+1))
                           DCF(ICFB( 1, 1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC+1)+W3*XCC(ICC+NXC)+W1*XCC(ICC+NXC+1))
                        ELSE IF (IZC==NZC) THEN
                           DCF(ICFB(-1,-1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC-1)+W3*XCC(ICC-NXC)+W1*XCC(ICC-NXC-1))
                           DCF(ICFB(-1, 1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC-1))
                           DCF(ICFB( 1,-1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC+1)+W3*XCC(ICC-NXC)+W1*XCC(ICC-NXC+1))
                           DCF(ICFB( 1, 1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC+1))
                        ELSE IF (IXC==1) THEN
                           DCF(ICFB(-1,-1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC-NXC))
                           DCF(ICFB(-1, 1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC+NXC))
                           DCF(ICFB( 1,-1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC+1)+W3*XCC(ICC-NXC)+W1*XCC(ICC-NXC+1))
                           DCF(ICFB( 1, 1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC+1)+W3*XCC(ICC+NXC)+W1*XCC(ICC+NXC+1))
                        ELSE IF (IXC==NXC) THEN
                           DCF(ICFB(-1,-1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC-1)+W3*XCC(ICC-NXC)+W1*XCC(ICC-NXC-1))
                           DCF(ICFB(-1, 1)) = SCAL*(W9 *XCC(ICC)+W3*XCC(ICC-1)+W3*XCC(ICC+NXC)+W1*XCC(ICC+NXC-1))
                           DCF(ICFB( 1,-1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC-NXC))
                           DCF(ICFB( 1, 1)) = SCAL*(W12*XCC(ICC)+W4*XCC(ICC+NXC))
                        ELSE
                           DCF(ICFB(-1,-1)) = SCAL*(W9*XCC(ICC)+W3*XCC(ICC-1)+W3*XCC(ICC-NXC)+W1*XCC(ICC-NXC-1))
                           DCF(ICFB(-1, 1)) = SCAL*(W9*XCC(ICC)+W3*XCC(ICC-1)+W3*XCC(ICC+NXC)+W1*XCC(ICC+NXC-1))
                           DCF(ICFB( 1,-1)) = SCAL*(W9*XCC(ICC)+W3*XCC(ICC+1)+W3*XCC(ICC-NXC)+W1*XCC(ICC-NXC+1))
                           DCF(ICFB( 1, 1)) = SCAL*(W9*XCC(ICC)+W3*XCC(ICC+1)+W3*XCC(ICC+NXC)+W1*XCC(ICC+NXC+1))
                        ENDIF
                     ENDDO
                  ENDDO

            END SELECT SELECT_INTERPOL

         ELSE

            ! Note: 3D-bilinear case is still missing
            DO IZC = 1, NZC
               DO IYC = 1, NYC
                  DO IXC = 1, NXC

                     IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.C%DISCRET(IXC,IYC,IZC,CGSC)/=IS_GASPHASE) CYCLE

                     IXF = 2*IXC
                     IYF = 2*IYC
                     IZF = 2*IZC

                     ICC = C%DISCRET(IXC, IYC, IZC, UNKH) 

                     ICF(1) = F%DISCRET(IXF-1, IYF-1, IZF-1, UNKH) 
                     ICF(2) = F%DISCRET(IXF-1, IYF-1, IZF  , UNKH) 
                     ICF(3) = F%DISCRET(IXF-1, IYF  , IZF-1, UNKH) 
                     ICF(4) = F%DISCRET(IXF-1, IYF  , IZF  , UNKH) 
                     ICF(5) = F%DISCRET(IXF  , IYF-1, IZF-1, UNKH) 
                     ICF(6) = F%DISCRET(IXF  , IYF-1, IZF  , UNKH) 
                     ICF(7) = F%DISCRET(IXF  , IYF  , IZF-1, UNKH) 
                     ICF(8) = F%DISCRET(IXF  , IYF  , IZF  , UNKH) 

                     DO I = 1, 8
                        DCF(ICF(I)) = XCC(ICC)
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
   IF (TYPE_COARSENING >= NSCARC_COARSENING_GMG) SCAL=2.0_EB

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      C => SCARC(NM)%LEVEL(NLC)
      F => SCARC(NM)%LEVEL(NLF)

      CALL POINT_TO_VECTOR(NVECTORC, NM, NLC, XCC)
      CALL POINT_TO_VECTOR(NVECTORF, NM, NLF, DCF)

      NCF  => F%NCS
      NCEF => F%NCE

      P   => F%P%VAL
      ROW => F%P%ROW
      COL => F%P%COL

      DO IC = 1, NCF
         AUX = 0.0_EB
         DO ICOL = ROW(IC), ROW(IC+1)-1
            ICC = COL(ICOL)
            AUX = XCC(ICC) * P(ICOL)
         ENDDO
         DCF(IC) = SCAL*AUX
      ENDDO
      DO IC = NCF+1, NCEF
         DCF(IC) = 0.0_EB
      ENDDO
   ENDDO

ENDIF

END SUBROUTINE SCARC_PROLONGATION

!> ------------------------------------------------------------------------------------------------
!> Copy final solution from GMG (as preconditioner) to corresponding vector of CG (as main solver)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRECONDITIONER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L => SCARC(NM)%LEVEL(NL)
   L%MAIN%G = L%PRECON%X
ENDDO
END SUBROUTINE SCARC_UPDATE_PRECONDITIONER

!> ------------------------------------------------------------------------------------------------
!> Finalize data
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRESSURE_MAINCELLS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, I, J, K, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   M => MESHES(NM)
   L => SCARC(NM)%LEVEL(NL)
   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF
   DO K = 1, M%KBAR
      DO J = 1, M%JBAR
         DO I = 1, M%IBAR
            IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. L%DISCRET(I,J,K,CGSC) /= IS_GASPHASE) CYCLE
            IC = L%DISCRET(I,J,K,UNKH)
            HP(I, J, K) = L%MAIN%X(IC)
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


DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   L => SCARC(NM)%LEVEL(NL)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

 !> compute ghost cell values
   WALL_CELLS_LOOP: DO IW = 1, L%N_EXTERNAL_WALL_CELLS

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
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) - DXI *M%BXS(IYW,IZW)
            ENDIF
         CASE (-1)
            IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXG,IYW,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXF(IYW,IZW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) + DXI *M%BXF(IYW,IZW)
            ENDIF
         CASE ( 2)
            IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYS(IXW,IZW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) - DETA *M%BYS(IXW,IZW)
            ENDIF
         CASE (-2)
            IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYF(IXW,IZW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) + DETA *M%BYF(IXW,IZW)
            ENDIF
         CASE ( 3)
            IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZS(IXW,IYW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) - DZETA *M%BZS(IXW,IYW)
            ENDIF
         CASE (-3)
         IF (L%WALL(IW)%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZF(IXW,IYW)
            ELSE IF (L%WALL(IW)%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) + DZETA *M%BZF(IXW,IYW)
            ENDIF
      END SELECT
   ENDDO WALL_CELLS_LOOP

ENDDO


!> -----------------------------------------------------------------------------------------------
!> Perform data exchange to achieve consistency of ghost values along internal boundaries
!> -----------------------------------------------------------------------------------------------
!IF (N_MPI_PROCESSES > 1) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_PRESSURE, NL)
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
TYPE (SCARC_LEVEL_TYPE), POINTER ::  OL

SCARC_ROUTINE = 'SCARC_EXCHANGE_RECEIVE'

RECEIVE_MESH_INDEX: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   RECEIVE_OMESH_INDEX: DO NOM = 1, NMESHES

      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_OMESH_INDEX

      S   => SCARC(NM)
      OS => SCARC(NM)%OSCARC(NOM)

      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE RECEIVE_OMESH_INDEX

      SELECT_EXCHANGE_TYPE: SELECT CASE (TYPE_EXCHANGE)


         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring step size along internal boundary
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT_BASIC(1),1,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring wall data
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALLINFO)

            OL => OS%LEVEL(NL)
            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT(1),SIZE(OS%RECV_INT),MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring wall data
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_DISCRET)

            OL => OS%LEVEL(NL)
            N_REQ = N_REQ+1
            OS%RECV_INT = 0
            CALL MPI_IRECV(OS%RECV_INT(1),SIZE(OS%RECV_INT),MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring grid dimensions
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MESHINFO)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT_BASIC(1),8,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring step size along internal boundary
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WIDTHINFO)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_REAL_BASIC(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange number of neighboring cells for AMG method (compact type only)
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_TRANSFER_SIZE)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT_BASIC(1),9,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange number of neighboring cells for AMG method (compact type only)
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT_BASIC(1),1,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange neighboring CELL_TYPEs
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_TYPE)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT(1),SIZE(OS%RECV_INT),MPI_INTEGER,&
                           SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)


         !> ---------------------------------------------------------------------------------------
         !> Perform exchanges for
         !>    - internal values for matrix-vector multiplication
         !>    - internal boundariy values
         !>    - internal subdiagonal matrix values
         !>    - internal subdiagonal or ghost matrix values
         !>    - internal measure/CELL_TYPE values
         !> ---------------------------------------------------------------------------------------
         CASE DEFAULT

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_REAL(1),SIZE(OS%RECV_REAL),MPI_DOUBLE_PRECISION,&
                           SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

      END SELECT SELECT_EXCHANGE_TYPE
   ENDDO RECEIVE_OMESH_INDEX
ENDDO RECEIVE_MESH_INDEX

END SUBROUTINE SCARC_EXCHANGE_RECEIVE


!> ------------------------------------------------------------------------------------------------
!> Send data to neighbors (corresponds to MESH_EXCHANGE)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_SEND (NL)
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, NOM, NCOL, NROW
INTEGER  :: IW, IWL, IWG, IWW, IW0, IROW, ICOL, IPTR, ICPL, ICELL_TYPE
INTEGER  :: IOR0, IC, JC, ICC, JCC, ICE, ICO, ICN, ICG, ICW
INTEGER  :: IX, IY, IZ
INTEGER  :: IXW, IYW, IZW
INTEGER  :: IXG, IYG, IZG
INTEGER  :: I, J, K, LL
REAL(EB) :: ZSUM, ZSUM1
INTEGER , POINTER, DIMENSION(:)     ::  RECV_INT, RECV_INT_BASIC
REAL(EB), POINTER, DIMENSION(:)     ::  RECV_REAL, RECV_REAL_BASIC
REAL(EB), POINTER, DIMENSION(:)     ::  VECTOR
REAL(EB), POINTER, DIMENSION(:,:,:) :: HVECTOR
TYPE (SCARC_TYPE)      ,  POINTER :: S
TYPE (OSCARC_TYPE)     ,  POINTER :: OS
TYPE (SCARC_LEVEL_TYPE),  POINTER :: L, OL
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, P, R, OP
TYPE (SCARC_MAPPING_TYPE), POINTER :: MAP, OMAP
TYPE (SCARC_AMG_TYPE),    POINTER :: AMG, OAMG

SCARC_ROUTINE = 'SCARC_EXCHANGE_SEND'

!> ------------------------------------------------------------------------------------------------
!> Collect data for sending corresponding to requested exchange type
!> ------------------------------------------------------------------------------------------------
MESH_PACK_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S  => SCARC(NM)
   L  => SCARC(NM)%LEVEL(NL)
   MAP => SCARC(NM)%LEVEL(NL)%MAP

   OMESH_PACK_LOOP: DO NOM = 1, NMESHES

      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      OS => SCARC(NM)%OSCARC(NOM)

      IF (OS%NICMAX_S == 0 .AND. OS%NICMAX_R == 0) CYCLE OMESH_PACK_LOOP

      OL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
      OMAP => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%MAP

      OMESH_PACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_BASIC: pack neighboring sizes for exchange
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC)

            OS%SEND_INT_BASIC(1)=OS%LEVEL(NL)%NWL

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT_BASIC(1),1,MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_WIDTHINFO: pack neighboring grid resolution information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WIDTHINFO)

            SELECT CASE(OL%IOR)
               CASE (1)
                  OL%DH = L%DXL(0)
               CASE (-1)
                  OL%DH = L%DXL(L%NX)
                CASE (2)
                  OL%DH = L%DYL(0)
               CASE (-2)
                  OL%DH = L%DYL(L%NY)
               CASE (3)
                  OL%DH = L%DZL(0)
               CASE (-3)
                  OL%DH = L%DZL(L%NZ)
            END SELECT
            OS%SEND_REAL_BASIC(1) = OL%DH

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL_BASIC(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_WALLINFO: pack neighboring wallinfo information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALLINFO)

            IPTR=1
            DO IWL = 1, OL%NWL
               IWG = OMAP%IWL_TO_IWG(IWL)
               OS%SEND_INT(IPTR   ) = L%WALL(IWG)%IXG
               OS%SEND_INT(IPTR+ 1) = L%WALL(IWG)%IYG
               OS%SEND_INT(IPTR+ 2) = L%WALL(IWG)%IZG
               OS%SEND_INT(IPTR+ 3) = L%WALL(IWG)%IXW
               OS%SEND_INT(IPTR+ 4) = L%WALL(IWG)%IYW
               OS%SEND_INT(IPTR+ 5) = L%WALL(IWG)%IZW
               OS%SEND_INT(IPTR+ 6) = L%WALL(IWG)%IXN(1)
               OS%SEND_INT(IPTR+ 7) = L%WALL(IWG)%IXN(2)
               OS%SEND_INT(IPTR+ 8) = L%WALL(IWG)%IYN(1)
               OS%SEND_INT(IPTR+ 9) = L%WALL(IWG)%IYN(2)
               OS%SEND_INT(IPTR+10) = L%WALL(IWG)%IZN(1)
               OS%SEND_INT(IPTR+11) = L%WALL(IWG)%IZN(2)
               OS%SEND_INT(IPTR+12) = L%WALL(IWG)%NOM
               IPTR = IPTR + 13
               DO ICPL=1,OL%NCPLS
                  OS%SEND_INT(IPTR)=L%WALL(IWG)%ICE(ICPL)
                  IPTR = IPTR + 1
               ENDDO
            ENDDO

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT(1),SIZE(OS%SEND_INT),MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_DISCRET: pack neighboring DISCRET information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_DISCRET)

            IPTR=1
            OS%SEND_INT=0
            DO IWL = 1, OL%NWL
               IWG = OMAP%IWL_TO_IWG(IWL)
               OS%SEND_INT(IPTR  ) = L%DISCRET(L%WALL(IWG)%IXW,L%WALL(IWG)%IYW,L%WALL(IWG)%IZW,CGSC)
               OS%SEND_INT(IPTR+1) = L%DISCRET(L%WALL(IWG)%IXW,L%WALL(IWG)%IYW,L%WALL(IWG)%IZW,UNKH)
               IPTR = IPTR + 2
            ENDDO

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT(1),SIZE(OS%SEND_INT),MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MESHINFO: pack neighboring mesh information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MESHINFO)

            L => S%LEVEL(NL)

            OS%SEND_INT_BASIC(1)=L%NX
            OS%SEND_INT_BASIC(2)=L%NY
            OS%SEND_INT_BASIC(3)=L%NZ
            OS%SEND_INT_BASIC(4)=L%NC
            OS%SEND_INT_BASIC(5)=L%NCS
            OS%SEND_INT_BASIC(6)=L%NW
            OS%SEND_INT_BASIC(7)=L%N_EXTERNAL_WALL_CELLS
            OS%SEND_INT_BASIC(8)=L%N_INTERNAL_WALL_CELLS

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT_BASIC(1),8,MPI_INTEGER,SNODE, &
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
               IWG = OMAP%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXW
               IY  = L%WALL(IWG)%IYW
               IZ  = L%WALL(IWG)%IZW
               OS%SEND_REAL(LL) = HVECTOR(IX,IY,IZ)
               LL = LL+1
            ENDDO

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_VECTOR: pack overlapping parts of a given vector
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_VECTOR)

            CALL POINT_TO_VECTOR(TYPE_VECTOR, NM, NL, VECTOR)

            LL = 1
            DO ICG= 1, OL%NCG
               !ZSUM = 0.0_EB
               IWG = OMAP%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXW
               IY  = L%WALL(IWG)%IYW
               IZ  = L%WALL(IWG)%IZW

               !ZSUM = VECTOR(L%DISCRET(IX, IY, IZ, UNKH))
               !OS%SEND_REAL(LL) = ZSUM/REAL(OL%NCPLR,EB)

               IF (PRES_ON_WHOLE_DOMAIN.OR.L%DISCRET(IX, IY, IZ, CGSC) == IS_GASPHASE) THEN
                  OS%SEND_REAL(LL) = VECTOR(L%DISCRET(IX, IY, IZ, UNKH))
               ELSE
                  OS%SEND_REAL(LL) = -987654321.0_EB
               ENDIF
               LL = LL + 1
            ENDDO 

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MEASURE_ADD: pack neighboring measure information (AMG only) with adding
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MEASURE_ADD)

            LL  = 1
            PACK_MEASURE_ADD: DO IWL=1, OL%NWL
               ICO = OMAP%IWL_TO_ICO(IWL)
               OS%SEND_REAL(LL) = L%AMG%MEASURE(ICO)
               LL = LL+1
            ENDDO PACK_MEASURE_ADD

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> Send data along internal boundaries corresponding to requested exchange type
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_TYPE)

            LL  = 1
            PACK_CELL_TYPE: DO IWL=1, OL%NWL
               ICW = OMAP%IWL_TO_ICW(IWL)
               OS%SEND_INT(LL) = L%AMG%CELL_TYPE(ICW)
               LL = LL+1
            ENDDO PACK_CELL_TYPE

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT(1), SIZE(OS%SEND_INT), MPI_INTEGER, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_SIZE: Send sizes of overlapping matrix areas
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            OL => OS%LEVEL(NL)
            OS%SEND_INT_BASIC(1) = OL%A%NAV

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT_BASIC(1),1,MPI_INTEGER,SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_SUBDIAG: pack subdiagonal entries of system matrix on overlapping parts
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG)

            LL = 1
            IWW = 0

            WRITE(*,*) '!>!>! EXCHANGE_MATRIX_SUBDIAG, PACK: Achtung, hier nochmal checken  AMG !>!>'

            DO ICG=1,OL%NCG
               IOR0 = ABS(OL%WALL(ICG)%IOR)
               IWG = OMAP%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXG
               IY  = L%WALL(IWG)%IYG
               IZ  = L%WALL(IWG)%IZG
            ENDDO 

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_STENCIL: pack matrix stencil information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_STENCIL)

            LL = 1
            OS%SEND_REAL = 0.0_EB
            A => SCARC(NM)%LEVEL(NL)%A

            DO ICG=1, OL%NCG
  
               WRITE(*,*) '!>!>! ACHTUNG: EXCHANGE_MATRIX_STENCIL, PACK:  nochmal checken  AMG !>!>'
 
               IWG = OMAP%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXG
               IY  = L%WALL(IWG)%IYG
               IZ  = L%WALL(IWG)%IZG
               DO ICPL = 1, OL%NCPLR
                  IC = L%DISCRET(IX, IY, IZ, UNKH)
                  NCOL = 0
                  DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
                     JC = A%COL(ICOL)
                     IF (JC > L%NC) THEN
                        IW0 = MAP%ICE_TO_IWG(JC)
                        IF (L%WALL(IW0)%NOM == NOM) NCOL = NCOL + 1
                     ELSE
                        NCOL = NCOL + 1
                     ENDIF
                  ENDDO
                  OS%SEND_REAL(LL) = REAL(NCOL)
                  LL = LL + 1
               ENDDO
            ENDDO 

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_SYSTEM: pack overlapping parts of system matrix
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SYSTEM)

            LL = 1
            OS%SEND_REAL = 0.0_EB
            A => SCARC(NM)%LEVEL(NL)%A

            WRITE(*,*) '!>!>! EXCHANGE_MATRIX_SYSTEM, PACK: Achtung, hier nochmal checken  AMG !>!>'
            !> Pack first cell layer
            PACK_MATRIX_SYSTEM: DO ICG=1, OL%NCG
               IF (OL%WALL(ICG)%NOM/=NM) CYCLE PACK_MATRIX_SYSTEM
               IWG = OMAP%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXG
               IY  = L%WALL(IWG)%IYG
               IZ  = L%WALL(IWG)%IZG

               DO ICPL = 1, OL%NCPLR
                  IC = L%DISCRET(IX, IY, IZ, UNKH)

                  MATRIX_COLUMN_LOOP: DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
                     JC = A%COL(ICOL)
                     IF (JC > L%NC) THEN
                        IWL = MAP%ICE_TO_IWL(JC)
                        IF (OL%WALL(IWL)%NOM /= NM) CYCLE MATRIX_COLUMN_LOOP
                        OS%SEND_REAL(LL) = - REAL(MAP%ICE_TO_IWG(JC),EB)
                        OS%SEND_REAL(LL) = - REAL(JC)
                     ELSE
                        OS%SEND_REAL(LL) =   REAL(JC,EB)
                     ENDIF
                     OS%SEND_REAL(LL+1) = A%VAL(ICOL)
                     LL = LL + 2
                  ENDDO MATRIX_COLUMN_LOOP

               ENDDO

            ENDDO PACK_MATRIX_SYSTEM

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_PROL: pack overlapping parts of prolongation matrix
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_PROL)

            LL = 1
            P   => SCARC(NM)%LEVEL(NL)%P
            AMG => SCARC(NM)%LEVEL(NL)%AMG

            WRITE(*,*) '!>!>! EXCHANGE_MATRIX_PROL, PACK: Achtung, hier nochmal checken  AMG !>!>'

            PACK_MATRIX_PROL: DO ICG=1, OL%NCG
               IF (OL%WALL(ICG)%NOM/=NM) CYCLE PACK_MATRIX_PROL
               IWG = OMAP%ICG_TO_IWG(ICG)
               IX  = L%WALL(IWG)%IXG
               IY  = L%WALL(IWG)%IYG
               IZ  = L%WALL(IWG)%IZG

               DO ICPL = 1, OL%NCPLR
                  IC = L%DISCRET(IX, IY, IZ, UNKH)

                  OS%SEND_REAL(LL)   = REAL(AMG%CELL_TYPE(IC),EB)
                  OS%SEND_REAL(LL+1) = REAL(P%ROW(IC+1)-P%ROW(IC),EB)

                  LL = LL + 2

                  DO ICOL = P%ROW(IC), P%ROW(IC+1)-1
                     JC = P%COL(ICOL)
                     IF (JC >= L%NC) THEN
                        OS%SEND_REAL(LL) = - REAL(JC,EB)
                     ELSE
                        OS%SEND_REAL(LL) =   REAL(JC,EB)
                     ENDIF
                     OS%SEND_REAL(LL+1) = P%VAL(ICOL)
                     LL = LL + 2
                  ENDDO

               ENDDO

            ENDDO PACK_MATRIX_PROL

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_MATRIX_REST: pack overlapping parts of restriction matrix
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_REST)

            WRITE(*,*) '!>!>! EXCHANGE_MATRIX_REST, PACK: Achtung, hier nochmal checken  AMG !>!>'
            LL = 1
            R   => SCARC(NM)%LEVEL(NL)%R
            AMG => SCARC(NM)%LEVEL(NL)%AMG

            PACK_MATRIX_REST: DO ICG=1, OL%NCG
               IF (OL%WALL(IWG)%NOM/=NM) CYCLE PACK_MATRIX_REST
               
                  WRITE(*,*) 'ACHTUNG: MATRIX_REST: Wegen DISCRET schauen!>'

                  IX = L%WALL(IWG)%IXG
                  IY = L%WALL(IWG)%IYG
                  IZ = L%WALL(IWG)%IZG
                  IC = L%DISCRET(IX, IY, IZ, UNKH)

                  ICC = AMG%CELL_TYPE(IC)
                  IF (ICC > 0) THEN
                     OS%SEND_REAL(LL)   = REAL(ICC)
                     OS%SEND_REAL(LL+1) = REAL(R%ROW(ICC+1)-R%ROW(ICC),EB)

                     LL = LL + 2

                     DO ICOL = R%ROW(ICC), R%ROW(ICC+1)-1
                        JCC = R%COL(ICOL)
                        IF (JCC >= AMG%NCC) THEN
                           OS%SEND_REAL(LL) = - REAL(JCC,EB)
                        ELSE
                           OS%SEND_REAL(LL) =   REAL(JCC,EB)
                        ENDIF
                        OS%SEND_REAL(LL+1) = R%VAL(ICOL)
                        LL = LL + 2
                     ENDDO

                  ENDIF

               !ENDDO REST_FIRST_LAYER_LOOP

            ENDDO PACK_MATRIX_REST

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> EXCHANGE_TRANSER_SIZE: pack sizes of transfer matrices (AMG only)
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_TRANSFER_SIZE)

            L    =>  S%LEVEL(NL)
            OL   => OS%LEVEL(NL)
            AMG  =>  S%LEVEL(NL)%AMG
            OAMG => OS%LEVEL(NL)%AMG

            OS%SEND_INT_BASIC(1) = L%NC
            OS%SEND_INT_BASIC(2) = L%NW
            OS%SEND_INT_BASIC(3) = L%NCE
            OS%SEND_INT_BASIC(4) = L%NW
            OS%SEND_INT_BASIC(5) = AMG%NCCI
            OS%SEND_INT_BASIC(6) = OAMG%NPS
            OS%SEND_INT_BASIC(7) = OAMG%NRS
            OS%SEND_INT_BASIC(8) = OAMG%NCCS
            OS%SEND_INT_BASIC(9) = OAMG%NCFS

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT_BASIC(1),9,MPI_INTEGER,SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
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
   OMESH_UNPACK_LOOP: DO NOM=1,NMESHES

      SNODE  = PROCESS(NM)
      RNODE  = PROCESS(NOM)

      OS => SCARC(NM)%OSCARC(NOM)

      OMESH_UNPACK_IF: IF (OS%NICMAX_S/=0 .AND. OS%NICMAX_R/=0) THEN

         L    => SCARC(NM)%LEVEL(NL)
         OL   => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
         OMAP => OL%MAP

         IF (RNODE/=SNODE) THEN
            RECV_INT_BASIC  => SCARC(NM)%OSCARC(NOM)%RECV_INT_BASIC
            RECV_REAL_BASIC => SCARC(NM)%OSCARC(NOM)%RECV_REAL_BASIC
            RECV_INT      => SCARC(NM)%OSCARC(NOM)%RECV_INT
            RECV_REAL     => SCARC(NM)%OSCARC(NOM)%RECV_REAL
         ELSE
            RECV_INT_BASIC  => SCARC(NOM)%OSCARC(NM)%SEND_INT_BASIC
            RECV_REAL_BASIC => SCARC(NOM)%OSCARC(NM)%SEND_REAL_BASIC
            RECV_INT      => SCARC(NOM)%OSCARC(NM)%SEND_INT
            RECV_REAL     => SCARC(NOM)%OSCARC(NM)%SEND_REAL
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

               OL%DH = RECV_REAL_BASIC(1)

               SELECT CASE (OL%IOR)
                  CASE ( 1)
                     L%DXL(0)    = 0.5_EB*(OL%DH + L%DXL(0))
                  CASE (-1)
                     L%DXL(L%NX) = 0.5_EB*(OL%DH + L%DXL(L%NX))
                  CASE ( 2)
                     L%DYL(0)    = 0.5_EB*(OL%DH + L%DYL(0))
                  CASE (-2)
                     L%DYL(L%NY) = 0.5_EB*(OL%DH + L%DYL(L%NY))
                  CASE ( 3)
                     L%DZL(0)    = 0.5_EB*(OL%DH + L%DZL(0))
                  CASE (-3)
                     L%DZL(L%NZ) = 0.5_EB*(OL%DH + L%DZL(L%NZ))
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
            !> EXCHANGE_DISCRET: unpack neighboring DISCRET information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_DISCRET)
               IPTR=1
               DO ICG = 1, OL%NCG

                  OL%WALL(ICG)%CGSC = RECV_INT(IPTR  )
                  OL%WALL(ICG)%UNKH = RECV_INT(IPTR+1)
                  IPTR = IPTR + 2

                  IWG = OMAP%ICG_TO_IWG(ICG)

                  IXG = L%WALL(IWG)%IXG
                  IYG = L%WALL(IWG)%IYG
                  IZG = L%WALL(IWG)%IZG

                  IXW= L%WALL(IWG)%IXW
                  IYW= L%WALL(IWG)%IYW
                  IZW= L%WALL(IWG)%IZW

                  !IF (OL%WALL(ICG)%CGSC == IS_GASPHASE.AND.L%DISCRET(IXW,IYW,IZW,CGSC)/=IS_SOLID) THEN
                  IF (OL%WALL(ICG)%CGSC == IS_GASPHASE) THEN
                     L%DISCRET(IXG, IYG, IZG, CGSC) = OL%WALL(ICG)%CGSC
                     L%DISCRET(IXG, IYG, IZG, UNKH) = OL%WALL(ICG)%UNKH
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
               OL%N_EXTERNAL_WALL_CELLS = RECV_INT_BASIC(7)
               OL%N_INTERNAL_WALL_CELLS = RECV_INT_BASIC(8)

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_WALLINFO: unpack neighboring matrix size information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

               OL%A%NAV = RECV_INT_BASIC(1)

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
               UNPACK_PRESSURE: DO IWL = 1, OL%NWL

                  ZSUM=0.0_EB
                  IWG = OMAP%IWL_TO_IWG(IWL)
                  DO ICPL = 1, OL%NCPL
                     ZSUM=ZSUM+RECV_REAL(LL)
                     LL = LL+1
                  ENDDO

                  I=L%WALL(IWG)%IXG
                  J=L%WALL(IWG)%IYG
                  K=L%WALL(IWG)%IZG

                  HVECTOR(I, J, K) = ZSUM/REAL(OL%NCPL,EB)

               ENDDO UNPACK_PRESSURE

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_VECTOR: unpack overlapping parts of a given vector
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_VECTOR)

               CALL POINT_TO_VECTOR (TYPE_VECTOR, NM, NL, VECTOR)

               LL = 1
               UNPACK_VECTOR: DO IWL = 1, OL%NWL
                  ZSUM = 0.0_EB
                  !ICO = OMAP%IWL_TO_ICO(IWL)

                  IWG = OMAP%IWL_TO_IWG(IWL)
                  IX = L%WALL(IWG)%IXG
                  IY = L%WALL(IWG)%IYG
                  IZ = L%WALL(IWG)%IZG

                  ICE = L%WALL(IWG)%ICE(1)

                  VECTOR(ICE) = RECV_REAL(LL)

                  !ICO = L%WALL(IWG)%ICO
                  !DO ICPL = 1, OL%NCPL
                  !   ZSUM = ZSUM + RECV_REAL(LL)
                  !   LL = LL + 1
                  !ENDDO
                  !VECTOR(ICO) = ZSUM/REAL(OL%NCPL,EB)

                  LL = LL + 1
               ENDDO UNPACK_VECTOR


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MEASURE_ADD: unpack neighboring measure information (AMG only) with adding
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MEASURE_ADD)

               LL = 1
               UNPACK_MEASURE_ADD: DO IWL = 1, OL%NWL
                  ICW = OMAP%IWL_TO_ICW(IWL)
                  L%AMG%MEASURE(ICW) = L%AMG%MEASURE(ICW)  + RECV_REAL(LL)
                  LL = LL + 1
               ENDDO UNPACK_MEASURE_ADD

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_CELL_TYPE: unpack neighboring CELL_TYPE information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_CELL_TYPE)

               LL = 1
               UNPACK_CELL_TYPE: DO IWL = 1, OL%NWL

                  ZSUM1=0.0_EB
                  DO ICPL = 1, OL%NCPL
                     ICG = OMAP%IWL_TO_ICG(IWL, ICPL)
                     OL%AMG%CELL_TYPE(ICG) = RECV_INT(LL)
                     LL = LL+1
                  ENDDO

               ENDDO UNPACK_CELL_TYPE


            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_SUBDIAG: unpack subdiagonal information from !neighboring system matrix
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG)

               LL = 1
               UNPACK_MATRIX_SUBDIAG: DO IWL = 1, OL%NWL

                  ZSUM=0.0_EB
                  ! ACHTUNG WEGEN NCPL SCHAUEN 
                  IWG = OMAP%IWL_TO_IWG(IWL)
                  IX = L%WALL(IWG)%IXG
                  IY = L%WALL(IWG)%IYG
                  IZ = L%WALL(IWG)%IZG
                  IC = L%DISCRET(IX, IY, IZ, UNKH)
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

               LL = 1
               UNPACK_MATRIX_STENCIL: DO IWL = 1, OL%NWL
                  !IWG = OMAP%IWL_TO_IWG(IWL)
                  !ICO = OMAP%IWL_TO_ICO(IWL)
                  !ICO = L%WALL(IWG)%ICO

                  OL%A%STENCIL(IWL) = NINT(RECV_REAL(LL))
                  LL = LL + 1

                ENDDO UNPACK_MATRIX_STENCIL

            !> ------------------------------------------------------------------------------------
            !> EXCHANGE_MATRIX_SYSTEM: extract neighboring system matrix information
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SYSTEM)

               LL = 1
               UNPACK_MATRIX_SYSTEM: DO IWL = 1, OL%NWL

                  IWG = OMAP%IWL_TO_IWG(IWL)

                  ICE = L%WALL(IWG)%ICE(1)
                  IX  = L%WALL(IWG)%IXG
                  IY  = L%WALL(IWG)%IYG
                  IZ  = L%WALL(IWG)%IZG
                  ICN = L%DISCRET(IX, IY, IZ, UNKH)
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

               AMG  => L%AMG
               OAMG => OL%AMG
               OP   => OL%P

               LL = 1
               NROW=1
               ICC = 1
               UNPACK_MATRIX_PROL: DO IWL = 1, OL%NWL

                  IWG = OMAP%IWL_TO_IWG(IWL)

                  ICG = L%WALL(IWG)%ICG(1)
                  ICE = L%WALL(IWG)%ICE(1)
                  IX  = L%WALL(IWG)%IXG
                  IY  = L%WALL(IWG)%IYG
                  IZ  = L%WALL(IWG)%IZG
                  ICN = L%DISCRET(IX, IY, IZ, UNKH)
                  ICELL_TYPE = NINT(RECV_REAL(LL))
                  NCOL       = NINT(RECV_REAL(LL+1))

                  OAMG%CELL_TYPE(ICG) = ICELL_TYPE
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
                 !> OL%AMG%CELL_TYPE(ICG) = ICELL_TYPE
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

IF (TYPE_DEBUG /= NSCARC_UNDEFINED ) CLOSE(LU_SCARC)

1001 FORMAT('| Scarc routine       |    Time (s)    |')
1002 FORMAT('|---------------------|----------------|')
2001 FORMAT('| Scarc routine       |  Time_min (s)  |  Time_max (s)  | Time_mean (s)  |')
2002 FORMAT('|---------------------|----------------|----------------|----------------|')
END SUBROUTINE SCARC_TIMINGS


!> -----------------------------------------------------------------------------
!> Setup DISCRET-array
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DISCRET

INTEGER :: NM
INTEGER :: I, J, K, IX, IY, IZ
INTEGER, PARAMETER :: IMPADD = 1
INTEGER, PARAMETER :: SHFTM(1:3,1:6) = RESHAPE((/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_LEVEL_TYPE), POINTER :: L


SCARC_ROUTINE = 'SCARC_SETUP_DISCRET'
!>
!> Initialize all cells as GASPHASE cells
!>
MESHES_LOOP1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   L => SCARC(NM)%LEVEL(NLEVEL_MIN)
   IF (.NOT. ALLOCATED(L%DISCRET)) ALLOCATE(L%DISCRET(0:L%NX+1,0:L%NY+1,0:L%NZ+1,NCVARS))
   L%DISCRET = 0
   L%DISCRET(:,:,:,CGSC) = IS_UNDEFINED
   L%DISCRET(1:L%NX, 1:L%NY, 1:L%NZ, CGSC) = IS_GASPHASE 
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
            IF (M%SOLID(M%CELL_INDEX(I, J, K))) L%DISCRET(I, J, K, CGSC) = IS_SOLID
         ENDDO
      ENDDO
   ENDDO
ENDDO MESHES_LOOP2

!>
!> Define local cell numbers for Poisson equation
!>
MESHES_LOOP3 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NLEVEL_MIN)
   L%DISCRET(:,:,:,UNKH) = IS_UNDEFINED

   !> consider gasphase and solid cells
   IF (PRES_ON_WHOLE_DOMAIN) THEN 
      DO K=1,L%NZ
         DO J=1,L%NY
            DO I=1,L%NX
               L%NUNKH_LOC(NM) = L%NUNKH_LOC(NM) + 1
               L%DISCRET(I,J,K,UNKH) = L%NUNKH_LOC(NM)
            ENDDO
         ENDDO
      ENDDO

   !> consider only gasphase cells
   ELSE 
      DO K=1,L%NZ
         DO J=1,L%NY
            DO I=1,L%NX
               IF (L%DISCRET(I,J,K,CGSC) == IS_GASPHASE ) THEN
                  L%NUNKH_LOC(NM) = L%NUNKH_LOC(NM) + 1
                  L%DISCRET(I,J,K,UNKH) = L%NUNKH_LOC(NM)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

   ENDIF 

   L%NCS = L%NUNKH_LOC(NM)

ENDDO MESHES_LOOP3

!> 
!> Print some debug messages if requested
!>
MESHES_LOOP4 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   L => SCARC(NM)%LEVEL(NLEVEL_MIN)
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) '     SET_DISCRET - 2 - ', L%NX, L%NY, L%NZ
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) 'TYPE_DISCRET=',TYPE_DISCRET
      WRITE(LU_SCARC,*) 'L%NCS=',L%NCS
      WRITE(LU_SCARC,*) 'DISCRET(.,.,.,CGSC)'
      DO IZ = L%NZ+1,0,-1
      IF (L%NX==11) THEN
         WRITE(LU_SCARC,'(13I6)') ((L%DISCRET(IX, IY, IZ, CGSC), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE IF (L%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((L%DISCRET(IX, IY, IZ, CGSC), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE IF (L%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((L%DISCRET(IX, IY, IZ, CGSC), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((L%DISCRET(IX, IY, IZ, CGSC), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ENDIF
      WRITE(LU_SCARC,*)
      ENDDO
      WRITE(LU_SCARC,*) 'DISCRET(.,.,.,UNKH)'
      DO IZ = L%NZ+1,0,-1
      IF (L%NX==11) THEN
         WRITE(LU_SCARC,'(13I6)') ((L%DISCRET(IX, IY, IZ, UNKH), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE IF (L%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((L%DISCRET(IX, IY, IZ, UNKH), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE IF (L%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((L%DISCRET(IX, IY, IZ, UNKH), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((L%DISCRET(IX, IY, IZ, UNKH), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ENDIF
      WRITE(LU_SCARC,*)
      ENDDO
   ENDIF

ENDDO MESHES_LOOP4

END SUBROUTINE SCARC_SETUP_DISCRET


!> -----------------------------------------------------------------------------
!> MV: Get information about global numbers of unknowns
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DISCRET_GLOBAL(NL)

INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2, IX, IY, IZ
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

SCARC_ROUTINE = 'SCARC_SETUP_DISCRETi_GLOBAL'

!> Get number of local cells and offset vector
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L => SCARC(NM)%LEVEL(NL)
   DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L%NUNKH_LOC(NM2) = SCARC(NM2)%LEVEL(NL)%NUNKH_LOC(NM2)
   ENDDO
ENDDO MESHES_LOOP1

!> Make visible local numbers of cells for all meshes
CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,L%NUNKH_LOC,COUNTS,DISPLS,&
                    MPI_INTEGER,MPI_COMM_WORLD,IERROR)

!> Get global number of cells and offset vector
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L => SCARC(NM)%LEVEL(NL)
   L%NUNKH_GLOB = SUM(L%NUNKH_LOC(1:NMESHES))
   NUNKH_GLOB(NL) = L%NUNKH_GLOB
   DO NM2=2,NMESHES
      L%UNKH_IND(NM2) = L%UNKH_IND(NM2-1) + L%NUNKH_LOC(NM2-1)
   ENDDO
ENDDO MESHES_LOOP2

!> Print some debug messages if requested
MESHES_LOOP3: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) '     SETUP_DISCRET_GLOBAL ', NL
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) 'DISCRET(.,.,.,CGSC)'
      DO IZ = L%NZ+1,0,-1
      IF (L%NX==11) THEN
         WRITE(LU_SCARC,'(13I6)') ((L%DISCRET(IX, IY, IZ, CGSC), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE IF (L%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((L%DISCRET(IX, IY, IZ, CGSC), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE IF (L%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((L%DISCRET(IX, IY, IZ, CGSC), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((L%DISCRET(IX, IY, IZ, CGSC), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ENDIF
      WRITE(LU_SCARC,*)
      ENDDO
      WRITE(LU_SCARC,*) 'DISCRET(.,.,.,UNKH)'
      DO IZ = L%NZ+1,0,-1
      IF (L%NX==11) THEN
         WRITE(LU_SCARC,'(13I6)') ((L%DISCRET(IX, IY, IZ, UNKH), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE IF (L%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((L%DISCRET(IX, IY, IZ, UNKH), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE IF (L%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((L%DISCRET(IX, IY, IZ, UNKH), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((L%DISCRET(IX, IY, IZ, UNKH), IX=0, L%NX+1), IY= L%NY+1,0,-1)
      ENDIF
      WRITE(LU_SCARC,*)
      ENDDO
      WRITE(LU_SCARC,*) 'NUNKH_LOC =', L%NUNKH_LOC(1:NMESHES)
      WRITE(LU_SCARC,*) 'NUNKH_GLOB=', L%NUNKH_GLOB
      WRITE(LU_SCARC,*) 'UNKH_IND  =', L%UNKH_IND(1:NMESHES)
      WRITE(LU_SCARC,*) 'NCS       =', L%NCS
   ENDIF
ENDDO MESHES_LOOP3

END SUBROUTINE SCARC_SETUP_DISCRET_GLOBAL


!> -----------------------------------------------------------------------------
!> MV: Setup DISCRET on coarser levels
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DISCRET_LEVEL(NL)

INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IXF, IYF, IZF, IX, IY, IZ, NSTEP
TYPE (SCARC_LEVEL_TYPE), POINTER :: F, C

SCARC_ROUTINE = 'SCARC_SETUP_DISCRET_LEVEL'

MESHES_LOOP1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   F => SCARC(NM)%LEVEL(NLEVEL_MIN)                  !> fine level
   C => SCARC(NM)%LEVEL(NL)                          !> coarse level

   !> Initialize and preset DISCRET-array on coarser levels
   IF (.NOT. ALLOCATED(C%DISCRET)) &
   ALLOCATE(C%DISCRET(0:C%NX+1,0:C%NY+1,0:C%NZ+1,NCVARS))
   C%DISCRET(:,:,:,CGSC) = IS_UNDEFINED
   C%DISCRET(:,:,:,UNKH) = IS_UNDEFINED

   NSTEP = 2**(NL - NLEVEL_MIN)

   !> In case of structured approach:
   IF (PRES_ON_WHOLE_DOMAIN) THEN
      DO IZ = 1, C%NZ
         IZF = (IZ-1)*NSTEP + 1
         DO IY = 1, C%NY
            IYF = (IY-1)*NSTEP + 1
            DO IX = 1, C%NX
               IXF = (IX-1)*NSTEP + 1
               C%DISCRET(IX,IY,IZ,CGSC) = F%DISCRET(IXF, IYF, IZF, CGSC)
               C%NUNKH_LOC(NM) = C%NUNKH_LOC(NM) + 1
               C%DISCRET(IX,IY,IZ,UNKH) = C%NUNKH_LOC(NM)
            ENDDO
         ENDDO
      ENDDO

   !> In case of unstructured approach:
   ELSE
      DO IZ = 1, C%NZ
         IZF = (IZ-1)*NSTEP + 1
         DO IY = 1, C%NY
            IYF = (IY-1)*NSTEP + 1
            DO IX = 1, C%NX
               IXF = (IX-1)*NSTEP + 1
               C%DISCRET(IX,IY,IZ,CGSC) = F%DISCRET(IXF, IYF, IZF, CGSC)
               IF (F%DISCRET(IXF, IYF, IZF, CGSC) == IS_GASPHASE) THEN
                  C%NUNKH_LOC(NM) = C%NUNKH_LOC(NM) + 1
                  C%DISCRET(IX,IY,IZ,UNKH) = C%NUNKH_LOC(NM)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   C%NCS = C%NUNKH_LOC(NM)

ENDDO MESHES_LOOP1

!> 
!> Print some debug messages if requested
!>
MESHES_LOOP2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   C => SCARC(NM)%LEVEL(NL)
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) '     SETUP_DISCRET_LEVEL ', NL, C%NCS
      WRITE(LU_SCARC,*) '============================================================='
      WRITE(LU_SCARC,*) 'NSTEP=',NSTEP
      WRITE(LU_SCARC,*) 'DISCRET(.,.,.,CGSC)'
      DO IZ = C%NZ+1,0,-1
      IF (C%NX==11) THEN
         WRITE(LU_SCARC,'(13I6)') ((C%DISCRET(IX, IY, IZ, CGSC), IX=0, C%NX+1), IY= C%NY+1,0,-1)
      ELSE IF (C%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((C%DISCRET(IX, IY, IZ, CGSC), IX=0, C%NX+1), IY= C%NY+1,0,-1)
      ELSE IF (C%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((C%DISCRET(IX, IY, IZ, CGSC), IX=0, C%NX+1), IY= C%NY+1,0,-1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((C%DISCRET(IX, IY, IZ, CGSC), IX=0, C%NX+1), IY= C%NY+1,0,-1)
      ENDIF
      WRITE(LU_SCARC,*)
      ENDDO
      WRITE(LU_SCARC,*) 'DISCRET(.,.,.,UNKH)'
      DO IZ = C%NZ+1,0,-1
      IF (C%NX==11) THEN
         WRITE(LU_SCARC,'(13I6)') ((C%DISCRET(IX, IY, IZ, UNKH), IX=0, C%NX+1), IY= C%NY+1,0,-1)
      ELSE IF (C%NX==16) THEN
         WRITE(LU_SCARC,'(18I6)') ((C%DISCRET(IX, IY, IZ, UNKH), IX=0, C%NX+1), IY= C%NY+1,0,-1)
      ELSE IF (C%NX==8) THEN
         WRITE(LU_SCARC,'(10I6)') ((C%DISCRET(IX, IY, IZ, UNKH), IX=0, C%NX+1), IY= C%NY+1,0,-1)
      ELSE
         WRITE(LU_SCARC,'(6I6)') ((C%DISCRET(IX, IY, IZ, UNKH), IX=0, C%NX+1), IY= C%NY+1,0,-1)
      ENDIF
      WRITE(LU_SCARC,*)
      ENDDO
      WRITE(LU_SCARC,*) 'NUNKH_LOC=', C%NUNKH_LOC(1:NMESHES)
      WRITE(LU_SCARC,*) 'NCS=', C%NCS
      WRITE(LU_SCARC,*) 'NX =', C%NX
      WRITE(LU_SCARC,*) 'NY =', C%NY
      WRITE(LU_SCARC,*) 'NZ =', C%NZ
  ENDIF
ENDDO MESHES_LOOP2


END SUBROUTINE SCARC_SETUP_DISCRET_LEVEL


!> ------------------------------------------------------------------------------------------------
!> Only for debugging reasons: print out vector information on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_LEVEL (NVECTOR, CNAME, NL)
INTEGER, INTENT(IN):: NVECTOR, NL
REAL (EB), POINTER, DIMENSION(:)     :: VC
REAL (EB):: VALUES(0:100)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CNAME
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

!IF (TYPE_SCOPE /= NSCARC_SCOPE_MAIN) RETURN
IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN


DO NM = 1, NMESHES

   IF (PROCESS(NM) /= MYID) CYCLE

   L => SCARC(NM)%LEVEL(NL)
   CALL POINT_TO_VECTOR (NVECTOR, NM, NL, VC)

   IF (L%NX>10.OR.L%NY>10.OR.L%NZ>10) CYCLE

   NX8=MIN(10,L%NX)
   NY8=MIN(10,L%NY)
   NZ8=MIN(10,L%NZ)

   WRITE(LU_SCARC,*) '=========================================================='
   WRITE(LU_SCARC,2001) SCARC_ROUTINE, CNAME, NM, NL
   WRITE(LU_SCARC,2002) L%NC, L%NCE, NX8, NY8, NZ8, NVECTOR, SIZE(VC)
   WRITE(LU_SCARC,*) '=========================================================='
   !IF (NL == NLEVEL_MIN) THEN
         DO KK = NZ8, 1, - 1
            DO JJ = NY8, 1, - 1
               DO II=1,NX8
                  IF (L%DISCRET(II,JJ,KK,CGSC) == IS_GASPHASE) THEN
                     IC=L%DISCRET(II,JJ,KK,UNKH)
                     IF (ABS(VC(IC))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VC(IC)
                     ENDIF
                  ELSE
                     VALUES(II)=-999999.0_EB
                  ENDIF
               ENDDO
               WRITE(LU_SCARC, '(10e11.3)') (VALUES(II), II=1, NX8)
            ENDDO
         ENDDO
   !ENDIF

ENDDO

!2000 FORMAT('=== ',A,' : ',A,' on mesh ',I8,' on level ',I8, ': NX, NY, NZ=',3I8,': NVECTOR=',I8)
2001 FORMAT('=== ',A,' : ',A,' on mesh ',I8,' on level ',I8)
2002 FORMAT('=== NC = ',I4,' NCE=',I4, ': NX, NY, NZ=',3I4,': NVECTOR=',I3,': Size=',I8)
END SUBROUTINE SCARC_DEBUG_LEVEL


!> ------------------------------------------------------------------------------------------------
!> Debug requested quantity
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_QUANTITY(NTYPE, NL, CNAME)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM, IP, IC, IW, I, J, K, IOR0, INBR
CHARACTER (*), INTENT(IN) :: CNAME
CHARACTER (20) :: LINE
TYPE (MESH_TYPE), POINTER :: M
TYPE (OSCARC_TYPE) , POINTER :: OS
TYPE (SCARC_LEVEL_TYPE), POINTER :: L, OL
TYPE (SCARC_MATRIX_TYPE), POINTER :: A, P, R, OP, OR, AS
TYPE (SCARC_AMG_TYPE), POINTER :: AMG, OAMG
TYPE (SCARC_MAPPING_TYPE), POINTER :: MAP, OMAP

IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN

SELECT CASE (NTYPE)

   !> ------------------------------------------------------------------------------------------------
   !> Debug system matrix A (corresponding to system type)
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIX)

      IF (NL > NLEVEL_MIN) RETURN
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%LEVEL(NL)%A
         WRITE(LU_SCARC,*) '----------- SHOWING FULL MATRIX ENTRIES'
         WRITE(LU_SCARC,*) 'NAV =',A%NAV
         WRITE(LU_SCARC,*) 'NAC =',A%NAC
         WRITE(LU_SCARC,*) 'NAR =',A%NAR
         WRITE(LU_SCARC,*) 'NCS =',L%NCS
         WRITE(LU_SCARC,*) 'SIZE(A%VAL) =',SIZE(A%VAL)
         WRITE(LU_SCARC,*) 'SIZE(A%COL) =',SIZE(A%COL)
         WRITE(LU_SCARC,*) 'SIZE(A$ROW) =',SIZE(A%ROW)
         WRITE(LU_SCARC,*) '---------------------- A%ROW:', L%NCS
         WRITE(LU_SCARC,'(4i9)') (A%ROW(IC), IC=1,L%NCS+1)
         WRITE(LU_SCARC,*) '---------------------- A%COL:'
         DO IC = 1, L%NCS
            WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(A%COL(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) '---------------------- A:', A%ROW(IC+1)- A%ROW(IC)
         DO IC = 1, L%NCS
            WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(A%VAL(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
         ENDDO
#ifdef WITH_MKL
         IF ((TYPE_METHOD == NSCARC_METHOD_MKL) .AND. (TYPE_MKL == NSCARC_MKL_GLOBAL)) THEN
            !WRITE(LU_SCARC,*) '---------------------- AG_COL:'
            !DO IC = 1, L%NCS
               !WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(A%COL_GLOBAL(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
               !WRITE(LU_SCARC,*)  IC,':',(A%COL_GLOBAL(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
            !ENDDO
         ENDIF
         DO IC=1,L%NCS
            DO IP=A%ROW(IC),A%ROW(IC+1)-1
                WRITE(LU_SCARC,'(2I8,F18.12)') IC,A%COL(IP),A%VAL(IP)
            ENDDO
            WRITE(LU_SCARC,*)
         ENDDO
         DO IC=1,L%NCS
            DO IP=A%ROW(IC),A%ROW(IC+1)-1
                IF (IC == A%COL(IP)) WRITE(LU_SCARC,'(2I8,F18.12)') IC,A%COL(IP),A%VAL(IP)
            ENDDO
         ENDDO
#endif
      ENDDO

         !CALL SCARC_PRINT_MATRIX (L%A, A%ROW, A%COL, L%NCS, L%NCS, NM, NL, 'A')
         !CALL SCARC_PRINT_MATRIX2(L%A, A%ROW, A%COL, L%NCS, L%NCS, L%NX, L%NY, L%NZ, NM, NL, 'A')
         !CALL SCARC_MATLAB_MATRIX(L%A, A%ROW, A%COL, L%NCS, L%NCS, NM, NL, 'A')

   !> ------------------------------------------------------------------------------------------------
   !> Debug symmetric system matrix AS
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIXS)

      IF (NL > NLEVEL_MIN) RETURN

#ifdef WITH_MKL
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         AS => SCARC(NM)%LEVEL(NL)%AS
         WRITE(LU_SCARC,*) '----------- SHOWING SYMMETRIC MATRIX ENTRIES'
         WRITE(LU_SCARC,*) 'SIZE(AS) =',    SIZE(AS%VAL)
         WRITE(LU_SCARC,*) 'SIZE(AS%COL) =',SIZE(AS%COL)
         WRITE(LU_SCARC,*) 'SIZE(AS%ROW) =',SIZE(AS%ROW)
         WRITE(LU_SCARC,*) '---------------------- AS%ROW:', L%NCS
         WRITE(LU_SCARC,'(4i9)') (AS%ROW(IC), IC=1,L%NCS+1)
         WRITE(LU_SCARC,*) '---------------------- AS%COL:'
         DO IC = 1, L%NC
            WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(AS%COL(IP),IP=AS%ROW(IC),AS%ROW(IC+1)-1)
         ENDDO
         !WRITE(LU_SCARC,*) '---------------------- AS:', AS%ROW(IC+1)- AS%ROW(IC)
         DO IC = 1, L%NC
            WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(AS%VAL(IP),IP=AS%ROW(IC),AS%ROW(IC+1)-1)
         ENDDO
         DO IC=1,L%NC
            DO IP=AS%ROW(IC),AS%ROW(IC+1)-1
                WRITE(LU_SCARC,'(2I8,F18.12)') IC,AS%COL(IP),AS%VAL(IP)
            ENDDO
            WRITE(LU_SCARC,*)
         ENDDO
            WRITE(LU_SCARC,*)
            WRITE(LU_SCARC,*)
         DO IC=1,L%NC
            DO IP=AS%ROW(IC),AS%ROW(IC+1)-1
                IF (IC == AS%COL(IP)) WRITE(LU_SCARC,'(2I8,F18.12)') IC,AS%COL(IP),AS%VAL(IP)
            ENDDO
            WRITE(LU_SCARC,*)
         ENDDO

         !CALL SCARC_PRINT_MATRIX (L%AS%VAL, AS%ROW, AS%COL, L%NC, L%NC, NM, NL, 'AS')
         !CALL SCARC_PRINT_MATRIX2(L%AS%VAL, AS%ROW, AS%COL, L%NC, L%NC, L%NX, L%NY, L%NZ, NM, NL, 'AS')
         !CALL SCARC_MATLAB_MATRIX(L%AS%VAL, AS%ROW, AS%COL, L%NC, L%NC, NM, NL, 'AS')
      ENDDO
#endif

   !> ------------------------------------------------------------------------------------------------
   !> Debug extended system matrix A (corresponding to system type)
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIXE)

      IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         A => SCARC(NM)%LEVEL(NL)%A
         WRITE(LU_SCARC,*) 'PRINTING EXTENDED MATRIX, NCE=',L%NC
         WRITE(LU_SCARC,*) '---------------------- A%ROW:', L%NC
         WRITE(LU_SCARC,'(4i9)') (A%ROW(IC), IC=1,L%NC+1)
         WRITE(LU_SCARC,*) '---------------------- A%COL:'
         DO IC = 1, L%NC
            WRITE(LU_SCARC,'(i5,a,20i9)') IC,':',(A%COL(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) '---------------------- A:'
         DO IC = 1, L%NC
            WRITE(LU_SCARC,'(i5,a,20f9.2)') IC,':',(A%VAL(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) 'SIZE(A%VAL) =',SIZE(A%VAL)
         WRITE(LU_SCARC,*) 'SIZE(A%COL) =',SIZE(A%COL)
         WRITE(LU_SCARC,*) 'SIZE(A%ROW) =',SIZE(A%ROW)
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug prolongation matrix P (only for compact system)
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_PROLONGATION)

      IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         AMG => SCARC(NM)%LEVEL(NL)%AMG
         P => SCARC(NM)%LEVEL(NL)%P
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         WRITE(LU_SCARC,*) '============= P%ROW:', NM
         WRITE(LU_SCARC,'(4i9)') (P%ROW(IC), IC=1, L%AMG%NCE+1)
         WRITE(LU_SCARC,*) '============= P%COL:', NM, L%AMG%NCE
         DO IC=1,L%NCE
            WRITE(LU_SCARC,'(i5,a,10i9)') IC,':',(P%COL(IP), IP=P%ROW(IC),P%ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) '============= P:', NM, L%NCE
         DO IC=1,AMG%NCE
            WRITE(LU_SCARC,'(i5,a,10f9.2)') IC,':',(P%VAL(IP), IP=P%ROW(IC),P%ROW(IC+1)-1)
         ENDDO
         DO NOM = 1, NMESHES
            IF (NOM == NM) CYCLE
            OS => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE
            OL => OS%LEVEL(NL)
            OP => OL%P
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P%ROW:', OL%NCG
            WRITE(LU_SCARC,'(4i9)') (OP%ROW(IC), IC=1, OL%NCG+1)
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P%COL:', OL%NCG
            DO IC=1,OL%NCG
               WRITE(LU_SCARC,'(i5,a,10i9)') IC,':',(OP%COL(IP), IP=OP%ROW(IC),OP%ROW(IC+1)-1)
            ENDDO
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P:'
            DO IC=1,OL%NCG
               WRITE(LU_SCARC,'(i5,a,10f9.2)') IC,':',(OP%VAL(IP), IP=OP%ROW(IC),OP%ROW(IC+1)-1)
            ENDDO
         ENDDO
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug restriction matrix R (only for compact system)
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_RESTRICTION)

      IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         R => SCARC(NM)%LEVEL(NL)%R
         AMG => SCARC(NM)%LEVEL(NL)%AMG
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         WRITE(LU_SCARC,*) '============= R%ROW:', NM, L%AMG%NCCE
         WRITE(LU_SCARC,'(4i9)') (R%ROW(IC), IC=1, L%AMG%NCCE)
         WRITE(LU_SCARC,*) '============= R%COL:', NM
         DO IC=1,AMG%NCCE
            WRITE(LU_SCARC,'(i5,a,10i9)') IC,':',(R%COL(IP), IP=R%ROW(IC),R%ROW(IC+1)-1)
         ENDDO
         WRITE(LU_SCARC,*) '============= R:', NM
         DO IC=1,AMG%NCCE
            WRITE(LU_SCARC,'(i5,a,10f9.2)') IC,':',(R%VAL(IP), IP=R%ROW(IC),R%ROW(IC+1)-1)
         ENDDO
         !DO NOM = 1, NMESHES
         DO NOM = 1, NMESHES
            IF (NOM == NM) CYCLE
            OS => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE
            OL => OS%LEVEL(NL)
            OAMG => OS%LEVEL(NL)%AMG
            OR => OL%R
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%R%ROW:', OAMG%NCC
            WRITE(LU_SCARC,'(4i9)') (OR%ROW(IC), IC=1, OAMG%NCC+1)
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%R%COL:', OAMG%NCC
            DO IC=1,OAMG%NCC
               WRITE(LU_SCARC,'(i5,a,10i9)') IC,':',(OR%COL(IP), IP=OR%ROW(IC),OR%ROW(IC+1)-1)
            ENDDO
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%P:', OAMG%NCC
            DO IC=1,OAMG%NCC
               WRITE(LU_SCARC,'(i5,a,10f9.2)') IC,':',(OR%VAL(IP), IP=OR%ROW(IC),OR%ROW(IC+1)-1)
            ENDDO
         ENDDO
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug FACEINFO
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_FACEINFO)

      IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         M => MESHES(NM)
         DO IOR0 = -3, 3
            IF (IOR0 == 0) CYCLE
            WRITE(LU_SCARC,*) '========================================='
            WRITE(LU_SCARC,*) '============= DEBUGGING FACE(',IOR0,'): FOR LEVEL ', NL
            WRITE(LU_SCARC,*) '========================================='
            WRITE(LU_SCARC,*) 'FACE(.)%NFX:', L%FACE(IOR0)%NFX
            WRITE(LU_SCARC,*) 'FACE(.)%NFY:', L%FACE(IOR0)%NFY
            WRITE(LU_SCARC,*) 'FACE(.)%NFZ:', L%FACE(IOR0)%NFZ
            WRITE(LU_SCARC,*) 'FACE(.)%NWL:', L%FACE(IOR0)%NFW
            WRITE(LU_SCARC,*) 'FACE(.)%IWG_PTR:', L%FACE(IOR0)%IWG_PTR
            WRITE(LU_SCARC,*) '----------------------------------------------'
            WRITE(LU_SCARC,*)
            !WRITE(LU_SCARC,*) 'FACE(.)%DH :', L%FACE(IOR0)%DH
            WRITE(LU_SCARC,*)
            IF (L%FACE(IOR0)%NUM_NEIGHBORS /= 0) THEN
               DO INBR=1,L%FACE(IOR0)%NUM_NEIGHBORS
                  NOM = L%FACE(IOR0)%NEIGHBORS(INBR)
                  OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  OMAP => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)%MAP
                  WRITE(LU_SCARC,*) 'NUM_NEIGHBORS:', L%FACE(IOR0)%NUM_NEIGHBORS
                  WRITE(LU_SCARC,*) 'NOM:', NOM
                  WRITE(LU_SCARC,*) 'SIZE(OMAP%ICG_TO_IWG)=',SIZE(OMAP%ICG_TO_IWG)
                  WRITE(LU_SCARC,*) 'SIZE(OMAP%ICG_TO_ICE)=',SIZE(OMAP%ICG_TO_ICE)
                  WRITE(LU_SCARC,*) 'SIZE(OMAP%ICG_TO_ICO)=',SIZE(OMAP%ICG_TO_ICO)
                  WRITE(LU_SCARC,*) 'SIZE(OMAP%IWL_TO_ICW)=',SIZE(OMAP%IWL_TO_ICW)
                  WRITE(LU_SCARC,*) 'SIZE(OMAP%IWL_TO_ICO)=',SIZE(OMAP%IWL_TO_ICO)
                  IF (BAMG) WRITE(LU_SCARC,*) 'SIZE(OMAP%IWL_TO_ICG)=',SIZE(OMAP%IWL_TO_ICG)
                  WRITE(LU_SCARC,'(a,i8,a,2f12.6)') '---OL(',NOM,')%DH :',OL%DH
                  WRITE(LU_SCARC,'(a,i8,a,2i8)') '---OL(',NOM,')%IOR   :',OL%IOR
                  WRITE(LU_SCARC,'(a,i8,a,2i8)') '---OL(',NOM,')%NWL(.):',OL%NWL
                  WRITE(LU_SCARC,'(a,i8,a,2i8)') '---OL(',NOM,')%NCG(.):',OL%NCG
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OL(',NOM,')%ICG_TO_IWG:'
                  DO IW = 1, OL%NCG
                     WRITE(LU_SCARC,'(32i8)') OMAP%ICG_TO_IWG(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OL(',NOM,')%ICG_TO_ICO:'
                  DO IW = 1, OL%NCG
                     WRITE(LU_SCARC,'(32i8)') OMAP%ICG_TO_ICO(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OL(',NOM,')%ICG_TO_ICE:'
                  DO IW = 1, OL%NCG
                     WRITE(LU_SCARC,'(32i8)') OMAP%ICG_TO_ICE(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OL(',NOM,')%IWL_TO_IWG:'
                  DO IW = 1, OL%NWL
                     WRITE(LU_SCARC,'(32i8)') OMAP%IWL_TO_IWG(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OL(',NOM,')%IWL_TO_ICW:'
                  DO IW = 1, OL%NWL
                     WRITE(LU_SCARC,'(32i8)') OMAP%IWL_TO_ICW(IW)
                  ENDDO
                  WRITE(LU_SCARC,*)
                  WRITE(LU_SCARC,'(a,i8,a)') '------OL(',NOM,')%IWL_TO_ICO:'
                  DO IW = 1, OL%NWL
                     WRITE(LU_SCARC,'(32i8)') OMAP%IWL_TO_ICO(IW)
                  ENDDO
                  !IF (BAMG) THEN
                  !   WRITE(LU_SCARC,*)
                  !   WRITE(LU_SCARC,'(a,i4,a)') '------OL(',NOM,')%IWL_TO_ICG:'
                  !   DO IW = 1, OMAP%NWL
                  !      WRITE(LU_SCARC,'(32i8)') OMAP%IWL_TO_ICG(IW, 1:OMAP%NCPL)
                  !   ENDDO
                  !ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug WALLINFO
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_WALLINFO)
      IF (TYPE_DEBUG < NSCARC_DEBUG_MUCH) RETURN

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         MAP => L%MAP
         M => MESHES(NM)
         WRITE(LU_SCARC,*) 'SIZE(MAP%ICE_TO_IWG)=',SIZE(MAP%ICE_TO_IWG)
         WRITE(LU_SCARC,*) 'SIZE(MAP%ICE_TO_IWL)=',SIZE(MAP%ICE_TO_IWL)
         WRITE(LU_SCARC,*) 'SIZE(MAP%ICE_TO_ICG)=',SIZE(MAP%ICE_TO_ICG)
         WRITE(LU_SCARC,*) 'NM  =',NM
         WRITE(LU_SCARC,*) 'NL  =',NL
         WRITE(LU_SCARC,*) 'NCE =',L%NCE
         WRITE(LU_SCARC,*) 'NCS =',L%NCS
         WRITE(LU_SCARC,*) 'NC  =',L%NC
         WRITE(LU_SCARC,*) 'NW  =',L%NW
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,'(a,i6,a)') '------MAP%ICE_TO_IWG:'
         WRITE(LU_SCARC,'(16i6)') (MAP%ICE_TO_IWG(IW), IW = L%NCS+1, L%NCE)
         WRITE(LU_SCARC,'(a,i6,a)') '------MAP%ICE_TO_IWL:'
         WRITE(LU_SCARC,'(16i6)') (MAP%ICE_TO_IWL(IW), IW = L%NCS+1, L%NCE)
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,'(a,i6,a)') '------MAP%ICE_TO_ICG:'
         WRITE(LU_SCARC,'(16i6)') (MAP%ICE_TO_ICG(IW), IW = L%NCS+1, L%NCE)
         WRITE(LU_SCARC,*)
         IF (NL == 1) THEN
         WRITE(LU_SCARC,*) '============= WALL(.)%IXG:', NM
         WRITE(LU_SCARC,'(16i6)') (L%WALL(IW)%IXG, IW=1,L%NW)
         WRITE(LU_SCARC,*) '============= WALL(.)%IYG:', NM
         WRITE(LU_SCARC,'(16i6)') (L%WALL(IW)%IYG, IW=1,L%NW)
         WRITE(LU_SCARC,*) '============= WALL(.)%IZG:', NM
         WRITE(LU_SCARC,'(16i6)') (L%WALL(IW)%IZG, IW=1,L%NW)
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%IXW:', NM
         WRITE(LU_SCARC,'(16i6)') (L%WALL(IW)%IXW, IW=1,L%NW)
         WRITE(LU_SCARC,*) '============= WALL(.)%IYW:', NM
         WRITE(LU_SCARC,'(16i6)') (L%WALL(IW)%IYW, IW=1,L%NW)
         WRITE(LU_SCARC,*) '============= WALL(.)%IZW:', NM
         WRITE(LU_SCARC,'(16i6)') (L%WALL(IW)%IZW, IW=1,L%NW)
         ENDIF
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%IOR:', NM
         WRITE(LU_SCARC,'(16i6)') (L%WALL(IW)%IOR, IW=1,L%NW)
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%NOM:', NM
         WRITE(LU_SCARC,'(16i6)') (L%WALL(IW)%NOM, IW=1,L%NW)
         WRITE(LU_SCARC,*)
         IF (NL == 1) THEN
         WRITE(LU_SCARC,*) '============= WALL(.)%NCPL:', NM
         !WRITE(LU_SCARC,'(16i10)') (L%WALL(IW)%NCPL, IW=1,L%NW)
         WRITE(LU_SCARC,*) (L%WALL(IW)%NCPL, IW=1,L%NW)
         ENDIF
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%ICW:', NM
         DO IW=1,L%NW
            WRITE(LU_SCARC,'(a,i6, a,i6)') 'IW=',IW,':',L%WALL(IW)%ICW
         ENDDO
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%ICO:', NM
         DO IW=1,L%NW
            IF (L%WALL(IW)%NOM /=0) &
               WRITE(LU_SCARC,'(a,i6, a,i6)') 'IW=',IW,':',L%WALL(IW)%ICO
         ENDDO
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%ICE:', NM
         DO IW=1,L%NW
            IF (L%WALL(IW)%NOM /=0) &
               WRITE(LU_SCARC,'(a,i6, a,i6)') 'IW=',IW,':',L%WALL(IW)%ICE
         ENDDO
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*)
         WRITE(LU_SCARC,*) '============= WALL(.)%ICG:', NM
         DO IW=1,L%NW
            IF (L%WALL(IW)%NOM /=0) &
               WRITE(LU_SCARC,'(a,i6, a,i6)') 'IW=',IW,':',L%WALL(IW)%ICG
         ENDDO
         WRITE(LU_SCARC,*) '====================================================='
         WRITE(LU_SCARC,*) ' Plotting out M%WALL-structure'
         WRITE(LU_SCARC,*) '====================================================='
         WRITE(LU_SCARC,*) ' M%WALL(.)%WALL_INDEX'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%WALL_INDEX, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%SURF_INDEX'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%SURF_INDEX, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%BACK_INDEX'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%BACK_INDEX, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%BOUNDARY_TYPE'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%BOUNDARY_TYPE, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%OBST_INDEX'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%OBST_INDEX, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%PRESSURE_BC_INDEX'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%PRESSURE_BC_INDEX, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%PRESSURE_ZONE'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%ONE_D%PRESSURE_ZONE, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%VENT_INDEX'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%VENT_INDEX, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%NOM'
         WRITE(LU_SCARC,'(16i6)') (M%EXTERNAL_WALL(IW)%NOM, IW=1,L%N_EXTERNAL_WALL_CELLS)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%II'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%ONE_D%II, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%JJ'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%ONE_D%JJ, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%KK'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%ONE_D%KK, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%IIG'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%ONE_D%IIG, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%JJG'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%ONE_D%JJG, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%KKG'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%ONE_D%KKG, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%N_LAYER_CELLS'
         WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%ONE_D%IOR, IW=1,L%NW)
         WRITE(LU_SCARC,*) ' M%WALL(.)%ONE_D%IOR'
         !WRITE(LU_SCARC,'(16i6)') (M%WALL(IW)%ONE_D%N_LAYER_CELLS, IW=1,L%NW)

      ENDDO
      !ENDIF


 !> ------------------------------------------------------------------------------------------------
 !> Debug PRESSURE_BC_INDEX
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_BCINDEX)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC, '(8i8)') (L%WALL(J)%BTYPE, J=1,L%NW)
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug WALL_CELL
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_ACELL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC, '(8i8)') (L%WALL(IW)%ICW, IW=1,L%NW)
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug GHOST_CELL
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_GCELL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,*) '-------- GHOST,  NGE=',L%NCG
         DO IW=1,L%NW
            !IF (L%WALL(IW)%NOM/=0) &
            !>   WRITE(LU_SCARC, '(8i8)') (L%WALL(IW)%ICG(I), I=1,L%WALL(IW)%NCPL)
         ENDDO
         !DO IW=1,L%NW
         !>   IF (L%WALL(IW)%NOM/=0.AND.NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) &
         !>      WRITE(LU_SCARC, '(8i8)') (L%WALL(IW)%ICG2(I), I=1,L%WALL(IW)%NCPL)
         !ENDDO
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug NOM_CELL
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_NCELL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,*) '-------- NOM, NCE= ', L%NCE
         WRITE(*,*) 'Achtung, hier wegen ICN schauen!'
         !DO IW=1,L%NW
         !   !IF (L%WALL(IW)%NOM/=0) &
         !   !>   WRITE(LU_SCARC, '(8i8)') (L%WALL(IW)%ICN(I), I=1,L%WALL(IW)%NCPL)
         !ENDDO
         !DO IW=1,L%NW
         !>   IF (L%WALL(IW)%NOM/=0.AND.NL==NLEVEL_MIN.AND.TYPE_LAYER == NSCARC_LAYER_TWO) &
         !>      WRITE(LU_SCARC, '(8i8)') (L%WALL(IW)%ICN2(I), I=1,L%WALL(IW)%NCPL)
         !ENDDO
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug SUBDIVISION
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_SUBDIVISION)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR= 1 '
         WRITE(LU_SCARC,'(3i8)') (L%SUBDIVISION(I, 1), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR=-1 '
         WRITE(LU_SCARC,'(3i8)') (L%SUBDIVISION(I,-1), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR= 2 '
         WRITE(LU_SCARC,'(3i8)') (L%SUBDIVISION(I, 2), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR=-2 '
         WRITE(LU_SCARC,'(3i8)') (L%SUBDIVISION(I,-2), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR= 3 '
         WRITE(LU_SCARC,'(3i8)') (L%SUBDIVISION(I, 3), I=1,3)
         WRITE(LU_SCARC,*) 'SUBDIVISION IOR=-3 '
         WRITE(LU_SCARC,'(3i8)') (L%SUBDIVISION(I,-3), I=1,3)
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug MEASURE
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MEASURE)

      WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, 1, NL
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         IF (NL == 1) THEN
            DO K = L%NZ, 1, -1
               WRITE(LU_SCARC, '(10f11.5)') (AMG%MEASURE((K-1)*L%NX*L%NY+I), I=1, L%NX)
            ENDDO
            WRITE(LU_SCARC,*) '--------------------'
            WRITE(LU_SCARC, '(4f11.5)') (AMG%MEASURE(IC), IC=L%NC+1, L%NCE)
         ELSE
            WRITE(LU_SCARC, '(4f11.5)') (AMG%MEASURE(IC), IC=1, L%NC)
         ENDIF
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug CELL_TYPE
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_CELL_TYPE)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
        WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
        L => SCARC(NM)%LEVEL(NL)
        IF (NL == 1) THEN
           DO K = L%NZ, 1, -1
              WRITE(LU_SCARC, '(8I4)') (AMG%CELL_TYPE((K-1)*L%NX*L%NY+I), I=1, L%NX)
           ENDDO
           WRITE(LU_SCARC,*) '--------------------'
           WRITE(LU_SCARC, '(4I4)') (AMG%CELL_TYPE(IC), IC=L%NC+1, L%NCE)
        ELSE
           WRITE(LU_SCARC, '(4I4)') (AMG%CELL_TYPE(IC), IC=1, L%NC)
        ENDIF
      ENDDO
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         DO NOM = 1, NMESHES
            IF (NOM == NM) CYCLE
            OS => SCARC(NM)%OSCARC(NOM)
            IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE
            OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            WRITE(LU_SCARC,'(a,i3,a,i3,a,i3)') '============= SCARC(',NM,')%OSCARC(',NOM,')%AMG%CELL_TYPE:', OL%NCG
            DO IC=1,OL%NCG
               WRITE(LU_SCARC,'(i8,a,i8)') IC,':',OL%AMG%CELL_TYPE(IC)
            ENDDO
         ENDDO
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug GRAPH
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_GRAPH)

      IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) THEN
       IF (NMESHES == 4 .OR. NMESHES==1) THEN                !> only temporarily
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
           WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
           L => SCARC(NM)%LEVEL(NL)
           IF (NL == 1) THEN
              DO K = L%NZ, 1, -1
                 WRITE(LU_SCARC, '(8I4)') (AMG%GRAPH((K-1)*L%NX*L%NY+I), I=1, L%NX)
              ENDDO
              WRITE(LU_SCARC,*) '--------------------'
              WRITE(LU_SCARC, '(4I4)') (AMG%GRAPH(IC), IC=L%NC+1, L%NCE)
           ELSE
              WRITE(LU_SCARC, '(4I4)') (AMG%GRAPH(IC), IC=1, L%NC)
           ENDIF
         ENDDO
       ELSEIF (NMESHES == 4.AND.NL==1) THEN
         L => SCARC(1)%LEVEL(NL)
         WRITE(LU_SCARC,'(10i8)') L%AMG%GRAPH(1:L%NC)
         L => SCARC(2)%LEVEL(NL)
         WRITE(LU_SCARC,'(10i8)') L%AMG%GRAPH(1:L%NC)
         L => SCARC(3)%LEVEL(NL)
         WRITE(LU_SCARC,'(9i8)') L%AMG%GRAPH(1:L%NC)
         L => SCARC(4)%LEVEL(NL)
         WRITE(LU_SCARC,'(9i8)') L%AMG%GRAPH(1:L%NC)
       ELSEIF (NMESHES == 4.AND.NL==2) THEN
         L => SCARC(1)%LEVEL(NL)
         WRITE(LU_SCARC,'(5i8)') L%AMG%GRAPH(1:L%NC)
         L => SCARC(2)%LEVEL(NL)
         WRITE(LU_SCARC,'(5i8)') L%AMG%GRAPH(1:L%NC)
         L => SCARC(3)%LEVEL(NL)
         WRITE(LU_SCARC,'(4i8)') L%AMG%GRAPH(1:L%NC)
         L => SCARC(4)%LEVEL(NL)
         WRITE(LU_SCARC,'(4i8)') L%AMG%GRAPH(1:L%NC)
       ELSEIF (NMESHES==1.AND.NL==1) THEN
         L => SCARC(1)%LEVEL(NL)
         WRITE(LU_SCARC,'(9i8)') L%AMG%GRAPH(1:L%NC)
       ELSEIF (NMESHES==1.AND.NL==2) THEN
         L => SCARC(1)%LEVEL(NL)
         WRITE(LU_SCARC,'(3i8)') L%AMG%GRAPH(1:L%NC)
       ENDIF
      ENDIF
      IF (TYPE_DEBUG>=NSCARC_DEBUG_LESS) THEN
      WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, 1, NL
      IF (NMESHES == 1.OR.NL>100) THEN                !> only temporarily
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         IF (NL == 1) THEN
            DO K = L%NZ, 1, -1
               WRITE(LU_SCARC, '(10i8)') (AMG%GRAPH((K-1)*L%NX*L%NY+I), I=1, L%NX)
            ENDDO
         ELSE
            WRITE(LU_SCARC, '(4i8)') (AMG%GRAPH(IC), IC=1, L%NC)
         ENDIF
      ENDDO
      ENDIF
      ENDIF

 !> ------------------------------------------------------------------------------------------------
 !> Debug COARSE
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_COARSE)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(LU_SCARC,1000) SCARC_ROUTINE, CNAME, NM, NL
         L => SCARC(NM)%LEVEL(NL)
         DO K = L%NZ,1,-1
            DO I = 1,L%NX
               IC = (K-1)*L%NX + I
               IF (AMG%CELL_TYPE(IC) < 0) THEN
                  LINE(I:I) = 'O'
               ELSE IF (AMG%CELL_TYPE(IC) > 0) THEN
                  LINE(I:I) = 'X'
               ENDIF
            ENDDO
            WRITE(LU_SCARC,'(12A4)') (LINE(I:I), I=1,L%NX)
         ENDDO
      ENDDO

END SELECT

1000 FORMAT('======================================================================================',/, &
            '=== ',A20,' : ', A25,' for mesh ',i3,' on level ', i3, /, &
            '======================================================================================')

END SUBROUTINE SCARC_DEBUG_QUANTITY

!> ------------------------------------------------------------------------------------------------
!> Set vector to a constant value - only for debugging reasons
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_VALUE(VAL, NVECTOR, NL)
INTEGER , INTENT(IN) :: NVECTOR, NL
REAL(EB), INTENT(IN) :: VAL
REAL(EB), POINTER, DIMENSION(:) :: V
INTEGER :: NM, IC
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L

IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
DO NM = 1, NMESHES
   IF (PROCESS(NM) /= MYID) CYCLE
   L => SCARC(NM)%LEVEL(NL)
   CALL POINT_TO_VECTOR(NVECTOR, NM, NL, V)
   DO IC = 1, L%NC
     V(IC) = VAL*NM
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_VALUE


!> ------------------------------------------------------------------------------------------------
!> Save dump of vector in dump-directory
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_QUANTITY (NVECTOR, CNAME, NL)
INTEGER, INTENT(IN):: NVECTOR, NL
CHARACTER(*), INTENT(IN):: CNAME
REAL (EB), POINTER, DIMENSION(:)     :: VC
CHARACTER(60) :: FN_DUMP
INTEGER :: LU_DUMP, IC, NM
TYPE (SCARC_LEVEL_TYPE), POINTER :: L

IF (TYPE_DEBUG <= NSCARC_DEBUG_EXTREME) RETURN

DO NM = 1, NMESHES
   L => SCARC(NM)%LEVEL(NL)
   CALL POINT_TO_VECTOR (NVECTOR, NM, NL, VC)
   WRITE (FN_DUMP, '(A,A,i3.3,A,I3.3,A,I3.3,A,I1,A,A,A,A)') &
      'dump/','CG',ITE_CG,'_MG',ITE_MG,'_SM',ITE_SM,'_L',NL,'_',TRIM(SCARC_ROUTINE),'_',TRIM(CNAME)
   LU_DUMP = GET_FILE_NUMBER()
   OPEN (LU_DUMP, FILE=FN_DUMP)
   DO IC = 1, L%NC
      !WRITE(LU_DUMP,'(F25.16)') VC(IC)
      WRITE(LU_DUMP,*) VC(IC)
   ENDDO
   CLOSE(LU_DUMP)
ENDDO

END SUBROUTINE SCARC_DUMP_QUANTITY

!> ------------------------------------------------------------------------------------------------
!> Preset right hand side in such a way that exact solution is known
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_RHS (NVECTOR, NL)
INTEGER, INTENT(IN):: NVECTOR, NL
REAL (EB), POINTER, DIMENSION(:) :: VC
INTEGER :: IC, NM, I, K
REAL (EB) :: X, Z
TYPE (SCARC_LEVEL_TYPE), POINTER :: L
IF (TYPE_DEBUG < NSCARC_DEBUG_EXTREME) RETURN
DO NM = 1, NMESHES
   L => SCARC(NM)%LEVEL(NL)
   CALL POINT_TO_VECTOR (NVECTOR, NM, NL, VC)
   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (.NOT.PRES_ON_WHOLE_DOMAIN.AND.L%DISCRET(I,1,K,CGSC) /= IS_GASPHASE) CYCLE
         IC = L%DISCRET(I,1,K,UNKH)
         X=L%XMID(I)
         Z=L%ZMID(K)
         VC(IC) = RHS(X,Z)
         IF (TYPE_DEBUG > NSCARC_DEBUG_MUCH) &
            WRITE(LU_SCARC,'(A,i3,a,e10.2,a,e10.2,a,e12.4)') 'IC=',IC,':X=',X,':Z=',Z,': RHS=',VC(IC)
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE SCARC_PRESET_RHS


!> ------------------------------------------------------------------------------------------------
!> Set right hand side
!> ------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION RHS(X,Z)
REAL (EB), INTENT(IN) :: X, Z
RHS = 2.0_EB*((1.0_EB - 6.0_EB*X**2)*Z**2*(1.0_EB-Z**2)+(1.0_EB-6.0_EB*Z**2)*X**2*(1.0_EB-X**2))
RETURN 
END FUNCTION RHS


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT1(WORKSPACE, NL1, NR1, CTEXT, BINIT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1
CHARACTER(*), INTENT(IN) :: CTEXT
LOGICAL, INTENT(IN) :: BINIT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR (SCARC_ROUTINE, CTEXT, IERROR)
   IF (BINIT) WORKSPACE = NSCARC_ZERO_INTEGER
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) THEN
      WRITE(SCARC_MESSAGE,'(2A)') 'SCARC_ALLOCATE_INT1: Inconsistent length for vector ',CTEXT
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_INT1

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 2
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT2(WORKSPACE, NL1, NR1, NL2, NR2, CTEXT, BINIT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2
CHARACTER(*), INTENT(IN) :: CTEXT
LOGICAL, INTENT(IN) :: BINIT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR (SCARC_ROUTINE, CTEXT, IERROR)
   IF (BINIT) WORKSPACE = NSCARC_ZERO_INTEGER
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
       SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      WRITE(SCARC_MESSAGE,'(2A)') 'SCARC_ALLOCATE_INT2: Inconsistent length for vector ',CTEXT
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_INT2

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 3
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, CTEXT, BINIT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3
CHARACTER(*), INTENT(IN) :: CTEXT
LOGICAL, INTENT(IN) :: BINIT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR (SCARC_ROUTINE, CTEXT, IERROR)
   IF (BINIT) WORKSPACE = NSCARC_ZERO_INTEGER
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
       SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
       SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      WRITE(SCARC_MESSAGE,'(2A)') 'SCARC_ALLOCATE_INT3: Inconsistent length for vector ',CTEXT
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_INT3

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL1(WORKSPACE, NL1, NR1, CTEXT, BINIT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1
CHARACTER(*), INTENT(IN) :: CTEXT
LOGICAL, INTENT(IN) :: BINIT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR (SCARC_ROUTINE, CTEXT, IERROR)
   IF (BINIT) WORKSPACE = NSCARC_ZERO_REAL
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) THEN
      WRITE(SCARC_MESSAGE,'(2A)') 'SCARC_ALLOCATE_REAL1: Inconsistent length for vector ',CTEXT
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_REAL1

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 2
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL2(WORKSPACE, NL1, NR1, NL2, NR2, CTEXT, BINIT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2
CHARACTER(*), INTENT(IN) :: CTEXT
LOGICAL, INTENT(IN) :: BINIT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR (SCARC_ROUTINE, CTEXT, IERROR)
   IF (BINIT) WORKSPACE = NSCARC_ZERO_REAL
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
       SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      WRITE(SCARC_MESSAGE,'(2A)') 'SCARC_ALLOCATE_REAL2: Inconsistent length for vector ',CTEXT
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_REAL2

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 3
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, CTEXT, BINIT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT):: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3
CHARACTER(*), INTENT(IN) :: CTEXT
LOGICAL, INTENT(IN) :: BINIT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR (SCARC_ROUTINE, CTEXT, IERROR)
   IF (BINIT) WORKSPACE = NSCARC_ZERO_REAL
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
       SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
       SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      WRITE(SCARC_MESSAGE,'(2A)') 'SCARC_ALLOCATE_REAL3: Inconsistent length for vector ',CTEXT
      CALL SHUTDOWN(SCARC_MESSAGE); RETURN
   ENDIF
ENDIF
END SUBROUTINE SCARC_ALLOCATE_REAL3

!> ------------------------------------------------------------------------------------------------
!> Allocate matrix with corresponding pointer and length structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_MATRIX(MAT, NAV, NAC, NAR, NSTENCIL, CMAT)
TYPE (SCARC_MATRIX_TYPE), INTENT(INOUT) :: MAT
INTEGER, INTENT(IN) :: NAV, NAC, NAR, NSTENCIL
CHARACTER(*), INTENT(IN) :: CMAT
CHARACTER(40) :: CINFO

WRITE(CINFO,'(A,A)') TRIM(CMAT),'.VAL'
CALL SCARC_ALLOCATE_REAL1(MAT%VAL, 1, NAV, CINFO, .TRUE.)

WRITE(CINFO,'(A,A)') TRIM(CMAT),'.COL'
CALL SCARC_ALLOCATE_INT1(MAT%COL, 1, NAC, CINFO, .TRUE.)

WRITE(CINFO,'(A,A)') TRIM(CMAT),'.ROW'
CALL SCARC_ALLOCATE_INT1(MAT%ROW, 1, NAR, CINFO, .TRUE.)

MAT%NAV = NAV
MAT%NAC = NAC
MAT%NAR = NAR
MAT%NSTENCIL = NSTENCIL

END SUBROUTINE SCARC_ALLOCATE_MATRIX

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

IF (MAT%NAV > NAV) THEN
   WRITE(CINFO,'(A,A)') TRIM(CMAT),'.VAL'
   CALL SCARC_COPY_MATRIX(MAT, AUX, 'AUX') 
   CALL SCARC_DEALLOCATE_MATRIX(MAT)
   CALL SCARC_COPY_MATRIX(AUX, MAT, 'MAT-copied')
ELSE
   WRITE(SCARC_MESSAGE,'(A)') 'SCARC_REDUCE_MATRIX: Length for resizing too big'
   CALL SHUTDOWN(SCARC_MESSAGE); RETURN
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
CALL SCARC_ALLOCATE_REAL1(MAT2%VAL, 1, MAT2%NAV, CINFO, .FALSE.)
CALL SCARC_COPY_REAL(MAT1%VAL, MAT2%VAL, 1.0_EB, MAT1%NAV)

WRITE(CINFO,'(A,A)') TRIM(CMAT2),'.COL'
CALL SCARC_ALLOCATE_INT1(MAT2%COL, 1, MAT2%NAC, CINFO, .FALSE.)
CALL SCARC_COPY_INT(MAT1%COL, MAT2%COL, 1.0_EB, MAT1%NAC)

WRITE(CINFO,'(A,A)') TRIM(CMAT2),'.ROW'
CALL SCARC_ALLOCATE_INT1(MAT2%ROW, 1, MAT2%NAR, CINFO, .FALSE.)
CALL SCARC_COPY_INT(MAT1%ROW, MAT2%ROW, 1.0_EB, MAT1%NAR)

END SUBROUTINE SCARC_COPY_MATRIX

END MODULE SCRC

!> Remember variable format
!> FORMAT( 3F<2*N+M>.1 )

