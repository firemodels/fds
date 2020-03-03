!> ================================================================================================================
!> Use of different directives possible
!>  - WITH_MKL                : use MKL routines PARDISO, CLUSTER_SPARSE_SOLVER, DDOT, DAXPY, DAXPBY, DCOPY, DSCAL
!>  - WITH_SCARC_VERBOSE      : print more detailed information about ScaRC iterations and workspace allocation
!>  - WITH_SCARC_DEBUG        : print detaild debugging info (only for debugging purposes)
!>  - WITH_SCARC_MGM          : experimental tests for McKeeney-Greengard-Mayo method
!>  - WITH_SCARC_BAROCLINIC   : experimental tests for  modified CG-routine which considers baroclinic effect
!>  - WITH_SCARC_STANDALONE   : dump environment for ScaRC standalone version
!> ================================================================================================================
!#undef WITH_MKL
#undef WITH_SCARC_VERBOSE
#undef WITH_SCARC_DEBUG
#undef WITH_SCARC_MGM
#undef WITH_SCARC_BAROCLINIC
#undef WITH_SCARC_STANDALONE


!> ================================================================================================================
!> MODULE 'SCARC_GLOBAL_CONSTANTS': 
!> Global SCARC parameters
!> ================================================================================================================
MODULE SCARC_GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS
IMPLICIT NONE

INTEGER, PARAMETER :: NSCARC_ACCURACY_ABSOLUTE       =  1, &    !> absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE       =  2       !> relative accuracy must be reached
                   
INTEGER, PARAMETER :: NSCARC_BUFFER_INT0             =  0, &    !> Point to RECV_INT0 buffer
                      NSCARC_BUFFER_INT1             =  1, &    !> Point to RECV_INT1 buffer
                      NSCARC_BUFFER_REAL0            =  2, &    !> Point to RECV_REAL0 buffer
                      NSCARC_BUFFER_REAL1            =  3, &    !> Point to RECV_REAL1 buffer
                      NSCARC_BUFFER_REAL2            =  4, &    !> Point to RECV_REAL2 buffer
                      NSCARC_BUFFER_REAL7            =  5       !> Point to RECV_REAL7 buffer

INTEGER, PARAMETER :: NSCARC_CELL_GASPHASE           =  1, &    !> gasphase cell
                      NSCARC_CELL_SOLID              =  2       !> solid cell

INTEGER, PARAMETER :: NSCARC_COARSE_ITERATIVE        =  1, &    !> iterative coarse grid solution
                      NSCARC_COARSE_DIRECT           =  2       !> direct coarse grid solution

INTEGER, PARAMETER :: NSCARC_COARSENING_BASIC        =  1, &    !> basic coarsening
                      NSCARC_COARSENING_FALGOUT      =  2, &    !> parallel Falgout
                      NSCARC_COARSENING_RS3          =  3, &    !> parallel RS3
                      NSCARC_COARSENING_A1           =  4, &    !> aggressive 1 (path =1, length =2)
                      NSCARC_COARSENING_A2           =  5, &    !> aggressive 2 (path =2, length =2)
                      NSCARC_COARSENING_BDRY         =  6       !> FDSA2  : FDS variant similar to A2

INTEGER, PARAMETER :: NSCARC_COUPLING_MAX            = 10       !> maximum of possible couplings in stencil

INTEGER, PARAMETER :: NSCARC_CYCLING_F               =  0, &    !> F-cycle for mg-method
                      NSCARC_CYCLING_V               =  1, &    !> V-cycle for mg-method
                      NSCARC_CYCLING_W               =  2, &    !> W-cycle for mg-method
                      NSCARC_CYCLING_FMG             =  3, &    !> Full multigrid grid cycle for mg-method
                      NSCARC_CYCLING_SETUP           =  4, &    !> initialize cycle
                      NSCARC_CYCLING_RESET           =  5, &    !> reset cycle counts
                      NSCARC_CYCLING_PROCEED         =  6, &    !> proceed cycle counts
                      NSCARC_CYCLING_NEXT            =  7, &    !> perform next cycling loop
                      NSCARC_CYCLING_EXIT            =  8, &    !> exit cycling loop
                      NSCARC_CYCLING_PRESMOOTH       = -1, &    !> presmoothing cycle
                      NSCARC_CYCLING_POSTSMOOTH      =  1       !> postsmoothing cycle

INTEGER, PARAMETER :: NSCARC_DEBUG_BDRY              =  1, &    !> show pressure_bc_index
                      NSCARC_DEBUG_FACE              =  2, &    !> show face information
                      NSCARC_DEBUG_GRID              =  3, &    !> show grid information
                      NSCARC_DEBUG_MATRIX            =  4, &    !> show matrix
                      NSCARC_DEBUG_MATRIX_SYM        =  5, &    !> show symmetric matrix
                      NSCARC_DEBUG_PRESSURE          =  6, &    !> show pressure quantities
                      NSCARC_DEBUG_STACK             =  7, &    !> show matrix
                      NSCARC_DEBUG_WALL              =  8       !> show wall information

#ifdef WITH_SCARC_STANDALONE
INTEGER, PARAMETER :: NSCARC_DUMP_AB                 =  1, &    !> dump banded AB for a convolutional neural network
                      NSCARC_DUMP_B_NEW              =  2, &    !> dump new B  ...
                      NSCARC_DUMP_B_OLD              =  3, &    !> dump old B  ...
                      NSCARC_DUMP_H_NEW              =  4, &    !> dump new H  ...
                      NSCARC_DUMP_H_OLD              =  5, &    !> dump old H  ...
                      NSCARC_DUMP_HS_NEW             =  6, &    !> dump new HS  ...
                      NSCARC_DUMP_HS_OLD             =  7       !> dump old HS ...
#endif

INTEGER, PARAMETER :: NSCARC_ERROR_AMG_DISABLED      =  1, &    !> AMG-method is currently missing
                      NSCARC_ERROR_BICG_DISABLED     =  2, &    !> BICG temporarily disabled
                      NSCARC_ERROR_BOUNDARY_SUM      =  3, &    !> wrong sum of elements along boundary
                      NSCARC_ERROR_BOUNDARY_TYPE     =  4, &    !> wrong boundary type
                      NSCARC_ERROR_DIRECT_NOMKL      =  5, &    !> wrong boundary type
                      NSCARC_ERROR_EXCHANGE_RECV     =  6, &    !> wrong receive exchange structure
                      NSCARC_ERROR_EXCHANGE_SEND     =  7, &    !> wrong send exchange structure
                      NSCARC_ERROR_FFT_GRID          =  8, &    !> wrong unstructured discretization for FFT
                      NSCARC_ERROR_GRID_INDEX        =  9, &    !> error with grid index
                      NSCARC_ERROR_GRID_NUMBER       = 10, &    !> cell number not divisable by two
                      NSCARC_ERROR_GRID_NUMBERX      = 11, &    !> cell number not divisable by two
                      NSCARC_ERROR_GRID_NUMBERY      = 12, &    !> cell number not divisable by two
                      NSCARC_ERROR_GRID_NUMBERZ      = 13, &    !> cell number not divisable by two
                      NSCARC_ERROR_GRID_RESOLUTION   = 14, &    !> error with grid resolution
                      NSCARC_ERROR_NEIGHBOR_NUMBER   = 15, &    !> wrong neighbor number
                      NSCARC_ERROR_NEIGHBOR_TYPE     = 16, &    !> wrong neighbor type
                      NSCARC_ERROR_MATRIX_ALLOCATION = 17, &    !> error in matrix allocation
                      NSCARC_ERROR_MATRIX_COPY       = 18, &    !> subdiagonal missing
                      NSCARC_ERROR_MATRIX_SETUP      = 19, &    !> error in matrix setup
                      NSCARC_ERROR_MATRIX_SIZE       = 20, &    !> error in matrix size
                      NSCARC_ERROR_MATRIX_SUBDIAG    = 21, &    !> subdiagonal missing
                      NSCARC_ERROR_MATRIX_SYMMETRY   = 22, &    !> matrix not symmetric
                      NSCARC_ERROR_MKL_CLUSTER       = 23, &    !> cluster_sparse_solver not available
                      NSCARC_ERROR_MKL_INTERNAL      = 24, &    !> pardiso solver not available
                      NSCARC_ERROR_MKL_PARDISO       = 25, &    !> pardiso solver not available
                      NSCARC_ERROR_MKL_STORAGE       = 26, &    !> wrong storage scheme for MKL-solvers
                      NSCARC_ERROR_MULTIGRID_LEVEL   = 27, &    !> wrong multigrid level
                      NSCARC_ERROR_PARSE_INPUT       = 28, &    !> wrong input parameter
                      NSCARC_ERROR_POC_STOP          = 29, &    !> temporary stop, only proof of concept
                      NSCARC_ERROR_STACK_MESSAGE     = 30, &    !> error with stack message
                      NSCARC_ERROR_STACK_SOLVER      = 31, &    !> error in solver stack
                      NSCARC_ERROR_STENCIL           = 32, &    !> error in matrix stencil
                      NSCARC_ERROR_RHS_SETUP         = 33, &    !> error with rhs setup
                      NSCARC_ERROR_VECTOR_LENGTH     = 34       !> error in vector length

INTEGER, PARAMETER :: NSCARC_EXCHANGE_FDS            =  1, &    !> use exchange information from surrounding FDS code
                      NSCARC_EXCHANGE_SCARC          =  2, &    !> use own ScaRC related exchange information
                      NSCARC_EXCHANGE_SETUP          =  3, &    !> initialize wall information
                      NSCARC_EXCHANGE_VECTOR         =  4, &    !> exchange vector for mat-vec multiplication
                      NSCARC_EXCHANGE_VECTOR_MEAN    =  5, &    !> exchange vector for mat-vec multiplication and build sum
                      NSCARC_EXCHANGE_PRESSURE       =  6, &    !> vector values along internal boundaries
                      NSCARC_EXCHANGE_CELL_NUMBER    =  7, &    !> neighboring cell numbers on overlap
                      NSCARC_EXCHANGE_CELL_WIDTH     =  8, &    !> neighboring grid widths on overlap
                      NSCARC_EXCHANGE_MATRIX_SIZE    =  9, &    !> size of neighboring matrix 
                      NSCARC_EXCHANGE_MATRIX_SUBDIAG = 10, &    !> subdiag value of neighboring matrix 
                      NSCARC_EXCHANGE_MATRIX_STENCIL = 11       !> stencil of neighboring matrix 

INTEGER, PARAMETER :: NSCARC_GRID_STRUCTURED         =  1, &    !> structured discretization
                      NSCARC_GRID_UNSTRUCTURED       =  2       !> unstructured discretization

INTEGER, PARAMETER :: NSCARC_HUGE_INT                = -999999999    !> undefined integer value
REAL(EB), PARAMETER:: NSCARC_HUGE_REAL               = -999999999_EB !> undefined integer value

INTEGER, PARAMETER :: NSCARC_INIT_UNDEF              = -99, &   !> initialize allocated array as undefined
                      NSCARC_INIT_NONE               =  -1, &   !> do not initialize allocated arrays
                      NSCARC_INIT_ZERO               =   0, &   !> initialize allocated array with zero
                      NSCARC_INIT_ONE                =   1, &   !> initialize allocated array with one
                      NSCARC_INIT_TRUE               =   2, &   !> initialize allocated array with one
                      NSCARC_INIT_FALSE              =   3      !> initialize allocated array with one

INTEGER, PARAMETER :: NSCARC_INTERPOL_STANDARD       =  1, &    !> standard interpolation
                      NSCARC_INTERPOL_CONSTANT       =  2, &    !> standard interpolation
                      NSCARC_INTERPOL_BILINEAR       =  3, &    !> standard interpolation
                      NSCARC_INTERPOL_CLASSICAL      =  4, &    !> classical interpolation
                      NSCARC_INTERPOL_CLASSICAL2     =  5, &    !> classical interpolation
                      NSCARC_INTERPOL_DIRECT         =  6, &    !> direct interpolation
                      NSCARC_INTERPOL_DIRECT_BDRY    =  7       !> direct interpolation with special boundary

INTEGER, PARAMETER :: NSCARC_KRYLOV_CG               =  1, &    !> CG   as Krylov solver
                      NSCARC_KRYLOV_CGBARO           =  2, &    !> CGBARO (Krylov with baroclinic effect), only testwise
                      NSCARC_KRYLOV_BICG             =  3, &    !> BICG as Krylov solver
                      NSCARC_KRYLOV_MAIN             =  1, &    !> Krylov solver as main solver
                      NSCARC_KRYLOV_COARSE           =  2       !> Krylov solver as coarse grid solver

INTEGER, PARAMETER :: NSCARC_LEVEL_MIN               =  0, &    !> minimum multigrid level
                      NSCARC_LEVEL_MAX               = 10, &    !> maximum multigrid level
                      NSCARC_LEVEL_SINGLE            =  1, &    !> only one grid level needed
                      NSCARC_LEVEL_MULTI             =  2, &    !> multiple grid levels needed for GMG meethod
                      NSCARC_LEVEL_AMG               =  3       !> number of grid levels specified by AMG method

INTEGER, PARAMETER :: NSCARC_MATRIX_COMPACT          =  1, &    !> matrix in compact  storage format
                      NSCARC_MATRIX_BANDED           =  2, &    !> matrix in bandwise storage format
                      NSCARC_MATRIX_CONDENSED        =  6, &    !> matrix condensing for Neumann problems
                      NSCARC_MATRIX_GENERAL          =  7, &    !> general matrix (not necessarily symmetric)
                      NSCARC_MATRIX_SYMMETRIC        =  8       !> symmetric matrix

INTEGER, PARAMETER :: NSCARC_MATVEC_GLOBAL           =  1, &    !> perform matrix vector product globally
                      NSCARC_MATVEC_LOCAL            =  2, &    !> perform matrix vector product locally
                      NSCARC_MATVEC_MEAN             =  3       !> perform matrix vector product with mean values of overlaps

INTEGER, PARAMETER :: NSCARC_MAX_FACE_NEIGHBORS      = 10       !> max number neighbors per mesh face
INTEGER, PARAMETER :: NSCARC_MAX_STENCIL             =  7

INTEGER, PARAMETER :: NSCARC_METHOD_KRYLOV           =  1, &    !> Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID        =  2, &    !> multigrid-method as global solver
                      NSCARC_METHOD_LU               =  3, &    !> LU-decomposition based on MKL solver
                      NSCARC_METHOD_FFT              =  4, &    !> FFT method based on CRAYFISHPAK
                      NSCARC_METHOD_MGM              =  5       !> MGM-method based on LU

INTEGER, PARAMETER :: NSCARC_MKL_LOCAL               =  1, &    !> local LU-decompositions (PARDISO)
                      NSCARC_MKL_GLOBAL              =  2, &    !> global LU-decomposition(CLUSTER_SPARSE_SOLVER)
                      NSCARC_MKL_COARSE              =  3       !> LU-decomposition on coarse grid level

INTEGER, PARAMETER :: NSCARC_MULTIGRID_GEOMETRIC     =  1, &    !> geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC     =  2, &    !> algebraic multigrid
                      NSCARC_MULTIGRID_MAIN          =  1, &    !> multigrid solver as main solver
                      NSCARC_MULTIGRID_PRECON        =  2       !> multigrid solver as preconditioner

INTEGER, PARAMETER :: NSCARC_PRECISION_SINGLE        =  1, &    !> single precision data 
                      NSCARC_PRECISION_DOUBLE        =  2       !> double precision data

INTEGER, PARAMETER :: NSCARC_RELAX_JAC               =  1, &    !> preconditioning by local JACOBI-methods
                      NSCARC_RELAX_SSOR              =  2, &    !> preconditioning by local SSOR-methods
                      NSCARC_RELAX_MJAC              =  3, &    !> preconditioning by local MJAC-methods (matrix form)
                      NSCARC_RELAX_MGS               =  4, &    !> preconditioning by local MGS-methods (matrix form)
                      NSCARC_RELAX_MSGS              =  5, &    !> preconditioning by local MSGS-methods (matrix form)
                      NSCARC_RELAX_MSOR              =  6, &    !> preconditioning by local MSOR-methods (matrix form)
                      NSCARC_RELAX_MSSOR             =  7, &    !> preconditioning by local MSSOR-methods (matrix form)
                      NSCARC_RELAX_ILU               =  8, &    !> preconditioning by local ILU-decompositions (own)
                      NSCARC_RELAX_FFT               =  9, &    !> preconditioning by local FFT-methods
                      NSCARC_RELAX_FFTO              = 10, &    !> preconditioning by local FFT-methods
                      NSCARC_RELAX_GMG               = 11, &    !> preconditioning by local GMG-methods
                      NSCARC_RELAX_MKL               = 12       !> preconditioning by local LU-decompositions (MKL)

INTEGER, PARAMETER :: NSCARC_RHS_HOMOGENEOUS         =  1, &    !> homogeneous boundary conditions for global problem
                      NSCARC_RHS_INHOMOGENEOUS       =  2, &    !> inhomogeneous boundary conditions for global problem
                      NSCARC_RHS_DEFECT              =  3       !> RHS of precond. is set to defect of main iteration

INTEGER, PARAMETER :: NSCARC_SCOPE_GLOBAL            =  0, &    !> scope of defect correction is global
                      NSCARC_SCOPE_LOCAL             =  1       !> scope of defect correction is local

INTEGER, PARAMETER :: NSCARC_SIZE_MATRIX             =  1, &    !> size of system matrix for compact system
                      NSCARC_SIZE_TRANSFER           =  2       !> size of transfer matrices for compact system

INTEGER, PARAMETER :: NSCARC_SOLVER_MAIN             =  1, &    !> Main solver
                      NSCARC_SOLVER_PRECON           =  2, &    !> Preconditioner
                      NSCARC_SOLVER_SMOOTH           =  3, &    !> Smoother
                      NSCARC_SOLVER_COARSE           =  4       !> Coarse grid solver

INTEGER, PARAMETER :: NSCARC_STACK_ZERO              =   0, &   !> zero stack
                      NSCARC_STACK_ROOT              =   1, &   !> root stack
                      NSCARC_STACK_MAX               =  10, &   !> maximum stack
                      NSCARC_STACK_NOPARENT          = -99      !> no stack information

INTEGER, PARAMETER :: NSCARC_STAGE_ONE               =  1, &    !> primary scope for solution vectors
                      NSCARC_STAGE_TWO               =  2       !> secondary scope for solution vectors

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED           =  0, &    !> proceed loop
                      NSCARC_STATE_CONV0             =  1, &    !> check convergence already for initial residual
                      NSCARC_STATE_CONV              =  2, &    !> check convergence for residual
                      NSCARC_STATE_DIVG              =  3       !> check divergence for residual

INTEGER, PARAMETER :: NSCARC_STENCIL_CONSTANT        =  1, &    !> constant matrix entries
                      NSCARC_STENCIL_VARIABLE        =  2, &    !> variable matrix entries, own version of scarc_daypyv
                      NSCARC_STENCIL_CONSTANT_LOCAL  =  3, &    !> constant-matrix entries with local corrections
                      NSCARC_STENCIL_VARIABLE_LOCAL  =  4       !> variable matrix entries with local matvec

REAL(EB), PARAMETER:: NSCARC_THRESHOLD_CONVERGENCE   = 1.0E-15_EB, & !> threshold for convergence
                      NSCARC_THRESHOLD_DIVGERGENCE   = 1.0E+15_EB    !> threshold for divergence

INTEGER, PARAMETER :: NSCARC_TWOLEVEL_NONE           =  0, &    !> no two levels, only one level
                      NSCARC_TWOLEVEL_ADD            =  1, &    !> additive 2-level method
                      NSCARC_TWOLEVEL_MUL            =  2, &    !> multiplicative 2-level method
                      NSCARC_TWOLEVEL_MUL2           =  3, &    !> multiplicative 2-level method
                      NSCARC_TWOLEVEL_COARSE         =  4       !> only coarse grid

INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_X            =  1, &    !> selection parameter for 1-stage vector X
                      NSCARC_VECTOR_ONE_B            =  2, &    !>       ...                              F
                      NSCARC_VECTOR_ONE_D            =  3, &    !>       ...                              G
                      NSCARC_VECTOR_ONE_R            =  4, &    !>       ...                              D
                      NSCARC_VECTOR_ONE_V            =  5, &    !>       ...                              R
                      NSCARC_VECTOR_ONE_Y            =  6, &    !>       ...                              Y
                      NSCARC_VECTOR_ONE_Z            =  7, &    !>       ...                              Z
                      NSCARC_VECTOR_ONE_E            =  8, &    !>       ...                              Z
                      NSCARC_VECTOR_TWO_X            =  9, &    !> selection parameter for 2-stage vector X
                      NSCARC_VECTOR_TWO_B            = 10, &    !>       ...                              F
                      NSCARC_VECTOR_TWO_D            = 11, &    !>       ...                              G
                      NSCARC_VECTOR_TWO_R            = 12, &    !>       ...                              D
                      NSCARC_VECTOR_TWO_V            = 13, &    !>       ...                              R
                      NSCARC_VECTOR_TWO_Y            = 14, &    !>       ...                              Y
                      NSCARC_VECTOR_TWO_Z            = 15, &    !>       ...                              Z
                      NSCARC_VECTOR_TWO_E            = 16, &    !>       ...                              Z
                      NSCARC_VECTOR_H                = 17, &    !> selection parameter for vector H
                      NSCARC_VECTOR_HS               = 18       !>       ...                      HS

INTEGER, PARAMETER :: NSCARC_UNDEF_INT               = -1, &    !> undefined integer value
                      NSCARC_ZERO_INT                =  0, &    !> zero integer value
                      NSCARC_ONE_INT                 =  1       !> one integer value

REAL(EB), PARAMETER:: NSCARC_UNDEF_REAL_EB           = -1.0_EB, &    !> undefined real value
                      NSCARC_ZERO_REAL_EB            =  0.0_EB, &    !> zero real value
                      NSCARC_ONE_REAL_EB             =  1.0_EB       !> one real value

REAL(EB), PARAMETER:: NSCARC_UNDEF_REAL_FB           = -1.0_FB, &    !> undefined real value
                      NSCARC_ZERO_REAL_FB            =  0.0_FB, &    !> zero real value
                      NSCARC_ONE_REAL_FB             =  1.0_FB       !> one real value

INTEGER, PARAMETER :: NSCARC_NONE                    = -123456789    !> dummy integer value for requests
CHARACTER(40), PARAMETER :: SCARC_NONE = 'NONE'                      !> dummy character value for requests

END MODULE SCARC_GLOBAL_CONSTANTS


!> ================================================================================================================
!> MODULE 'SCARC_VARIABLES': 
!> Global SCARC variables 
!> ================================================================================================================
MODULE SCARC_VARIABLES
USE PRECISION_PARAMETERS
USE SCARC_GLOBAL_CONSTANTS
IMPLICIT NONE

!>
!> ---------- Basic definitions
!>
CHARACTER(40) :: SCARC_GRID              = 'STRUCTURED'     !> Type of discretization (STRUCTURED/UNSTRUCTURED)
CHARACTER(40) :: SCARC_METHOD            = 'NONE'           !> Requested solver method (KRYLOV/MULTIGRID)
CHARACTER(40) :: SCARC_MATRIX            = 'NONE'           !> Type of matrix storage (COMPACT/BANDED)
CHARACTER(40) :: SCARC_STENCIL           = 'VARIABLE'       !> Type of matrix stencil (CONSTANT/VARIABLE)
CHARACTER(40) :: SCARC_TWOLEVEL          = 'NONE'           !> Type of two-level method (NONE/ADDITIVE/MULTIPLICATIVE)

!>
!> ---------- General iteration parameters
!>
CHARACTER(40) :: SCARC_ACCURACY          = 'ABSOLUTE'       !> Accuracy type (ABSOLUTE/RELATIVE)
REAL (EB)     :: SCARC_CAPPA             =  0.0_EB          !> Convergence rate of selected ScarC solver
INTEGER       :: SCARC_ITERATIONS        =  0               !> Number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_RESIDUAL          =  0.0_EB          !> Residual of global selected solver

!>
!> ---------- Parameters for coarse grid method
!>
CHARACTER(40) :: SCARC_COARSE            = 'DIRECT'        !> Type of coarse grid solver (ITERATIVE/DIRECT)
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-14_EB        !> Requested accuracy for iterative solver
INTEGER       :: SCARC_COARSE_ITERATIONS = 100              !> Max number of iterations for iterative solver
INTEGER       :: SCARC_COARSE_LEVEL      =  1               !> Coarse grid level for twolevel-Krylov method (default minimum level)
REAL (EB)     :: SCARC_COARSE_OMEGA      = 0.80E+0_EB       !> Relaxation parameter

!>
!> ---------- Parameters for Krylov-type methods
!>
CHARACTER(40) :: SCARC_KRYLOV            = 'CG'             !> Type of Krylov-method (CG/BICG)
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-8_EB         !> Requested accuracy for convergence
CHARACTER(40) :: SCARC_KRYLOV_INTERPOL   = 'CONSTANT'       !> twolevel-interpolation (CONSTANT/BILINEAR)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 1000             !> Max number of iterations

!>
!> ---------- Parameters for multigrid-type methods
!>
CHARACTER(40) :: SCARC_MULTIGRID            = 'GEOMETRIC'   !> Type of MG-method (GEOMETRIC/ALGEBRAIC)
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY   = 1.E-8_EB      !> Requested accuracy for convergence
CHARACTER(40) :: SCARC_MULTIGRID_COARSENING = 'FALGOUT'     !> Coarsening strategy  (Falgout/RS3/A1/A2/...)
CHARACTER(3)  :: SCARC_MULTIGRID_CYCLE      = 'V'           !> Cycling type  (F/V/W)
CHARACTER(40) :: SCARC_MULTIGRID_INTERPOL   = 'CONSTANT'    !> Interpolation strategy (CONSTANT/BILINEAR)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS = 100           !> Max number of iterations
INTEGER       :: SCARC_MULTIGRID_PRESMOOTH  = 8             !> Number of presmoothing iterations
INTEGER       :: SCARC_MULTIGRID_POSTSMOOTH = 8             !> Number of postsmoothing iterations
INTEGER       :: SCARC_MULTIGRID_LEVEL      = -1            !> User defined number of MG-levels (optionally, otherwise maximum)

!>
!> ---------- Parameters for smoothing method (used in multigrids-methods)
!>
CHARACTER(40) :: SCARC_SMOOTH            = 'SSOR'           !> Smoother for MG (JACOBI/SSOR)
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-8_EB         !> Requested accuracy for convergence
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 5                !> Max number of iterations
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.80E+0_EB       !> Relaxation parameter
CHARACTER(40) :: SCARC_SMOOTH_SCOPE      = 'GLOBAL'         !> Scope of action (LOCAL/GLOBAL)

!>
!> ---------- Parameters for preconditioning method (used in Krylov-methods)
!>
CHARACTER(40) :: SCARC_PRECON            = 'NONE'           !> Preconditioner for CG/BICG (JACOBI/SSOR/FFT/PARDISO/MG)
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-10_EB        !> Requested accuracy for convergence
INTEGER       :: SCARC_PRECON_ITERATIONS = 100              !> Max number of iterations
REAL (EB)     :: SCARC_PRECON_OMEGA      = 1.50E+0_EB       !> Relaxation parameter
CHARACTER(40) :: SCARC_PRECON_SCOPE      = 'LOCAL'          !> Scope of action (LOCAL/GLOBAL)

!>
!> ---------- Parameter for MKL solver
!>
CHARACTER(40) :: SCARC_MKL            = 'GLOBAL'            !> Type of MKL solver (LOCAL:Pardiso/GLOBAL:Cluster_Sparse_solver)
CHARACTER(40) :: SCARC_MKL_MTYPE      = 'SYMMETRIC'         !> Type of MKL matrix (SYMMETRIC/UNSYMMETRIC)
CHARACTER(6)  :: SCARC_MKL_PRECISION  = 'DOUBLE'            !> Single/double precision for MKL solver

!>
!> ---------- Dump out of error information and error handling
!>
LOGICAL :: SCARC_ERROR_FILE = .FALSE.                       !> Print ScaRC statistics into chid_scarc.csv (TRUE/FALSE)
INTEGER :: IERROR = 0                                       !> general error flag - used at different positions


!> 
!> ---------- Logical indicators for different methods and mechanisms
!> 
LOGICAL :: IS_STRUCTURED      = .FALSE.                     !> is structured discretization?
LOGICAL :: IS_UNSTRUCTURED    = .FALSE.                     !> is unstructured discretization?
LOGICAL :: IS_PURE_NEUMANN    = .FALSE.                     !> pure Neumann system
LOGICAL :: IS_CG              = .FALSE.                     !> is Krylov-method?
LOGICAL :: IS_CG_GMG          = .FALSE.                     !> is Krylov-method with GMG-preconditioning?
LOGICAL :: IS_CG_ADD          = .FALSE.                     !> is additive twolevel-Krylov-method?
LOGICAL :: IS_CG_MUL          = .FALSE.                     !> is multiplicative Twolevel-Krylov-method ?
LOGICAL :: IS_CG_COARSE       = .FALSE.                     !> is only coarse grid preconditiner?
LOGICAL :: IS_MG              = .FALSE.                     !> is Multigrid-method?
LOGICAL :: IS_GMG             = .FALSE.                     !> is Geometric Multigrid-method?
LOGICAL :: IS_FFT             = .FALSE.                     !> is FFT-method?
LOGICAL :: IS_FFTO            = .FALSE.                     !> is FFTO-method?
LOGICAL :: IS_MGM             = .FALSE.                     !> is McKeeney-Greengard-Mayo method?
LOGICAL :: IS_MKL             = .FALSE.                     !> is MKL-method?
LOGICAL :: IS_MKL_LEVEL(10)   = .FALSE.                     !> is level-dependent MKL method?
 
LOGICAL :: HAS_CSV_DUMP       = .FALSE.                     !> has CSV-file to be dumped out
LOGICAL :: HAS_TWO_LEVELS     = .FALSE.                     !> has two grid levels?
LOGICAL :: HAS_MULTI_LEVELS   = .FALSE.                     !> has multiple grid levels?
LOGICAL :: HAS_MULTI_GRIDS    = .FALSE.                     !> has multiple discretization types?

!> 
!> ---------- Globally used types for description of different solvers
!> 
INTEGER :: TYPE_ACCURACY      = NSCARC_ACCURACY_ABSOLUTE    !> default type of requested accuracy
INTEGER :: TYPE_COARSE        = NSCARC_COARSE_ITERATIVE     !> default type of coarse grid solver for multigrid method
INTEGER :: TYPE_COARSENING    = NSCARC_UNDEF_INT            !> no default type of coarsening algorithm for AMG
INTEGER :: TYPE_CYCLING       = NSCARC_CYCLING_V            !> default type of cycling for multigrid method
INTEGER :: TYPE_GRID          = NSCARC_GRID_STRUCTURED      !> default type of discretization (structured/unstructured)
INTEGER :: TYPE_EXCHANGE      = NSCARC_UNDEF_INT            !> no default type of data exchange
INTEGER :: TYPE_INTERPOL      = NSCARC_INTERPOL_CONSTANT    !> default type of interpolation method
INTEGER :: TYPE_KRYLOV        = NSCARC_KRYLOV_CG            !> default type of Krylov method (CG/BICG)
INTEGER :: TYPE_LEVEL(0:2)    = NSCARC_UNDEF_INT            !> default type of levels
INTEGER :: TYPE_MATRIX        = NSCARC_MATRIX_COMPACT       !> default type of storage for matrix
INTEGER :: TYPE_MATVEC        = NSCARC_MATVEC_GLOBAL        !> default type of matrix-vector multiplication
INTEGER :: TYPE_METHOD        = NSCARC_METHOD_KRYLOV        !> default type of ScaRC method
INTEGER :: TYPE_MKL(0:10)     = NSCARC_UNDEF_INT            !> no default type of MKL for single levels
INTEGER :: TYPE_MKL_PRECISION = NSCARC_PRECISION_DOUBLE     !> default double precision MKL solver
INTEGER :: TYPE_MULTIGRID     = NSCARC_MULTIGRID_GEOMETRIC  !> default type of multigrid method (GMG/AMG)
INTEGER :: TYPE_PARENT        = NSCARC_UNDEF_INT            !> no default type of parent (calling) solver
INTEGER :: TYPE_PRECON        = NSCARC_UNDEF_INT            !> no default type of preconditioner for iterative solver
INTEGER :: TYPE_RELAX         = NSCARC_UNDEF_INT            !> no default type of preconditioner for iterative solver
INTEGER :: TYPE_SCOPE(0:2)    = NSCARC_SCOPE_LOCAL          !> default scope types are local
INTEGER :: TYPE_SMOOTH        = NSCARC_UNDEF_INT            !> no default type of smoother for multigrid method
INTEGER :: TYPE_SOLVER        = NSCARC_SOLVER_MAIN          !> default type of surrounding solver stage
INTEGER :: TYPE_STAGE         = NSCARC_STAGE_ONE            !> default type of surrounding solver stage
INTEGER :: TYPE_STENCIL       = NSCARC_STENCIL_VARIABLE     !> default type of storage for matrix
INTEGER :: TYPE_TWOLEVEL      = NSCARC_TWOLEVEL_NONE        !> default type of two-level method
INTEGER :: TYPE_VECTOR        = NSCARC_UNDEF_INT            !> no default type of vector to point to

!> 
!> ---------- Globally used parameters
!> 
INTEGER :: NLEVEL_MIN, NLEVEL_MAX                           !> minimum and maximum number of multigrid levels
INTEGER :: NC_GLOBAL(20) = 0                                !> number of global cells
INTEGER :: N_DIRIC_GLOBAL(20) = 0                           !> global number of Dirichlet BCs
INTEGER :: N_STACK_TOTAL                                    !> maximum number of used solvers in stack

INTEGER :: N_REQ, N_EXCHANGES, TAG                           !> Information for data exchange
INTEGER :: SNODE, RNODE                                     !> Process identifier for data exchange

INTEGER,  ALLOCATABLE, DIMENSION (:)  :: REQ                !> Request array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: COUNTS             !> Counter array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: DISPLS             !> Displacement array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: MESH_INT           !> Local integer data array for data exchange
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: MESH_REAL          !> Local real data array for data exchange

INTEGER  :: GLOBAL_INT,  RANK_INT
REAL(EB) :: GLOBAL_REAL, RANK_REAL

INTEGER :: FACE_ORIENTATION(6) = (/1,-1,2,-2,3,-3/)         !> Coordinate direction related order of mesh faces

!> 
!> ---------- Iteration parameters vectors and iteration parameters
!> 
REAL (EB) :: DT, DTI
REAL (EB) :: OMEGA, EPS, RES, RESIN = -1.0_EB, CAPPA = -1.0_EB, ERR
INTEGER :: X, B, D, R, V, Y, Z
#ifdef WITH_SCARC_DEBUG
INTEGER :: E
#endif
INTEGER :: NIT, ITE
INTEGER :: ITE_PRES=0, ITE_TOTAL=0, ITE_CG=0, ITE_MG=0, ITE_LU=0, ITE_SMOOTH=0, ITE_COARSE=0, ITE_GLOBAL = 0
CHARACTER(60) :: CNAME

!> 
!> ---------- Grid transfer information
!> 
REAL(EB), PARAMETER :: SCALR  = 0.015625_EB                 !> some scaling parameters
REAL(EB), PARAMETER :: SCALP  = 0.0625_EB                    
REAL(EB), PARAMETER :: W1     =  1.0_EB                     !> some weighting parameters
REAL(EB), PARAMETER :: W3     =  3.0_EB                      
REAL(EB), PARAMETER :: W4     =  4.0_EB                     
REAL(EB), PARAMETER :: W9     =  9.0_EB                    
REAL(EB), PARAMETER :: W12    = 12.0_EB                   
REAL(EB), PARAMETER :: W16    = 16.0_EB                  


END MODULE SCARC_VARIABLES


!> ================================================================================================================
!> MODULE 'SCARC_TYPES': 
!> Collection of different types for ScaRC
!> ================================================================================================================
MODULE SCARC_TYPES
USE PRECISION_PARAMETERS
USE SCARC_GLOBAL_CONSTANTS
#ifdef WITH_MKL
USE MKL_PARDISO
USE MKL_CLUSTER_SPARSE_SOLVER
#endif
IMPLICIT NONE

!>
!> ---------- Messaging and debugging mechanisms
!>
TYPE SCARC_MESSAGE_TYPE

   CHARACTER(100) :: TEXT
   CHARACTER(60)  :: FILE_DEBUG, FILE_STAT, FILE_DUMP, FILE_VERBOSE, FILE_CPU, FILE_CNN
   CHARACTER(120) :: FILE_SCARC

   INTEGER :: LU_DEBUG = 0, LU_STAT = 0, LU_DUMP = 0, LU_VERBOSE = 0, LU_SCARC = 0, LU_CPU = 0, LU_CNN = 0

END TYPE SCARC_MESSAGE_TYPE


!>
!> ---------- Face information related to wall cells and neighbors
!>
TYPE SCARC_FACE_TYPE

   REAL(EB), POINTER, DIMENSION(:) :: DH              !> step size vector between adjacent faces
   REAL(EB) :: SCAL_FACE                              !> increments for matrix subdiagonal for cells right face
   REAL(EB) :: SCAL_INSIDE                            !> increments for matrix subdiagonal for cells between opposite faces
   REAL(EB) :: SCAL_BOUNDARY                          !> increments for boundary conditions
   REAL(EB) :: SCAL_DIRICHLET                         !> scaling factors for Dirichlet BC's
   REAL(EB) :: SCAL_NEUMANN                           !> scaling factors for Neumann BC's

   INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS    !> adjacent neighbors at that face
   INTEGER  :: N_NEIGHBORS = 0                        !> number of adjacent neighbors 

   INTEGER  :: NOP = 0                                !> number of cells between opposite faces
   INTEGER  :: NX, NY, NZ                             !> cells in different directions on that face
   INTEGER  :: NCW0 = 0, NCW = 0                      !> number of first wall cell and total number of wall cells

END TYPE SCARC_FACE_TYPE


!>
!> ---------- Wall information related to neighbors and BC's
!>
TYPE SCARC_WALL_TYPE

   !> different properties of wall cell
   INTEGER :: BTYPE = 0                               !> BTYPE of wall cell (Dirichlet/Neumann/Internal)
   INTEGER :: BOUNDARY_TYPE = 0                       !> boundary type of wall cell (Solid/Interpolated/Open))
   INTEGER :: IOR = 0                                 !> orientation of wall cell
   INTEGER :: NOM = 0                                 !> adjacent neighbor at wall cell

   INTEGER :: ICW = NSCARC_UNDEF_INT                  !> internal wall cell for IW
   INTEGER :: IXG, IYG, IZG                           !> x-, y- and z-indices of ghost cells
   INTEGER :: IXW, IYW, IZW                           !> x-, y- and z-indices of (internal) wall cells
   INTEGER :: IXN(2), IYN(2), IZN(2)                  !> x-, y- and z-indices of neighboring cells
   INTEGER :: ICE, ICG

END TYPE SCARC_WALL_TYPE


!>
!> ---------- Obstruction information
!>
TYPE SCARC_OBST_TYPE
   INTEGER :: I1, I2, J1, J2, K1, K2                         !> cell indices of obstructions
END TYPE SCARC_OBST_TYPE


!>
!> ---------- Compact matrix entries which will be exchanged during generation of condensed system
!>
TYPE SCARC_MATRIX_COMPACT_CONDENSED_TYPE

   REAL(EB) :: VAL1(NSCARC_MAX_STENCIL) = 0.0_EB         !> original values (double precision)
   REAL(EB) :: VAL2(NSCARC_MAX_STENCIL) = 0.0_EB         !> condensed values (double precision)

   INTEGER :: COL(NSCARC_MAX_STENCIL) = 0                !> column pointers
   INTEGER :: PTR(NSCARC_MAX_STENCIL) = 0                !> storage pointer
   INTEGER :: N_COL

END TYPE SCARC_MATRIX_COMPACT_CONDENSED_TYPE


!>
!> ---------- Banded matrix entries which will exchanged during generation of condensed system
!>
TYPE SCARC_MATRIX_BANDED_CONDENSED_TYPE

   REAL(EB) :: VAL1(NSCARC_MAX_STENCIL) = 0.0_EB         !> original values (double precision)
   REAL(EB) :: VAL2(NSCARC_MAX_STENCIL) = 0.0_EB         !> condensed values (double precision)
   INTEGER :: IOR0 = 0                                   !> position pointer
   INTEGER :: ICO = 0                                    !> cell pointer

END TYPE SCARC_MATRIX_BANDED_CONDENSED_TYPE


!>
!> ---------- Compact sparse row (COMPACT) storage technique for matrices
!> Is based on three arrays:
!>    - non-zero matrix values
!>    - corresponding columns pointers
!>    - row pointers
!>
TYPE SCARC_MATRIX_COMPACT_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:) :: VAL                !> values of matrix (real precision)
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: ILU                !> ILU-decomposition
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: RELAX              !> workspace for relaxation
   REAL(EB), DIMENSION (-3:3)           :: STENCIL            !> store basic stencil information in single precision

   REAL(FB), ALLOCATABLE, DIMENSION (:) :: VAL_FB             !> values of matrix (single precision)
   REAL(FB), ALLOCATABLE, DIMENSION (:) :: RELAX_FB           !> workspace for relaxation
   REAL(FB), DIMENSION (-3:3)           :: STENCIL_FB         !> store basic stencil information in single precision

   TYPE (SCARC_MATRIX_COMPACT_CONDENSED_TYPE) :: CONDENSED(NSCARC_MAX_STENCIL)

   INTEGER,  ALLOCATABLE, DIMENSION (:) :: ROW                !> row pointer
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: COL                !> column pointers
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: COL_GLOBAL         !> column pointer for global numbering

   INTEGER :: POS(-3:3) = 0                                   !> Position of IOR's in STENCIL
   INTEGER :: N_STENCIL = 0                                   !> number of points in matrix stencil
   INTEGER :: N_CONDENSED = 0                                 !> number of condensed entries in matrix
   INTEGER :: N_VAL = 0                                       !> number of matrix values
   INTEGER :: N_ROW = 0                                       !> number of matrix rows

   CHARACTER(40) :: CNAME                                     !> Name of matrix

END TYPE SCARC_MATRIX_COMPACT_TYPE

!> 
!> ---------- Banded storage technique for matrices
!> The entries are stored one diagonal after the other
!> Missing entries of subdiagonals are filled with zero
!> Is based on two arrays:
!>      - non-zero matrix entries diagonal-wise
!>      - the offsets from the main diagonal
!>
TYPE SCARC_MATRIX_BANDED_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: AUX           !> auxiliary vector (double precision)
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: VAL           !> values of matrix (double precision)
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: RELAX         !> workspace for relaxation
   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: RELAXD        !> workspace for relaxation - only for diagonal scaling
   REAL(EB), DIMENSION (-3:3)             :: STENCIL       !> store basic stencil information (double precision)

   REAL(FB), ALLOCATABLE, DIMENSION (:)   :: AUX_FB        !> auxiliary vector (single precision)
   REAL(FB), ALLOCATABLE, DIMENSION (:,:) :: VAL_FB        !> values of matrix (single precision)
   REAL(FB), ALLOCATABLE, DIMENSION (:,:) :: RELAX_FB      !> workspace for relaxation
   REAL(FB), DIMENSION (-3:3)             :: STENCIL_FB    !> store basic stencil information (single precision)

   TYPE (SCARC_MATRIX_BANDED_CONDENSED_TYPE) :: CONDENSED(NSCARC_MAX_STENCIL)

   CHARACTER(40) :: CNAME                                  !> Name of matrix

   INTEGER,  DIMENSION (-3:3) :: OFFSET                    !> offset pointers
   INTEGER,  DIMENSION (-3:3) :: LENGTH                    !> Relevant diagonal length 
   INTEGER,  DIMENSION (-3:3) :: SOURCE                    !> Source address in corresponding diagonal
   INTEGER,  DIMENSION (-3:3) :: TARGET                    !> Target address in corresponding diagonal

   INTEGER :: POS(-3:3) = 0                                !> position of IOR's in STENCIL and in matrix storage array
   INTEGER :: N_STENCIL = 0                                !> number of points in matrix stencil
   INTEGER :: N_CONDENSED = 0                              !> number of condensed entries in matrix
   INTEGER :: N_VAL = 0                                    !> number of matrix values in general and symmetric cass
   INTEGER :: N_DIAG = 0                                   !> length of main diagonal

END TYPE SCARC_MATRIX_BANDED_TYPE

!> 
!> ---------- Pressure information
!> 
TYPE SCARC_PRESSURE_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)       :: B_OLD, B_NEW
   REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: H_NEW, HS_NEW
   REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: H_OLD, HS_OLD
   REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: H_DIFF
   REAL(EB) :: DIFF_H = 0.0_EB, DIFF_HS = 0.0_EB

END TYPE SCARC_PRESSURE_TYPE

!> 
!> ---------- Workspace for FFT preconditioners
!> 
TYPE SCARC_FFT_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)       :: SAVE1, WORK, HX
   REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BXS, BXF, BYS, BYF, BZS, BZF
   REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: PRHS
   REAL(EB) :: XS, XF, YS, YF, ZS, ZF
   REAL(EB) :: POIS_PTB = 0.0_EB, XLM = 0.0_EB

   INTEGER :: LSAVE, LWORK
   INTEGER :: LBC, MBC, NBC
   INTEGER :: ITRN, JTRN, KTRN
   INTEGER :: IBAR, JBAR, KBAR

END TYPE SCARC_FFT_TYPE

#ifdef WITH_MKL
!> 
!> ---------- MKL information
!> 
TYPE SCARC_MKL_TYPE

   CHARACTER(40) :: CNAME                                  !> Name of matrix

   INTEGER, ALLOCATABLE :: IPARM(:)
   INTEGER :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR, MSGLVL
   INTEGER :: PERM(1)

   TYPE (MKL_PARDISO_HANDLE),               ALLOCATABLE :: PT_H(:), PT(:)
   TYPE (MKL_CLUSTER_SPARSE_SOLVER_HANDLE), ALLOCATABLE :: CT_H(:), CT(:)

END TYPE SCARC_MKL_TYPE
#endif

!> 
!> ---------- Different scopes for solution, rhs and auxiliary vectors of different solvers
!> 
TYPE SCARC_STAGE_TYPE

   REAL (EB), ALLOCATABLE, DIMENSION (:) :: X          !> solution vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: B          !> right hand side vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: R          !> residual vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: D          !> auxiliary vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: V          !> auxiliary vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: Y          !> auxiliary vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: Z          !> auxiliary vector double precision

   REAL (FB), ALLOCATABLE, DIMENSION (:) :: X_FB       !> solution vector vector single precision
   REAL (FB), ALLOCATABLE, DIMENSION (:) :: B_FB       !> right hand side vector single precision
   REAL (FB), ALLOCATABLE, DIMENSION (:) :: R_FB       !> residual vector single precision
   REAL (FB), ALLOCATABLE, DIMENSION (:) :: V_FB       !> auxiliary vector single precision

#ifdef WITH_SCARC_DEBUG
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: E          !> error vector double precision
#endif

END TYPE SCARC_STAGE_TYPE

!> 
!> ---------- Multigrid type - to be extended for algebraic multigrid
!> 
TYPE SCARC_MULTIGRID_TYPE
   INTEGER :: CYCLING(2) = 0                               !> Counter for multigrid cycling
   INTEGER :: N_PRESMOOTH, N_POSTSMOOTH
END TYPE SCARC_MULTIGRID_TYPE

!> 
!> ---------- Information related to discretization type (structured/unstructured)
!> 
TYPE SCARC_GRID_TYPE

   !> Matrices in different storage techniques
   TYPE (SCARC_MATRIX_BANDED_TYPE)  :: AB                      !> Poisson matrix in banded storage technique
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: AC                      !> Poisson matrix in compact storage technique
#ifdef WITH_MKL
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: ACS                     !> symmetric part of Poisson matrix (only for MKL)
#endif
   TYPE (SCARC_WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: WALL   !> wall information

   !> Pointer arrays for data exchange with neighbors
   REAL(EB), ALLOCATABLE, DIMENSION (:)  :: ICE_TO_VAL         !> mapping from ICE to VAL

   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_IWG         !> mapping from ICE to IWG
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_ICG         !> mapping from ICE to ICG
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_ICN         !> mapping from ICE to ICN

   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_IWG         !> mapping from ICG to IWG
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_ICE         !> mapping from ICG to ICE
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_ICW         !> mapping from ICG to ICW

   !> Assignment of cell coordinates
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICX                !> I-coordinate of cell IC
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICY                !> J-coordinate of cell IC
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICZ                !> J-coordinate of cell IC

   !> Number and state of single cells
   INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: CELL_NUMBER      !> numbering of single cells

   !> Cell numbers of all meshes and offsets between meshes
   INTEGER, ALLOCATABLE, DIMENSION (:) :: NC_LOCAL             !> number of cells in local meshes
   INTEGER, ALLOCATABLE, DIMENSION (:) :: NC_OFFSET            !> offset in cell numbering between meshes
   INTEGER :: NC_GLOBAL = NSCARC_ZERO_INT                      !> global number of cells in all meshes

   !> Local numbers of internal, extended and ghost cells
   INTEGER :: NC  = NSCARC_ZERO_INT                            !> number of cells needed for matrix
   INTEGER :: NW  = NSCARC_ZERO_INT                            !> number of wall cells
   INTEGER :: NCE = NSCARC_ZERO_INT                            !> number of extended cells
   INTEGER :: NCG = NSCARC_ZERO_INT                            !> number of ghost cells 

   !> Number of Dirichlet and Neumann boundary cells
   INTEGER :: N_DIRIC   = NSCARC_ZERO_INT                      !> number of Dirichlet BCs
   INTEGER :: N_NEUMANN = NSCARC_ZERO_INT                      !> number of Neumann BCs

   !> Pointer variables and arrays for data exchange with neighbors
   INTEGER :: ICG = 0, ICG2 = 0, ICE = 0                       !> ghost cell and extended cell pointers

END TYPE SCARC_GRID_TYPE


#ifdef WITH_SCARC_MGM
!> 
!> ---------- McKenney-Greengard-Mayo method
!> 
TYPE SCARC_MGM_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: H1, H2       !> H pressure vectors of different parts
   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: US, VS, WS   !> Velocity components along internal BC's
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: A, IA        !> Lower part of LU-decomposition for unstructured AC
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: L, IL        !> Lower part of LU-decomposition for unstructured AC
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: U, IU        !> Upper part of LU-decomposition for unstructured AC
   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: B, X, Y      !> Right hand side, solution vectors

   INTEGER,  ALLOCATABLE, DIMENSION (:)   :: PERM         !> Permutation vector for reordering of matrix rows
   INTEGER :: NW1, NW2, NWI, NWE                          !> Range of IW's with non-zero B-values

END TYPE SCARC_MGM_TYPE
#endif


!> 
!> ---------- Collection of grid level related information on single mesh
!> 
TYPE SCARC_LEVEL_TYPE

   !> Administrative structures for different components based on given grid level
   TYPE (SCARC_FACE_TYPE), ALLOCATABLE, DIMENSION(:) :: FACE   !> Face information
   TYPE (SCARC_OBST_TYPE), ALLOCATABLE, DIMENSION(:) :: OBST   !> Obstruction information
   TYPE (SCARC_STAGE_TYPE), DIMENSION(2) :: STAGE              !> Hierarchy of solvers and related working vectors
   TYPE (SCARC_GRID_TYPE)      :: STRUCTURED, UNSTRUCTURED     !> Structured and unstructured grid information
   TYPE (SCARC_MULTIGRID_TYPE) :: MG                           !> Multigrid method information
   TYPE (SCARC_FFT_TYPE)       :: FFT                          !> FFT preconditioner based on CRAYFISHPAK
   TYPE (SCARC_PRESSURE_TYPE)  :: PRES                         !> Multigrid method information
#ifdef WITH_SCARC_MGM
   TYPE (SCARC_MGM_TYPE)       :: MGM                          !> McKenney-Greengar-Mayo method 
#endif
#ifdef WITH_MKL
   TYPE (SCARC_MKL_TYPE)       :: MKL                          !> MKL preconditioner based on Intel MKL
#endif

   !> Coordinate information
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: XCOR, YCOR, ZCOR    !> Coordinate vectors in x-, y- and z-direction
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: XMID, YMID, ZMID    !> Midpoint vectors in x-, y- and z-direction
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: DXL, DYL, DZL       !> Step size vectors in x-, y- and z-direction
   REAL(EB) :: DX , DY , DZ                                    !> Step sizes in x-, y- and z-direction
   REAL(EB) :: DXI, DYI, DZI                                   !> Inversed of step sizes in x-, y- and z-direction
   REAL(EB) :: DXI2, DYI2, DZI2                                !> Squared and inversed step sizes in x-, y- and z-direction

   !> Cell and wall index information:
   !>    - on the finest level, the original arrays from FDS are used
   !>    - separate arrays will only be allocated for coarser levels
   !>    - to address them on all levels, corresponding pointers are used
   INTEGER, POINTER, DIMENSION (:,:,:)     :: CELL_INDEX_PTR   !> Pointer to CELL_INDEX
   INTEGER, POINTER, DIMENSION (:,:)       :: WALL_INDEX_PTR   !> Pointer to WALL_INDEX
   INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: CELL_INDEX       !> Cell index list (only allocated for coarser levels)
   INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: WALL_INDEX       !> Wall index list (only allocated for coarser levels)

   LOGICAL, ALLOCATABLE, DIMENSION (:,:,:) :: IS_SOLID         !> State of single cells (.TRUE. if solid/.FALSE. otherwise)

   !> Orientation and different cell related lengths
   INTEGER :: NCG0(-3:3), NCG(-3:3)                            !> First and total wall cell numbers for all faces
   INTEGER :: NX, NY, NZ                                       !> Number of grid cells in x-, y- and z-direction

   !> Number of discretizations and obstructions
   INTEGER :: N_DISCRET                                        !> Number of discretization types used
   INTEGER :: N_OBST                                           !> Number of obstructions
   INTEGER :: N_CELL_INDEX                                     !> Number of entries in CELL_INDEX array
   INTEGER :: N_CELLS                                          !> Number of cells in structured discretization

   !> Different wall related lengths
   INTEGER :: N_WALL_CELLS                                     !> Number of wall cells
   INTEGER :: N_WALL_CELLS_EXT                                 !> Number of external wall cells
   INTEGER :: N_WALL_CELLS_INT                                 !> Number of internal wall cells
   INTEGER :: N_WALL_CELLS_LOCAL                               !> Number of local wall cells

END TYPE SCARC_LEVEL_TYPE

!> 
!> ---------- Sample sequence of used solvers in stack
!> 
TYPE SCARC_SOLVER_TYPE

   CHARACTER(30) :: CNAME = 'NONE'                          !> Name of current solver

   !> Types of different solver components
   INTEGER :: TYPE_ACCURACY      = NSCARC_ACCURACY_ABSOLUTE    !> Default type of requested accuracy
   INTEGER :: TYPE_COARSE        = NSCARC_COARSE_ITERATIVE     !> Default type of coarse grid solver for multigrid method
   INTEGER :: TYPE_COARSENING    = NSCARC_UNDEF_INT            !> No default type of coarsening algorithm for AMG
   INTEGER :: TYPE_CYCLING       = NSCARC_CYCLING_V            !> Default type of cycling for multigrid method
   INTEGER :: TYPE_GRID          = NSCARC_GRID_STRUCTURED      !> Default type of discretization (structured/unstructured)
   INTEGER :: TYPE_EXCHANGE      = NSCARC_UNDEF_INT            !> No default type of data exchange
   INTEGER :: TYPE_INTERPOL      = NSCARC_INTERPOL_CONSTANT    !> Default type of interpolation method
   INTEGER :: TYPE_KRYLOV        = NSCARC_KRYLOV_CG            !> Default type of Krylov method (CG/BICG)
   INTEGER :: TYPE_LEVEL(0:2)    = NSCARC_UNDEF_INT            !> Default type of levels
   INTEGER :: TYPE_MATRIX        = NSCARC_MATRIX_COMPACT       !> Default type of storage for matrix
   INTEGER :: TYPE_METHOD        = NSCARC_METHOD_KRYLOV        !> Default type of ScaRC method
   INTEGER :: TYPE_MKL(0:10)     = NSCARC_UNDEF_INT            !> No default type of MKL for single levels
   INTEGER :: TYPE_MKL_PRECISION = NSCARC_PRECISION_DOUBLE     !> Default double precision MKL solver
   INTEGER :: TYPE_MULTIGRID     = NSCARC_MULTIGRID_GEOMETRIC  !> Default type of multigrid method (GMG/AMG)
   INTEGER :: TYPE_PARENT        = NSCARC_UNDEF_INT            !> No default type of parent (calling) solver
   INTEGER :: TYPE_PRECON        = NSCARC_UNDEF_INT            !> No default type of preconditioner for iterative solver
   INTEGER :: TYPE_RELAX         = NSCARC_UNDEF_INT            !> No default type of preconditioner for iterative solver
   INTEGER :: TYPE_SCOPE(0:2)    = NSCARC_SCOPE_LOCAL          !> Default scope types are local
   INTEGER :: TYPE_SMOOTH        = NSCARC_UNDEF_INT            !> No default type of smoother for multigrid method
   INTEGER :: TYPE_SOLVER        = NSCARC_SOLVER_MAIN          !> Default type of surrounding solver stage
   INTEGER :: TYPE_STAGE         = NSCARC_STAGE_ONE            !> Default type of surrounding solver stage
   INTEGER :: TYPE_STENCIL       = NSCARC_STENCIL_CONSTANT     !> Default type of storage for matrix
   INTEGER :: TYPE_TWOLEVEL      = NSCARC_TWOLEVEL_NONE        !> Default type of two-level method
   INTEGER :: TYPE_VECTOR        = NSCARC_UNDEF_INT            !> No default type of vector to point to


   !> References to different vectors which are needed for the current solver
   INTEGER :: X = NSCARC_UNDEF_INT                 !> Reference to local X-vector, double precision
   INTEGER :: B = NSCARC_UNDEF_INT                 !> Reference to local B-vector, double precision
   INTEGER :: D = NSCARC_UNDEF_INT                 !> Reference to local D-vector, double precision
   INTEGER :: R = NSCARC_UNDEF_INT                 !> Reference to local R-vector, double precision
   INTEGER :: V = NSCARC_UNDEF_INT                 !> Reference to local V-vector, double precision
   INTEGER :: Y = NSCARC_UNDEF_INT                 !> Reference to local Y-vector, double precision
   INTEGER :: Z = NSCARC_UNDEF_INT                 !> Reference to local Z-vector, double precision

   INTEGER :: X_FB = NSCARC_UNDEF_INT              !> Reference to local X-vector, single precision
   INTEGER :: B_FB = NSCARC_UNDEF_INT              !> Reference to local B-vector, single precision
   INTEGER :: R_FB = NSCARC_UNDEF_INT              !> Reference to local R-vector, single precision
   INTEGER :: V_FB = NSCARC_UNDEF_INT              !> Reference to local V-vector, single precision

#ifdef WITH_SCARC_DEBUG
   INTEGER :: E = NSCARC_UNDEF_INT                 !> Reference to local E-vector, double precision
#endif

   !> Converegence requirements for current solver
   INTEGER  :: NIT   = NSCARC_UNDEF_INT            !> Maximum iteration number
   INTEGER  :: ITE   = NSCARC_UNDEF_INT            !> Current iteration number
   REAL(EB) :: EPS   = NSCARC_UNDEF_REAL_EB        !> Required accuracy
   REAL(EB) :: RES   = NSCARC_UNDEF_REAL_EB        !> Current residual
   REAL(EB) :: RESIN = NSCARC_UNDEF_REAL_EB        !> Initial residual
   REAL(EB) :: ERR   = NSCARC_UNDEF_REAL_EB        !> Initial residual
   REAL(EB) :: OMEGA = NSCARC_UNDEF_REAL_EB        !> Relaxation parameter
   REAL(EB) :: CAPPA = NSCARC_UNDEF_REAL_EB        !> Convergence rate

END TYPE SCARC_SOLVER_TYPE

!> 
!> ---------- Stack type
!> 
TYPE SCARC_STACK_TYPE
   TYPE (SCARC_SOLVER_TYPE), POINTER :: SOLVER     !> Type of current solver
END TYPE SCARC_STACK_TYPE

!> 
!> ---------- Administration other mesh data needed for the coupling of adjacent neighbors
!> 
TYPE SCARC_NEIGHBOR_TYPE

   TYPE (SCARC_LEVEL_TYPE), ALLOCATABLE, DIMENSION(:) :: LEVEL          !> Description of level related information

   REAL(EB) :: SEND_REAL0(50), RECV_REAL0(50)                           !> Exchange buffers for setup 
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: SEND_REAL1, RECV_REAL1       !> Exchange buffers (usual length)
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: SEND_REAL2, RECV_REAL2       !> Exchange buffers (double length)
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: SEND_REAL7, RECV_REAL7       !> Exchange buffers (sevenfold length)

   INTEGER :: SEND_INT0(50), RECV_INT0(50)                              !> Integer exchange buffers for setup
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: SEND_INT1, RECV_INT1         !> Integer exchange buffers

END TYPE SCARC_NEIGHBOR_TYPE

!> 
!> ---------- Measurement of CPU times
!> 
TYPE SCARC_CPU_TYPE
   REAL(EB) :: OVERALL          = 0.0_EB                      !> complete time for ScaRC
   REAL(EB) :: SETUP            = 0.0_EB                      !> time for setup of requested ScaRC solver
   REAL(EB) :: SOLVER           = 0.0_EB                      !> time for solver 
   REAL(EB) :: ITERATION        = 0.0_EB                      !> time for krylov solver
   REAL(EB) :: MATVEC_PRODUCT   = 0.0_EB                      !> time for matrix vector multiplication
   REAL(EB) :: SCALAR_PRODUCT   = 0.0_EB                      !> time for scalar product
   REAL(EB) :: L2NORM           = 0.0_EB                      !> time for l2-norm
   REAL(EB) :: RELAXATION       = 0.0_EB                      !> time for relaxation
   REAL(EB) :: SMOOTHER         = 0.0_EB                      !> time for smoothing
   REAL(EB) :: COARSE           = 0.0_EB                      !> time for coarse grid solver
   REAL(EB) :: EXCHANGE         = 0.0_EB                      !> time for data exchange
   REAL(EB) :: BUFFER_PACKING   = 0.0_EB                      !> time for data exchange
   REAL(EB) :: BUFFER_UNPACKING = 0.0_EB                      !> time for data exchange
   INTEGER  :: N_TIMER          = 12                          !> total number of timers
END TYPE SCARC_CPU_TYPE

!>
!> ---------- Basic administration type for ScaRC-method
!>
TYPE SCARC_TYPE

   TYPE (SCARC_NEIGHBOR_TYPE) , ALLOCATABLE, DIMENSION(:)   :: OSCARC    !> ScaRC type on other mesh
   TYPE (SCARC_LEVEL_TYPE)    , ALLOCATABLE, DIMENSION(:)   :: LEVEL     !> Level related information

   REAL(EB) :: XS, XF, YS, YF, ZS, ZF                                    !> x-, y- and z-bounds (corresponding to main prg)
   REAL(EB) :: RHS_END = 0.0_EB                                          !> Very last RHS entry, needed for matrix condensing

   INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS                       !> List of adjacent neighbors of whole mesh
   INTEGER :: N_NEIGHBORS = 0                                            !> Number of adjacent neighbors of whole mesh
   INTEGER :: NC = 0                                                     !> Total number of cells on that mesh
   INTEGER :: IBAR, JBAR, KBAR                                           !> Number of cells (corresponding to main prg)

END TYPE SCARC_TYPE

END MODULE SCARC_TYPES


!> ================================================================================================================
!> MODULE 'SCARC_POINTERS' : 
!> Collection of different pointers to specify the different meshes, grid levels, discretizations and matrices
!> ================================================================================================================
MODULE SCARC_POINTERS
USE MESH_VARIABLES
USE SCARC_TYPES
IMPLICIT NONE

TYPE (MESH_TYPE),  POINTER :: M=>NULL()
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()
TYPE (WALL_TYPE),  POINTER :: MWC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()

TYPE (SCARC_TYPE),          POINTER :: S=>NULL()
TYPE (SCARC_NEIGHBOR_TYPE), POINTER :: OS=>NULL()
TYPE (SCARC_LEVEL_TYPE),    POINTER :: L=>NULL(), LF=>NULL(), LC=>NULL(), OL=>NULL(), OLF=>NULL(), OLC=>NULL()
TYPE (SCARC_GRID_TYPE),     POINTER :: G=>NULL(), GF=>NULL(), GC=>NULL(), OG=>NULL(), OGF=>NULL(), OGC=>NULL()
TYPE (SCARC_FACE_TYPE),     POINTER :: F=>NULL(), FF=>NULL(), FC=>NULL()
TYPE (SCARC_PRESSURE_TYPE), POINTER :: P=>NULL()
TYPE (SCARC_OBST_TYPE),     POINTER :: OB=>NULL()

TYPE (SCARC_WALL_TYPE), POINTER :: GWC=>NULL()
TYPE (SCARC_WALL_TYPE), DIMENSION(:), POINTER :: W=>NULL(), WF=>NULL(), WC=>NULL()

TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL(), SVP=>NULL()
TYPE (SCARC_STAGE_TYPE),  POINTER :: ST=>NULL(), STP=>NULL()

TYPE (SCARC_FFT_TYPE), POINTER :: FFT=>NULL()
TYPE (SCARC_MULTIGRID_TYPE), POINTER :: MG=>NULL()
#ifdef WITH_SCARC_MGM
TYPE (SCARC_MGM_TYPE), POINTER :: MGM=>NULL()
#endif

TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()
TYPE (SCARC_MATRIX_BANDED_TYPE),  POINTER :: AB=>NULL()

TYPE (SCARC_MATRIX_COMPACT_CONDENSED_TYPE), POINTER :: ACC=>NULL()
TYPE (SCARC_MATRIX_BANDED_CONDENSED_TYPE),  POINTER :: ABC=>NULL()

REAL(EB), DIMENSION(:), POINTER :: XCOR=>NULL(), YCOR=>NULL(), ZCOR=>NULL()
REAL(EB), DIMENSION(:), POINTER :: XMID=>NULL(), YMID=>NULL(), ZMID=>NULL()

REAL(EB), DIMENSION(:), POINTER :: VC=>NULL(), VF=>NULL()
REAL(EB), DIMENSION(:), POINTER :: V1=>NULL(), V2=>NULL()
REAL(FB), DIMENSION(:), POINTER :: V1_FB=>NULL(), V2_FB=>NULL()

REAL(EB), DIMENSION(:,:,:), POINTER :: HP=>NULL(), HVC=>NULL(), PRHS=>NULL()
REAL(EB), DIMENSION(:,:,:), POINTER :: UU=>NULL(), VV=>NULL(), WW=>NULL()

REAL(EB), DIMENSION(:), POINTER ::  RECV_REAL0, RECV_REAL1, RECV_REAL2, RECV_REAL7
INTEGER,  DIMENSION(:), POINTER ::  RECV_INT0, RECV_INT1

INTEGER, POINTER :: I_MIN_R, I_MAX_R, I_MIN_S, I_MAX_S
INTEGER, POINTER :: J_MIN_R, J_MAX_R, J_MIN_S, J_MAX_S
INTEGER, POINTER :: K_MIN_R, K_MAX_R, K_MIN_S, K_MAX_S

INTEGER, POINTER, DIMENSION(:) :: IOR_R, IOR_S
INTEGER, POINTER, DIMENSION(:) :: IIO_R, IIO_S
INTEGER, POINTER, DIMENSION(:) :: JJO_R, JJO_S
INTEGER, POINTER, DIMENSION(:) :: KKO_R, KKO_S

#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: MKL=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: ACS=>NULL()
#endif

END MODULE SCARC_POINTERS


!> ================================================================================================================
!> ROUTINE 'SCRC'
!> ================================================================================================================
MODULE SCRC

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME, GET_FILE_NUMBER, SHUTDOWN
USE MPI
USE SCARC_GLOBAL_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_TYPES
USE SCARC_POINTERS

#ifdef WITH_MKL
USE MKL_PARDISO
USE MKL_CLUSTER_SPARSE_SOLVER
#endif

IMPLICIT NONE

!> 
!> ---------- Public variables
!> 
PUBLIC :: SCARC_SETUP, SCARC_METHOD, SCARC_SOLVER, SCARC_GRID, SCARC_TWOLEVEL      
PUBLIC :: SCARC_ACCURACY, SCARC_CAPPA, SCARC_ITERATIONS, SCARC_RESIDUAL, SCARC_MATRIX
PUBLIC :: SCARC_DUMP_TIMERS, SCARC_ERROR_FILE, SCARC_MKL_PRECISION

PUBLIC :: SCARC_KRYLOV, SCARC_KRYLOV_ACCURACY, SCARC_KRYLOV_ITERATIONS, SCARC_KRYLOV_INTERPOL
PUBLIC :: SCARC_COARSE, SCARC_COARSE_ACCURACY, SCARC_COARSE_ITERATIONS, SCARC_COARSE_OMEGA, SCARC_COARSE_LEVEL     
PUBLIC :: SCARC_PRECON, SCARC_PRECON_ACCURACY, SCARC_PRECON_ITERATIONS, SCARC_PRECON_OMEGA, SCARC_PRECON_SCOPE 
PUBLIC :: SCARC_SMOOTH, SCARC_SMOOTH_ACCURACY, SCARC_SMOOTH_ITERATIONS, SCARC_SMOOTH_OMEGA, SCARC_SMOOTH_SCOPE 

PUBLIC :: SCARC_MULTIGRID, SCARC_MULTIGRID_ACCURACY, SCARC_MULTIGRID_ITERATIONS, SCARC_MULTIGRID_INTERPOL
PUBLIC :: SCARC_MULTIGRID_CYCLE, SCARC_MULTIGRID_COARSENING , SCARC_MULTIGRID_LEVEL
PUBLIC :: SCARC_MULTIGRID_PRESMOOTH, SCARC_MULTIGRID_POSTSMOOTH


!> 
!> ---------- Type declarations
!> 
TYPE (SCARC_TYPE)        , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: SCARC
TYPE (SCARC_STACK_TYPE)  , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: STACK
TYPE (SCARC_CPU_TYPE)    , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: CPU

TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_CG, MAIN_CG_STRUCTURED, MAIN_CG_UNSTRUCTURED
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_GMG
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: COARSE_KRYLOV
#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_LU
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: COARSE_CLUSTER, COARSE_PARDISO
#endif
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_JAC, PRECON_SSOR, PRECON_ILU, PRECON_FFT, PRECON_FFTO, PRECON_GMG
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_MJAC, PRECON_MGS, PRECON_MSGS, PRECON_MSOR, PRECON_MSSOR
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: SMOOTH_JACOBI, SMOOTH_SSOR, SMOOTH_ILU, SMOOTH_FFT, SMOOTH_FFTO
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: SMOOTH_MJAC, SMOOTH_MGS, SMOOTH_MSGS, SMOOTH_MSOR, SMOOTH_MSSOR
#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_MKL, SMOOTH_MKL
#endif
TYPE (SCARC_MESSAGE_TYPE), SAVE :: MSG


CONTAINS


!> ------------------------------------------------------------------------------------------------
!> Initialize ScaRC structures based on SCARC-input parameters from &PRES namelist
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

!>
!> Setup messaging/debugging mechanisms and time measurements
!>
CALL SCARC_SETUP_MESSAGING
CALL SCARC_SETUP_TIMING

!>
!> Parse all ScaRC parameters which have been read in read.f90
!>
CALL SCARC_PARSE_INPUT                      ; IF (STOP_STATUS==SETUP_STOP) RETURN

!>
!> Setup different basic components of ScaRC solver
!>
CALL SCARC_SETUP_LEVELS                     ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_TYPES                      ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_GRIDS                      ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_NEIGHBORS                  ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_GLOBALS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_FACES                      ; IF (STOP_STATUS==SETUP_STOP) RETURN

!>
!> Setup wall information according to specified discretization type/method
!>
IF (HAS_MULTI_GRIDS) THEN
   CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_STRUCTURED)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
   CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
ELSE
   CALL SCARC_ASSIGN_GRID_TYPE (TYPE_GRID)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
ENDIF

!>
!> Setup information for data exchanges, matrix systems, used solvers and vectors
!>
CALL SCARC_SETUP_EXCHANGES                  ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_SYSTEMS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_METHODS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_VECTORS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN
#ifdef WITH_SCARC_STANDALONE
CALL SCARC_SETUP_PRESSURE                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
#endif


CPU(MYID)%SETUP   = CPU(MYID)%SETUP   + CURRENT_TIME() - TNOW
CPU(MYID)%OVERALL = CPU(MYID)%OVERALL + CURRENT_TIME() - TNOW
END SUBROUTINE SCARC_SETUP


!> ------------------------------------------------------------------------------------------------
!> Setup time measurements
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TIMING
ALLOCATE (CPU(0:N_MPI_PROCESSES-1), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_TIMING', 'CPU', IERROR)
END SUBROUTINE SCARC_SETUP_TIMING


!> ------------------------------------------------------------------------------------------------
!> Setup debug file if requested
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MESSAGING
#if defined(WITH_SCARC_VERBOSE) || defined(WITH_SCARC_DEBUG)
INTEGER :: NM, LASTID
#endif

IF (SCARC_ERROR_FILE) HAS_CSV_DUMP = .TRUE.

!>
!> If requested, open file for CSV-information about convergence of different solvers
!>
IF (HAS_CSV_DUMP) THEN
   IF (MYID == 0) THEN
      WRITE (MSG%FILE_STAT, '(A,A)') TRIM(CHID),'_scarc.csv'
      MSG%LU_STAT = GET_FILE_NUMBER()
      OPEN (MSG%LU_STAT, FILE=MSG%FILE_STAT)
      WRITE(MSG%LU_STAT,*) '  #Pres,   Stack,  #ScaRC,     #CG,     #MG,   Level, #Smooth, SmoType, ', &
                           '#Coarse,     #LU,    Residual,   Cappa'
   ENDIF
ENDIF

#ifdef WITH_SCARC_VERBOSE
!>
!> If requested, open file for log-information
!>
LASTID = -99999
DO NM=LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (MYID == LASTID) CYCLE
   WRITE (MSG%FILE_VERBOSE, '(A,A,i3.3)') TRIM(CHID),'.log',MYID+1
   MSG%LU_VERBOSE = GET_FILE_NUMBER()
   OPEN (MSG%LU_VERBOSE, FILE=MSG%FILE_VERBOSE, ACTION = 'readwrite')
   LASTID = MYID
ENDDO
#endif

#ifdef WITH_SCARC_DEBUG
!>
!> If requested, open file for debug messages
!>
LASTID = -99999
DO NM=LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (MYID == LASTID) CYCLE
   WRITE (MSG%FILE_DEBUG, '(A,A,i3.3)') TRIM(CHID),'.debug',MYID+1
   MSG%LU_DEBUG = GET_FILE_NUMBER()
   OPEN (MSG%LU_DEBUG, FILE=MSG%FILE_DEBUG, ACTION = 'readwrite')
   LASTID = MYID
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_MESSAGING


!> ------------------------------------------------------------------------------------------------
!> Shutdown ScaRC with error message
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SHUTDOWN(NERROR, CPARAM, NPARAM)
CHARACTER(*), INTENT(IN) :: CPARAM
INTEGER, INTENT(IN) :: NERROR, NPARAM
CHARACTER(80) :: CERROR

!>
!> Assign error message according to specified error
!>
SELECT CASE (NERROR)
   CASE (NSCARC_ERROR_PARSE_INPUT)
      CERROR = 'Wrong input parameter'
   CASE (NSCARC_ERROR_MKL_INTERNAL)
      CERROR = 'The following MKL error was detected'
   CASE (NSCARC_ERROR_MKL_PARDISO)
      CERROR = 'MKL Library compile flag not defined, Pardiso solver not available'
   CASE (NSCARC_ERROR_MKL_CLUSTER)
      CERROR = 'MKL Library compile flag not defined, Cluster_Sparse_Solver not available'
   CASE (NSCARC_ERROR_MKL_STORAGE)
      CERROR = 'Wrong matrix storage scheme for MKL solvers, only COMPACT storage available'
   CASE (NSCARC_ERROR_GRID_NUMBER)
      CERROR = 'Number not divisable by 2'
   CASE (NSCARC_ERROR_GRID_NUMBERX)
      CERROR = 'Number of cells not divisable by 2 in x-direction, NC'
   CASE (NSCARC_ERROR_GRID_NUMBERY)
      CERROR = 'Number of cells not divisable by 2 in y-direction, NC'
   CASE (NSCARC_ERROR_GRID_NUMBERZ)
      CERROR = 'Number of cells not divisable by 2 in z-direction, NC'
   CASE (NSCARC_ERROR_BOUNDARY_SUM)
      CERROR = 'Wrong boundary sum for IOR'
   CASE (NSCARC_ERROR_DIRECT_NOMKL)
      CERROR = 'Direct coarse grid solver is only working in combination with MKL'
   CASE (NSCARC_ERROR_BOUNDARY_TYPE)
      CERROR = 'Wrong boundary type'
   CASE (NSCARC_ERROR_NEIGHBOR_TYPE)
      CERROR = 'Wrong neighbor'
   CASE (NSCARC_ERROR_NEIGHBOR_NUMBER)
      CERROR = 'More than 20 neighbors along one face not allowed'
   CASE (NSCARC_ERROR_MATRIX_ALLOCATION)
      CERROR = 'Wrong specifier during allocation or deallocation of  matrix'
   CASE (NSCARC_ERROR_MATRIX_SYMMETRY)
      CERROR = 'Matrix not symmetric for mesh'
   CASE (NSCARC_ERROR_MATRIX_SUBDIAG)
      CERROR = 'Subdiagonal missing for system matrix'
   CASE (NSCARC_ERROR_MATRIX_SETUP)
      CERROR = 'Matrix setup failed for level type'
   CASE (NSCARC_ERROR_MATRIX_SIZE)
      CERROR = 'Matrix resized failed because new length is too big'
   CASE (NSCARC_ERROR_MATRIX_COPY)
      CERROR = 'Matrix copy failed due to too already existing array'
   CASE (NSCARC_ERROR_STENCIL)
      CERROR = 'Wrong type for matrix stencil - only constant or variable allowed'
   CASE (NSCARC_ERROR_STACK_SOLVER)
      CERROR = 'Wrong number of solvers in stack'
   CASE (NSCARC_ERROR_STACK_MESSAGE)
      CERROR = 'Too many messages in calling stack'
   CASE (NSCARC_ERROR_MULTIGRID_LEVEL)
      CERROR = 'Wrong level for multigrid method'
   CASE (NSCARC_ERROR_GRID_RESOLUTION)
      CERROR = 'Wrong grid resolution at IOR'
   CASE (NSCARC_ERROR_GRID_INDEX)
      CERROR = 'Wrong index for J'
   CASE (NSCARC_ERROR_RHS_SETUP)
      CERROR = 'Wrong level for presetting RHS'
   CASE (NSCARC_ERROR_AMG_DISABLED)
      CERROR = 'Algebraic multigrid is currently missing'
   CASE (NSCARC_ERROR_VECTOR_LENGTH)
      CERROR = 'Inconsistent length for vector allocation'
   CASE (NSCARC_ERROR_POC_STOP)
      CERROR = 'Only one call of solver due to proof of concept'
   CASE (NSCARC_ERROR_BICG_DISABLED)
      CERROR = 'Krylov solver BICG temporarily disabled'
   CASE (NSCARC_ERROR_EXCHANGE_RECV)
      CERROR = 'Wrong receive exchange structure'
   CASE (NSCARC_ERROR_EXCHANGE_SEND)
      CERROR = 'Wrong send exchange structure'
END SELECT

!>
!> Specify more detailed information if available
!>
IF (CPARAM /= SCARC_NONE) THEN
   IF (MYID == 0) WRITE(LU_ERR,1000)  CERROR, CPARAM, TRIM(CHID)
ELSE IF (NPARAM /= NSCARC_NONE) THEN
   IF (MYID == 0) WRITE(LU_ERR,2000)  CERROR, NPARAM, TRIM(CHID)
ELSE
   IF (MYID == 0) WRITE(LU_ERR,3000)  CERROR, TRIM(CHID)
ENDIF

!>
!> Also print verbose message if enabled
!>
#ifdef WITH_SCARC_VERBOSE
IF (CPARAM /= SCARC_NONE) THEN
   WRITE(MSG%LU_VERBOSE,1000)  CERROR, CPARAM, TRIM(CHID)
ELSE IF (NPARAM /= NSCARC_NONE) THEN
   WRITE(MSG%LU_VERBOSE,2000)  CERROR, NPARAM, TRIM(CHID)
ELSE
   WRITE(MSG%LU_VERBOSE,3000)  CERROR, TRIM(CHID)
ENDIF
CLOSE(MSG%LU_VERBOSE)
#endif

#ifdef WITH_SCARC_DEBUG
CLOSE(MSG%LU_DEBUG)
#endif

STOP_STATUS = SETUP_STOP
RETURN

1000 FORMAT('Stop in ScaRC-solver: ', A,' : ',   A, ' (CHID: ',A,')' )
2000 FORMAT('Stop in ScaRC-solver: ', A,' : ', I12, ' (CHID: ',A,')' )
3000 FORMAT('Stop in ScaRC-solver: ', A, ' (CHID: ',A,')' )
END SUBROUTINE SCARC_SHUTDOWN


!> ----------------------------------------------------------------------------------------------------
!> Determine types of input parameters
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PARSE_INPUT

ITERATE_PRESSURE = .TRUE.  !> Although there is no need to do pressure iterations to drive down 
                           !> velocity error leave it .TRUE. to write out velocity error diagnostics

!>
!> ------------- Set type of discretization
!>
SELECT CASE (TRIM(PRES_METHOD))
   CASE ('SCARC')
      TYPE_GRID = NSCARC_GRID_STRUCTURED
      IS_STRUCTURED   = .TRUE.
      IS_UNSTRUCTURED = .FALSE.
   CASE ('USCARC')
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED
      IS_STRUCTURED   = .FALSE.
      IS_UNSTRUCTURED = .TRUE.
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_GRID, NSCARC_NONE)
END SELECT

!>
!> ------------ Set type of matrix storage (COMPACT/BANDED)
!>
SELECT CASE (TRIM(SCARC_MATRIX))
   CASE ('COMPACT')
      TYPE_MATRIX = NSCARC_MATRIX_COMPACT
   CASE ('BANDED')
      TYPE_MATRIX = NSCARC_MATRIX_BANDED
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MATRIX, NSCARC_NONE)
END SELECT

!>
!> ------------ Set type of matrix stencil (CONSTANT/VARIABLE)
!>
SELECT CASE (TRIM(SCARC_STENCIL))
   CASE ('CONSTANT')
      TYPE_STENCIL = NSCARC_STENCIL_CONSTANT
   CASE ('VARIABLE')
      TYPE_STENCIL = NSCARC_STENCIL_VARIABLE
   CASE ('CONSTANT_LOCAL')
      TYPE_STENCIL = NSCARC_STENCIL_CONSTANT_LOCAL
   CASE ('VARIABLE_LOCAL')
      TYPE_STENCIL = NSCARC_STENCIL_VARIABLE_LOCAL
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_STENCIL, NSCARC_NONE)
END SELECT

!>
!> ------------ Set type of global solver
!>
SELECT CASE (TRIM(SCARC_METHOD))

   !> ------------------------- Global Krylov solver ----------------------------------
   CASE ('KRYLOV')

      TYPE_METHOD = NSCARC_METHOD_KRYLOV

      !> Set type of Krylov-method (CG/BICG)
      SELECT CASE (TRIM(SCARC_KRYLOV))
         CASE ('CG')
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
#ifdef WITH_SCARC_BAROCLINIC
         CASE ('BAROCLINIC')
            TYPE_KRYLOV = NSCARC_KRYLOV_CGBARO
#endif
         CASE ('BICG')
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_BICG_DISABLED, SCARC_KRYLOV, NSCARC_NONE)    !> only temporarily
            TYPE_KRYLOV = NSCARC_KRYLOV_BICG
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_KRYLOV, NSCARC_NONE)
      END SELECT

      !> Set type of two-level method
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
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_TWOLEVEL, NSCARC_NONE)
      END SELECT

      !> Set type of interpolation for two-level Krylov method
      SELECT CASE (TRIM(SCARC_KRYLOV_INTERPOL))
         CASE ('NONE')
            TYPE_INTERPOL = NSCARC_UNDEF_INT
         CASE ('CONSTANT')
            TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
         CASE ('BILINEAR')
            TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_KRYLOV_INTERPOL, NSCARC_NONE)
      END SELECT

      !> Set type of preconditioner (JACOBI/SSOR/MGS/MSGS/MSOR/MSSOR/ILU/FFT/GMG/PARDISO/CLUSTER)
      SELECT CASE (TRIM(SCARC_PRECON))
         CASE ('JACOBI')                                    !> Jacobi preconditioner
            TYPE_PRECON = NSCARC_RELAX_JAC
         CASE ('SSOR')                                      !> Symmetric SOR preconditioner
            TYPE_PRECON = NSCARC_RELAX_SSOR
         CASE ('MJAC')                                      !> Jacobi preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MJAC
         CASE ('MGS')                                       !> Gauss-Seidel preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MGS
         CASE ('MSGS')                                      !> Symmetric Gauss-Seidel preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MSGS
         CASE ('MSOR')                                      !> SOR preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MSOR
         CASE ('MSSOR')                                     !> Symmetric SOR preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MSSOR
         CASE ('ILU')                                       !> ILU preconditioner
            TYPE_PRECON = NSCARC_RELAX_ILU
         CASE ('GMG')                                       !> Geometric multigrid preconditioner
            TYPE_PRECON = NSCARC_RELAX_GMG
            SELECT CASE (TRIM(SCARC_SMOOTH))
               CASE ('JACOBI')
                  TYPE_SMOOTH = NSCARC_RELAX_JAC
               CASE ('SSOR')
                  TYPE_SMOOTH = NSCARC_RELAX_SSOR
               CASE ('MJAC')
                  TYPE_SMOOTH = NSCARC_RELAX_MJAC
               CASE ('MGS')
                  TYPE_SMOOTH = NSCARC_RELAX_MGS
               CASE ('MSGS')
                  TYPE_SMOOTH = NSCARC_RELAX_MSGS
               CASE ('MSOR')
                  TYPE_SMOOTH = NSCARC_RELAX_MSOR
               CASE ('MSSOR')
                  TYPE_SMOOTH = NSCARC_RELAX_MSSOR
               CASE ('ILU')
                  TYPE_SMOOTH = NSCARC_RELAX_ILU
               CASE ('FFT')
                  IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
                  TYPE_SMOOTH = NSCARC_RELAX_FFT
               CASE ('FFTO')
                  IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
                  IF (NMESHES == 1) THEN
                     TYPE_SMOOTH = NSCARC_RELAX_FFT
                  ELSE
                     TYPE_SMOOTH = NSCARC_RELAX_FFTO
                  ENDIF
               CASE ('PARDISO')
#ifdef WITH_MKL
                  TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
                  CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif

               CASE ('CLUSTER')
#ifdef WITH_MKL
                  TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
                  CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
#endif
            END SELECT
         CASE ('FFT')                                                !> FFT preconditioner
            IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            TYPE_PRECON = NSCARC_RELAX_FFT
         CASE ('FFTO')                                               !> FFT with overlap preconditioner
            IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            IF (NMESHES == 1) THEN
               TYPE_PRECON = NSCARC_RELAX_FFT
            ELSE
               TYPE_PRECON = NSCARC_RELAX_FFTO
            ENDIF
         CASE ('PARDISO')                                            !> LU preconditioner based on MKL-PARDISO
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_MKL
            TYPE_MKL(0)   = NSCARC_MKL_LOCAL
            TYPE_SCOPE(0) = NSCARC_SCOPE_LOCAL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('CLUSTER')                            !>  LU-preconditioner based on MKL Cluster_Sparse_Solver
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_MKL
            TYPE_MKL(0)   = NSCARC_MKL_GLOBAL
            TYPE_SCOPE(0) = NSCARC_SCOPE_GLOBAL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_PRECON, NSCARC_NONE)
      END SELECT

      !> set type scope for preconditioner (GLOBAL/LOCAL)
      SELECT CASE (TRIM(SCARC_PRECON_SCOPE))
         CASE ('GLOBAL')
            TYPE_SCOPE(1) = NSCARC_SCOPE_GLOBAL
         CASE ('LOCAL')
            TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_PRECON_SCOPE, NSCARC_NONE)
      END SELECT

   !> ------------------------- Global geometric multigrid solver -------------------------------
   CASE ('MULTIGRID')

      TYPE_METHOD = NSCARC_METHOD_MULTIGRID

      !> Set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_MULTIGRID))
         CASE ('GEOMETRIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
         CASE ('ALGEBRAIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID, NSCARC_NONE)
      END SELECT

      !> Set type of smoother (JACOBI/SGS/SSOR/MSSOR/ILU/PARDISO/CLUSTER)
      SELECT CASE (TRIM(SCARC_SMOOTH))                        !> use same parameters as for preconditioner
         CASE ('JACOBI')
            TYPE_SMOOTH = NSCARC_RELAX_JAC
         CASE ('SSOR')
            TYPE_SMOOTH = NSCARC_RELAX_SSOR
         CASE ('MJAC')                                      !> Jacobi preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MJAC
         CASE ('MGS')                                       !> Gauss-Seidel preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MGS
         CASE ('MSGS')                                      !> Symmetric Gauss-Seidel preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MSGS
         CASE ('MSOR')                                      !> SOR preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MSOR
         CASE ('MSSOR')                                     !> Symmetric SOR preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MSSOR
         CASE ('ILU')
            TYPE_SMOOTH = NSCARC_RELAX_ILU
         CASE ('FFT')
            IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            TYPE_SMOOTH = NSCARC_RELAX_FFT
         CASE ('FFTO')
            IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            IF (NMESHES == 1) THEN
               TYPE_SMOOTH = NSCARC_RELAX_FFT
            ELSE
               TYPE_SMOOTH = NSCARC_RELAX_FFTO
            ENDIF
         CASE ('PARDISO')
#ifdef WITH_MKL
            TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('CLUSTER')
#ifdef WITH_MKL
            TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_SMOOTH, NSCARC_NONE)
      END SELECT

      !> set type scope for smoother (GLOBAL/LOCAL)
      SELECT CASE (TRIM(SCARC_SMOOTH_SCOPE))
         CASE ('GLOBAL')
            TYPE_SCOPE(2) = NSCARC_SCOPE_GLOBAL
         CASE ('LOCAL')
            TYPE_SCOPE(2) = NSCARC_SCOPE_LOCAL
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_SMOOTH_SCOPE, NSCARC_NONE)
      END SELECT

   !> ------------------------- Global LU-decomposition solver -------------------------------
#ifdef WITH_MKL
   CASE ('MKL')

      TYPE_METHOD  = NSCARC_METHOD_LU

      !> Set type of MKL method (global/local)
      SELECT CASE (TRIM(SCARC_MKL))                  
         CASE ('GLOBAL')
#ifdef WITH_MKL
            TYPE_MKL(0)   = NSCARC_MKL_GLOBAL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('LOCAL')
#ifdef WITH_MKL
            TYPE_MKL(0)   = NSCARC_MKL_LOCAL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MKL, NSCARC_NONE)
      END SELECT
#endif

#ifdef WITH_SCARC_MGM
   !> ------------------------- McKenny-Greengard-Mayo solver -------------------------
   CASE ('MGM')

      !> Just preset some values for proof of concept
      HAS_MULTI_GRIDS = .TRUE.

      TYPE_METHOD   = NSCARC_METHOD_MGM
      TYPE_KRYLOV   = NSCARC_KRYLOV_CG
      TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE
      TYPE_PRECON   = NSCARC_RELAX_FFT
      TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
#endif

   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_METHOD, NSCARC_NONE)

END SELECT

!>
!> If a multigrid solver is used (either as main solver or as preconditioner)
!> set types for multigrid, coarse grid solver and cycling pattern
!>
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_RELAX_GMG) THEN

   !> Set type of multigrid (GEOMETRIC/ALGEBRAIC with corresponding coarsening strategy)
   SELECT CASE (TRIM(SCARC_MULTIGRID))

      CASE ('GEOMETRIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC

      CASE ('ALGEBRAIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC

         !> Set type of coarsening strategy (STANDARD/AGGRESSIVE)
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
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_COARSENING, NSCARC_NONE)
         END SELECT

      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID, NSCARC_NONE)
   END SELECT

   !> Set type of cycling pattern (F/V/W)
   SELECT CASE (TRIM(SCARC_MULTIGRID_CYCLE))
      CASE ('F')
         TYPE_CYCLING = NSCARC_CYCLING_F
      CASE ('V')
         TYPE_CYCLING = NSCARC_CYCLING_V
      CASE ('W')
         TYPE_CYCLING = NSCARC_CYCLING_W
      CASE ('FMG')
         TYPE_CYCLING = NSCARC_CYCLING_FMG
      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_CYCLE, NSCARC_NONE)
   END SELECT

   !> Set type of interpolation (STANDARD/DIRECT/MULTIPASS)
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
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_INTERPOL, NSCARC_NONE)
   END SELECT

ENDIF

!> Set type of coarse grid solver
SELECT CASE (TRIM(SCARC_COARSE))
   CASE ('ITERATIVE')
      TYPE_COARSE = NSCARC_COARSE_ITERATIVE
      TYPE_KRYLOV = NSCARC_KRYLOV_CG
   CASE ('DIRECT')
#ifdef WITH_MKL
      TYPE_COARSE   = NSCARC_COARSE_DIRECT
      TYPE_MKL(0)   = NSCARC_MKL_COARSE
#else
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_COARSE, NSCARC_NONE)
END SELECT

!>
!> Set type of accuracy (ABSOLUTE/RELATIVE)
!>
SELECT CASE (TRIM(SCARC_ACCURACY))
   CASE ('ABSOLUTE')
      TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
   CASE ('RELATIVE')
      TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_ACCURACY, NSCARC_NONE)
END SELECT

!>
!> Set type of precision for MKL solver (SINGLE/DOUBLE)
!>
SELECT CASE (TRIM(SCARC_MKL_PRECISION))
   CASE ('SINGLE')
      TYPE_MKL_PRECISION = NSCARC_PRECISION_SINGLE
   CASE ('DOUBLE')
      TYPE_MKL_PRECISION = NSCARC_PRECISION_DOUBLE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MKL_PRECISION, NSCARC_NONE)
END SELECT

!>
!> -------- define some logical variables - just for notational convenience
!>
IS_STRUCTURED   = TYPE_GRID == NSCARC_GRID_STRUCTURED
IS_UNSTRUCTURED = TYPE_GRID == NSCARC_GRID_UNSTRUCTURED

IS_CG     = TYPE_METHOD == NSCARC_METHOD_KRYLOV
IS_CG_GMG = IS_CG .AND. TYPE_PRECON == NSCARC_RELAX_GMG .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC

IS_MG  = TYPE_METHOD == NSCARC_METHOD_MULTIGRID
IS_GMG = IS_MG .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC

IS_CG_ADD    = HAS_TWO_LEVELS .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_ADD
IS_CG_MUL    = HAS_TWO_LEVELS .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_MUL
IS_CG_COARSE = HAS_TWO_LEVELS .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_COARSE

HAS_TWO_LEVELS   = IS_CG .AND. TYPE_PRECON /= NSCARC_RELAX_GMG .AND. TYPE_TWOLEVEL > NSCARC_TWOLEVEL_NONE
HAS_MULTI_LEVELS = IS_GMG .OR. IS_CG_GMG .OR. HAS_TWO_LEVELS

IS_FFT =  TYPE_PRECON == NSCARC_RELAX_FFT  .OR.  TYPE_SMOOTH == NSCARC_RELAX_FFT
IS_FFTO=  TYPE_PRECON == NSCARC_RELAX_FFTO .OR.  TYPE_SMOOTH == NSCARC_RELAX_FFTO
IS_MKL = (TYPE_PRECON >= NSCARC_RELAX_MKL) .OR. (TYPE_SMOOTH >= NSCARC_RELAX_MKL) .OR. &
   (HAS_MULTI_LEVELS .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT)

#ifdef WITH_SCARC_MGM
IS_MGM = TYPE_METHOD == NSCARC_METHOD_MGM
#endif

END SUBROUTINE SCARC_PARSE_INPUT


!> ------------------------------------------------------------------------------------------------
!> Determine number of grid levels  (1 for CG/BICG-method, NLEVEL for MG-method)
!> Note: NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LEVELS
#ifdef WITH_MKL
INTEGER :: NL
#endif

SELECT_SCARC_METHOD: SELECT CASE (TYPE_METHOD)

   !>
   !> ----------------------- Global data-parallel Krylov method --------------------------------------------
   !>
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

#ifdef WITH_MKL
         !>
         !> Preconditioning by defect correction based on LU-decomposition
         !> if two-level method, also use coarse grid level, otherwise only use single (finest) grid level
         !>
         CASE (NSCARC_RELAX_MKL)

            IF (HAS_TWO_LEVELS) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF

            !> Either using a global CLUSTER_SPARSE_SOLVER or local PARDISO solvers from MKL
            SELECT CASE (TYPE_SCOPE(1))
               CASE(NSCARC_SCOPE_GLOBAL)
                  TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
               CASE(NSCARC_SCOPE_LOCAL)
                  TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
            END SELECT

#endif

         !>
         !> Preconditioning by defect correction based on geometric multigrid method,
         !> use specified hierarchy of grid levels
         !>
         CASE (NSCARC_RELAX_GMG)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
            IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif

         !>
         !> Preconditioning by defect correction based on local basic iterations (JACOBI/SSOR),
         !> if two-level method, also use coarse grid, otherwise only use single (finest) grid level
         !>
         CASE DEFAULT
            IF (HAS_TWO_LEVELS) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF

      END SELECT SELECT_KRYLOV_PRECON

   !>
   !> ----------------------- Global data-parallel Multigrid method -------------------------------------
   !>
   CASE (NSCARC_METHOD_MULTIGRID)

      SELECT_MULTIGRID_TYPE: SELECT CASE (TYPE_MULTIGRID)

         !>
         !> Use of geometric multigrid method
         !> If not specified by user, determine number of possible grid levels
         !>
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)

#ifdef WITH_MKL
            !> in case of smoothing by different MKL solvers, mark levels for the use of MKL,
            !> either by locally acting PARDISO solvers or Globally acting CLUSTER_SPARSE_SOLVER
            IF (TYPE_SMOOTH == NSCARC_RELAX_MKL) THEN

               MKL_SCOPE_SELECT: SELECT CASE (TYPE_SCOPE(2))
               CASE (NSCARC_SCOPE_LOCAL)
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     TYPE_MKL(NL) = NSCARC_MKL_LOCAL
                  ENDDO
               CASE (NSCARC_SCOPE_GLOBAL)
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     TYPE_MKL(NL) = NSCARC_MKL_GLOBAL
                  ENDDO
               END SELECT MKL_SCOPE_SELECT

            ENDIF

            IF (TYPE_MKL(0) == NSCARC_MKL_COARSE) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif

         !>
         !> Use of algebraic multigrid method - currently disabled, still in progres
         !> first, only finest level is set, further levels are defined during coarsening process
         !>
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)

            CALL SCARC_SHUTDOWN(NSCARC_ERROR_AMG_DISABLED, SCARC_NONE, NSCARC_NONE)
            !CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_AMG)

      END SELECT SELECT_MULTIGRID_TYPE

   !>
   !> ----------------------- Global LU-decomposition -----------------------------------------
   !>
   CASE (NSCARC_METHOD_LU)

      !> Only use single (finest) grid level
      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)

      !> and mark this level for the use of MKL methods
      SELECT_MKL_SCOPE: SELECT CASE (TYPE_MKL(0))
         CASE (NSCARC_MKL_LOCAL)
            TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
         CASE (NSCARC_MKL_GLOBAL)
            TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
      END SELECT SELECT_MKL_SCOPE

#ifdef WITH_SCARC_MGM
   !>
   !> ----------------------- Global McKenney-Greengard-Mayo method - only finest level --------------------
   !>
   CASE (NSCARC_METHOD_MGM)

      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)

#endif

END SELECT SELECT_SCARC_METHOD


#ifdef WITH_MKL
!>
!> Define MKL related logical short names based on number of levels
!>
DO NL = NLEVEL_MIN, NLEVEL_MAX
   IS_MKL_LEVEL(NL) = (TYPE_MKL(0)  == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
                      (TYPE_MKL(0)  == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
                      (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL)
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_LEVELS


!> ------------------------------------------------------------------------------------------------
!> Setup single level in case of default Krylov method
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS(NTYPE)
USE SCARC_POINTERS, ONLY: M
INTEGER, INTENT(IN) :: NTYPE
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM, NLEVEL

SELECT_LEVEL_TYPE: SELECT CASE (NTYPE)

   !> only use finest grid level
   CASE(NSCARC_LEVEL_SINGLE)
   
      NLEVEL     = 1
      NLEVEL_MIN = 1
      NLEVEL_MAX = 1
   
   !>
   !> determine maximum number of possible levels based on number of grid cells (based on doubling)
   !>
   CASE(NSCARC_LEVEL_MULTI)
   
      NLEVEL = NSCARC_LEVEL_MAX
      KLEVEL = NSCARC_LEVEL_MAX
   
      DO NM=1,NMESHES

         M => MESHES(NM)

         KLEVEL(1)=SCARC_GET_MAX_LEVEL(M%IBAR,1)
         IF (.NOT.TWO_D) KLEVEL(2)=SCARC_GET_MAX_LEVEL(M%JBAR,2)
         KLEVEL(3)=SCARC_GET_MAX_LEVEL(M%KBAR,3)

         KLEVEL_MIN = MINVAL(KLEVEL)
         IF (KLEVEL_MIN<NLEVEL) NLEVEL=KLEVEL_MIN

      ENDDO
   
      NLEVEL_MIN  = 1
      IF (IS_GMG.OR.IS_CG_GMG) THEN

         IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
            NLEVEL_MAX  = NLEVEL_MIN + SCARC_MULTIGRID_LEVEL - 1
         ELSE
            NLEVEL_MAX  = NLEVEL
         ENDIF

      ELSE IF (HAS_TWO_LEVELS) THEN

         IF (SCARC_COARSE_LEVEL /= -1) THEN
            NLEVEL_MAX  = NLEVEL_MIN + SCARC_COARSE_LEVEL - 1
         ELSE
            NLEVEL_MAX  = NLEVEL
         ENDIF

      ENDIF
   
   !>
   !> use user specified number of grid levels - currently disabled, still in progress
   !>
   CASE(NSCARC_LEVEL_AMG)
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_AMG_DISABLED, SCARC_NONE, NSCARC_NONE)
      !   NLEVEL_MIN = 1
      !   IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
      !      NLEVEL_MAX  = SCARC_MULTIGRID_LEVEL
      !   ELSE
      !      NLEVEL_MAX  = NSCARC_LEVEL_MAX
      !   ENDIF
      !   NLEVEL = NLEVEL_MAX
   
END SELECT SELECT_LEVEL_TYPE

END SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS


!> ------------------------------------------------------------------------------------------------
!> Determine maximum number of possible levels on direction IOR0 of mesh NM
!> In case of GMG- or 2-Level-method, NC must be divisable by 2 at least one time
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_GET_MAX_LEVEL(NC, IOR0)
INTEGER, INTENT(IN) :: NC, IOR0
INTEGER :: NC0, NL

!> Print error message if not divisable by 2
IF (HAS_MULTI_LEVELS .AND. SCARC_MULTIGRID_LEVEL > 1 .AND. MOD(NC,2)/=0) THEN
   SELECT CASE (ABS(IOR0))
      CASE (1)
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBERX, SCARC_NONE, NC)
      CASE (2)
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBERY, SCARC_NONE, NC)
      CASE (3)
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBERZ, SCARC_NONE, NC)
   END SELECT
ENDIF

!> Divide by 2 as often as possible or until user defined max-level is reached
IF (SCARC_MULTIGRID_LEVEL > 1) THEN
   NC0=NC
   DO NL=1,NSCARC_LEVEL_MAX
      NC0=NC0/2
      IF (MOD(NC0,2)/=0) EXIT                !> if no longer divisable by two, leave loop ...
      IF (NL==SCARC_MULTIGRID_LEVEL) EXIT    !> if max possible number of levels reached, leave loop ...
      IF (NC0==1) EXIT                       !> if corresponding power of two has been found, leave loop ...
   ENDDO
   SCARC_GET_MAX_LEVEL=NL
ELSE
   SCARC_GET_MAX_LEVEL=NLEVEL_MIN
ENDIF

RETURN
END FUNCTION SCARC_GET_MAX_LEVEL


!> ------------------------------------------------------------------------------------------------
!> Allocate basic ScaRC-structures for all needed levels
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TYPES
USE SCARC_POINTERS, ONLY: S
INTEGER :: NM

!> Basic information for all requested grid levels
ALLOCATE (SCARC(NMESHES), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERROR)

!> Basic solver stack
ALLOCATE (STACK(NSCARC_STACK_MAX), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'STACK', IERROR)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MESH(NM)

   !> Needed information about other meshes
   ALLOCATE (S%OSCARC(NMESHES), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSCARC', IERROR)

   !> Information for single grid levels
   ALLOCATE (S%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'LEVEL', IERROR)

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_TYPES


!> ----------------------------------------------------------------------------------------------------
!> This routine assumes, that L already points to the correct level NL of mesh NL and
!> additionally sets the requested discretization type
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ASSIGN_GRID_TYPE(NTYPE)
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE (NTYPE)
   CASE (NSCARC_GRID_STRUCTURED)
      TYPE_GRID = NSCARC_GRID_STRUCTURED
      IS_STRUCTURED   = .TRUE.
      IS_UNSTRUCTURED = .FALSE.
   CASE (NSCARC_GRID_UNSTRUCTURED)
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED
      IS_STRUCTURED   = .FALSE.
      IS_UNSTRUCTURED = .TRUE.
END SELECT

END SUBROUTINE SCARC_ASSIGN_GRID_TYPE


!> -----------------------------------------------------------------------------
!> Setup discretization information
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRIDS
USE SCARC_POINTERS, ONLY: M, S, L, G, XCOR, YCOR, ZCOR, XMID, YMID, ZMID
INTEGER :: NL, NM, NC, IX, IY, IZ, IO
INTEGER :: IBAR, JBAR, KBAR

!>
!> -------------------- On all grid levels --------------------------------------
!> Specify general mesh related geometry information
!>

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MESH(NM)

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

   IBAR = M%IBAR
   JBAR = M%JBAR
   KBAR = M%KBAR

   LEVEL_LOOP1: DO NL = NLEVEL_MIN, NLEVEL_MAX

      CALL SCARC_POINT_TO_GRID(NM, NL)

      L%NX = IBAR
      L%NY = JBAR
      L%NZ = KBAR

      L%N_CELLS  = L%NX * L%NY * L%NZ

      L%N_WALL_CELLS_EXT = M%N_EXTERNAL_WALL_CELLS
      L%N_WALL_CELLS_INT = M%N_INTERNAL_WALL_CELLS

      L%N_WALL_CELLS = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT

      !> get coordination information
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
      IBAR=IBAR/2
      IF (.NOT.TWO_D) JBAR=JBAR/2
      KBAR=KBAR/2

      IF (NL == NLEVEL_MIN) THEN

         !> On finest level store information about obstructions
         L%N_OBST = M%N_OBST
         ALLOCATE(L%OBST(L%N_OBST), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_GRIDS','OBST',IERROR)

         DO IO = 1, L%N_OBST
            L%OBST(IO)%I1  = M%OBSTRUCTION(IO)%I1
            L%OBST(IO)%I2  = M%OBSTRUCTION(IO)%I2
            L%OBST(IO)%J1  = M%OBSTRUCTION(IO)%J1
            L%OBST(IO)%J2  = M%OBSTRUCTION(IO)%J2
            L%OBST(IO)%K1  = M%OBSTRUCTION(IO)%K1
            L%OBST(IO)%K2  = M%OBSTRUCTION(IO)%K2
         ENDDO

         !> Point to already existing arrays from main FDS program
         XCOR => M%X ; YCOR => M%Y ; ZCOR => M%Z
         XMID => M%XC; YMID => M%YC; ZMID => M%ZC

      ELSE

         !> Allocate and compute coordinate information for coarser levels
         CALL SCARC_ALLOCATE_REAL1(L%XCOR, 0, L%NX, NSCARC_INIT_NONE, 'XCOR')
         CALL SCARC_ALLOCATE_REAL1(L%YCOR, 0, L%NY, NSCARC_INIT_NONE, 'YCOR')
         CALL SCARC_ALLOCATE_REAL1(L%ZCOR, 0, L%NZ, NSCARC_INIT_NONE, 'ZCOR')

         DO IX = 0, L%NX
            L%XCOR(IX) = S%XS + IX*L%DX
         ENDDO
         DO IY = 0, L%NY
            L%YCOR(IY) = S%YS + IY*L%DY
         ENDDO
         DO IZ = 0, L%NZ
            L%ZCOR(IZ) = S%ZS + IZ*L%DZ
         ENDDO

         !> Allocate and compute midpoint information for coarser levels
         CALL SCARC_ALLOCATE_REAL1(L%XMID, 0, L%NX+1, NSCARC_INIT_NONE, 'XMID')
         CALL SCARC_ALLOCATE_REAL1(L%YMID, 0, L%NY+1, NSCARC_INIT_NONE, 'YMID')
         CALL SCARC_ALLOCATE_REAL1(L%ZMID, 0, L%NZ+1, NSCARC_INIT_NONE, 'ZMID')

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

         XCOR => L%XCOR
         YCOR => L%YCOR
         ZCOR => L%ZCOR

         XMID => L%XMID
         YMID => L%YMID
         ZMID => L%ZMID

      ENDIF

      !> Allocate vectors for step sizes in different directions
      CALL SCARC_ALLOCATE_REAL1(L%DXL, 0, L%NX, NSCARC_INIT_ZERO, 'DXL')
      CALL SCARC_ALLOCATE_REAL1(L%DYL, 0, L%NY, NSCARC_INIT_ZERO, 'DYL')
      CALL SCARC_ALLOCATE_REAL1(L%DZL, 0, L%NZ, NSCARC_INIT_ZERO, 'DZL')

      !> set step sizes between cell midpoints, use interior step sizes for ghost cells as initial values
      !> correct sizes for ghost cells are exchanged later
      DO IX = 1, L%NX-1
         L%DXL(IX) = XMID(IX+1) - XMID(IX)
      ENDDO
      L%DXL(0)    = L%DXL(1)
      L%DXL(L%NX) = L%DXL(L%NX-1)

      DO IY = 1, L%NY-1
         L%DYL(IY) = YMID(IY+1) - YMID(IY)
      ENDDO
      L%DYL(0)    = L%DYL(1)
      L%DYL(L%NY) = L%DYL(L%NY-1)

      DO IZ = 1, L%NZ-1
         L%DZL(IZ) = ZMID(IZ+1) - ZMID(IZ)
      ENDDO
      L%DZL(0)    = L%DZL(1)
      L%DZL(L%NZ) = L%DZL(L%NZ-1)

   ENDDO LEVEL_LOOP1
ENDDO MESHES_LOOP1

!>
!> ---------------------- On finest grid level -------------------------------------------------
!> Allocate several arrays for the administration of discretization related data
!>
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_LEVEL(NM, NLEVEL_MIN)

   !> Set pointers to already existing CELL_INDEX WALL_INDEX arrays from main program (on finest level)
   CALL SCARC_SETUP_CELL_INDEX(NM, NLEVEL_MIN)
   CALL SCARC_SETUP_WALL_INDEX(NM, NLEVEL_MIN)

   !> Allocate and initialize IS_SOLID array which indicates the state of a cell (gasphase/solid)
   CALL SCARC_ALLOCATE_LOG3(L%IS_SOLID, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_TRUE, 'IS_SOLID')
   L%IS_SOLID (1:L%NX, 1:L%NY, 1:L%NZ) = .FALSE.

   !> Identify and mark solid obstruction cells in IS_SOLID-part of the discretization
   DO IZ = 1, L%NZ
      DO IY = 1, L%NY
         DO IX = 1, L%NX
            IF (M%SOLID(M%CELL_INDEX(IX, IY, IZ))) L%IS_SOLID(IX, IY, IZ) = .TRUE.
         ENDDO
      ENDDO
   ENDDO

   !>
   !> If both discretization types (structured/unstructured) must be administrated (MGM method only):
   !> Allocate all arrays which are related to a specific discretization type
   !>
   IF (HAS_MULTI_GRIDS) THEN

      !>
      !> ---------------- First process structured discretization
      !>
      CALL SCARC_ASSIGN_GRID_TYPE(NSCARC_GRID_STRUCTURED)
      CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)

      !> Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET')

      !> Allocate wall information array
      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      !> Allocate and preset cell numbers array
      CALL SCARC_ALLOCATE_INT3(G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'CELL_NUMBER')

      !> Define local cell numbers for Poisson equation
      DO IZ=1,L%NZ
         DO IY=1,L%NY
            DO IX=1,L%NX
               G%NC_LOCAL(NM) = G%NC_LOCAL(NM) + 1
               G%CELL_NUMBER(IX,IY,IZ) = G%NC_LOCAL(NM)
            ENDDO
         ENDDO
      ENDDO
      G%NC  = G%NC_LOCAL(NM)
      G%NCE = G%NC_LOCAL(NM)

      !>
      !> ---------------- Then process unstructured discretization
      !>
      CALL SCARC_ASSIGN_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)
      CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)

      !> Also allocate and preset cell numbers and state arrays for unstructured discretization
      CALL SCARC_ALLOCATE_INT3(G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'CELL_NUMBER')

      !> Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET')

      !> Allocate wall information array
      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      !> Define local cell numbers for Poisson equation
      DO IZ=1,L%NZ
         DO IY=1,L%NY
            DO IX=1,L%NX
               IF (.NOT.L%IS_SOLID(IX,IY,IZ)) THEN
                  G%NC_LOCAL(NM) = G%NC_LOCAL(NM) + 1
                  G%CELL_NUMBER(IX,IY,IZ) = G%NC_LOCAL(NM)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      G%NC  = G%NC_LOCAL(NM)
      G%NCE = G%NC_LOCAL(NM)

   !>
   !> If only one specified type of discretization must be admistrated:
   !> Allocate and preset cell numbers and state arrays for requested type of discretization
   !>
   ELSE

      !> ---------------- Only process specified type of discretization
      CALL SCARC_ASSIGN_GRID_TYPE(TYPE_GRID)
      CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)

      !> Also allocate and preset cell numbers and state arrays for unstructured discretization
      CALL SCARC_ALLOCATE_INT3(G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'CELL_NUMBER')

      !> Allocate index array which specifies I, J, K components for all degrees of freedom
      NC = L%NX * L%NY * L%NZ
      CALL SCARC_ALLOCATE_INT1(G%ICX, 1, NC, NSCARC_INIT_UNDEF, 'ICX')
      CALL SCARC_ALLOCATE_INT1(G%ICY, 1, NC, NSCARC_INIT_UNDEF, 'ICY')
      CALL SCARC_ALLOCATE_INT1(G%ICZ, 1, NC, NSCARC_INIT_UNDEF, 'ICZ')

      !> Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET')

      !> Allocate wall information array
      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      !> Define local cell numbers for Poisson equation
      DO IZ=1,L%NZ
         DO IY=1,L%NY
            DO IX=1,L%NX
               IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX,IY,IZ)) THEN

                  G%NC_LOCAL(NM) = G%NC_LOCAL(NM) + 1
                  G%CELL_NUMBER(IX, IY, IZ) = G%NC_LOCAL(NM)

                  G%ICX(G%NC_LOCAL(NM)) = IX
                  G%ICY(G%NC_LOCAL(NM)) = IY
                  G%ICZ(G%NC_LOCAL(NM)) = IZ

               ENDIF
            ENDDO
         ENDDO
      ENDDO
      G%NC  = G%NC_LOCAL(NM)
      G%NCE = G%NC_LOCAL(NM)

   ENDIF

ENDDO MESHES_LOOP2

END SUBROUTINE SCARC_SETUP_GRIDS


!> -----------------------------------------------------------------------------
!> Setup discretization information on coarser levels
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRID_LEVEL(NL)
USE SCARC_POINTERS, ONLY: LF, LC, GC
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NC, IXF, IYF, IZF, IX, IY, IZ, NSTEP

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LF => SCARC(NM)%LEVEL(NLEVEL_MIN)
   LC => SCARC(NM)%LEVEL(NL)

   CALL SCARC_ALLOCATE_LOG3(LC%IS_SOLID, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_NONE, 'IS_SOLID')
   LC%IS_SOLID (1:LC%NX, 1:LC%NY, 1:LC%NZ)  = .FALSE.

   NSTEP = 2**(NL - NLEVEL_MIN)

   SELECT CASE(TYPE_GRID)

      !>
      !> Get cell numberings for coarser grid in case of structured discretization
      !>
      CASE (NSCARC_GRID_STRUCTURED)

         GC => LC%STRUCTURED

         CALL SCARC_ALLOCATE_INT1(GC%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL')
         CALL SCARC_ALLOCATE_INT1(GC%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET')

         CALL SCARC_ALLOCATE_INT3(GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, &
                                  NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER')

         NC = LC%NX * LC%NY * LC%NZ

         CALL SCARC_ALLOCATE_INT1(GC%ICX , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICX')
         CALL SCARC_ALLOCATE_INT1(GC%ICY , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICY')
         CALL SCARC_ALLOCATE_INT1(GC%ICZ , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICZ')

         GC%CELL_NUMBER = NSCARC_UNDEF_INT

         DO IZ = 1, LC%NZ
            IZF = (IZ-1)*NSTEP + 1
            DO IY = 1, LC%NY
               IYF = (IY-1)*NSTEP + 1
               DO IX = 1, LC%NX

                  IXF = (IX-1)*NSTEP + 1
                  LC%IS_SOLID(IX,IY,IZ) = LF%IS_SOLID(IXF, IYF, IZF)

                  GC%NC_LOCAL(NM)  = GC%NC_LOCAL(NM) + 1
                  GC%CELL_NUMBER(IX,IY,IZ) = GC%NC_LOCAL(NM)

                  GC%ICX(GC%NC_LOCAL(NM)) = IX
                  GC%ICY(GC%NC_LOCAL(NM)) = IY
                  GC%ICZ(GC%NC_LOCAL(NM)) = IZ

               ENDDO
            ENDDO
         ENDDO

         GC%NC = GC%NC_LOCAL(NM)

      !>
      !> Get cell numberings for coarser grid in case of unstructured discretization
      !>
      CASE (NSCARC_GRID_UNSTRUCTURED)

         GC => LC%UNSTRUCTURED

         CALL SCARC_ALLOCATE_INT1(GC%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL')
         CALL SCARC_ALLOCATE_INT1(GC%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET')

         CALL SCARC_ALLOCATE_INT3(GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, &
                                  NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER')

         NC = LC%NX * LC%NY * LC%NZ

         CALL SCARC_ALLOCATE_INT1(GC%ICX , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICX')
         CALL SCARC_ALLOCATE_INT1(GC%ICY , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICY')
         CALL SCARC_ALLOCATE_INT1(GC%ICZ , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICZ')

         GC%CELL_NUMBER = NSCARC_UNDEF_INT

         DO IZ = 1, LC%NZ
            IZF = (IZ-1)*NSTEP + 1
            DO IY = 1, LC%NY
               IYF = (IY-1)*NSTEP + 1
               DO IX = 1, LC%NX

                  IXF = (IX-1)*NSTEP + 1
                  LC%IS_SOLID(IX,IY,IZ) = LF%IS_SOLID(IXF, IYF, IZF)

                  IF (.NOT.LF%IS_SOLID(IXF, IYF, IZF)) THEN

                     GC%NC_LOCAL(NM) = GC%NC_LOCAL(NM) + 1
                     GC%CELL_NUMBER(IX,IY,IZ) = GC%NC_LOCAL(NM)

                     GC%ICX(GC%NC_LOCAL(NM)) = IX
                     GC%ICY(GC%NC_LOCAL(NM)) = IY
                     GC%ICZ(GC%NC_LOCAL(NM)) = IZ

                  ENDIF

               ENDDO
            ENDDO
         ENDDO

         GC%NC = GC%NC_LOCAL(NM)

   END SELECT

ENDDO MESHES_LOOP1

END SUBROUTINE SCARC_SETUP_GRID_LEVEL


!> ------------------------------------------------------------------------------------------------
!> Setup communication structure for data exchange along mesh interfaces
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_NEIGHBORS
USE SCARC_POINTERS, ONLY: OS, OLF, OLC
INTEGER :: NM, NOM, NL

!> Initialize communication counter for ScaRC, use same TAG for all communications
TAG   = 99
N_REQ =  0
N_EXCHANGES = 0

!>
!> Initialize level structures on neighboring meshes
!>
LEVEL_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LEVEL_NEIGHBOR_LOOP: DO NOM = 1, NMESHES

      !>
      !> On finest level point to exchange structures from surrounding FDS 
      !>
      IF (.NOT. ARE_NEIGHBORS(NM, NOM)) CYCLE LEVEL_NEIGHBOR_LOOP

      N_EXCHANGES = N_EXCHANGES+1                                         !> count number of exchanges

      OS => SCARC(NM)%OSCARC(NOM)
      ALLOCATE (OS%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)             !> allocate neighboring structures
      CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OS%LEVEL', IERROR)

      OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)                      !> point to neighbor on finest grid level

      OLF%NX = MESHES(NOM)%IBAR                                           !> number of cells in x-direction on neighbor
      OLF%NY = MESHES(NOM)%JBAR                                           !> number of cells in y-direction on neighbor
      OLF%NZ = MESHES(NOM)%KBAR                                           !> number of cells in z-direction on neighbor

      OLF%N_WALL_CELLS_EXT = MESHES(NOM)%N_EXTERNAL_WALL_CELLS            !> number of external wall cells on neighbor
      OLF%N_WALL_CELLS_INT = MESHES(NOM)%N_INTERNAL_WALL_CELLS            !> number of external wall cells on neighbor
      OLF%N_WALL_CELLS     = OLF%N_WALL_CELLS_EXT + OLF%N_WALL_CELLS_INT  !> number of walls cell on neighbor

      OLF%N_CELLS = OLF%NX*OLF%NY*OLF%NZ                                  !> number of cells on neighbor (structured)

      !>
      !> In case of GMG with a predefined grid hierarchy define corresponding level-structures
      !>
      IF (NLEVEL_MAX > NLEVEL_MIN) THEN                                   

         DO NL=NLEVEL_MIN+1,NLEVEL_MAX

            OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)                        !> OLF points to finer, OLC to coarser level

            OLC%NX = OLF%NX/2                                             !> use double grid width
            IF (TWO_D) THEN
               OLC%NY = 1
            ELSE
               OLC%NY = OLF%NY/2
            ENDIF
            OLC%NZ = OLF%NZ/2

            OLC%N_CELLS          = OLC%NX * OLC%NY * OLC%NZ               !> set new number of cells
            OLC%N_WALL_CELLS     = OLC%N_WALL_CELLS_EXT                   !> set new number of wall cells
            OLC%N_WALL_CELLS_EXT = 2 * (OLC%NX*OLC%NZ + OLC%NX*OLC%NY + OLC%NY*OLC%NZ)    !> TODO: CHECK!

         ENDDO
      ENDIF

   ENDDO LEVEL_NEIGHBOR_LOOP
ENDDO LEVEL_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_NEIGHBORS


!> ----------------------------------------------------------------------------------------------------
!> Setup FACE related structures on finest grid level:
!>   - get dimensions for each of the 6 faces of a mesh
!>   - get grid width vector along face
!>   - get information for adjacent neighbors
!>   - allocate pointer arrays for data exchanges with neighbors
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACES
USE SCARC_POINTERS, ONLY: M, S, L, LC, F
INTEGER :: NL, NM, NOM
INTEGER :: IFACE, IOR0, JOR0, INBR, IWG, ICW
LOGICAL :: IS_KNOWN(-3:3)
INTEGER :: FACE_NEIGHBORS(-3:3, NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: MESH_NEIGHBORS(6*NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: N_FACE_NEIGHBORS(-3:3)
INTEGER :: N_MESH_NEIGHBORS

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_LEVEL(NM, NLEVEL_MIN)             !> consider only finest grid level

   !> Allocate FACE arrays on different grid levels
   ALLOCATE(L%FACE(-3:3), STAT=IERROR)
   CALL ChkMemErr('SCARC_SETUP_FACES','FACE',IERROR)

   IF (NLEVEL_MAX > NLEVEL_MIN) THEN
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX
         LC => SCARC(NM)%LEVEL(NL)
         ALLOCATE(LC%FACE(-3:3), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_FACES','FACE',IERROR)
      ENDDO
   ENDIF

   FACE_NEIGHBORS = -1
   MESH_NEIGHBORS = -1

   N_FACE_NEIGHBORS = 0
   N_MESH_NEIGHBORS = 0

   CALL SCARC_SETUP_FACE_BASICS(NM, NLEVEL_MIN)

   !>
   !> Store first wall cell number for each face
   !>
   ICW = 1
   FACE_ORDER_LOOP: DO IFACE = 1, 6
      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)
      F%NCW0 = ICW
      ICW = ICW + F%NCW
   ENDDO FACE_ORDER_LOOP

   !>
   !> loop over external wall cells:
   !> store basic data and determine number of adajacent neighbors to each face
   !>
   EXTERNAL_WALL_CELLS_LOOP: DO IWG = 1, L%N_WALL_CELLS_EXT

      NOM  = M%EXTERNAL_WALL(IWG)%NOM
      IOR0 = M%WALL(IWG)%ONE_D%IOR

      IF (NOM /= 0) THEN
         IS_KNOWN = .FALSE.
         DO JOR0 = -3, 3                                              !> neighbor already known?
            IF (JOR0 == 0) CYCLE
            DO INBR = 1, N_FACE_NEIGHBORS(JOR0)
               IF (FACE_NEIGHBORS(JOR0, INBR) == NOM) THEN
                  IS_KNOWN(JOR0) = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDDO
         IF (.NOT.IS_KNOWN(IOR0)) THEN
            N_FACE_NEIGHBORS(IOR0) = N_FACE_NEIGHBORS(IOR0) + 1       !> increase neighbor counter for face
            FACE_NEIGHBORS(IOR0, N_FACE_NEIGHBORS(IOR0)) = NOM        !> store number of neighbor for face
         ENDIF
         IF (.NOT.ANY(IS_KNOWN)) THEN
            N_MESH_NEIGHBORS = N_MESH_NEIGHBORS + 1                   !> increase neighbor counter for mesh
            MESH_NEIGHBORS(N_FACE_NEIGHBORS(IOR0)) = NOM              !> store number of neighbor for mesh
         ENDIF
      ENDIF

   ENDDO EXTERNAL_WALL_CELLS_LOOP

   !>
   !> Allocate array which stores numbers of all neighboring meshes
   !>
   IF (N_MESH_NEIGHBORS /= 0) &
      CALL SCARC_ALLOCATE_INT1(S%NEIGHBORS, 1, N_MESH_NEIGHBORS, NSCARC_INIT_UNDEF, 'NEIGHBORS')
   S%N_NEIGHBORS = N_MESH_NEIGHBORS

   NEIGHBORS_OF_FACE_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE NEIGHBORS_OF_FACE_LOOP

      !> if there are neighbors at face IOR0 store information about them
      F => L%FACE(IOR0)
      IF (N_FACE_NEIGHBORS(IOR0) /= 0) THEN

         !> allocate array for storing the numbers of the single neighbors
         F%N_NEIGHBORS = N_FACE_NEIGHBORS(IOR0)
         CALL SCARC_ALLOCATE_INT1(F%NEIGHBORS, 1, N_FACE_NEIGHBORS(IOR0), NSCARC_INIT_NONE, 'NEIGHBORS')

         !> store every neighbor and allocate corresponding administration arrays on finest level
         DO INBR = 1, N_FACE_NEIGHBORS(IOR0)

            NOM = FACE_NEIGHBORS(IOR0, INBR)
            F%NEIGHBORS(INBR) = NOM                          !> store NOM as a neighbor of that face and if
            CALL SCARC_STORE_NEIGHBOR(NM, NOM)               !> not already done also as mesh neighbor itself

            CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NLEVEL_MIN)
            IF (.NOT.ALLOCATED(OL%FACE)) THEN
               ALLOCATE(OL%FACE(-3:3), STAT=IERROR)
               CALL ChkMemErr('SCARC_SETUP_FACES','OL%FACE',IERROR)
            ENDIF

         ENDDO

      ENDIF
   ENDDO NEIGHBORS_OF_FACE_LOOP

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FACES

!> ----------------------------------------------------------------------------------------------------
!> Determine basic data for single faces (orientation, dimensions, numbers)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE_BASICS(NM, NL)
USE SCARC_POINTERS, ONLY: L, F
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: IOR0

L => SCARC(NM)%LEVEL(NL)

FACES_OF_MESH_LOOP: DO IOR0 = -3, 3

   IF (IOR0 == 0) CYCLE FACES_OF_MESH_LOOP

   
   F => L%FACE(IOR0)
   
   SELECT CASE (ABS(IOR0))

      !> 
      !> ---------- Faces in x-direction
      !> 
      CASE (1)

         F%NOP =  L%NX                           !> number of cells between opposite mesh faces

         F%NX  = 1                               !> number of cells in x-direction
         F%NY  = L%NY                            !> number of cells in y-direction
         F%NZ  = L%NZ                            !> number of cells in z-direction

         F%NCW = L%NY*L%NZ                       !> number of wall cells at that face

         F%DH  => L%DXL                          !> step size vector between opposite mesh faces

         F%SCAL_BOUNDARY  = L%DXI2               !> contribution due to boundary condition 
         F%SCAL_DIRICHLET = -2.0_EB * L%DXI2

         IF (IOR0 > 0) THEN
            F%SCAL_NEUMANN = L%DXI
            F%SCAL_FACE    = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
            F%SCAL_INSIDE  = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
         ELSE
            F%SCAL_NEUMANN = -L%DXI
            F%SCAL_FACE    = 2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
            F%SCAL_INSIDE  = 2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
         ENDIF     

      !> 
      !> ---------- Faces in y-direction
      !> 
      CASE (2)

         F%NOP =  L%NY                   !> dito

         F%NX  = L%NX
         F%NY  = 1
         F%NZ  = L%NZ

         F%NCW = L%NX*L%NZ

         F%DH  => L%DYL

         F%SCAL_BOUNDARY  = L%DYI2
         F%SCAL_DIRICHLET = -2.0_EB * L%DYI2

         IF (IOR0>0) THEN
            F%SCAL_NEUMANN = L%DYI
            IF (.NOT.TWO_D) THEN
               F%SCAL_FACE   = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
               F%SCAL_INSIDE = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
            ENDIF
         ELSE
            F%SCAL_NEUMANN = -L%DYI
            IF (.NOT.TWO_D) THEN
               F%SCAL_FACE   =  2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
               F%SCAL_INSIDE =  2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
            ENDIF
         ENDIF

      !> 
      !> ---------- Faces in z-direction
      !> 
      CASE (3)

         F%NOP =  L%NZ                   !> dito

         F%NX  = L%NX
         F%NY  = L%NY
         F%NZ  = 1

         F%NCW = L%NX*L%NY

         F%DH  => L%DZL

         F%SCAL_BOUNDARY  = L%DZI2
         F%SCAL_DIRICHLET = -2.0_EB * L%DZI2

         IF (IOR0>0) THEN
            F%SCAL_NEUMANN = L%DZI
            F%SCAL_FACE    = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
            F%SCAL_INSIDE  = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
         ELSE
            F%SCAL_NEUMANN = -L%DZI
            F%SCAL_FACE    =  2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
            F%SCAL_INSIDE  =  2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
         ENDIF
   END SELECT

ENDDO FACES_OF_MESH_LOOP

END SUBROUTINE SCARC_SETUP_FACE_BASICS


!> ----------------------------------------------------------------------------------------------------
!> Setup WALL related structures and boundary conditions
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLS
USE SCARC_POINTERS, ONLY: M, L, LF, LC, FF, FC, OL, OLF, OLC, G, GC, GF, WC, WF, &
                          OGC, OGF, GWC, MWC, EWC
INTEGER :: NL, NM, NOM
INTEGER :: IREFINE, IFACE, IOR0, JOR0, INBR, IWG, IWC, ICW
LOGICAL :: IS_KNOWN(-3:3), IS_DIRIC, IS_OPEN

!>
!> -------- Get dimensionings for wall cells
!>
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)

   !>
   !> First loop over external wall cells:
   !> Determine number of adajacent neighbors to each face with corresponding number of IW's
   !>
   EXTERNAL_WALL_CELLS_LOOP1: DO IWG = 1, L%N_WALL_CELLS_EXT

      MWC => M%WALL(IWG)
      EWC => M%EXTERNAL_WALL(IWG)

      !> Determine and store neighbors, orientation and number of couplings for a single wall cell
      NOM  =  EWC%NOM
      IOR0 =  MWC%ONE_D%IOR

      GWC => G%WALL(IWG)

      GWC%NOM  = NOM                                    !> store number of neighbor in wall cell
      GWC%IOR  = IOR0                                   !> store orientation of that cell

      IF (NOM /= 0) THEN

         IS_KNOWN = .FALSE.

         DO JOR0 = -3, 3
            IF (JOR0 == 0) CYCLE
            DO INBR = 1, L%FACE(JOR0)%N_NEIGHBORS
               IF (L%FACE(JOR0)%NEIGHBORS(INBR) == NOM) THEN
                  IS_KNOWN(JOR0) = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         G%NCE = G%NCE + 1                                            !> increase number of extended grid cells

         CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NLEVEL_MIN)

         IF (ANY(IS_KNOWN)) OG%NCG = OG%NCG + 1                       !> increase counter for local ghost cells
         IF (OL%NCG0(IOR0) == 0) OL%NCG0(IOR0) = OG%NCG               !> save first ghost cell for -IOR0
         OL%NCG(IOR0) = OG%NCG                                     

      ENDIF

   ENDDO EXTERNAL_WALL_CELLS_LOOP1

   !>
   !> Then process internal wall cells
   !>
   INTERNAL_WALL_CELLS_LOOP1: DO IWG = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT

      MWC => M%WALL(IWG)
      GWC => G%WALL(IWG)

      GWC%IOR  = MWC%ONE_D%IOR
      GWC%NOM  = 0

      GWC%BTYPE = NEUMANN
      GWC%BOUNDARY_TYPE = M%WALL(IWG)%BOUNDARY_TYPE

      GWC%IXG =  MWC%ONE_D%II                        !> ghost cell indices
      GWC%IYG =  MWC%ONE_D%JJ
      GWC%IZG =  MWC%ONE_D%KK

      GWC%IXW =  MWC%ONE_D%IIG                       !> (internal) wall cell indices
      GWC%IYW =  MWC%ONE_D%JJG
      GWC%IZW =  MWC%ONE_D%KKG

   ENDDO INTERNAL_WALL_CELLS_LOOP1

   !>
   !> Allocate corresponding pointer arrays for data exchanges with neighbors
   !>
   IF (G%NCE > G%NC) THEN
      CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_IWG, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'ICE_TO_IWG')
      CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_ICG, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'ICE_TO_ICG')
      CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_ICN, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'ICE_TO_ICN')
      CALL SCARC_ALLOCATE_REAL1(G%ICE_TO_VAL, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'ICE_TO_VAL')
   ENDIF

   FACE_NEIGHBORS_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE FACE_NEIGHBORS_LOOP

      DO INBR = 1, L%FACE(IOR0)%N_NEIGHBORS

         NOM = L%FACE(IOR0)%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NLEVEL_MIN)

         CALL SCARC_ALLOCATE_INT1(OG%ICG_TO_IWG, 1, OG%NCG, NSCARC_INIT_ZERO, 'ICG_TO_IWG')
         CALL SCARC_ALLOCATE_INT1(OG%ICG_TO_ICE, 1, OG%NCG, NSCARC_INIT_ZERO, 'ICG_TO_ICE')
         CALL SCARC_ALLOCATE_INT1(OG%ICG_TO_ICW, 1, OG%NCG, NSCARC_INIT_ZERO, 'ICG_TO_ICW')
         CALL SCARC_ALLOCATE_INT1(OG%ICG_TO_ICW, 1, OG%NCG, NSCARC_INIT_ZERO, 'ICG_TO_ICW')

      ENDDO

   ENDDO FACE_NEIGHBORS_LOOP

   !>
   !> Second loop over external wall cells:
   !> Store detailed coordinate and cell data and get type of boundary condition
   !>
   G%ICE = G%NC
   WALL_CELLS_LOOP2: DO IWG = 1, L%N_WALL_CELLS_EXT

      NOM  = G%WALL(IWG)%NOM
      IOR0 = G%WALL(IWG)%IOR

      MWC => M%WALL(IWG)
      EWC => M%EXTERNAL_WALL(IWG)

      !>
      !> Preset ScaRC's boundary type indicator BTYPE
      !> INTERNAL  : the global Poisson problem is solved, so no BC's along mesh interfaces are needed
      !> DIRICHLET : - in the structured case face-wise BC-settings are used ccording to original FFT-solver
      !>               (this also allows to use FFT as local preconditioner)
      !>             - in the unstructured case Dirichlet BC's are only used for open boundary cells
      !> NEUMANN   : is used for the rest
      !>
      IS_DIRIC = MWC%PRESSURE_BC_INDEX == DIRICHLET
      IS_OPEN  = MWC%BOUNDARY_TYPE     == OPEN_BOUNDARY

      GWC => G%WALL(IWG)

      IF (EWC%NOM /= 0) THEN
         GWC%BTYPE = INTERNAL
      ELSE IF ((IS_STRUCTURED .AND. IS_DIRIC) .OR. (IS_UNSTRUCTURED .AND. IS_OPEN)) THEN
         GWC%BTYPE = DIRICHLET
         G%N_DIRIC = G%N_DIRIC + 1
      ELSE
         GWC%BTYPE = NEUMANN
         G%N_NEUMANN = G%N_NEUMANN + 1
      ENDIF

      GWC%BOUNDARY_TYPE = MWC%BOUNDARY_TYPE

      GWC%IXG = MWC%ONE_D%II                                 !> ghost cell indices
      GWC%IYG = MWC%ONE_D%JJ
      GWC%IZG = MWC%ONE_D%KK

      GWC%IXW = MWC%ONE_D%IIG                                !> (internal) wall cell indices
      GWC%IYW = MWC%ONE_D%JJG
      GWC%IZW = MWC%ONE_D%KKG

      !> If there exists a neighbor for that wall cell, setup corresponding neighborship information
      IF (NOM /= 0) THEN
         CALL SCARC_SETUP_WALL_NEIGHBOR(EWC%IIO_MIN, EWC%IIO_MAX, &
                                        EWC%JJO_MIN, EWC%JJO_MAX, &
                                        EWC%KKO_MIN, EWC%KKO_MAX, &
                                        IWG, NM, NOM, NLEVEL_MIN)
      ENDIF

   ENDDO WALL_CELLS_LOOP2

ENDDO MESHES_LOOP1

!>
!> Set dimensions on finest level for requested type(s) of discretization
!>
CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

!>
!> -------- For multi-level variants get discretization information and dimensions on coarser levels
!>
DO NL = NLEVEL_MIN+1, NLEVEL_MAX
   CALL SCARC_SETUP_GRID_LEVEL(NL)
   CALL SCARC_SETUP_DIMENSIONS(NL)
ENDDO

!>
!> -------- Check whether there are no Dirichlet BC's available - TODO: Check !!!
!>
MESH_INT = 0                            
RANK_INT = 0

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)
   MESH_INT(NM) = G%N_DIRIC   
   RANK_INT = RANK_INT + MESH_INT(NM)
ENDDO

IF (N_MPI_PROCESSES>1) &
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_INT, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERROR)
N_DIRIC_GLOBAL(NLEVEL_MIN) = RANK_INT

IS_PURE_NEUMANN = N_DIRIC_GLOBAL(NLEVEL_MIN) == 0 .AND. &
                  (TYPE_PRECON /= NSCARC_RELAX_FFT .OR. TYPE_PRECON /= NSCARC_RELAX_FFTO)


!>
!> -------- Only for multi-level variants 
!> (twolevel-CG or GMG-method as main solver or preconditioner):
!> Determine WALL, FACE and OSCARC types for coarser levels
!>
MULTI_LEVEL_IF: IF (HAS_MULTI_LEVELS) THEN

   MESHES_LOOP3: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      CALL SCARC_POINT_TO_MESH(NM)

      IREFINE=1
      LEVEL_GMG_LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX

         LF => SCARC(NM)%LEVEL(NL-1)                 !> LF points to finer level
         LC => SCARC(NM)%LEVEL(NL)                   !> LC points to coarser level

         SELECT CASE(TYPE_GRID)
            CASE (NSCARC_GRID_STRUCTURED)
               GF => LF%STRUCTURED
               GC => LC%STRUCTURED
            CASE (NSCARC_GRID_UNSTRUCTURED)
               GF => LF%UNSTRUCTURED
               GC => LC%UNSTRUCTURED
         END SELECT

         IREFINE=IREFINE*2
         CALL SCARC_CHECK_DIVISIBILITY(GF%NCE-GF%NC, 'GF%NCE')

         !>
         !> Initialize counts for overlapping and wall cells
         !>
         GC%NCE = GC%NC + (GF%NCE-GF%NC)/2
         GC%ICE = GC%NC

         LC%N_WALL_CELLS_EXT = SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
         LC%N_WALL_CELLS_INT = SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)
         LC%N_WALL_CELLS     = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT

         SELECT CASE(TYPE_GRID)
            CASE (NSCARC_GRID_STRUCTURED)
               GC%NW = LC%N_WALL_CELLS_EXT
            CASE (NSCARC_GRID_UNSTRUCTURED)
               GC%NW = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT
         END SELECT

         ALLOCATE(GC%WALL(LC%N_WALL_CELLS), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_NEIGHBORS','WALL',IERROR)

         !>
         !> First allocate administrative mapping arrays for own mesh
         !>
         IF (GC%NCE > GC%NC) THEN
            CALL SCARC_ALLOCATE_INT1 (GC%ICE_TO_IWG, GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'ICE_TO_IWG')
            CALL SCARC_ALLOCATE_INT1 (GC%ICE_TO_ICG, GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'ICE_TO_ICG')
            CALL SCARC_ALLOCATE_INT1 (GC%ICE_TO_ICN, GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'ICE_TO_ICN')
            CALL SCARC_ALLOCATE_REAL1(GC%ICE_TO_VAL, GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'ICE_TO_VAL')
         ENDIF

         !>
         !> set basic and wall information for all faces of coarser level
         !>
         CALL SCARC_SETUP_FACE_BASICS(NM, NL)

         WC => GC%WALL
         WF => GF%WALL

         IWC = 1
         IWG = 1
         ICW = 1
         DO IFACE = 1, 6

            IOR0 = FACE_ORIENTATION(IFACE)

            FF => LF%FACE(IOR0)
            FC => LC%FACE(IOR0)

            !> initialize FACE type for coarser mesh
            FC%NCW0 = ICW
            FC%N_NEIGHBORS = FF%N_NEIGHBORS

            IF (FC%N_NEIGHBORS /= 0) &
               CALL SCARC_ALLOCATE_INT1(FC%NEIGHBORS, 1, FC%N_NEIGHBORS, NSCARC_INIT_NONE, 'FACE_NEIGHBORS')
            DO INBR= 1, FC%N_NEIGHBORS
               FC%NEIGHBORS(INBR) = FF%NEIGHBORS(INBR)
            ENDDO
            
            FC%NCW = FC%NX * FC%NY * FC%NZ                                !> get number of wall cells for that face
            ICW = ICW + FC%NCW                                            !> increase global wall cell counter

            !> get related data and pointer structures for every mesh neighbor
            IF (LF%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
               DO INBR = 1, LF%FACE(IOR0)%N_NEIGHBORS

                  NOM = LF%FACE(IOR0)%NEIGHBORS(INBR)

                  OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL-1)
                  OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

                  SELECT CASE(TYPE_GRID)
                     CASE (NSCARC_GRID_STRUCTURED)
                        OGF => OLF%STRUCTURED
                        OGC => OLC%STRUCTURED
                     CASE (NSCARC_GRID_UNSTRUCTURED)
                        OGF => OLF%UNSTRUCTURED
                        OGC => OLC%UNSTRUCTURED
                  END SELECT

                  CALL SCARC_CHECK_DIVISIBILITY(OLF%N_WALL_CELLS_LOCAL, 'OLF%N_WALL_CELLS_LOCAL')
                  CALL SCARC_CHECK_DIVISIBILITY(OGF%NCG, 'OGF%NCG')

                  IF (.NOT.TWO_D) THEN
                     OLC%N_WALL_CELLS_LOCAL = OLF%N_WALL_CELLS_LOCAL/4
                     OGC%NCG = OGF%NCG/4
                  ELSE
                     OLC%N_WALL_CELLS_LOCAL = OLF%N_WALL_CELLS_LOCAL/2
                     OGC%NCG = OGF%NCG/2
                  ENDIF

                  CALL SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)

                  CALL SCARC_ALLOCATE_INT1(OGC%ICG_TO_IWG, 1, OGC%NCG, NSCARC_INIT_ZERO, 'ICG_TO_IWG')
                  CALL SCARC_ALLOCATE_INT1(OGC%ICG_TO_ICE, 1, OGC%NCG, NSCARC_INIT_ZERO, 'ICG_TO_ICE')
                  CALL SCARC_ALLOCATE_INT1(OGC%ICG_TO_ICW, 1, OGC%NCG, NSCARC_INIT_ZERO, 'ICG_TO_ICW')

               ENDDO
            ENDIF

            !> setup complete face information for coarser mesh
            CALL SCARC_SETUP_WALL_LEVEL(IOR0, IWC, IREFINE, NM, NL)

         ENDDO

         CALL SCARC_POINT_TO_GRID(NM, NL)
         CALL SCARC_SETUP_CELL_INDEX (NM, NL)
         CALL SCARC_SETUP_WALL_COORDS(NM, NL)
         CALL SCARC_SETUP_WALL_INDEX (NM, NL)
         CALL SCARC_SETUP_WALL_GHOSTS(NM, NL)

      ENDDO LEVEL_GMG_LEVEL_LOOP
   ENDDO MESHES_LOOP3
ENDIF MULTI_LEVEL_IF

!> -------------------------------------------------------------------------------------------
!> Debug FACE, WALL and DISCRET structures - only if directive SCARC_DEBUG is set
!> -------------------------------------------------------------------------------------------
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_STACK, NLEVEL_MIN, 'STACK')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACE , NLEVEL_MIN, 'FACE')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALL , NLEVEL_MIN, 'WALL')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_GRID , NLEVEL_MIN, 'DISCRET')
#endif

END SUBROUTINE SCARC_SETUP_WALLS


!> ------------------------------------------------------------------------------------------------
!> Initialize arrays for data exchange
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGES
USE SCARC_POINTERS, ONLY:  OS, OG
INTEGER :: NL, NM, NOM, NLEN
INTEGER :: INBR

!>
!> Allocate request array for data exchanges
!> Exchange basic information about wall sizes (needed for the dimensioning of the exchange buffers)
!>
IF (N_MPI_PROCESSES>1) THEN
   ALLOCATE (REQ(N_EXCHANGES*40), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_EXCHANGES', 'REQ', IERROR)
   REQ = MPI_REQUEST_NULL
ENDIF
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_SETUP, NLEVEL_MIN)

!>
!> Allocate send and receive buffers (real and integer) in correct lengths
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      CALL SCARC_POINT_TO_LEVEL(NM, NLEVEL_MIN)                          

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NLEVEL_MIN)

      !> allocate send and receive buffers for neighbor NOM in correct lengths
      NLEN = OG%NCG
      CALL SCARC_ALLOCATE_INT1 (OS%SEND_INT1 , 1, NLEN, NSCARC_INIT_ZERO, 'SEND_INT1')
      CALL SCARC_ALLOCATE_INT1 (OS%RECV_INT1 , 1, NLEN, NSCARC_INIT_ZERO, 'RECV_INT1')
      CALL SCARC_ALLOCATE_REAL1(OS%SEND_REAL1, 1, NLEN, NSCARC_INIT_ZERO, 'SEND_REAL1')
      CALL SCARC_ALLOCATE_REAL1(OS%RECV_REAL1, 1, NLEN, NSCARC_INIT_ZERO, 'RECV_REAL1')

      NLEN = 2 * OG%NCG
      calL SCARC_ALLOCATE_REAL1(OS%SEND_REAL2, 1, NLEN, NSCARC_INIT_ZERO, 'SEND_REAL2')
      CALL SCARC_ALLOCATE_REAL1(OS%RECV_REAL2, 1, NLEN, NSCARC_INIT_ZERO, 'RECV_REAL2')

      NLEN = NSCARC_MAX_STENCIL * OG%NCG
      calL SCARC_ALLOCATE_REAL1(OS%SEND_REAL7, 1, NLEN, NSCARC_INIT_ZERO, 'SEND_REAL7')
      CALL SCARC_ALLOCATE_REAL1(OS%RECV_REAL7, 1, NLEN, NSCARC_INIT_ZERO, 'RECV_REAL7')

      !> neighboring wall structures for common wall cells
      ALLOCATE (OG%WALL(OG%NCG), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OG%WALL', IERROR)

      !> In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN) THEN
         DO NL=NLEVEL_MIN+1,NLEVEL_MAX
            CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)
            ALLOCATE (OG%WALL(OG%NCG), STAT=IERROR)
            CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OG%WALL', IERROR)
         ENDDO
      ENDIF

   ENDDO
ENDDO

!>
!> If there is more than 1 MPI-process, initialize communication structures on finest level 
!>
!IF (N_MPI_PROCESSES > 1) THEN
IF (NMESHES > 1) THEN
   DO NL = NLEVEL_MIN, NLEVEL_MAX
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_NUMBER, NL)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_WIDTH, NL)
   ENDDO
ENDIF

!>
!> Correct boundary types for cells adjacent to obstructions on ghost cells
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (IS_UNSTRUCTURED) THEN
      CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NLEVEL_MIN)
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX
         CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
      ENDDO
   ENDIF
ENDDO

END SUBROUTINE SCARC_SETUP_EXCHANGES


!> -----------------------------------------------------------------------------------------
!> --- Store neighbors of mesh
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_STORE_NEIGHBOR(NM, NOM)
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: INBR
DO INBR = 1, SCARC(NM)%N_NEIGHBORS
   IF (SCARC(NM)%NEIGHBORS(INBR) == NSCARC_UNDEF_INT) EXIT      !> not found, to be stored
   IF (SCARC(NM)%NEIGHBORS(INBR) == NOM) RETURN                 !> nothing to do, already stored
ENDDO
SCARC(NM)%NEIGHBORS(INBR) = NOM
RETURN
END SUBROUTINE SCARC_STORE_NEIGHBOR


!> -----------------------------------------------------------------------------------------
!> --- Setup CELL_INDEX array on coarser grid levels in case of MG-method
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CELL_INDEX(NM, NL)
USE SCARC_POINTERS, ONLY: M, L, OB
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, NOBST

CALL SCARC_POINT_TO_GRID(NM, NL)

!> if finest level, the corresponding CELL_INDEX array is already available by surrounding routines
IF (NL == NLEVEL_MIN) THEN
   L%CELL_INDEX_PTR => M%CELL_INDEX

   !> on coarser levels, it must still be computed
ELSE

   CALL SCARC_ALLOCATE_INT3(L%CELL_INDEX, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'CELL_INDEX')
   L%CELL_INDEX_PTR => L%CELL_INDEX
   L%N_CELL_INDEX = 0

   !>
   !> Preset it for all grid cells
   !>
   DO K=0,L%NZ+1
      DO J=0,L%NY+1
         DO I=0,1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
         DO I=L%NX,L%NX+1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=0,L%NZ+1
      DO I=0,L%NX+1
         DO J=0,1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
         DO J=L%NY,L%NY+1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO J=0,L%NY+1
      DO I=0,L%NX+1
         DO K=0,1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
         DO K=L%NZ,L%NZ+1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   !>
   !> Consider cells in obstructions
   !>
   DO NOBST=1,L%N_OBST
      OB => SCARC(NM)%LEVEL(NL)%OBST(NOBST)
      DO K=OB%K1,OB%K2+1
         DO J=OB%J1,OB%J2+1
            DO I=OB%I1,OB%I2+1
               IF (L%CELL_INDEX(I,J,K)==0) THEN
                  L%N_CELL_INDEX = L%N_CELL_INDEX + 1
                  L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

ENDIF

END SUBROUTINE SCARC_SETUP_CELL_INDEX


!> -----------------------------------------------------------------------------------------
!> In case of an MG-method:
!> Setup WALL_INDEX array on coarser grid levels
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_INDEX(NM, NL)
USE SCARC_POINTERS, ONLY: M, L, G
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, ICG, IW, IOR0


CALL SCARC_POINT_TO_GRID(NM, NL)

!> if on finest level, the array WALL_INDEX is already available by surrounding routines
IF (NL == NLEVEL_MIN) THEN

   L%WALL_INDEX_PTR => M%WALL_INDEX

   !> if on coarser levels, it must still be computed
ELSE

   CALL SCARC_ALLOCATE_INT2(L%WALL_INDEX, 1, L%N_CELL_INDEX, -3, 3, NSCARC_INIT_ZERO, 'WALL_INDEX')
   L%WALL_INDEX_PTR => L%WALL_INDEX

   DO IW = 1, L%N_WALL_CELLS_EXT

      I = G%WALL(IW)%IXW
      J = G%WALL(IW)%IYW
      K = G%WALL(IW)%IZW

      IOR0 = G%WALL(IW)%IOR
      ICG  = L%CELL_INDEX(I,J,K)

      L%WALL_INDEX(ICG,-IOR0) = IW

   ENDDO

ENDIF

END SUBROUTINE SCARC_SETUP_WALL_INDEX

!> -----------------------------------------------------------------------------------------
!> In case of an MG-method:
!> Setup ghost cells order 
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_GHOSTS(NM, NL)
USE SCARC_POINTERS, ONLY: OL, OG, GWC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IW, IOR0, NOM

CALL SCARC_POINT_TO_GRID(NM, NL)
DO IW = 1, L%N_WALL_CELLS

  GWC => G%WALL(IW)
  NOM  = GWC%NOM
  IF (NOM == 0) CYCLE

  IOR0 = GWC%IOR
  CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)

  OG%ICG2 = OG%ICG2 + 1
  IF (OL%NCG0(IOR0) == 0) OL%NCG0(IOR0) = OG%ICG2
  OL%NCG(IOR0) = OG%ICG2

ENDDO

END SUBROUTINE SCARC_SETUP_WALL_GHOSTS

!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor in case of MG-method
!> Number of obstructions on coarse level is the same as on fine level
!> TODO: Only works for special cases which run for GMG, must still be extended!!
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_COORDS(NM, NL)
USE SCARC_POINTERS, ONLY: L, G, OB
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, IO, IWC
INTEGER :: I, J, K

CALL SCARC_POINT_TO_GRID(NM, NL)
IWC = L%N_WALL_CELLS_EXT + 1

DO IO = 1, L%N_OBST

   OB => L%OBST(IO)

   !> Analyze IOR = 1
   I = OB%I1
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = G%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I+1; GWC%IYW = J; GWC%IZW = K
         GWC%IXG = I  ; GWC%IYG = J; GWC%IZG = K
         GWC%IOR = 1
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = G%CELL_NUMBER(I+1, J, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I  ; GWC%IYW = J; GWC%IZW = K
         GWC%IXG = I+1; GWC%IYG = J; GWC%IZG = K
         GWC%IOR =-1
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 2
   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = G%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I; GWC%IYW = J+1; GWC%IZW = K
         GWC%IXG = I; GWC%IYG = J  ; GWC%IZG = K
         GWC%IOR = 2
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -2
   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = G%CELL_NUMBER(I, J+1, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I; GWC%IYW = J  ; GWC%IZW = K
         GWC%IXG = I; GWC%IYG = J+1; GWC%IZG = K
         GWC%IOR =-2
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 3
   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = G%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I; GWC%IYW = J; GWC%IZW = K+1
         GWC%IXG = I; GWC%IYG = J; GWC%IZG = K
         GWC%IOR = 3
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -3
   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = G%CELL_NUMBER(I, J, K+1)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I; GWC%IYW = J; GWC%IZW = K
         GWC%IXG = I; GWC%IYG = J; GWC%IZG = K+1
         GWC%IOR =-3
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

ENDDO

END SUBROUTINE SCARC_SETUP_WALL_COORDS


!> -------------------------------------------------------------------------------------------------
!> Correct BTYPE related to internal obstructions on ghost cells
!> TODO: to check again !!!
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
USE SCARC_POINTERS, ONLY: L, G, GWC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWG
INTEGER :: IX, IY, IZ, IOR0, BTYPE0

CALL SCARC_POINT_TO_GRID(NM, NL)

DO IWG = 1, L%N_WALL_CELLS_EXT

   GWC => G%WALL(IWG)

   IF (GWC%NOM == 0) CYCLE                    !> TODO: equal or not equal ??

   IX = GWC%IXW
   IY = GWC%IYW
   IZ = GWC%IZW

   IOR0   = GWC%IOR
   BTYPE0 = GWC%BTYPE

   ! ICG = L%CELL_INDEX_PTR(IX, IY, IZ)
   ! IWG = L%WALL_INDEX_PTR(ICG, IOR0)
   ! IF (GWC%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY) GWC%BTYPE=NEUMANN
   ! IF (L%IS_SOLID(IX, IY, IZ)) GWC%BTYPE=NEUMANN

   IF (GWC%BOUNDARY_TYPE == SOLID_BOUNDARY) GWC%BTYPE=NEUMANN

ENDDO

END SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
USE SCARC_POINTERS, ONLY: LF, LC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IXC, IYC, IZC
INTEGER :: IXF, IYF, IZF
INTEGER :: NLF, NLC
INTEGER :: IWC, ICF(4)=0, IWF(4)=0, IOR0

NLF = NL-1
NLC = NL

LF => SCARC(NM)%LEVEL(NLF)
LC => SCARC(NM)%LEVEL(NLC)

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
      ICF(1) = LF%CELL_INDEX_PTR(1  , IYF  , IZF  )
      ICF(2) = LF%CELL_INDEX_PTR(1  , IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

   !> IOR = -1
   IOR0 = -1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      ICF(1) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF  )
      ICF(2) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

   !> IOR = 2
   IOR0 = 2
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF  , IYF, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1, IYF, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = -2
   IOR0 = -2
   IXF  = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF    , IYF, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1  , IYF, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = 3
   IOR0 = 3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF  , 1)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF  , 1)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

   !> IOR = -3
   IOR0 = -3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF  , LF%NZ)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF  , LF%NZ)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

ELSE

   !> IOR = 1
   IOR0 = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IYC = 1, LC%NY
         IYF = 2*IYC - 1
         ICF(1) = LF%CELL_INDEX_PTR(1  , IYF  , IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(1  , IYF+1, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(1  , IYF  , IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(1  , IYF+1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = -1
   IOR0 = -1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IYC = 1, LC%NY
         IYF = 2*IYC - 1
         ICF(1) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(LF%NX, IYF+1, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(LF%NX, IYF+1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = 2
   IOR0 = 2
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF  , 1, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1, 1, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF  , 1, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1, 1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = -2
   IOR0 = -2
   IXF  = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF    , LF%NY, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , LF%NY, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF    , LF%NY, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1  , LF%NY, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = 3
   IOR0 = 3
   DO IYC = 1, LC%NY
      IYF = 2*IYC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF  , 1)
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF  , 1)
         ICF(3) = LF%CELL_INDEX_PTR(IXF    , IYF+1, 1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1  , IYF+1, 1)
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

   !> IOR = -3
   IOR0 = -3
   DO IYC = 1, LC%NY
      IYF = 2*IYC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF  , LF%NZ)
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF  , LF%NZ)
         ICF(3) = LF%CELL_INDEX_PTR(IXF  , IYF+1, LF%NZ)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1, IYF+1, LF%NZ)
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

ENDIF

SCARC_COUNT_EXTERNAL_WALL_CELLS = IWC
END FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)
USE SCARC_POINTERS, ONLY: LF, LC, GF, GC, OB
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWC, NW_INT
INTEGER :: IC, IO
INTEGER :: I, J, K

LC => SCARC(NM)%LEVEL(NL)
LF => SCARC(NM)%LEVEL(NL-1)

SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      GC => LC%STRUCTURED
      GF => LF%STRUCTURED
   CASE (NSCARC_GRID_UNSTRUCTURED)
      GC => LC%UNSTRUCTURED
      GF => LF%UNSTRUCTURED
END SELECT

LC%N_OBST = LF%N_OBST                   !> Number of obstructions is the same on all levels

ALLOCATE(LC%OBST(LC%N_OBST), STAT=IERROR)
CALL ChkMemErr('SCARC_COUNT_INTERNAL_WALL_CELLS','OBST',IERROR)

NW_INT = 0
IWC = LC%N_WALL_CELLS_EXT + 1

DO IO = 1, LF%N_OBST

   OB => SCARC(NM)%LEVEL(NL)%OBST(IO)

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
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = GC%CELL_NUMBER(I+1, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 2
   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -2
   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = GC%CELL_NUMBER(I, J+1, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 3
   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -3
   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = GC%CELL_NUMBER(I, J, K+1)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

ENDDO

SCARC_COUNT_INTERNAL_WALL_CELLS = NW_INT
END FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS


!> -------------------------------------------------------------------------------------------------
!> Count external wall cells on face IOR
!> -------------------------------------------------------------------------------------------------
LOGICAL FUNCTION IS_EXTERNAL_WALLCELL(IOR0, ICF, NCNT, NM, NL)
USE SCARC_POINTERS, ONLY: L, G
INTEGER, INTENT(IN) :: IOR0, NCNT, NM, NL
INTEGER, DIMENSION(:), INTENT(IN) :: ICF
INTEGER :: I, IWF_LAST, IWF(4)=0
REAL(EB) :: BSUM

IS_EXTERNAL_WALLCELL = .FALSE.

CALL SCARC_POINT_TO_GRID(NM, NL)

DO I = 1, NCNT
   IWF(I) = L%WALL_INDEX_PTR(ICF(I), -IOR0)
ENDDO

BSUM = 0.0_EB
IWF_LAST = 0

DO I = 1, NCNT
   IF (IWF(I)>0) THEN
      BSUM = BSUM + REAL(G%WALL(IWF(I))%BTYPE,EB)
      IWF_LAST = IWF(I)
   ENDIF
ENDDO

IF (IWF_LAST == 0) RETURN
IF (ABS(BSUM/REAL(NCNT,EB) - REAL(G%WALL(IWF_LAST)%BTYPE,EB)) < 1E-12) THEN
   IS_EXTERNAL_WALLCELL = .TRUE.
   RETURN
ELSE
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_BOUNDARY_SUM, SCARC_NONE, IOR0)
ENDIF

END FUNCTION IS_EXTERNAL_WALLCELL


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_NEIGHBOR(NX1, NX2, NY1, NY2, NZ1, NZ2, IWG, NM, NOM, NL)
USE SCARC_POINTERS, ONLY: G, OG
INTEGER, INTENT(IN) :: NX1, NX2, NY1, NY2, NZ1, NZ2
INTEGER, INTENT(IN) :: IWG, NM, NOM, NL
INTEGER :: NOMX, NOMY, NOMZ
INTEGER :: ICG, ICE, IX, IY, IZ, JL, IXW, IYW, IZW, IXG, IYG, IZG

CALL SCARC_POINT_TO_GRID(NM, NL)
CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)

ICE  = G%ICE
ICG  = OG%ICG

!> set neighboring coordinates
GWC => G%WALL(IWG)

GWC%IXN(1) = NX1
GWC%IXN(2) = NX2
GWC%IYN(1) = NY1
GWC%IYN(2) = NY2
GWC%IZN(1) = NZ1
GWC%IZN(2) = NZ2

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
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ICG  = ICG  + 1
         ICE  = ICE  + 1

         GWC%ICE = ICE                                         !> number of extended grid cell
         GWC%ICG = ICG                                         !> number of ghost grid cell

         G%ICE_TO_IWG(ICE) = IWG                               !> map extended cell to global wall cell
         G%ICE_TO_ICG(ICE) = ICG                               !> map extended cell to ghost cell
         
         IXG = G%WALL(IWG)%IXG
         IYG = G%WALL(IWG)%IYG
         IZG = G%WALL(IWG)%IZG

         IXW = G%WALL(IWG)%IXW
         IYW = G%WALL(IWG)%IYW
         IZW = G%WALL(IWG)%IZW

         G%CELL_NUMBER(IXG, IYG, IZG) = ICE

         OG%ICG_TO_IWG(ICG) = IWG                              !> map ghost cell to global wall cell
         OG%ICG_TO_ICW(ICG) = G%CELL_NUMBER(IXW, IYW, IZW)     !> get cell number of adjacent internal cell

      ENDDO
   ENDDO
ENDDO

G%ICE  = ICE                                                   !> store extended cell counter
OG%ICG = ICG                                                   !> store ghost cell counter
OG%ICG_TO_IWG(ICG) = IWG                                       !> map local wall cell to global wall cell
OG%ICG_TO_ICE(ICG) = ICE                                       !> map local wall cell to global wall cell

END SUBROUTINE SCARC_SETUP_WALL_NEIGHBOR


!> ----------------------------------------------------------------------------------------------------
!> Check divisibility by 2 of a given number of elements (in one grid direction)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_DIVISIBILITY(NN, CDIR)
INTEGER, INTENT(IN) :: NN
CHARACTER(*) , INTENT(IN) :: CDIR
IF (MOD(NN,2) /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBER, CDIR, NSCARC_NONE)
END SUBROUTINE SCARC_CHECK_DIVISIBILITY


!> ----------------------------------------------------------------------------------------------------
!> Set wall cell information on coarse level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_LEVEL(IOR0, IWC, IREFINE, NM, NL)
USE SCARC_POINTERS, ONLY: LF, LC, FF, WF, WC, OLC, OGC
INTEGER, INTENT(INOUT) :: IWC
INTEGER, INTENT(IN) :: NM, NL
INTEGER, INTENT(IN) :: IOR0, IREFINE
INTEGER :: IWF(4) , IBCF(4), NOMF(4)
INTEGER :: IX,  IY,  IZ, I
INTEGER :: NX1, NY1, NZ1
INTEGER :: NX2, NY2, NZ2
INTEGER :: IX1, IY1, IZ1
INTEGER :: IX2, IY2, IZ2
INTEGER :: IDIFF, JDIFF, KDIFF


!> set coordinate dimensions for correspoding face
SELECT CASE (ABS(IOR0))
   CASE (1)
      IF (IOR0 > 0) THEN                                           !> set dimensions for wall cell counting
         NX1 = 0; NX2 = 0
      ELSE
         NX1 = LC%NX+1; NX2 = LC%NX+1
      ENDIF
      NY1 = 1;  NY2 = LC%NY
      NZ1 = 1;  NZ2 = LC%NZ
   CASE (2)
      NX1 = 1; NX2 = LC%NX
      IF (IOR0 > 0) THEN
         NY1 = 0; NY2 = 0
      ELSE
         NY1 = LC%NY+1; NY2 = LC%NY+1
      ENDIF
      NZ1 = 1; NZ2 = LC%NZ
   CASE (3)
      NX1 = 1; NX2 = LC%NX
      NY1 = 1; NY2 = LC%NY
      IF (IOR0 > 0) THEN
         NZ1 = 0; NZ2 = 0
      ELSE
         NZ1 =LC%NZ+1; NZ2 =LC%NZ+1
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
                  IWF(1) = FF%NCW0 + 2*(IZ-1)
               CASE ( 2)
                  IWF(1) = FF%NCW0 + 2*(IZ-1)*LF%NX + 2*(IX - 1)
               CASE ( 3)
                  IWF(1) = FF%NCW0 + 2*(IX-1)
            END SELECT
            IWF(2) = IWF(1)+1

            !> set fine cell neighbors (they must be the same for all fine IW's)
            NOMF(1) = WF(IWF(1))%NOM
            NOMF(2) = WF(IWF(2))%NOM
            IF (NOMF(1) /= NOMF(2)) CALL SCARC_SHUTDOWN(NSCARC_ERROR_NEIGHBOR_TYPE, SCARC_NONE, NOMF(1))

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
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_BOUNDARY_TYPE, SCARC_NONE, IBCF(1))
            ENDIF

            !> in case of an internal boundary set neighboring WALL cells
            IF (NOMF(1) > 0) THEN

               OLC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)
               SELECT CASE(TYPE_GRID)
                  CASE (NSCARC_GRID_STRUCTURED)
                     OGC => OLC%STRUCTURED
                  CASE (NSCARC_GRID_UNSTRUCTURED)
                     OGC => OLC%UNSTRUCTURED
               END SELECT

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
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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
               CALL SCARC_SETUP_WALL_NEIGHBOR(IX1, IX2, 1, 1, IZ1, IZ2, IWC, NM, NOMF(1), NL)

            ENDIF

         !> ------------------------------------------------------------
         !> 3D-version
         !> ------------------------------------------------------------
         ELSE

            !> determine fine IW's, which must be merged to one coarse IW
            SELECT CASE (ABS(IOR0))
               CASE (1)
                  IWF(1) = FF%NCW0 + (2*IZ-2)*LF%NY + 2*IY - 2
                  IWF(3) = FF%NCW0 + (2*IZ-1)*LF%NY + 2*IY - 2
               CASE (2)
                  IWF(1) = FF%NCW0 + (2*IZ-2)*LF%NX + 2*IX - 2
                  IWF(3) = FF%NCW0 + (2*IZ-1)*LF%NX + 2*IX - 2
               CASE (3)
                  IWF(1) = FF%NCW0 + (2*IY-2)*LF%NX + 2*IX - 2
                  IWF(3) = FF%NCW0 + (2*IY-1)*LF%NX + 2*IX - 2
            END SELECT
            IWF(2) = IWF(1)+1
            IWF(4) = IWF(3)+1

            !> set fine cell neighbors (they must be the same for all fine IW's)
            DO I=1,4
               NOMF(I) = WF(IWF(I))%NOM
            ENDDO

            IF (NOMF(1) /= NOMF(2) .OR. NOMF(1) /= NOMF(3) .OR. NOMF(1) /= NOMF(4)) &
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_NEIGHBOR_TYPE, SCARC_NONE, IOR0)
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
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_BOUNDARY_TYPE, SCARC_NONE, -999)
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
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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

               CALL SCARC_SETUP_WALL_NEIGHBOR(IX1, IX2, IY1, IY2, IZ1, IZ2, IWC, NM, NOMF(1), NL)

            ENDIF
         ENDIF
         IWC = IWC + 1
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_WALL_LEVEL


!> -------------------------------------------------------------------------------------------------
!> Set analogues to I_MIN, IMAX, J_MIN, J_MAX, K_MIN, K_MAX on coarse level (only GMG!)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)
USE SCARC_POINTERS, ONLY: G, OG
INTEGER, INTENT(IN) :: IREFINE, NM, NOM, NL
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IW
LOGICAL :: FOUND

CALL SCARC_POINT_TO_GRID(NM, NL)
CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)

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

FOUND = .FALSE.
SEARCH_LOOP: DO IW=1, OG%NCG

   !> neighborship structure already known from finest level
   IF (G%WALL(IW)%NOM/=NOM) CYCLE SEARCH_LOOP
   FOUND = .TRUE.

   SELECT CASE (G%WALL(IW)%IOR)
      CASE ( 1)
         IMIN=MAX(IMIN,G%WALL(NM)%IXN(1)-1)
      CASE (-1)
         IMAX=MIN(IMAX,G%WALL(NM)%IXN(2)+1)
      CASE ( 2)
         JMIN=MAX(JMIN,G%WALL(NM)%IYN(1)-1)
      CASE (-2)
         JMAX=MIN(JMAX,G%WALL(NM)%IYN(2)+1)
      CASE ( 3)
         KMIN=MAX(KMIN,G%WALL(NM)%IZN(1)-1)
      CASE (-3)
         KMAX=MIN(KMAX,G%WALL(NM)%IZN(2)+1)
   END SELECT
ENDDO SEARCH_LOOP

N_EXCHANGES = N_EXCHANGES+1

END SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS


!> ----------------------------------------------------------------------------------------------------
!> Allocate several global structures for data exchange
!> --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NP

IF (N_MPI_PROCESSES > 1) THEN

   !> Allocate and preset counter and displacement vector for global data exchanges
   CALL SCARC_ALLOCATE_INT1 (COUNTS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'COUNTS')
   CALL SCARC_ALLOCATE_INT1 (DISPLS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'DISPLS')

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

ENDIF

CALL SCARC_ALLOCATE_INT1 (MESH_INT , 1, NMESHES, NSCARC_INIT_ZERO, 'MESH_INT')
CALL SCARC_ALLOCATE_REAL1(MESH_REAL, 1, NMESHES, NSCARC_INIT_ZERO, 'MESH_REAL')

CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

END SUBROUTINE SCARC_SETUP_GLOBALS


!> -----------------------------------------------------------------------------
!> Get information about global numbers of unknowns for unstructured case
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DIMENSIONS(NL)
USE SCARC_POINTERS, ONLY: G
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2

!> Preset communication array MESH_INT with local numbers of cells for all meshes depending on type of discretization
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      MESH_INT(NM2) = G%NC_LOCAL(NM2)
   ENDDO

ENDDO MESHES_LOOP1

!> Broadcast number of local mesh cells on level NL to all and build global sum
IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,MESH_INT,COUNTS,DISPLS, MPI_INTEGER,MPI_COMM_WORLD,IERROR)
NC_GLOBAL(NL) = SUM(MESH_INT(1:NMESHES))

!> Store information on local and global cells numbers on data structure of corresponding discretization type
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   G%NC_LOCAL(1:NMESHES) = MESH_INT(1:NMESHES)
   G%NC_GLOBAL = SUM(MESH_INT(1:NMESHES))

   !> compute offset between local grid numberings
   IF (NMESHES > 1) THEN
      DO NM2=2,NMESHES
         G%NC_OFFSET(NM2) = G%NC_OFFSET(NM2-1) + G%NC_LOCAL(NM2-1)
      ENDDO
   ENDIF

ENDDO MESHES_LOOP2

IF (NL == NLEVEL_MIN) THEN
   DO NM = 1, NMESHES
      SCARC(NM)%NC = MESH_INT(NM)
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_DIMENSIONS


!> ----------------------------------------------------------------------------------------------------
!> Setup system of equation:
!> Define matrix stencils and initialize matrices and boundary conditions on all needed levels
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEMS
INTEGER :: NM, NL

!> 
!> ------ Setup sizes for system matrices
!> 
SELECT_SCARC_METHOD_SIZES: SELECT CASE (TYPE_METHOD)

   !>
   !> -------- Global Krylov method
   !>
   CASE (NSCARC_METHOD_KRYLOV)
   
      CALL SCARC_ASSIGN_GRID_TYPE (TYPE_GRID)                                !> process specified discretization type
      CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)         !> use specified discretization
   
      IF (HAS_TWO_LEVELS) &                                                  !> two grid levels used?
         CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MAX)      !> setup size for coarse grid
   
      IF (IS_CG_GMG) THEN                                                    !> GMG-preconditioning?
         DO NL=NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX , NL)          !> setup size for all levels
         ENDDO
      ENDIF
   
   !>
   !> -------- Global Multigrid method
   !>
   CASE (NSCARC_METHOD_MULTIGRID)
   
      CALL SCARC_ASSIGN_GRID_TYPE (TYPE_GRID)                                !> process specified discretization type
      SELECT CASE (TYPE_MULTIGRID)
         CASE (NSCARC_MULTIGRID_GEOMETRIC)                                   !> Geometric MG?
            DO NL=NLEVEL_MIN, NLEVEL_MAX
               CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NL)        !> setup size for all levels
            ENDDO
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)                                   !> Algebraic MG?
            CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)   !> setup sizes for AMG
      END SELECT
   
#ifdef WITH_SCARC_MGM
   !>
   !> -------- Global MGM method - currently just proof of concept
   !>
   CASE (NSCARC_METHOD_MGM)
   
      !>
      !> first process structured discretization
      !>
      CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_STRUCTURED)              
      CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)        
   
      !>
      !> then process unstructured discretization
      !>
      CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)            
      CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)         
   
#endif

END SELECT SELECT_SCARC_METHOD_SIZES

!> 
!> ------ Assemble system matrices on requested grid levels and set boundary conditions
!> 
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   SELECT_SCARC_METHOD: SELECT CASE (TYPE_METHOD)

      !>
      !> ---------- Krylov method (CG/BICG) as main solver, different preconditioners possible
      !> 
      CASE (NSCARC_METHOD_KRYLOV)

         SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

            !> 
            !> in case of multigrid as preconditioner:
            !> 
            CASE (NSCARC_RELAX_GMG)

               SELECT_KRYLOV_PRECON_MG: SELECT CASE (TYPE_MULTIGRID)

                  !>
                  !> Geometric multigrid:
                  !>    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
                  !>
                  CASE (NSCARC_MULTIGRID_GEOMETRIC)

                     DO NL = NLEVEL_MIN, NLEVEL_MAX
                        CALL SCARC_SETUP_MATRIX  (NM, NL)
                        CALL SCARC_SETUP_BOUNDARY(NM, NL)
#ifdef WITH_MKL
                        IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL .OR. TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) THEN
                           IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
                              CALL SCARC_SETUP_MATRIX_MKL_SINGLE(NM, NL)
                           ELSE
                              CALL SCARC_SETUP_MATRIX_MKL_DOUBLE(NM, NL)
                           ENDIF
                        ENDIF
#endif
                     ENDDO

                  !>
                  !> Algebraic multigrid:
                  !>    -  use compact storage technique on all levels (no other choise possible!)
                  !>    -  assemble standard n-point-matrix only on finest level
                  !>    -  construct all coarser levels by requested coarsening strategy
                  !>
                  CASE (NSCARC_MULTIGRID_ALGEBRAIC)
   
                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               END SELECT SELECT_KRYLOV_PRECON_MG

#ifdef WITH_MKL
            !> 
            !> in case of LU-decomposition as preconditioner
            !> 
            CASE (NSCARC_RELAX_MKL)

               SELECT_KRYLOV_PRECON_MKL: SELECT CASE(TYPE_SCOPE(1))

                  !> Locally acting: PARDISO from MKL as preconditioners with possible coarse grid correction
                  CASE (NSCARC_SCOPE_LOCAL)
   
                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
                     IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
                        CALL SCARC_SETUP_MATRIX_MKL_SINGLE(NM, NLEVEL_MIN)
                     ELSE
                        CALL SCARC_SETUP_MATRIX_MKL_DOUBLE(NM, NLEVEL_MIN)
                     ENDIF
   
                     IF (HAS_TWO_LEVELS) THEN
                        CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MAX)
                        CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
                        IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
                           IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
                              CALL SCARC_SETUP_MATRIX_MKL_SINGLE(NM, NLEVEL_MAX)
                           ELSE
                              CALL SCARC_SETUP_MATRIX_MKL_DOUBLE(NM, NLEVEL_MAX)
                           ENDIF
                        ENDIF
                     ENDIF

                  !> Globally acting: Cluster_Sparse_Solver from MKL as preconditioner
                  CASE (NSCARC_SCOPE_GLOBAL)
   
                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
                     IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
                        CALL SCARC_SETUP_MATRIX_MKL_SINGLE(NM, NLEVEL_MIN)
                     ELSE
                        CALL SCARC_SETUP_MATRIX_MKL_DOUBLE(NM, NLEVEL_MIN)
                     ENDIF
   
               END SELECT SELECT_KRYLOV_PRECON_MKL
#endif

            !> 
            !> in case of one-level preconditioners (JACOBI/SSOR/FFT):
            !> assemble standard n-point-matrix on finest level with possible coarse grid correction
            !> 
            CASE DEFAULT
   
               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               IF (HAS_TWO_LEVELS) THEN
                  CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
#ifdef WITH_MKL
                  IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
                     IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
                        CALL SCARC_SETUP_MATRIX_MKL_SINGLE(NM, NLEVEL_MAX)
                     ELSE
                        CALL SCARC_SETUP_MATRIX_MKL_DOUBLE(NM, NLEVEL_MAX)
                     ENDIF
                  ENDIF
#endif
               ENDIF

         END SELECT SELECT_KRYLOV_PRECON

      !> 
      !> ---------- Multigrid as main solver
      !> 
      CASE (NSCARC_METHOD_MULTIGRID)
   
         SELECT_MULTIGRID_TYPE: SELECT CASE (TYPE_MULTIGRID)
   
            !> 
            !>  Geometric multigrid:
            !>    -  assemble standard n-point-matrix hierarchy on all levels
            !>    -  use MKL coarse grid solver if requested
            !> 
            CASE (NSCARC_MULTIGRID_GEOMETRIC)

               DO NL = NLEVEL_MIN, NLEVEL_MAX
                  CALL SCARC_SETUP_MATRIX  (NM, NL)
                  CALL SCARC_SETUP_BOUNDARY(NM, NL)
               ENDDO

#ifdef WITH_MKL
               IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL .OR. TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) &
                        CALL SCARC_SETUP_MATRIX_MKL_SINGLE(NM, NL)
                  ENDDO
                  IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_MATRIX_MKL_SINGLE(NM, NLEVEL_MAX)
               ELSE
                  DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                     IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL .OR. TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) &
                        CALL SCARC_SETUP_MATRIX_MKL_DOUBLE(NM, NL)
                  ENDDO
                  IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_MATRIX_MKL_DOUBLE(NM, NLEVEL_MAX)
               ENDIF
#endif

            !> 
            !>  Algebraic multigrid:
            !>    -  use compact storage technique (no other choice possible!)
            !>    -  assemble standard n-point-matrix only on finest level
            !>    -  construct all coarser levels later by requested coarsening strategy
            !> 
            CASE (NSCARC_MULTIGRID_ALGEBRAIC)

               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         END SELECT SELECT_MULTIGRID_TYPE

#ifdef WITH_MKL
      !> 
      !> ---------- MKL-LU decomposition as main solver
      !> 
      CASE (NSCARC_METHOD_LU)
  
         CALL SCARC_SETUP_MATRIX(NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
         IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
            CALL SCARC_SETUP_MATRIX_MKL_SINGLE(NM, NLEVEL_MIN)
         ELSE
            CALL SCARC_SETUP_MATRIX_MKL_DOUBLE(NM, NLEVEL_MIN)
         ENDIF
#endif

#ifdef WITH_SCARC_MGM
      !> 
      !> ---------- McKenny-Greengard-Mayo method:
      !> Solving for the structured and unstructured Poisson matrix
      !> Assemble both, the structured and unstructured Poisson matrix
      !> temporarily they will be stored separately in matrices AC and ACU due to the different
      !> settings along internal boundary cells,
      !> in the medium term, a toggle mechanism will be implemented which only switches the corresponding
      !> entries while keeping the entries which are the same for both discretization types
      !> 
      CASE (NSCARC_METHOD_MGM)
   
         !>
         !> First assemble unstructured matrix
         !>
         CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
         CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         !>
         !> Then assemble structured matrix
         !>
         CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_STRUCTURED)
         CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
#endif

   END SELECT SELECT_SCARC_METHOD

ENDDO MESHES_LOOP

!>
!> If there is more than 1 MPI-process, exchange matrix values in overlapping parts
!>
IF (NMESHES > 1) THEN
   DO NL = NLEVEL_MIN, NLEVEL_MAX
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SUBDIAG, NL)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_STENCIL, NL)
   ENDDO
ENDIF

!> -------------------------------------------------------------------------------------------
!> Debug matrix and wall structures - only if directive SCARC_DEBUG is set
!> -------------------------------------------------------------------------------------------
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALL, NLEVEL_MIN, 'WALL')
DO NL = NLEVEL_MIN, NLEVEL_MAX
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NL, 'SYSTEM-MATRIX')
ENDDO
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACE, NLEVEL_MIN, 'FACE_AFTER_SYSTEM')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALL, NLEVEL_MIN, 'WALL_AFTER_SYSTEM')
#endif

END SUBROUTINE SCARC_SETUP_SYSTEMS


!> ------------------------------------------------------------------------------------------------
!> Allocate matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
!> Use compact storage technique:
!> Compression technique to store sparse matrices, non-zero entries are stored
!> in a 1D-vector B(.), row after row,
!> Each row starts with its diagonal entry followed by the other non-zero entries
!> In order to identify each element, pointer arrays ROW and COL are needed,
!> ROW points to the several diagonal entries in vector B(.),
!> COL points to the columns which non-zero entries in the matrix stencil
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX (NM, NL)
USE SCARC_POINTERS, ONLY: L, G
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IX, IY, IZ, IC, IP

!>
!> Compute single matrix entries and corresponding row and column pointers
!> Along internal boundaries use placeholders for the neighboring matrix entries
!> which will be communicated in a following step
!>
SELECT_STORAGE_TYPE: SELECT CASE (SCARC_MATRIX_LEVEL(NL))

   !>
   !> ---------- COMPACT Storage technique
   !>
   CASE (NSCARC_MATRIX_COMPACT)
   
      CALL SCARC_POINT_TO_GRID(NM, NL)
      CALL SCARC_ALLOCATE_MATRIX_COMPACT (G%AC, NL, NSCARC_PRECISION_DOUBLE, NSCARC_INIT_ZERO, 'AC')

      IP  = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
   
               !> Main diagonal entry
               IC = G%CELL_NUMBER(IX, IY, IZ)
               CALL SCARC_SETUP_MATRIX_MAINDIAG_COMPACT (IC, IX, IY, IZ, IP, NM, NL)
   
               !> Lower subdiagonals
               IF (VALID_SUBDIAG(IX, IY, IZ,  3)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX, IY, IZ-1, IP,  3, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ,  2)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX, IY-1, IZ, IP,  2, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ,  1)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX-1, IY, IZ, IP,  1, NM, NL)
   
               !> Upper subdiagonals
               IF (VALID_SUBDIAG(IX, IY, IZ, -1)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX+1, IY, IZ, IP, -1, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ, -2)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX, IY+1, IZ, IP, -2, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ, -3)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX, IY, IZ+1, IP, -3, NM, NL)
   
            ENDDO
         ENDDO
      ENDDO
   
      G%AC%ROW(G%AC%N_ROW) = IP
      G%AC%N_VAL           = IP-1                         !> set correct number of matrix entries
   
      CALL SCARC_RESIZE_MATRIX_COMPACT(G%AC, G%AC%N_VAL, 'AC')
   
   !>
   !> ---------- BANDED Storage technique
   !>
   CASE (NSCARC_MATRIX_BANDED)
   
      CALL SCARC_POINT_TO_GRID(NM, NL)
      CALL SCARC_ALLOCATE_MATRIX_BANDED(G%AB, NL, NSCARC_INIT_ZERO, 'AB')
   
      IP  = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
   
               IC = G%CELL_NUMBER(IX, IY, IZ)
   
               !> Lower diagonal entries
               IF (VALID_SUBDIAG(IX, IY, IZ, 3)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX  , IY  , IZ-1, 3)
               IF (VALID_SUBDIAG(IX, IY, IZ, 2)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX  , IY-1, IZ  , 2)
               IF (VALID_SUBDIAG(IX, IY, IZ, 1)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX-1, IY  , IZ  , 1)
   
               !> Main diagonal entry
               CALL SCARC_SETUP_MATRIX_MAINDIAG_BANDED (IC, IX, IY, IZ)
   
               !> Upper diagonal entries
               IF (VALID_SUBDIAG(IX, IY, IZ, -1)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX+1, IY  , IZ  , -1)
               IF (VALID_SUBDIAG(IX, IY, IZ, -2)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX  , IY+1, IZ  , -2)
               IF (VALID_SUBDIAG(IX, IY, IZ, -3)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX  , IY  , IZ+1, -3)
   
            ENDDO
         ENDDO
      ENDDO
   
END SELECT SELECT_STORAGE_TYPE

END SUBROUTINE SCARC_SETUP_MATRIX


!> ------------------------------------------------------------------------------------------------
!> Check if a subdiagonal entry must be computed in direction IOR0
!> ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION VALID_SUBDIAG(IX, IY, IZ, IOR0)
USE SCARC_POINTERS, ONLY: L, G
INTEGER, INTENT(IN)  :: IX, IY, IZ, IOR0
INTEGER :: IC_INDEX, IW_INDEX

VALID_SUBDIAG = .FALSE.
IF (TWO_D .AND. ABS(IOR0) == 2) RETURN

SELECT CASE (TYPE_GRID)

   !>
   !> If a structured discretization is used, then subdiagonals are built in every direction
   !>
   CASE (NSCARC_GRID_STRUCTURED)

      VALID_SUBDIAG = .TRUE.
      RETURN

   !>
   !> Else check the type of the neighboring cell in direction IOR0
   !>
   CASE (NSCARC_GRID_UNSTRUCTURED)

      !> get cell index of corresponding cell and check its wall index
      IC_INDEX = L%CELL_INDEX_PTR(IX, IY, IZ)
      IW_INDEX = 0
      IF (IC_INDEX /= 0) IW_INDEX  = L%WALL_INDEX_PTR(IC_INDEX, -IOR0)

      !> If this wall index is zero, build a subdiagonal in this direction
      IF (IW_INDEX == 0) THEN
         VALID_SUBDIAG = .TRUE.
         RETURN
   
      !> Else build the subdiagonal entry only if the wall cell is of interpolated type
      ELSE
         IF (G%WALL(IW_INDEX)%BOUNDARY_TYPE== INTERPOLATED_BOUNDARY) THEN
            VALID_SUBDIAG = .TRUE.
            RETURN
         ENDIF
      ENDIF

END SELECT
RETURN

END FUNCTION VALID_SUBDIAG


!> ------------------------------------------------------------------------------------------------
!> Set main diagonal entry for matrix - full matrix of the global problem
!> In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
!> If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG_COMPACT (IC, IX, IY, IZ, IP, NM, NL)
USE SCARC_POINTERS, ONLY: L, G
INTEGER, INTENT(IN) :: IC, IX, IY, IZ, NM, NL
INTEGER, INTENT(INOUT) :: IP
#ifndef WITH_MKL
INTEGER :: NDUMMY
#endif

G%AC%VAL(IP) = - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D) G%AC%VAL(IP) = G%AC%VAL(IP) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
G%AC%VAL(IP) = G%AC%VAL(IP) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

G%AC%ROW(IC) = IP
G%AC%COL(IP) = IC

#ifdef WITH_MKL
IF (IS_MKL_LEVEL(NL)) THEN
   G%AC%COL_GLOBAL(IP) = G%AC%COL(IP) + G%NC_OFFSET(NM)
ENDIF
#else
NDUMMY = NL                                      !> set these variables for simplicity, so that they can be used 
NDUMMY = NM                                      !> as an input parameter even for the non-MKL case
#endif

G%AC%STENCIL(0) = G%AC%VAL(IP)

IP = IP + 1
END SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG_COMPACT


!> ------------------------------------------------------------------------------------------------
!> Compute subdiagonal contribution in direction IOR0 in case of a COMPACT matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG_COMPACT (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IP, IOR0, NM, NL)
USE SCARC_POINTERS, ONLY: L, F, G, AC
INTEGER, INTENT(IN) :: IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0, NM, NL
INTEGER, INTENT(INOUT) :: IP
INTEGER :: IW
#ifdef WITH_MKL
INTEGER :: IX, IY, IZ, ICE
#else
INTEGER :: NDUMMY
#endif
LOGICAL  :: IS_INTERNAL_CELL

F => L%FACE(IOR0)
AC => G%AC

!> Decide wheter cell is interior or exterior cell
SELECT CASE (IOR0)
   CASE ( 1)
      IS_INTERNAL_CELL = IX1 > 1
   CASE (-1)
      IS_INTERNAL_CELL = IX1 < F%NOP
   CASE ( 2)
      IS_INTERNAL_CELL = IY1 > 1
   CASE (-2)
      IS_INTERNAL_CELL = IY1 < F%NOP
   CASE ( 3)
      IS_INTERNAL_CELL = IZ1 > 1
   CASE (-3)
      IS_INTERNAL_CELL = IZ1 < F%NOP
END SELECT

!>
!> if IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
!>
IF (IS_INTERNAL_CELL) THEN

   IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX2, IY2, IZ2)) THEN

      AC%VAL(IP) = AC%VAL(IP) + F%SCAL_INSIDE
      AC%COL(IP) = G%CELL_NUMBER(IX2, IY2, IZ2)

      AC%STENCIL(-IOR0) = AC%VAL(IP)

#ifdef WITH_MKL
      IF (IS_MKL_LEVEL(NL)) THEN
         AC%COL_GLOBAL(IP)= AC%COL(IP) + G%NC_OFFSET(NM)
      ENDIF
#endif
      IP = IP + 1
   ELSE
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

!>
!> if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
!>
ELSE IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN

   IW = SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)           !> get IW of a possibly suitable neighbor at face IOR0
   IF (IW > 0) then                                    !> if available, build corresponding subdiagonal entry

      AC%VAL(IP) = AC%VAL(IP) + F%SCAL_FACE
      AC%COL(IP) = G%WALL(IW)%ICE                      !> store its extended number in matrix column pointers

      AC%STENCIL(-IOR0) = AC%VAL(IP)

#ifdef WITH_MKL
      IF (IS_MKL_LEVEL(NL)) THEN                       !> if MKL method used, also store its global number
         IX  = G%WALL(IW)%IXG
         IY  = G%WALL(IW)%IYG
         IZ  = G%WALL(IW)%IZG
         ICE = G%CELL_NUMBER(IX, IY, IZ)
         AC%COL_GLOBAL(IP) = G%ICE_TO_ICN(ICE) + G%NC_OFFSET(G%WALL(IW)%NOM)
      ENDIF
#else
      NDUMMY = NL                                      !> set these variables for simplicity, so that they can be used 
      NDUMMY = NM                                      !> as an input parameter even for the non-MKL case
#endif

      IP = IP + 1
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG_COMPACT


!> ------------------------------------------------------------------------------------------------
!> Determine if cell IC has a neighbor and, if yes, return corresponding IW-value
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)
USE SCARC_POINTERS, ONLY: L, G, F, GWC
INTEGER, INTENT(IN) :: IC, IOR0
INTEGER :: IXW, IYW, IZW
INTEGER :: IXG, IYG, IZG
INTEGER :: IW

F => L%FACE(IOR0)
SCARC_ASSIGN_SUBDIAG_TYPE = -1
SEARCH_WALL_CELLS_LOOP: DO IW = F%NCW0, F%NCW0 + F%NCW - 1

   GWC => G%WALL(IW)

   IF (GWC%NOM == 0) CYCLE

   IXW = GWC%IXW
   IYW = GWC%IYW
   IZW = GWC%IZW

   IF (G%CELL_NUMBER(IXW, IYW, IZW) /= IC) CYCLE

   IXG = GWC%IXG
   IYG = GWC%IYG
   IZG = GWC%IZG

   IF (IS_UNSTRUCTURED.AND.L%IS_SOLID(IXG, IYG, IZG)) RETURN

   SCARC_ASSIGN_SUBDIAG_TYPE = IW
   RETURN

ENDDO SEARCH_WALL_CELLS_LOOP

END FUNCTION SCARC_ASSIGN_SUBDIAG_TYPE


!> ------------------------------------------------------------------------------------------------
!> Get type of matrix storage scheme for level NL
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_MATRIX_LEVEL(NL)
INTEGER, INTENT(IN) :: NL

IF (NL == NLEVEL_MAX .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
   SCARC_MATRIX_LEVEL = NSCARC_MATRIX_COMPACT
ELSE
   SCARC_MATRIX_LEVEL = TYPE_MATRIX
ENDIF

END FUNCTION SCARC_MATRIX_LEVEL


!> ------------------------------------------------------------------------------------------------
!> Set main diagonal entry for matrix - full matrix of the global problem
!> In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
!> If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG_BANDED (IC, IX, IY, IZ)
USE SCARC_POINTERS, ONLY: L, G
INTEGER, INTENT(IN)  :: IC, IX, IY, IZ
INTEGER :: ID

ID = G%AB%POS(0)               !> get column vector corresponding to matrix diagonal

G%AB%VAL(IC, ID) = G%AB%VAL(IC, ID) - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D)    G%AB%VAL(IC, ID) = G%AB%VAL(IC, ID) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
G%AB%VAL(IC, ID) = G%AB%VAL(IC, ID) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

G%AB%STENCIL(0) = G%AB%VAL(IC, ID)

END SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG_BANDED


!> ------------------------------------------------------------------------------------------------
!> Compute subdiagonal contribution in direction IOR0 in case of a COMPACT matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG_BANDED (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0)
USE SCARC_POINTERS, ONLY: L, F, G
INTEGER, INTENT(IN) :: IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0
INTEGER :: IW, ID
LOGICAL  :: IS_INTERNAL_CELL

F => L%FACE(IOR0)

!> Decide wheter cell is interior or exterior cell
ID = G%AB%POS(IOR0)                                
SELECT CASE (IOR0)
   CASE ( 1)
      IS_INTERNAL_CELL = IX1 > 1
   CASE (-1)
      IS_INTERNAL_CELL = IX1 < F%NOP
   CASE ( 2)
      IS_INTERNAL_CELL = IY1 > 1
   CASE (-2)
      IS_INTERNAL_CELL = IY1 < F%NOP
   CASE ( 3)
      IS_INTERNAL_CELL = IZ1 > 1
   CASE (-3)
      IS_INTERNAL_CELL = IZ1 < F%NOP
END SELECT

!>
!> if IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
!>
IF (IS_INTERNAL_CELL) THEN

   IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX2, IY2, IZ2)) THEN
      G%AB%VAL(IC, ID)   = G%AB%VAL(IC, ID) + F%SCAL_INSIDE
      G%AB%STENCIL(IOR0) = G%AB%VAL(IC, ID)
   ELSE
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

!>
!> if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
!>
!ELSE IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
ELSE IF (L%FACE(IOR0)%N_NEIGHBORS == 123456) THEN       !> CAUTION: TO FIX AGAIN, ONLY FOR TESTING, IMPOSSIBLE CONDITION

   IW = SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)            !> get IW of a possibly suitable neighbor at face IOR0
   IF (IW /= 0) THEN
      G%AB%VAL(IC, ID)   = G%AB%VAL(IC, ID) + F%SCAL_FACE
      G%AB%STENCIL(IOR0) = G%AB%VAL(IC, ID)
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG_BANDED


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Build system matrix for MKL solver in double precision
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MKL_DOUBLE (NM, NL)
USE SCARC_POINTERS, ONLY: G, AC, ACS
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, JC, JC0
INTEGER :: ICS, JCS
INTEGER :: ICOL, JCOL, IAS
INTEGER :: ISYM, JSYM, NSYM
REAL(EB) :: VAL = 0.0_EB, VALS = 0.0_EB, DIFF
LOGICAL  :: BSYM
INTEGER, DIMENSION(:), ALLOCATABLE :: KCOL_AUX, KC_AUX

CALL SCARC_POINT_TO_GRID(NM, NL)

AC  => G%AC
ACS => G%ACS

!> 
!> ---------- Store only symmetric parts of matrix (diagonal and upper part)
!> 
IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN

   !> First check whether symmetry of system matrix is guaranteed
   DO IC = 1, G%NC

      COLUMN_LOOP: DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1
         ICS = AC%COL(ICOL)
         VAL = AC%VAL(ICOL)
         IF (ICS > IC .AND. ICS <= G%NC) THEN
            BSYM = .FALSE.
            DO JCOL = AC%ROW(ICS)+1, AC%ROW(ICS+1)-1
               JCS = AC%COL(JCOL)
               IF (JCS == IC) THEN
                  VALS = AC%VAL(JCOL)
                  DIFF = ABS(VAL-VALS)
                  IF (ABS(VAL - VALS) < 1E-6) THEN
                     BSYM=.TRUE.
                     CYCLE COLUMN_LOOP
                  ENDIF
               ENDIF
            ENDDO
            IF (.NOT.BSYM) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SYMMETRY, SCARC_NONE, NM)
         ENDIF
      ENDDO COLUMN_LOOP

   ENDDO

   !>
   !> Compute number of entries in symmetric matrix
   !>
   ACS%N_VAL = 0
   DO IC = 1, G%NC
      DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
         IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL) THEN
            JC = AC%COL(ICOL)
            IF (JC >= IC .AND. JC <= G%NC) ACS%N_VAL = ACS%N_VAL+1
         ELSE IF (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) THEN
            JC = AC%COL_GLOBAL(ICOL)
            IF (JC >= IC + G%NC_OFFSET(NM)) ACS%N_VAL = ACS%N_VAL+1
         ELSE
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SETUP, SCARC_NONE, TYPE_MKL(NL))
         ENDIF
      ENDDO
   ENDDO

ELSE
   ACS%N_VAL = AC%N_VAL
ENDIF

!>
!> allocate storage for symmetric matrix and its column and row pointers
!> 
ACS%N_ROW = AC%N_ROW
CALL SCARC_ALLOCATE_MATRIX_COMPACT(ACS, NL, NSCARC_PRECISION_DOUBLE, NSCARC_INIT_ZERO, 'ACS')

!> if global MKL method is used, also allocate auxiliary space for computation of global numbering
IF (IS_MKL_LEVEL(NL)) THEN
   CALL SCARC_ALLOCATE_INT1(KCOL_AUX, 1, AC%N_STENCIL, NSCARC_INIT_NONE, 'KCOL_AUX')
   CALL SCARC_ALLOCATE_INT1(KC_AUX  , 1, AC%N_STENCIL, NSCARC_INIT_NONE, 'KC_AUX')
ENDIF

!> 
!> extract symmetric matrix part from usual system matrix
!> 
IAS = 1
DO IC = 1, G%NC
   ACS%ROW(IC) = IAS


   !> blockwise use of local MKL solvers - no global numbering required
   IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL) THEN

      DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
         JC = AC%COL(ICOL)

         IF (JC >= IC .AND. JC <= G%NC) THEN
            ACS%COL(IAS) = AC%COL(ICOL)
            ACS%VAL(IAS) = AC%VAL(ICOL)
            IAS = IAS + 1
         ENDIF
      ENDDO

   !> global use of MKL solver - get global numbering of matrix elements
   ELSE IF (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) THEN

      !> store indices of all diagonal and upper-diagonal entries
      KCOL_AUX = 0
      KC_AUX   = 99999999
      ISYM = 1
      JC0 = AC%COL_GLOBAL(AC%ROW(IC))
      DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
         !JC = AC%COL(ICOL)
         JC = AC%COL_GLOBAL(ICOL)
         IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
            IF (JC >= JC0) THEN
               KCOL_AUX(ISYM) = ICOL
               KC_AUX(ISYM) = JC
               ISYM  = ISYM  + 1
            ENDIF
         ELSE
            KCOL_AUX(ISYM) = ICOL
            KC_AUX(ISYM) = JC
            ISYM  = ISYM  + 1
         ENDIF
      ENDDO
      NSYM = ISYM - 1

      !> sort them in increasing order (for the use of Cluster_Sparse_Solver and PARDISO functionality)
      JSYM = 1
      SORT_LOOP: DO WHILE (JSYM <= NSYM)
         DO ISYM = 1, NSYM
            JC = KC_AUX(ISYM)
            IF (JC == 99999999) CYCLE
            IF (JC <= MINVAL(KC_AUX)) THEN
               ICOL = KCOL_AUX(ISYM)
               ACS%VAL(IAS) = AC%VAL(ICOL)
               !ACS%COL(IAS) = AC%COL(ICOL)
               ACS%COL(IAS) = AC%COL_GLOBAL(ICOL)
               KC_AUX(ISYM) = 99999999            ! mark entry as already used
               IAS  = IAS  + 1
            ENDIF
         ENDDO
         JSYM = JSYM + 1
      ENDDO SORT_LOOP
   ENDIF
ENDDO

ACS%ROW(G%NC+1) = IAS

IF (IS_MKL_LEVEL(NL)) THEN
   DEALLOCATE (KCOL_AUX, STAT=IERROR)
   DEALLOCATE (KC_AUX,   STAT=IERROR)
   DEALLOCATE (AC%COL_GLOBAL, STAT=IERROR)
ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_MKL_DOUBLE


!> ------------------------------------------------------------------------------------------------
!> Build system matrix for MKL solver in single precision
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MKL_SINGLE (NM, NL)
USE SCARC_POINTERS, ONLY: G, AC, ACS
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, JC, JC0
INTEGER :: ICS, JCS
INTEGER :: ICOL, JCOL, IAS
INTEGER :: ISYM, JSYM, NSYM
REAL(EB) :: VAL = 0.0_EB, VALS = 0.0_EB, DIFF
LOGICAL  :: BSYM
INTEGER, DIMENSION(:), ALLOCATABLE :: KCOL_AUX, KC_AUX

CALL SCARC_POINT_TO_GRID(NM, NL)

AC  => G%AC
ACS => G%ACS

!> 
!> ---------- Store only symmetric parts of matrix (diagonal and upper part)
!> 
IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN

   !> First check whether symmetry of system matrix is guaranteed
   DO IC = 1, G%NC

      COLUMN_LOOP: DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1
         ICS = AC%COL(ICOL)
         VAL = AC%VAL(ICOL)
         IF (ICS > IC .AND. ICS <= G%NC) THEN
            BSYM = .FALSE.
            DO JCOL = AC%ROW(ICS)+1, AC%ROW(ICS+1)-1
               JCS = AC%COL(JCOL)
               IF (JCS == IC) THEN
                  VALS = AC%VAL(JCOL)
                  DIFF = ABS(VAL-VALS)
                  IF (ABS(VAL - VALS) < 1E-6) THEN
                     BSYM=.TRUE.
                     CYCLE COLUMN_LOOP
                  ENDIF
               ENDIF
            ENDDO
            IF (.NOT.BSYM) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SYMMETRY, SCARC_NONE, NM)
         ENDIF
      ENDDO COLUMN_LOOP

   ENDDO

   !>
   !> Compute number of entries in symmetric matrix
   !>
   ACS%N_VAL = 0
   DO IC = 1, G%NC
      DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
         IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL) THEN
            JC = AC%COL(ICOL)
            IF (JC >= IC .AND. JC <= G%NC) ACS%N_VAL = ACS%N_VAL+1
         ELSE IF (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) THEN
            JC = AC%COL_GLOBAL(ICOL)
            IF (JC >= IC + G%NC_OFFSET(NM)) ACS%N_VAL = ACS%N_VAL+1
         ELSE
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SETUP, SCARC_NONE, TYPE_MKL(NL))
         ENDIF
      ENDDO
   ENDDO

ELSE
   ACS%N_VAL = AC%N_VAL
ENDIF

!>
!> allocate storage for symmetric matrix and its column and row pointers
!> 
ACS%N_ROW = AC%N_ROW
CALL SCARC_ALLOCATE_MATRIX_COMPACT(ACS, NL, NSCARC_PRECISION_SINGLE, NSCARC_INIT_ZERO, 'ACS')

!> if global MKL method is used, also allocate auxiliary space for computation of global numbering
IF (IS_MKL_LEVEL(NL)) THEN
   CALL SCARC_ALLOCATE_INT1(KCOL_AUX, 1, AC%N_STENCIL, NSCARC_INIT_NONE, 'KCOL_AUX')
   CALL SCARC_ALLOCATE_INT1(KC_AUX  , 1, AC%N_STENCIL, NSCARC_INIT_NONE, 'KC_AUX')
ENDIF

!> 
!> extract symmetric matrix part from usual system matrix
!> 
IAS = 1
DO IC = 1, G%NC
   ACS%ROW(IC) = IAS

   !> blockwise use of local MKL solvers - no global numbering required
   IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL) THEN

      DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
         JC = AC%COL(ICOL)

         IF (JC >= IC .AND. JC <= G%NC) THEN
            ACS%COL(IAS) = AC%COL(ICOL)
            ACS%VAL_FB(IAS) = REAL(AC%VAL(ICOL),FB)
            IAS = IAS + 1
         ENDIF
      ENDDO

   !> global use of MKL solver - get global numbering of matrix elements
   ELSE IF (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) THEN

      !> store indices of all diagonal and upper-diagonal entries
      KCOL_AUX = 0
      KC_AUX   = 99999999
      ISYM = 1
      JC0 = AC%COL_GLOBAL(AC%ROW(IC))
      DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
         !JC = AC%COL(ICOL)
         JC = AC%COL_GLOBAL(ICOL)
         IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
            IF (JC >= JC0) THEN
               KCOL_AUX(ISYM) = ICOL
               KC_AUX(ISYM) = JC
               ISYM  = ISYM  + 1
            ENDIF
         ELSE
            KCOL_AUX(ISYM) = ICOL
            KC_AUX(ISYM) = JC
            ISYM  = ISYM  + 1
         ENDIF
      ENDDO
      NSYM = ISYM - 1

      !> sort them in increasing order (for the use of Cluster_Sparse_Solver and PARDISO functionality)
      JSYM = 1
      SORT_LOOP: DO WHILE (JSYM <= NSYM)
         DO ISYM = 1, NSYM
            JC = KC_AUX(ISYM)
            IF (JC == 99999999) CYCLE
            IF (JC <= MINVAL(KC_AUX)) THEN
               ICOL = KCOL_AUX(ISYM)
               ACS%VAL_FB(IAS) = REAL(AC%VAL(ICOL),FB)
               !ACS%COL(IAS) = AC%COL(ICOL)
               ACS%COL(IAS) = AC%COL_GLOBAL(ICOL)
               KC_AUX(ISYM) = 99999999            ! mark entry as already used
               IAS  = IAS  + 1
            ENDIF
         ENDDO
         JSYM = JSYM + 1
      ENDDO SORT_LOOP
   ENDIF
ENDDO

ACS%ROW(G%NC+1) = IAS

IF (IS_MKL_LEVEL(NL)) THEN
   DEALLOCATE (KCOL_AUX, STAT=IERROR)
   DEALLOCATE (KC_AUX,   STAT=IERROR)
   DEALLOCATE (AC%COL_GLOBAL, STAT=IERROR)
ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_MKL_SINGLE
#endif


!> ------------------------------------------------------------------------------------------------
!> Insert correct boundary conditions into system matrix
!>
!> If A is a pure Neumann matrix, get neighboring cell indices of communicated stencil legs for 
!> condensed system, also save values and column indices of last matrix row of last mesh
!>
!> Set correct boundary conditions for system matrix
!> Take care of whether the structured or unstructured discretization is used
!>
!> If there are no Dirichlet BC's transform sytem into condensed one by replacing the
!> matrix entries in last column and row by the stored ones (zeros and one at diaonal position)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
USE SCARC_POINTERS, ONLY: L, G, F, GWC, AC, AB, ACC, ABC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP, ICO, ICOL

CALL SCARC_POINT_TO_GRID(NM, NL)

SELECT CASE (SCARC_MATRIX_LEVEL(NL))

   !>
   !> ---------- Matrix in COMPACT storage technique
   !>
   CASE (NSCARC_MATRIX_COMPACT)

      AC => G%AC

      !> Setup condensing if there are no Dirichlet BC's 
      IF (IS_PURE_NEUMANN) CALL SCARC_SETUP_MATRIX_COMPACT_CONDENSED(NM)

      !> Set correct boundary conditions 
      DO IW = 1, G%NW

         GWC => G%WALL(IW)
         IOR0 = GWC%IOR
         IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE       

         F  => L%FACE(IOR0)

         I = GWC%IXW
         J = GWC%IYW
         K = GWC%IZW

         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

         NOM = GWC%NOM
         IC  = G%CELL_NUMBER(I, J, K)
         GWC%ICW = IC

         !> SPD-matrix with mixture of Dirichlet and Neumann BC's according to BTYPE
         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN

            IP = AC%ROW(IC)
            SELECT CASE (GWC%BTYPE)
               CASE (DIRICHLET)
                  AC%VAL(IP) = AC%VAL(IP) - F%SCAL_BOUNDARY
               CASE (NEUMANN)
                  AC%VAL(IP) = AC%VAL(IP) + F%SCAL_BOUNDARY
            END SELECT

         !> purely Neumann matrix
         ELSE IF (GWC%BTYPE == NEUMANN) THEN
            IP = AC%ROW(IC)
            AC%VAL(IP) = AC%VAL(IP) + F%SCAL_BOUNDARY
         ENDIF

      ENDDO 

      !> Transform into condensed system, if there are no Dirichlet BC's 
      IF (IS_PURE_NEUMANN) THEN
         DO ICO = 1, AC%N_CONDENSED
            ACC => AC%CONDENSED(ICO)
            DO ICOL = 1, ACC%N_COL
               IP = ACC%PTR(ICOL)
               AC%VAL(IP) = ACC%VAL2(ICOL)
            ENDDO
         ENDDO
      ENDIF 

   !>
   !> ---------- Matrix in Banded storage technique
   !>
   CASE (NSCARC_MATRIX_BANDED)

      !> Preset matrix switch if no Dirichlet BC's available
      AB => G%AB
      IF (IS_PURE_NEUMANN) CALL SCARC_SETUP_MATRIX_BANDED_CONDENSED(NM)

      !> Set right boundary conditions 
      DO IW = 1, G%NW

         GWC => G%WALL(IW)
         IOR0 = GWC%IOR
         IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE     

         F  => L%FACE(IOR0)

         I = GWC%IXW
         J = GWC%IYW
         K = GWC%IZW

         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

         NOM  = GWC%NOM
         GWC%ICW =G%CELL_NUMBER(I, J, K)
         IC = G%CELL_NUMBER(I, J, K)

         !> SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE
         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN

            SELECT CASE (GWC%BTYPE)
               CASE (DIRICHLET)
                  AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) - F%SCAL_BOUNDARY
               CASE (NEUMANN)
                  AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) + F%SCAL_BOUNDARY
            END SELECT

         !> Purely Neumann matrix
         ELSE
            IF (GWC%BTYPE == NEUMANN) AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) + F%SCAL_BOUNDARY
         ENDIF

      ENDDO 
   
      !> Transform into condensed system, if there are no Dirichlet BC's 
      IF (IS_PURE_NEUMANN) THEN
         DO ICO = 1, AB%N_CONDENSED
            ABC => AB%CONDENSED(ICO)
            IF (ICO == 1) THEN
               AB%VAL(ABC%ICO, 1:AB%N_STENCIL) = ABC%VAL2(1:AB%N_STENCIL)
            ELSE
               IP = AB%POS(ABC%IOR0)
               AB%VAL(ABC%ICO, IP) = ABC%VAL2(IP)
            ENDIF
         ENDDO
      ENDIF 
 
END SELECT 

END SUBROUTINE SCARC_SETUP_BOUNDARY


!> ------------------------------------------------------------------------------------------------
!> Compact matrix storage technique:
!> Define switch entries for toggle between original and condensed values
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_COMPACT_CONDENSED (NM)
USE SCARC_POINTERS, ONLY: L, G, AC, ACC
INTEGER, INTENT(IN) :: NM
INTEGER :: ICO = 0, NC, NOM, IP, IC, JC, ICE, ICN, ICOL, IOR0, IW, I, J, K

AC => G%AC
LAST_CELL_IN_LAST_MESH_IF: IF (NM == NMESHES) THEN

   NC = G%NC_LOCAL(NMESHES)
   IP = AC%ROW(NC)

   !> store column indices and values of diagonal and all off-diagonal entries in last row
   !> index '1' corresponds to main diagonal entry
   ICO = ICO + 1
   ACC => AC%CONDENSED(ICO)

   ICOL = 1
   ACC%PTR(ICOL)  = IP
   ACC%COL(ICOL)  = AC%COL(IP)
   ACC%VAL1(ICOL) = AC%VAL(IP)
   ACC%VAL2(ICOL) = 1.0_EB

   DO IP = AC%ROW(NC)+1, AC%ROW(NC+1)-1
      ICOL = ICOL + 1
      ACC%PTR(ICOL)  = IP
      ACC%COL(ICOL)  = AC%COL(IP)
      ACC%VAL1(ICOL) = AC%VAL(IP)
      ACC%VAL2(ICOL) = 0.0_EB
   ENDDO
   ACC%N_COL = ICOL                                !> number of stored columns

   !>
   !> within last mesh: check which other cells have a connection to the last cell;
   !> in each corresponding matrix row store the column index and value of just that matrix entry
   !> for each direction only one value has to be stored
   !>
   JC = NC - 1
   DO IP = AC%ROW(JC)+1, AC%ROW(JC+1)-1
      IF (AC%COL(IP) == NC) THEN
         ICO = ICO + 1
         ACC => AC%CONDENSED(ICO)
         ACC%PTR(1)  = IP
         ACC%COL(1)  = JC
         ACC%VAL1(1) = AC%VAL(IP)                  ! store original value of system matrix
         ACC%VAL2(1) = 0.0_EB                      ! store new value of condensed system matrix
         ACC%N_COL   = 1
         EXIT
      ENDIF
   ENDDO

   JC = NC - L%NX
   DO IP = AC%ROW(JC)+1, AC%ROW(JC+1)-1
      IF (AC%COL(IP) == NC) THEN
         ICO = ICO + 1
         ACC => AC%CONDENSED(ICO)
         ACC%PTR(1)  = IP
         ACC%COL(1)  = JC
         ACC%VAL1(1) = AC%VAL(IP)                  ! store original value of system matrix
         ACC%VAL2(1) = 0.0_EB                      ! store new value of condensed system matrix
         ACC%N_COL   = 1
         EXIT
      ENDIF
   ENDDO

   IF (.NOT.TWO_D) THEN
      JC = NC - L%NX * L%NY
      DO IP = AC%ROW(JC)+1, AC%ROW(JC+1)-1
         IF (AC%COL(IP) == NC) THEN
            ICO = ICO + 1
            ACC => AC%CONDENSED(ICO)
            ACC%PTR(1)  = IP
            ACC%COL(1)  = JC
            ACC%VAL1(1) = AC%VAL(IP)                  ! store original value of system matrix
            ACC%VAL2(1) = 0.0_EB                      ! store new value of condensed system matrix
            ACC%N_COL   = 1
            EXIT
         ENDIF
      ENDDO
   ENDIF

ENDIF LAST_CELL_IN_LAST_MESH_IF

!>
!> cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled
!> with the last cell of last mesh;
!> this can be a cell on the opposite side of the own mesh or on a different mesh
!> if such a cell exists, store corresponding matrix entry
!>
DO IW = 1, G%NW

   GWC => G%WALL(IW)

   IOR0 = GWC%IOR
   IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE

   I    = GWC%IXW
   J    = GWC%IYW
   K    = GWC%IZW

   IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

   NOM = GWC%NOM
   IC  = G%CELL_NUMBER(I, J, K)
   GWC%ICW = IC

   IF (NOM == NMESHES) THEN

      ICE = GWC%ICE                               !> adjacent ghost cell number
      ICN = G%ICE_TO_ICN(ICE)                     !> get column index of neighboring offdiagonal matrix entry
      IF (ICN /= SCARC(NMESHES)%NC) CYCLE         !> if no relation to last cell in last mesh, cycle

      DO IP = AC%ROW(IC)+1, AC%ROW(IC+1)-1
         IF (AC%COL(IP) == ICE) THEN
            ICO = ICO + 1
            ACC => AC%CONDENSED(ICO)
            ACC%PTR(1)  = IP
            ACC%COL(1)  = ICN
            ACC%VAL1(1) = AC%VAL(IP)
            ACC%VAL2(1) = 0.0_EB
            ACC%N_COL   = 1
            EXIT
         ENDIF
      ENDDO

   ENDIF 
ENDDO 

AC%N_CONDENSED = ICO

END SUBROUTINE SCARC_SETUP_MATRIX_COMPACT_CONDENSED


!> ------------------------------------------------------------------------------------------------
!> Banded matrix storage technique:
!> Define switch entries for toggle between original and condensed values
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_BANDED_CONDENSED (NM)
USE SCARC_POINTERS, ONLY: L, G, AB, ABC
INTEGER, INTENT(IN) :: NM
INTEGER :: ICO = 0, NC, NOM, IOR0, IC, JC, ICE, ICN, IW, I, J, K

AB => G%AB
LAST_CELL_IN_LAST_MESH_BANDED_IF: IF (NM == NMESHES) THEN

   NC = G%NC_LOCAL(NMESHES)

   !> store column indices and values of diagonal and all off-diagonal entries in last row
   !> index '1' corresponds to main diagonal entry
   ICO = ICO + 1
   ABC => AB%CONDENSED(ICO)

   ABC%IOR0 = 0
   ABC%ICO  = NC
   ABC%VAL1(1:AB%N_STENCIL) = AB%VAL(NC, 1:AB%N_STENCIL)
   ABC%VAL2(1:AB%N_STENCIL) = 0.0_EB
   ABC%VAL2(AB%POS(0)) = 1.0_EB

   !>
   !> within last mesh: check which other cells have a connection to the last cell;
   !> in each corresponding matrix row store the column index and value of just that matrix entry
   !> for each direction only one value has to be stored
   !>
   JC = NC - 1
   DO IOR0 = -3, 3
      IF (JC + AB%OFFSET(IOR0) == NC) THEN
         ICO = ICO + 1
         ABC => AB%CONDENSED(ICO)
         ABC%IOR0 = IOR0
         ABC%ICO  = JC
         ABC%VAL1(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABC%VAL2(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABC%VAL2(AB%POS(ABC%IOR0)) = 0.0_EB
         EXIT
      ENDIF
   ENDDO

   JC = NC - L%NX
   DO IOR0 = -3, 3
      IF (JC + AB%OFFSET(IOR0) == NC) THEN
         ICO = ICO + 1
         ABC => AB%CONDENSED(ICO)
         ABC%IOR0 = IOR0
         ABC%ICO  = JC
         ABC%VAL1(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABC%VAL2(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABC%VAL2(AB%POS(ABC%IOR0)) = 0.0_EB
         EXIT
      ENDIF
   ENDDO

   IF (.NOT.TWO_D) THEN
      JC = NC - L%NX * L%NY
      DO IOR0 = -3, 3
         IF (JC + AB%OFFSET(IOR0) == NC) THEN
            ICO = ICO + 1
            ABC => AB%CONDENSED(ICO)
            ABC%IOR0 = IOR0
            ABC%ICO  = JC
            ABC%VAL1(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
            ABC%VAL2(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
            ABC%VAL2(AB%POS(ABC%IOR0)) = 0.0_EB
            EXIT
         ENDIF
      ENDDO
   ENDIF

ENDIF LAST_CELL_IN_LAST_MESH_BANDED_IF

!>
!> cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled
!> with the last cell of last mesh;
!> this can be a cell on the opposite side of the own mesh or a cell on a different mesh
!> if such a cell exists, store corresponding matrix entry
!>
DO IW = 1, G%NW

   GWC => G%WALL(IW)

   IOR0 = GWC%IOR
   IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE

   I    = GWC%IXW
   J    = GWC%IYW
   K    = GWC%IZW

   IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

   NOM = GWC%NOM
   IC  = G%CELL_NUMBER(I, J, K)
   GWC%ICW = IC

   IF (NOM == NMESHES) THEN
      ICE = GWC%ICE                               !> adjacent ghost cell number
      ICN = G%ICE_TO_ICN(ICE)                     !> get column index of neighboring offdiagonal matrix entry
      IF (ICN /= SCARC(NMESHES)%NC) CYCLE         !> if no relation to last cell in last mesh, cycle
      ICO = ICO + 1
      ABC => AB%CONDENSED(ICO)
      ABC%IOR0 = IOR0
      ABC%ICO  = IC
      ABC%VAL1(1:AB%N_STENCIL) = AB%VAL(IC, 1:AB%N_STENCIL)
      ABC%VAL2(1:AB%N_STENCIL) = AB%VAL(IC, 1:AB%N_STENCIL)
      ABC%VAL2(AB%POS(ABC%IOR0)) = 0.0_EB
      EXIT
   ENDIF 
ENDDO 

AB%N_CONDENSED = ICO

END SUBROUTINE SCARC_SETUP_MATRIX_BANDED_CONDENSED



!> ------------------------------------------------------------------------------------------------
!> Setup condensed system in case of periodic or pure Neumann boundary conditions
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM_CONDENSED (NV, NL, ITYPE)
USE SCARC_POINTERS, ONLY: G, F, OL, AC, AB, ACC, ABC, VC
INTEGER, INTENT(IN) :: NV, NL, ITYPE
INTEGER :: NM, NOM, IFACE, ICN, ICE, ICW, JC, NC, ICO, IOR0, IP, ICG, INBR

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0 .OR. &
    TYPE_PRECON == NSCARC_RELAX_FFT .OR. TYPE_PRECON == NSCARC_RELAX_FFTO) RETURN

!>
!> In last mesh:
!> Subtract B*RHS(end) for internal legs of stencil
!>
MESH_REAL = 0.0_EB
IF (UPPER_MESH_INDEX == NMESHES) THEN

   CALL SCARC_POINT_TO_GRID(NMESHES, NL)

   NC =  G%NC_LOCAL(NMESHES)
   VC => SCARC_POINT_TO_VECTOR(NMESHES, NL, NV)

   !>
   !> process last column entries of all rows except of last one
   !> for those rows only one matrix entry was stored, namely that one which connects to the last cell
   !>
   SELECT CASE (SCARC_MATRIX_LEVEL(NL))

      CASE (NSCARC_MATRIX_COMPACT)
         AC => G%AC
         DO ICO = 2, AC%N_CONDENSED
            ACC => AC%CONDENSED(ICO)
            JC = ACC%COL(1)
            IF (JC < NC) VC(JC) = VC(JC) - ACC%VAL1(1)*VC(NC)
         ENDDO

      CASE (NSCARC_MATRIX_BANDED)
         AB => G%AB
         DO ICO = 2, AB%N_CONDENSED
            ABC => AB%CONDENSED(ICO)
            IP = AB%POS(ABC%IOR0)
            JC = ABC%ICO
            IF (JC < NC) VC(JC) = VC(JC) - ABC%VAL1(IP)*VC(NC)
        ENDDO

   END SELECT

   MESH_REAL(NMESHES) = VC(NC)     !> store last entry of RHS
   VC(NC) = 0.0_EB                 !> set last entry of last mesh to zero

ENDIF

IF (ITYPE == 0) RETURN

!>
!> broadcast last RHS-value of last cell in last mesh to all meshes
!>
IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MESH_REAL, 1, MPI_DOUBLE_PRECISION,&
                      MPI_COMM_WORLD, IERROR)

DO NM = 1, NMESHES
   SCARC(NM)%RHS_END = MESH_REAL(NMESHES)
ENDDO

!>
!> Only in case of periodic BC's:
!> Subtract B*RHS(end) for corresponding entries of all periodic communication partners
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   SNODE = PROCESS(NM)
   RNODE = PROCESS(NMESHES)

   IF (.NOT. ARE_NEIGHBORS(NM, NMESHES)) CYCLE

   CALL SCARC_POINT_TO_NEIGHBOR(NM, NMESHES, NL)
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   !>
   !> subtract B*RHS(end) at corresponding positions
   !>
   DO IFACE = 1, 6                                         !> check if this face has connection to last cell

      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)

      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
         IF (NOM /= NMESHES) CYCLE                      !> only check for common matrix entries with last mesh
         CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)

         DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)

            ICW = OG%ICG_TO_ICW(ICG)
            ICE = OG%ICG_TO_ICE(ICG)
            ICN = G%ICE_TO_ICN(ICE)                        !> get column index of neighboring offdiagonal matrix entry

            IF (ICN /= SCARC(NMESHES)%NC) CYCLE            !> if no relation to last cell in last mesh, cycle

            !VC(ICW) = VC(ICW) - G%ICE_TO_VAL(ICE) * SCARC(NM)%RHS_END
            VC(ICW) = VC(ICW) - F%SCAL_FACE * SCARC(NM)%RHS_END

         ENDDO

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_SYSTEM_CONDENSED


!> ----------------------------------------------------------------------------------------------------
!> Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_METHODS
INTEGER :: NSTACK

SELECT_METHOD: SELECT CASE(TYPE_METHOD)


   !>
   !> ------------------ Global Krylov method -------------------------------------
   !>
   CASE (NSCARC_METHOD_KRYLOV)

      !>
      !> Setup basic CG solver
      !>
      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_CG
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !>
      !> Setup preconditioner for Krylov solver
      !>
      NSTACK = NSTACK + 1
      SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

         !> Jacobi-preconditioning (acting locally by default)
         CASE (NSCARC_RELAX_JAC)
            STACK(NSTACK)%SOLVER => PRECON_JAC
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

         !> SSOR-preconditioning (acting locally by default)
         CASE (NSCARC_RELAX_SSOR)
            STACK(NSTACK)%SOLVER => PRECON_SSOR
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

         !> JACOBI-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_MJAC)
            STACK(NSTACK)%SOLVER => PRECON_MJAC
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MJAC(NLEVEL_MIN, NLEVEL_MAX)

         !> GS-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_MGS)
            STACK(NSTACK)%SOLVER => PRECON_MGS
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MGS(NLEVEL_MIN, NLEVEL_MAX)

         !> SGS-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_MSGS)
            STACK(NSTACK)%SOLVER => PRECON_MSGS
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MSGS(NLEVEL_MIN, NLEVEL_MAX)

         !> SOR-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_MSOR)
            STACK(NSTACK)%SOLVER => PRECON_MSOR
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

         !> SSOR-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_MSSOR)
            STACK(NSTACK)%SOLVER => PRECON_MSSOR
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MSSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

         !> ILU(0)-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_ILU)
            STACK(NSTACK)%SOLVER => PRECON_ILU
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_ILU(NLEVEL_MIN, NLEVEL_MAX)

         !> FFT-preconditioning (acting locally by default)
         CASE (NSCARC_RELAX_FFT)
            STACK(NSTACK)%SOLVER => PRECON_FFT
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MIN)

         !> FFT-preconditioning (acting locally by default)
         CASE (NSCARC_RELAX_FFTO)
            STACK(NSTACK)%SOLVER => PRECON_FFT
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_FFTO(NLEVEL_MIN, NLEVEL_MIN)

#ifdef WITH_MKL
         !> LU-preconditioning based on MKL (either locally or Globally acting depending on user specification)
         CASE (NSCARC_RELAX_MKL)
            STACK(NSTACK)%SOLVER => PRECON_MKL

            SELECT CASE(TYPE_SCOPE(1))

               !> Globally acting - call global CLUSTER_SPARSE_SOLVER from MKL
               CASE (NSCARC_SCOPE_GLOBAL)
                  CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_GLOBAL)
                  CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)

               !> locally acting - call global PARDISO solver from MKL
               CASE (NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)

            END SELECT
#endif

         !>
         !> Preconditioning by Geometric multigrid,
         !> either locally or Globally acting, depending on user specification stored in TYPE_SCOPE(1)
         !>
         CASE (NSCARC_RELAX_GMG)

            STACK(NSTACK)%SOLVER => PRECON_GMG
            CALL SCARC_SETUP_PRECON(NSTACK, TYPE_SCOPE(1))
            CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_PRECON, TYPE_SCOPE(1), NSCARC_STAGE_TWO, NSTACK, &
                                       NLEVEL_MIN, NLEVEL_MAX)

            NSTACK = NSTACK + 1
            SELECT CASE (TYPE_SMOOTH)

               !> Jacobi-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_JAC)
                  STACK(NSTACK)%SOLVER => SMOOTH_JACOBI
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

               !> SSOR-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_SSOR)
                  STACK(NSTACK)%SOLVER => SMOOTH_SSOR
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

               !> JACOBI-preconditioning in matrix form (acting locally by default)
               CASE (NSCARC_RELAX_MJAC)
                  STACK(NSTACK)%SOLVER => SMOOTH_MJAC
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_MJAC(NLEVEL_MIN, NLEVEL_MAX-1)

               !> GS-preconditioning in matrix form (acting locally by default)
               CASE (NSCARC_RELAX_MGS)
                  STACK(NSTACK)%SOLVER => SMOOTH_MGS
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_MGS(NLEVEL_MIN, NLEVEL_MAX-1)

               !> SGS-preconditioning in matrix form (acting locally by default)
               CASE (NSCARC_RELAX_MSGS)
                  STACK(NSTACK)%SOLVER => SMOOTH_MSGS
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_MSGS(NLEVEL_MIN, NLEVEL_MAX-1)

               !> SOR-preconditioning in matrix form (acting locally by default)
               CASE (NSCARC_RELAX_MSOR)
                  STACK(NSTACK)%SOLVER => SMOOTH_MSOR
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_MSOR(NLEVEL_MIN, NLEVEL_MAX-1, NSTACK)

               !> SSOR-preconditioning in matrix form (acting locally by default)
               CASE (NSCARC_RELAX_MSSOR)
                  STACK(NSTACK)%SOLVER => SMOOTH_MSSOR
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_MSSOR(NLEVEL_MIN, NLEVEL_MAX-1, NSTACK)

               !> FFT-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_FFT)
                  STACK(NSTACK)%SOLVER => SMOOTH_FFT
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MAX-1)

               !> FFTO-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_FFTO)
                  STACK(NSTACK)%SOLVER => SMOOTH_FFTO
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_FFTO(NLEVEL_MIN, NLEVEL_MAX-1)
#ifdef WITH_MKL
               !> LU-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_MKL)
                  STACK(NSTACK)%SOLVER => SMOOTH_MKL
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
#endif
            END SELECT

            !> Coarse grid solver (same scope of action as calling GMG)
            NSTACK = NSTACK + 1
            CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_TWO, TYPE_SCOPE(1), NSTACK, NLEVEL_MAX, NLEVEL_MAX)

      END SELECT SELECT_KRYLOV_PRECON

      !>
      !> If two-level Krylov, allocate intermediate structures for interpolation and workspace for global coarse solver
      !>
      IF (HAS_TWO_LEVELS) THEN

         CALL SCARC_SETUP_INTERPOLATION(NSCARC_STAGE_ONE, NLEVEL_MIN+1, NLEVEL_MAX)

         NSTACK = NSTACK + 1
         CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)

      ENDIF

    !>
    !> ------------------ Global Multigrid method -------------------------------------
    !>
    CASE (NSCARC_METHOD_MULTIGRID)

       NSTACK = NSCARC_STACK_ROOT
       STACK(NSTACK)%SOLVER => MAIN_GMG
       CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MAX)

       NSTACK = NSTACK + 1
       SELECT CASE(TYPE_SMOOTH)

          !> Jacobi-smoothing (acting locally by default)
          CASE (NSCARC_RELAX_JAC)
             STACK(NSTACK)%SOLVER => SMOOTH_JACOBI
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

          !> SSOR-smoothing (acting locally by default)
          CASE (NSCARC_RELAX_SSOR)
             STACK(NSTACK)%SOLVER => SMOOTH_SSOR
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

          !> JACOBI-preconditioning in matrix form (acting locally by default)
          CASE (NSCARC_RELAX_MJAC)
             STACK(NSTACK)%SOLVER => SMOOTH_MJAC
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
             CALL SCARC_SETUP_MJAC(NLEVEL_MIN, NLEVEL_MAX)

          !> GS-preconditioning in matrix form (acting locally by default)
          CASE (NSCARC_RELAX_MGS)
             STACK(NSTACK)%SOLVER => SMOOTH_MGS
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
             CALL SCARC_SETUP_MGS(NLEVEL_MIN, NLEVEL_MAX)

          !> SGS-preconditioning in matrix form (acting locally by default)
          CASE (NSCARC_RELAX_MSGS)
             STACK(NSTACK)%SOLVER => SMOOTH_MSGS
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
             CALL SCARC_SETUP_MSGS(NLEVEL_MIN, NLEVEL_MAX)

          !> SOR-preconditioning in matrix form (acting locally by default)
          CASE (NSCARC_RELAX_MSOR)
             STACK(NSTACK)%SOLVER => SMOOTH_MSOR
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
             CALL SCARC_SETUP_MSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

          !> SSOR-preconditioning in matrix form (acting locally by default)
          CASE (NSCARC_RELAX_MSSOR)
             STACK(NSTACK)%SOLVER => SMOOTH_MSSOR
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
             CALL SCARC_SETUP_MSSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

          !> FFT-smoothing (acting locally by default)
          CASE (NSCARC_RELAX_FFT)
             STACK(NSTACK)%SOLVER => SMOOTH_FFT
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
             CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MAX-1)

          !> FFTO-smoothing (acting locally by default)
          CASE (NSCARC_RELAX_FFTO)
             STACK(NSTACK)%SOLVER => SMOOTH_FFTO
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
             CALL SCARC_SETUP_FFTO(NLEVEL_MIN, NLEVEL_MAX-1)

#ifdef WITH_MKL
          !> smoothing by LU-decomposition
          CASE (NSCARC_RELAX_MKL)
             CALL SCARC_SETUP_SMOOTH(NSTACK, TYPE_SCOPE(2))

             SELECT CASE(TYPE_SCOPE(2))

                !> Globally acting - call global CLUSTER_SPARSE_SOLVER on MKL
                CASE (NSCARC_SCOPE_GLOBAL)
                   !WRITE(MSG%LU_VERBOSE,*) 'B: SELECT GLOBAL SMOOTHING'
                   STACK(NSTACK)%SOLVER => SMOOTH_MKL
                   CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)

                !> locally acting - call local PARDISO solvers based on MKL
                CASE (NSCARC_SCOPE_LOCAL)
                   !WRITE(MSG%LU_VERBOSE,*) 'B: SELECT LOCAL SMOOTHING'
                   STACK(NSTACK)%SOLVER => SMOOTH_MKL
                   CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)

             END SELECT
#endif

       END SELECT

       !> Globally acting coarse grid solver
       NSTACK = NSTACK + 1
       CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)

#ifdef WITH_SCARC_MGM
   !>
   !> ---------------- McKenny-Greengard-Mayo method (MGM) --------------------
   !>
   CASE (NSCARC_METHOD_MGM)

      !> Allocate velocity vectors along internal obstructions for the setting of internal BC's
      CALL SCARC_SETUP_MGM(NLEVEL_MIN, NLEVEL_MIN)

      !> ------- First part of method: Setup CG solver for inhomogeneous problem on structured discretization
      TYPE_GRID = NSCARC_GRID_STRUCTURED

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_CG_STRUCTURED
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !> use FFT-preconditioning
      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => PRECON_FFT
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MIN)

      !> ------- Second part of method: Setup CG solver for homogeneous problem on unstructured discretization
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED

      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => MAIN_CG_UNSTRUCTURED
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !> for a first proof of concept only use SSOR-preconditioning (may be extended later to other preconditioners)
      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => PRECON_ILU
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

#endif

   !>
   !> ------------------ MKL method -------------------------------------
   !>
#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LU)

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_LU

      CALL SCARC_SETUP_MKL(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !> In the multi-mesh case use CLUSTER_SPARSE_SOLVER, else PARDISO solver (only on finest grid level)
      IF (NMESHES > 1) THEN
         CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)
      ELSE
         CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
      ENDIF
#endif

END SELECT SELECT_METHOD

!> Store total number of stack entries (used solvers)
N_STACK_TOTAL = NSTACK

END SUBROUTINE SCARC_SETUP_METHODS


!> ----------------------------------------------------------------------------------------------------
!> Setup references to solution vectors related to used scope (main/relax)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_REFERENCES(BX, BB, BD, BR, BV, BY, BZ, NSTACK)
USE SCARC_POINTERS, ONLY: SV
LOGICAL, INTENT(IN) :: BX, BB, BD, BR, BV, BY, BZ
INTEGER, INTENT(IN) :: NSTACK

SV  => STACK(NSTACK)%SOLVER

SELECT CASE (SV%TYPE_STAGE)

   !>
   !> Solver from working stage ONE, e.g. main Krylov or MG solver
   !>
   CASE (NSCARC_STAGE_ONE)
      IF (BX) SV%X = NSCARC_VECTOR_ONE_X
      IF (BB) SV%B = NSCARC_VECTOR_ONE_B
      IF (BD) SV%D = NSCARC_VECTOR_ONE_D
      IF (BR) SV%R = NSCARC_VECTOR_ONE_R
      IF (BV) SV%V = NSCARC_VECTOR_ONE_V
      IF (BY) SV%Y = NSCARC_VECTOR_ONE_Y
      IF (BZ) SV%Z = NSCARC_VECTOR_ONE_Z

#ifdef WITH_MKL
     IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
        IF (BX) SV%X_FB = NSCARC_VECTOR_ONE_X
        IF (BB) SV%B_FB = NSCARC_VECTOR_ONE_B
        IF (BR) SV%R_FB = NSCARC_VECTOR_ONE_R
        IF (BV) SV%V_FB = NSCARC_VECTOR_ONE_V
     ENDIF
#endif

#ifdef WITH_SCARC_DEBUG
      SV%E = NSCARC_VECTOR_ONE_E
#endif

   !>
   !> Solver from working stage TWO, e.g. MG solver as preconditioner
   !>
   CASE (NSCARC_STAGE_TWO)
      IF (BX) SV%X = NSCARC_VECTOR_TWO_X
      IF (BB) SV%B = NSCARC_VECTOR_TWO_B
      IF (BD) SV%D = NSCARC_VECTOR_TWO_D
      IF (BR) SV%R = NSCARC_VECTOR_TWO_R
      IF (BV) SV%V = NSCARC_VECTOR_TWO_V
      IF (BY) SV%Y = NSCARC_VECTOR_TWO_Y
      IF (BZ) SV%Z = NSCARC_VECTOR_TWO_Z
#ifdef WITH_SCARC_DEBUG
      SV%E = NSCARC_VECTOR_TWO_E
#endif

#ifdef WITH_MKL
     IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
        IF (BX) SV%X_FB = NSCARC_VECTOR_TWO_X
        IF (BB) SV%B_FB = NSCARC_VECTOR_TWO_B
        IF (BR) SV%R_FB = NSCARC_VECTOR_TWO_R
        IF (BV) SV%V_FB = NSCARC_VECTOR_TWO_V
     ENDIF
#endif

END SELECT

END SUBROUTINE SCARC_SETUP_REFERENCES


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_VECTORS()
USE SCARC_POINTERS, ONLY: G, SV, ST
INTEGER :: NM, NSTACK, NL

DO NSTACK = 1, N_STACK_TOTAL

   SV  => STACK(NSTACK)%SOLVER

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      DO NL = SV%TYPE_LEVEL(1), SV%TYPE_LEVEL(2)

         CALL SCARC_POINT_TO_GRID(NM, NL)
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)

         IF (SV%X /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%X, 1, G%NCE, NSCARC_INIT_ZERO, 'X')
         IF (SV%B /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%B, 1, G%NCE, NSCARC_INIT_ZERO, 'B')
         IF (SV%D /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%D, 1, G%NCE, NSCARC_INIT_ZERO, 'D')
         IF (SV%R /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%R, 1, G%NCE, NSCARC_INIT_ZERO, 'R')
         IF (SV%V /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%V, 1, G%NCE, NSCARC_INIT_ZERO, 'V')
         IF (SV%Y /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%Y, 1, G%NCE, NSCARC_INIT_ZERO, 'Y')
         IF (SV%Z /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%Z, 1, G%NCE, NSCARC_INIT_ZERO, 'Z')

#ifdef WITH_SCARC_DEBUG
         IF (SV%E /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%E, 1, G%NCE, NSCARC_INIT_ZERO, 'E')
#endif

         IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
            IF (SV%X_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%X_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'X_FB')
            IF (SV%B_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%B_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'B_FB')
            IF (SV%R_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%R_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'R_FB')
            IF (SV%V_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%V_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'V_FB')
         ENDIF

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_VECTORS




!> ----------------------------------------------------------------------------------------------------
!> Basic setup of environent on specified stack level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STACK(NSTACK)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN):: NSTACK

SV => STACK(NSTACK)%SOLVER

SV%TYPE_ACCURACY      = TYPE_ACCURACY 
SV%TYPE_COARSE        = TYPE_COARSE
SV%TYPE_COARSENING    = TYPE_COARSENING
SV%TYPE_CYCLING       = TYPE_CYCLING
SV%TYPE_EXCHANGE      = TYPE_EXCHANGE
SV%TYPE_GRID          = TYPE_GRID
SV%TYPE_INTERPOL      = TYPE_INTERPOL
SV%TYPE_KRYLOV        = TYPE_KRYLOV
SV%TYPE_LEVEL         = TYPE_LEVEL 
SV%TYPE_MATRIX        = TYPE_MATRIX
SV%TYPE_METHOD        = TYPE_METHOD
SV%TYPE_MKL           = TYPE_MKL
SV%TYPE_MKL_PRECISION = TYPE_MKL_PRECISION
SV%TYPE_MULTIGRID     = TYPE_MULTIGRID
SV%TYPE_PARENT        = TYPE_PARENT
SV%TYPE_PRECON        = TYPE_PRECON
SV%TYPE_RELAX         = TYPE_RELAX
SV%TYPE_SCOPE         = TYPE_SCOPE
SV%TYPE_SMOOTH        = TYPE_SMOOTH
SV%TYPE_SOLVER        = TYPE_SOLVER
SV%TYPE_STAGE         = TYPE_STAGE
SV%TYPE_STENCIL       = TYPE_STENCIL
SV%TYPE_TWOLEVEL      = TYPE_TWOLEVEL
SV%TYPE_VECTOR        = TYPE_VECTOR

END SUBROUTINE SCARC_SETUP_STACK


!> ----------------------------------------------------------------------------------------------------
!> Basic setup of environent on specified stack level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_STACK(NSTACK)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSTACK

SV => STACK(NSTACK)%SOLVER

TYPE_ACCURACY      = SV%TYPE_ACCURACY 
TYPE_COARSE        = SV%TYPE_COARSE
TYPE_COARSENING    = SV%TYPE_COARSENING
TYPE_CYCLING       = SV%TYPE_CYCLING
TYPE_EXCHANGE      = SV%TYPE_EXCHANGE
TYPE_GRID          = SV%TYPE_GRID
TYPE_INTERPOL      = SV%TYPE_INTERPOL
TYPE_KRYLOV        = SV%TYPE_KRYLOV
TYPE_LEVEL         = SV%TYPE_LEVEL 
TYPE_MATRIX        = SV%TYPE_MATRIX
TYPE_METHOD        = SV%TYPE_METHOD
TYPE_MKL           = SV%TYPE_MKL
TYPE_MKL_PRECISION = SV%TYPE_MKL_PRECISION
TYPE_MULTIGRID     = SV%TYPE_MULTIGRID
TYPE_PARENT        = SV%TYPE_PARENT
TYPE_PRECON        = SV%TYPE_PRECON
TYPE_RELAX         = SV%TYPE_RELAX
TYPE_SCOPE         = SV%TYPE_SCOPE
TYPE_SMOOTH        = SV%TYPE_SMOOTH
TYPE_SOLVER        = SV%TYPE_SOLVER
TYPE_STAGE         = SV%TYPE_STAGE
TYPE_STENCIL       = SV%TYPE_STENCIL
TYPE_TWOLEVEL      = SV%TYPE_TWOLEVEL
TYPE_VECTOR        = SV%TYPE_VECTOR

END SUBROUTINE SCARC_RESTORE_STACK


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_KRYLOV(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX

!>
!> Basic setup of stack information and types for Krylov method
!>
CALL SCARC_SETUP_STACK(NSTACK)

SV  => STACK(NSTACK)%SOLVER
SV%TYPE_METHOD   = NSCARC_METHOD_KRYLOV
SV%TYPE_SOLVER   = NSOLVER
SV%TYPE_SCOPE(0) = NSCOPE
SV%TYPE_STAGE    = NSTAGE
SV%TYPE_LEVEL(1) = NLMIN
SV%TYPE_LEVEL(2) = NLMAX
SV%TYPE_MATRIX   = TYPE_MATRIX

!>
!> Preset iteration parameters for Krylov method
!>
SELECT CASE(NSOLVER)

   !> -------------- Krylov method is used as main solver
   CASE (NSCARC_SOLVER_MAIN)
   
      SV%CNAME = 'SCARC_MAIN_KRYLOV'
   
      SV%EPS = SCARC_KRYLOV_ACCURACY
      SV%NIT = SCARC_KRYLOV_ITERATIONS
   
      SV%TYPE_RELAX    = TYPE_PRECON
      SV%TYPE_TWOLEVEL = TYPE_TWOLEVEL
   
   !> -------------- Krylov method is used as coarse grid solver solver
   CASE (NSCARC_SOLVER_COARSE)
   
      SV%CNAME = 'SCARC_COARSE_KRYLOV'
   
      SV%EPS = SCARC_COARSE_ACCURACY
      SV%NIT = SCARC_COARSE_ITERATIONS
   
      SV%TYPE_RELAX  = NSCARC_RELAX_SSOR               !> only use SSOR-preconditioning for coarse solver
      SV%TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE          !> only use one level for coarse solver
   
   !> -------------- Otherwise: print error message
   CASE DEFAULT
   
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)
   
END SELECT

!>
!> Point to solution vectors (in corresponding scope)
!>
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_KRYLOV


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Geometric Multigrid method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MULTIGRID(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX

!>
!> Basic setup of stack information and types for multigrid method
!>
CALL SCARC_SETUP_STACK(NSTACK)

SV  => STACK(NSTACK)%SOLVER
SV%TYPE_METHOD        = NSCARC_METHOD_MULTIGRID
SV%TYPE_SOLVER        = NSOLVER
SV%TYPE_SCOPE(0)      = NSCOPE
SV%TYPE_STAGE         = NSTAGE
SV%TYPE_LEVEL(0)      = NLMAX-NLMIN+1
SV%TYPE_LEVEL(1)      = NLMIN
SV%TYPE_LEVEL(2)      = NLMAX
SV%TYPE_RELAX         = TYPE_SMOOTH
SV%TYPE_GRID          = TYPE_GRID
SV%TYPE_INTERPOL      = TYPE_INTERPOL
SV%TYPE_ACCURACY      = TYPE_ACCURACY
SV%TYPE_CYCLING       = TYPE_CYCLING
SV%TYPE_MATRIX        = TYPE_MATRIX
SV%TYPE_MKL_PRECISION = TYPE_MKL_PRECISION

SELECT CASE(NSOLVER)

   !> -------------- Multigrid method is used as main solver
   CASE (NSCARC_SOLVER_MAIN)
      SV%CNAME = 'SCARC_MAIN_GMG'

   !> -------------- Multigrid method is only used as preconditioner for global Krylov-method
   CASE (NSCARC_SOLVER_PRECON)
      SV%CNAME = 'SCARC_PRECON_GMG'

   !> -------------- Otherwise: print error message
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)

END SELECT

!>
!> Preset iteration parameters for Multigrid method
!>
SV%EPS = SCARC_MULTIGRID_ACCURACY
SV%NIT = SCARC_MULTIGRID_ITERATIONS

!>
!> Point to solution vectors (in corresponding scope)
!>
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_MULTIGRID


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for MKL-methods
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSE_SOLVER(NSTAGE, NSCOPE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN)    :: NSCOPE, NSTAGE, NLMIN, NLMAX
INTEGER, INTENT(INOUT) :: NSTACK

SELECT_COARSE: SELECT CASE (TYPE_COARSE)

   !> -------------- CG-method is used as iterative coarse grid solver
   CASE (NSCARC_COARSE_ITERATIVE)

      !> initialize current stack position as CG-method
      STACK(NSTACK)%SOLVER => COARSE_KRYLOV
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)

      !> and next stack position as its SSOR-preconditioner
      NSTACK = NSTACK + 1
      TYPE_PRECON = NSCARC_RELAX_SSOR
      STACK(NSTACK)%SOLVER => PRECON_SSOR
      CALL SCARC_SETUP_PRECON(NSTACK, NSCOPE)

   !> -------------- LU-decomposition (from MKL) is used as direct coarse grid solver
#ifdef WITH_MKL 
   CASE (NSCARC_COARSE_DIRECT)

      !> Global scope in the multi-mesh case:
      !> initialize current stack position as global CLUSTER_SPARSE_SOLVER
      IF (NSCOPE == NSCARC_SCOPE_GLOBAL .AND. NMESHES > 1) THEN
         STACK(NSTACK)%SOLVER => COARSE_CLUSTER
         CALL SCARC_SETUP_MKL(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
         CALL SCARC_SETUP_CLUSTER(NLMIN, NLMAX)

      !> Local scope:
      !> initialize current stack position as PARDISO solver
      ELSE
         STACK(NSTACK)%SOLVER => COARSE_PARDISO
         CALL SCARC_SETUP_MKL(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
         CALL SCARC_SETUP_PARDISO(NLMIN, NLMAX)
      ENDIF
#endif

   !> -------------- Otherwise: print error message
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_COARSE)

END SELECT SELECT_COARSE
END SUBROUTINE SCARC_SETUP_COARSE_SOLVER


#ifdef WITH_MKL
!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for LU-solvers (based on MKL)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MKL(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX

!>
!> Basic setup of stack information and types for MKL
!>
CALL SCARC_SETUP_STACK(NSTACK)

SV  => STACK(NSTACK)%SOLVER
SELECT CASE (NSOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      SV%CNAME = 'SCARC_MAIN_LU'
   CASE (NSCARC_SOLVER_COARSE)
      SV%CNAME = 'SCARC_COARSE_LU'
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)
END SELECT

!>
!> Reset types for LU-decomposition method
!>
SV%TYPE_METHOD        = NSCARC_METHOD_LU
SV%TYPE_SOLVER        = NSOLVER
SV%TYPE_SCOPE(0)      = NSCOPE
SV%TYPE_STAGE         = NSTAGE
SV%TYPE_GRID          = TYPE_GRID
SV%TYPE_LEVEL(0)      = NLMAX-NLMIN+1
SV%TYPE_LEVEL(1)      = NLMIN
SV%TYPE_LEVEL(2)      = NLMAX
SV%TYPE_MATRIX        = NSCARC_MATRIX_COMPACT    
SV%TYPE_MKL_PRECISION = TYPE_MKL_PRECISION

!>
!> Point to solution vectors (in corresponding scope)
!>
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., NSTACK)

END SUBROUTINE SCARC_SETUP_MKL
#endif


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PRECON(NSTACK, NSCOPE)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NSTACK, NSCOPE

!>
!> Point to stack entry of called preconditioner and types for preconditioner
!>
SV => STACK(NSTACK)%SOLVER

SV%TYPE_SCOPE(0) = NSCOPE

SV%EPS   =  SCARC_PRECON_ACCURACY
SV%NIT   =  SCARC_PRECON_ITERATIONS
SV%OMEGA =  SCARC_PRECON_OMEGA

!>
!> Preset name of preconditioning method
!>
SELECT CASE(TYPE_PRECON)
   CASE (NSCARC_RELAX_JAC)
      SV%CNAME = 'SCARC_PRECON_JAC'
   CASE (NSCARC_RELAX_SSOR)
      SV%CNAME = 'SCARC_PRECON_SSOR'
   CASE (NSCARC_RELAX_MJAC)
      SV%CNAME = 'SCARC_PRECON_MJAC'
   CASE (NSCARC_RELAX_MGS)
      SV%CNAME = 'SCARC_PRECON_MGS'
   CASE (NSCARC_RELAX_MSGS)
      SV%CNAME = 'SCARC_PRECON_MSGS'
   CASE (NSCARC_RELAX_MSOR)
      SV%CNAME = 'SCARC_PRECON_MSOR'
   CASE (NSCARC_RELAX_MSSOR)
      SV%CNAME = 'SCARC_PRECON_MSSOR'
   CASE (NSCARC_RELAX_ILU)
      SV%CNAME = 'SCARC_PRECON_ILU'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_FFT)
      SV%CNAME = 'SCARC_PRECON_FFT'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_FFTO)
      SV%CNAME = 'SCARC_PRECON_FFTO'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_GMG)
      SV%CNAME = 'SCARC_PRECON_GMG'
#ifdef WITH_MKL
   CASE (NSCARC_RELAX_MKL)
      SV%OMEGA = 1.0_EB
      IF (NSCOPE == NSCARC_SCOPE_LOCAL) THEN
         SV%CNAME = 'SCARC_PRECON_PARDISO'
      ELSE
         SV%CNAME = 'SCARC_PRECON_CLUSTER'
         SV%NIT   = 1
      ENDIF
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_PRECON)
END SELECT

!>
!> Preset types for preconditioner (use same as for calling solver)
!>
SVP => STACK(NSTACK-1)%SOLVER

SV%TYPE_ACCURACY      = SVP%TYPE_ACCURACY
SV%TYPE_GRID          = SVP%TYPE_GRID
SV%TYPE_LEVEL(0:2)    = SVP%TYPE_LEVEL(0:2)
SV%TYPE_INTERPOL      = SVP%TYPE_INTERPOL
SV%TYPE_MATRIX        = SVP%TYPE_MATRIX
SV%TYPE_MKL_PRECISION = SVP%TYPE_MKL_PRECISION
SV%TYPE_RELAX         = SVP%TYPE_RELAX
SV%TYPE_SOLVER        = SVP%TYPE_SOLVER
SV%TYPE_STAGE         = SVP%TYPE_STAGE

!>
!> Preset pointers for preconditioner (use same as for alling solver)
!>
SV%X = SVP%X
SV%B = SVP%B
SV%D = SVP%D
SV%R = SVP%R
SV%V = SVP%V
SV%Y = SVP%Y
SV%Z = SVP%Z

IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
   SV%X_FB = SVP%X_FB
   SV%B_FB = SVP%B_FB
   SV%R_FB = SVP%R_FB
   SV%V_FB = SVP%V_FB
ENDIF

#ifdef WITH_SCARC_DEBUG
SV%E = SVP%E
#endif

END SUBROUTINE SCARC_SETUP_PRECON


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SMOOTH(NSTACK, NSCOPE)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NSTACK, NSCOPE

!>
!> Basic setup of stack information and types/names for smoother
!>
CALL SCARC_SETUP_STACK(NSTACK)

SV => STACK(NSTACK)%SOLVER
SV%TYPE_SCOPE(0) = NSCOPE
SV%OMEGA =  SCARC_SMOOTH_OMEGA

SELECT CASE(TYPE_SMOOTH)
   CASE (NSCARC_RELAX_JAC)
      SV%CNAME = 'SCARC_SMOOTH_JACOBI'
   CASE (NSCARC_RELAX_SSOR)
      SV%CNAME = 'SCARC_SMOOTH_SSOR'
   CASE (NSCARC_RELAX_MJAC)
      SV%CNAME = 'SCARC_SMOOTH_MJAC'
   CASE (NSCARC_RELAX_MGS)
      SV%CNAME = 'SCARC_SMOOTH_MGS'
   CASE (NSCARC_RELAX_MSGS)
      SV%CNAME = 'SCARC_SMOOTH_MSGS'
   CASE (NSCARC_RELAX_MSOR)
      SV%CNAME = 'SCARC_SMOOTH_MSOR'
   CASE (NSCARC_RELAX_MSSOR)
      SV%CNAME = 'SCARC_SMOOTH_MSSOR'
   CASE (NSCARC_RELAX_ILU)
      SV%CNAME = 'SCARC_SMOOTH_ILU'
   CASE (NSCARC_RELAX_FFT)
      SV%CNAME = 'SCARC_SMOOTH_FFT'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_FFTO)
      SV%CNAME = 'SCARC_SMOOTH_FFTO'
      SV%OMEGA = 1.0_EB
#ifdef WITH_MKL
   CASE (NSCARC_RELAX_MKL)
      SV%OMEGA = 1.0_EB
      IF (NSCOPE == NSCARC_SCOPE_GLOBAL) THEN
         SV%CNAME = 'SCARC_SMOOTH_CLUSTER'
      ELSE
         SV%CNAME = 'SCARC_SMOOTH_PARDISO'
      ENDIF
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_SMOOTH)
END SELECT

!>
!> Preset iteration parameters for Multigrid method
!>
SV%EPS = SCARC_SMOOTH_ACCURACY
SV%NIT = SCARC_SMOOTH_ITERATIONS

!>
!> Preset types for preconditioner (use same descriptors as calling solver)
!>
SVP => STACK(NSTACK-1)%SOLVER

SV%TYPE_SOLVER    = NSCARC_SOLVER_SMOOTH
SV%TYPE_STAGE     = SVP%TYPE_STAGE
SV%TYPE_LEVEL(0:2)= SVP%TYPE_LEVEL(0:2)
SV%TYPE_RELAX     = SVP%TYPE_RELAX
SV%TYPE_GRID      = SVP%TYPE_GRID
SV%TYPE_INTERPOL  = SVP%TYPE_INTERPOL
SV%TYPE_ACCURACY  = SVP%TYPE_ACCURACY
SV%TYPE_MKL_PRECISION = SVP%TYPE_MKL_PRECISION
SV%TYPE_MATRIX    = SVP%TYPE_MATRIX

!>
!>
!> Preset references for preconditioner (use same pointers as calling solver)
SV%X = SVP%X
SV%B = SVP%B
SV%D = SVP%D
SV%R = SVP%R
SV%V = SVP%V
SV%Y = SVP%Y
SV%Z = SVP%Z

IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
   SV%X_FB = SVP%X_FB
   SV%B_FB = SVP%B_FB
   SV%R_FB = SVP%R_FB
   SV%V_FB = SVP%V_FB
ENDIF

#ifdef WITH_SCARC_DEBUG
SV%E = SVP%E
SV%R = SVP%R
#endif

END SUBROUTINE SCARC_SETUP_SMOOTH


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for additive or multiplicative coarse grid
!> (corresponding to Schwarz domain decomposition method)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERPOLATION(NSTAGE, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, ST
INTEGER, INTENT(IN) :: NSTAGE, NLMIN, NLMAX
INTEGER :: NM, NL

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)
      ST => SCARC(NM)%LEVEL(NL)%STAGE(NSTAGE)

      CALL SCARC_ALLOCATE_REAL1(ST%X, 1, G%NCE, NSCARC_INIT_ZERO, 'X')
      CALL SCARC_ALLOCATE_REAL1(ST%B, 1, G%NCE, NSCARC_INIT_ZERO, 'B')
      CALL SCARC_ALLOCATE_REAL1(ST%V, 1, G%NCE, NSCARC_INIT_ZERO, 'Q')
      CALL SCARC_ALLOCATE_REAL1(ST%R, 1, G%NCE, NSCARC_INIT_ZERO, 'W')
      CALL SCARC_ALLOCATE_REAL1(ST%Y, 1, G%NCE, NSCARC_INIT_ZERO, 'Y')
      CALL SCARC_ALLOCATE_REAL1(ST%Z, 1, G%NCE, NSCARC_INIT_ZERO, 'Z')

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_INTERPOLATION


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for blockwise FFT methods
!> New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFT(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: M, S, L, FFT
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IERR = 0

!>
!> Allocate working space for FFT routine
!>
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)
      FFT => SCARC(NM)%LEVEL(NL)%FFT

      FFT%LBC = M%LBC
      FFT%MBC = M%MBC
      FFT%NBC = M%NBC

      FFT%XS = S%XS
      FFT%XF = S%XF
      FFT%YS = S%YS
      FFT%YF = S%YF
      FFT%ZS = S%ZS
      FFT%ZF = S%ZF

      FFT%IBAR = L%NX
      FFT%JBAR = L%NY
      FFT%KBAR = L%NZ

      FFT%ITRN = L%NX+1
      IF (TWO_D) THEN
         FFT%JTRN = 1
      ELSE
         FFT%JTRN = L%NY+1
      ENDIF
      FFT%KTRN = L%NZ+1

      FFT%LSAVE = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN+7*FFT%ITRN+5*FFT%JTRN+6*FFT%KTRN+56
      FFT%LWORK = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN

      CALL SCARC_ALLOCATE_REAL1(FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT')
      CALL SCARC_ALLOCATE_REAL1(FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT')

      !> Allocate stretching vector (set to 1)
      CALL SCARC_ALLOCATE_REAL1(FFT%HX, 0, FFT%ITRN, NSCARC_INIT_ONE, 'FFT')

      !> Allocate RHS vector for FFT routine
      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS')
      ELSE
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS')
      ENDIF

      !> Allocate boundary data vector for XS
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXS')
      ENDIF

      !> Allocate boundary data vector for XF
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXF')
      ENDIF

      !> Allocate boundary data vector for YS
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYS')
      ENDIF

      !> Allocate boundary data vector for YF
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYF')
      ENDIF

      !> Allocate boundary data vector for ZS
      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZS')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZS')

      !> Allocate boundary data vector for ZF
      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZF')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZF')

      IF (TWO_D) THEN
         CALL H2CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,IERR,FFT%SAVE1)
      ELSE
         CALL H3CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%YS,FFT%YF,FFT%JBAR,FFT%MBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERR,FFT%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FFT

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for blockwise FFT methods with overlap
!> New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFTO(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: M, S, L, FFT
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IERR = 0

!>
!> Allocate working space for FFT routine
!>
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)
      FFT => SCARC(NM)%LEVEL(NL)%FFT

      FFT%LBC = M%LBC
      FFT%MBC = M%MBC
      FFT%NBC = M%NBC

      IF (NM == 1) THEN
         FFT%IBAR = S%IBAR+1
         FFT%XS = S%XS
         FFT%XF = S%XF + L%DX
      ELSE IF (NM == NMESHES) THEN
         FFT%IBAR = S%IBAR+1
         FFT%XS = S%XS - L%DX
         FFT%XF = S%XF
      ELSE 
         FFT%IBAR = S%IBAR+2
         FFT%XS = S%XS - L%DX
         FFT%XF = S%XF + L%DX
      ENDIF

      FFT%JBAR = S%JBAR
      FFT%KBAR = S%KBAR

      FFT%YS = S%YS
      FFT%YF = S%YF
      FFT%ZS = S%ZS
      FFT%ZF = S%ZF

      FFT%ITRN = FFT%IBAR+1
      IF (TWO_D) THEN
         FFT%JTRN = 1
      ELSE
         FFT%JTRN = L%NY+1
      ENDIF
      FFT%KTRN = L%NZ+1

      FFT%LSAVE = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN+7*FFT%ITRN+5*FFT%JTRN+6*FFT%KTRN+56
      FFT%LWORK = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN

      CALL SCARC_ALLOCATE_REAL1(FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT')
      CALL SCARC_ALLOCATE_REAL1(FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT')

      !> Allocate stretching vector (set to 1)
      CALL SCARC_ALLOCATE_REAL1(FFT%HX, 0, FFT%ITRN, NSCARC_INIT_ONE, 'FFT')

      !> Allocate RHS vector for FFT routine
      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT')
      ELSE
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT')
      ENDIF

      !> Allocate boundary data vector for XS
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'BXS')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'BXS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'BXS')
      ENDIF

      !> Allocate boundary data vector for XF
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'BXF')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'BXF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'BXF')
      ENDIF

      !> Allocate boundary data vector for YS
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'BYS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'BYS')
      ENDIF

      !> Allocate boundary data vector for YF
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'BYF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'BYF')
      ENDIF

      !> Allocate boundary data vector for ZS
      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'BZS')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'BZS')

      !> Allocate boundary data vector for ZF
      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'BZF')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'BZF')

      IF (TWO_D) THEN
         CALL H2CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,IERR,FFT%SAVE1)
      ELSE
         CALL H3CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%YS,FFT%YF,FFT%JBAR,FFT%MBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERR,FFT%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FFTO


!> ----------------------------------------------------------------------------------------------------
!> Store JACOBI preconditioner in matrix form
!> Based on the following splitting of A = D - E - F
!> where :   D is the diagonal part
!> the MJAC-preconditioner in matrix form is defined
!>           M_MJAC = D 
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MJAC(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, AB, AC
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)

      SELECT CASE(TYPE_MATRIX)
         CASE (NSCARC_MATRIX_COMPACT)
            AC => G%AC
            CALL SCARC_ALLOCATE_REAL1(AC%RELAX, 1, G%NC, NSCARC_INIT_ZERO, 'MJAC')
            DO IC = 1, G%NC
               AC%RELAX(IC) = 1.0_EB/AC%VAL(AC%ROW(IC))
            ENDDO 
         CASE (NSCARC_MATRIX_BANDED)
            AB => G%AB
            CALL SCARC_ALLOCATE_REAL1(AB%RELAXD, 1, AB%N_DIAG,  NSCARC_INIT_ZERO, 'MJAC')
            DO IC = 1, G%NC
               AB%RELAXD(IC) = 1.0_EB/AB%VAL(IC, AB%POS(0))
            ENDDO 
       END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MJAC


!> ----------------------------------------------------------------------------------------------------
!> Store GS preconditioner in matrix form
!> Based on the following splitting of A = D - E - F
!> where :   D is the diagonal part
!>          -E is the strictly lower part
!>          -F is the strictly upper part
!> the SGS-preconditioner in matrix form is defined
!>           M_MGS = (D - E) = (I - E D^{-1}) D
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGS(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, AB, AC
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, IPTR

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)

      SELECT CASE (SCARC_MATRIX_LEVEL(NL))

         !>
         !> matrix in COMPACT storage technique
         !>
         CASE (NSCARC_MATRIX_COMPACT)

            AC => G%AC
            CALL SCARC_ALLOCATE_REAL1(AC%RELAX, 1, AC%N_VAL, NSCARC_INIT_ZERO, 'MGS')

            DO IC = 1, G%NC
               DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
                  JC = AC%COL(IPTR)
                  IF (JC <  IC) AC%RELAX(IPTR) = AC%VAL(IPTR) / AC%VAL(AC%ROW(JC))
                  IF (JC == IC) AC%RELAX(IPTR) = AC%VAL(IPTR)
               ENDDO
            ENDDO 

         !>
         !> matrix in BANDED storage technique
         !>
         CASE (NSCARC_MATRIX_BANDED)

            AB => G%AB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'MGS')
            WRITE(*,*) 'SCARC_SETUP_MGS: BANDED: NOT FINISHED YET'

      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MGS


!> ----------------------------------------------------------------------------------------------------
!> Store symmetric Gauss-Seidel preconditioner in matrix form
!> Based on the following splitting of A = D - E - F
!> where :   D is the diagonal part
!>          -E is the strictly lower part
!>          -F is the strictly upper part
!> the SGS-preconditioner in matrix form is defined
!>           M_MSGS = (D - E) D^{-1} (D - F)  =  (I - E D^{-1}) (D - F)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MSGS(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, AB, AC
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, IPTR, I, IS, IL, IOR0

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)

      SELECT CASE(TYPE_MATRIX)

         !>
         !> matrix in COMPACT storage technique
         !>
         CASE (NSCARC_MATRIX_COMPACT)

            AC => G%AC
            CALL SCARC_ALLOCATE_REAL1(AC%RELAX, 1, AC%N_VAL, NSCARC_INIT_ZERO, 'MSGS')
            AC%RELAX = AC%VAL

            DO IC = 1, G%NC
               !>  l(i,j) = a(i,j)/a(j,j)
               DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
                  JC = AC%COL(IPTR)
                  IF (JC < IC)  AC%RELAX(IPTR) = AC%VAL(IPTR) / AC%VAL(AC%ROW(JC))
               ENDDO
            ENDDO 

         !>
         !> matrix in BANDED storage technique
         !>
         CASE (NSCARC_MATRIX_BANDED)

            AB => G%AB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'MSGS')
            AB%RELAX = AB%VAL
            DO IOR0 = 3, 1, -1
               IS = AB%TARGET(IOR0)
               IL = AB%LENGTH(IOR0)
               DO I = 1, AB%LENGTH(IOR0)
                  AB%RELAX(IS+I-1, AB%POS(IOR0)) = AB%VAL(IS+I-1, AB%POS(IOR0))/AB%VAL(I,AB%POS(0))
               ENDDO
            ENDDO

      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP
      
END SUBROUTINE SCARC_SETUP_MSGS


!> ----------------------------------------------------------------------------------------------------
!> Store SOR preconditioner in matrix form
!> Based on the following splitting of A = D - E - F
!> where :   D is the diagonal part
!>          -E is the strictly lower part
!>          -F is the strictly upper part
!> the SOR-preconditioner in matrix form is defined
!>           M_MSOR = (D−ωE) = (I−ωE D^{-1}) D
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MSOR(NLMIN, NLMAX, NSTACK)
USE SCARC_POINTERS, ONLY: G, AB, AC
INTEGER, INTENT(IN) :: NLMIN, NLMAX, NSTACK
REAL (EB) :: OMEGA
INTEGER :: NM, NL, IC, JC, IPTR

OMEGA = STACK(NSTACK)%SOLVER%OMEGA

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)

      SELECT CASE(TYPE_MATRIX)

         !>
         !> matrix in COMPACT storage technique
         !>
         CASE (NSCARC_MATRIX_COMPACT)

            AC => G%AC
            CALL SCARC_ALLOCATE_REAL1(AC%RELAX, 1, AC%N_VAL, NSCARC_INIT_ZERO, 'MSOR')
      
            DO IC = 1, G%NC
      
               !DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
               !   JC = AC%COL(IPTR)
               !   IF (JC <  IC) AC%RELAX(IPTR) = OMEGA * AC%VAL(IPTR) / AC%VAL(AC%ROW(JC))
               !   IF (JC == IC) AC%RELAX(IPTR) = AC%VAL(IPTR)
               !ENDDO
               DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
                  JC = AC%COL(IPTR)
                  IF (JC < IC) THEN
                     AC%RELAX(IPTR) = OMEGA * AC%VAL(IPTR) / AC%VAL(AC%ROW(JC))
                  ELSE IF (JC == IC) THEN
                     AC%RELAX(IPTR) = AC%VAL(IPTR)
                  ELSE IF (JC > IC) THEN
                     AC%RELAX(IPTR) = OMEGA * AC%VAL(IPTR) 
                  ENDIF
               ENDDO
      
            ENDDO 

         !>
         !> matrix in BANDED storage technique
         !>
         CASE (NSCARC_MATRIX_BANDED)

            AB => G%AB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'MSOR')
            WRITE(*,*) 'SCARC_SETUP_MSOR: BANDED: NOT FINISHED YET'
      
      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MSOR


!> ----------------------------------------------------------------------------------------------------
!> Store SSOR preconditioner in matrix form
!> Based on the following splitting of A = D - E - F
!> where :   D is the diagonal part
!>          -E is the strictly lower part
!>          -F is the strictly upper part
!> the SSOR-preconditioner in matrix form is defined
!>           B_SSOR = 1/(omega * (2-omega)) * (D - omega * E) D^{-1} (D - omega * F)
!>                  = (I - omega E D^{-1}) * [1/(omega * (2-omega)) * D -  1/(2-omega) * F]
!> Defining the triangular matrices
!>               L  = I - omega E D^{-1}
!>               U  = 1/(omega * (2-omega)) * D -  1/(2-omega) * F
!> the SSOR-preconditioning can be thought as the solution of two triangular systems
!> Both matrices can be stored as a single matrix that occupies the same amount of storage as A
!> where the same row and column pointers can be used as for A (identical pattern)
!> Note that the diagonal elements of L are 1 (and are omitted, only the diagonal of U is stored there)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MSSOR(NLMIN, NLMAX, NSTACK)
USE SCARC_POINTERS, ONLY: G, AB, AC
INTEGER, INTENT(IN) :: NLMIN, NLMAX, NSTACK
INTEGER :: NM, NL, IC, JC, IPTR, INCR, IS, IL, IOR0, I
REAL(EB) :: OMEGA, SCAL1, SCAL2

OMEGA = STACK(NSTACK)%SOLVER%OMEGA
SCAL1  = 1.0_EB / (OMEGA * (2.0_EB - OMEGA))
SCAL2  = 1.0_EB / (2.0_EB - OMEGA)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)

      SELECT CASE(TYPE_MATRIX)

         !>
         !> matrix in COMPACT storage technique
         !>
         CASE (NSCARC_MATRIX_COMPACT)

            AC => G%AC
            CALL SCARC_ALLOCATE_REAL1(AC%RELAX, 1, AC%N_VAL, NSCARC_INIT_ZERO, 'MSSOR')
      
            DO IC = 1, G%NC
      
               DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
                  JC = AC%COL(IPTR)
                  IF (JC < IC) THEN
                     AC%RELAX(IPTR) = OMEGA * AC%VAL(IPTR) / AC%VAL(AC%ROW(JC))
                  ELSE IF (JC == IC) THEN
                     AC%RELAX(IPTR) = SCAL1 * AC%VAL(IPTR)
                  ELSE IF (JC > IC) THEN
                     AC%RELAX(IPTR) = SCAL2 * AC%VAL(IPTR)
                  ENDIF
               ENDDO

            ENDDO 

         !>
         !> matrix in BANDED storage technique
         !>
         CASE (NSCARC_MATRIX_BANDED)

            AB => G%AB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'MSSOR')

            IF (TWO_D) THEN
               INCR = -2
            ELSE
               INCR = -1
            ENDIF
            
            DO IOR0 = 3, 1, INCR
               IS = AB%TARGET(IOR0)
               IL = AB%LENGTH(IOR0)
               DO I = 1, AB%LENGTH(IOR0)
                  AB%RELAX(IS+I-1, AB%POS(IOR0)) = OMEGA * AB%VAL(IS+I-1, AB%POS(IOR0))/AB%VAL(I,AB%POS(0))
               ENDDO
            ENDDO
            DO I = 1, AB%LENGTH(0)
               AB%RELAX(I, AB%POS(0)) = SCAL1 * AB%VAL(I, AB%POS(0))
            ENDDO
            DO IOR0 = -1, -3, INCR
               IS = AB%TARGET(IOR0)
               IL = AB%LENGTH(IOR0)
               DO I = IS, IS + AB%LENGTH(IOR0) -1
                  AB%RELAX(I, AB%POS(IOR0)) = SCAL2 * AB%VAL(I, AB%POS(IOR0))
               ENDDO
            ENDDO

!>               L  = I - omega E D^{-1}
!>               U  = 1/(omega * (2-omega)) * D -  1/(2-omega) * F

      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

!STACK(NSTACK)%SOLVER%OMEGA = 1.0_EB

END SUBROUTINE SCARC_SETUP_MSSOR


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize ILU(0) decomposition of Poisson matrix
!> L- and U-parts are stored in the same array, diagonal elements of L are supposed to be 1
!> Based on Saad-algorithm 10.4 from 'Iterative Methods for Sparse Linear Systems':
!>   for i = 2 , ... , n do
!>      for k = 1 , ... , i-1 and for (i,k) in NZ(A) do
!>         compute a_ik = a_ik / a_kk
!>         for j = k+1 , ... , n and for (i,j) in NZ(A) do
!>            compute a_ij = a_ij - a_ik a_kj
!>         enddo
!>      enddo
!>   enddo
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_ILU(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, AB, AC
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, KC, IPTR, JPTR, KPTR, KPTR0, IOR0, JOR0, KOR0
LOGICAL :: BFOUND


MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)

      SELECT CASE(TYPE_MATRIX)

         !>
         !> matrix in COMPACT storage technique
         !>
         CASE (NSCARC_MATRIX_COMPACT)

            AC => G%AC
            CALL SCARC_ALLOCATE_REAL1(AC%RELAX, 1, AC%N_VAL, NSCARC_INIT_ZERO, 'ILU')
            AC%RELAX = AC%VAL
      
            CELL_LOOP: DO IC = 2, G%NC
      
               COLUMN_LOOP: DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
      
                  KC = AC%COL(IPTR)                        !> get number of neighboring cell
                  IF (KC >= IC) CYCLE                      !> only consider neighbors with lower cell numbers than IC
                  IF (AC%RELAX(IPTR) == 0) CYCLE
      
                  KPTR = AC%ROW(KC)                        !> get diagonal entry of neighbor
                  AC%RELAX(IPTR) = AC%RELAX(IPTR)/AC%RELAX(KPTR)
      
                  DO JPTR = AC%ROW(IC), AC%ROW(IC+1)-1
      
                     JC = AC%COL(JPTR)
                     IF (JC<=KC) CYCLE                     !> only consider neighbors with higher cell numbers than IC
                     IF (AC%RELAX(JPTR) == 0) CYCLE
      
                     KPTR = -1
                     DO KPTR0 = AC%ROW(KC), AC%ROW(KC+1)-1
                        IF (AC%COL(KPTR0) == JC) THEN
                          KPTR = KPTR0
                        ENDIF
                     ENDDO
                     IF (KPTR>0) AC%RELAX(JPTR) = AC%RELAX(JPTR) - AC%RELAX(IPTR) * AC%RELAX(KPTR)
      
                  ENDDO
      
               ENDDO COLUMN_LOOP
            ENDDO CELL_LOOP
      
         !>
         !> matrix in BANDED storage technique
         !>
         CASE (NSCARC_MATRIX_BANDED)

            AB => G%AB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'ILU')
            AB%RELAX = AB%VAL
      
            CELL_BANDED_LOOP: DO IC = 2, G%NC
      
               COLUMN_BANDED_LOOP: DO IOR0 = 3, -3, -1 
      
                  IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE

                  KC = IC + AB%OFFSET(IOR0)                           !> get number of neighboring cell
                  IF (KC<=0  .OR. KC >= IC) CYCLE                     !> only consider neighbors with lower cell numbers than IC
                  IF (AB%RELAX(IC, AB%POS(IOR0)) == 0) CYCLE
      
                  AB%RELAX(IC,AB%POS(IOR0)) = AB%RELAX(IC, AB%POS(IOR0))/AB%RELAX(KC,AB%POS(0))
      
                  DO JOR0 = 3, -3, -1

                     IF (TWO_D .AND. ABS(JOR0) == 2) CYCLE
      
                     JC = IC + AB%OFFSET(JOR0)                       !> get number of neighboring cell

                     IF (JC<=KC .OR. JC >G%NC) CYCLE                 !> only consider neighbors with higher cell numbers than IC
                     IF (AB%RELAX(IC, AB%POS(JOR0)) == 0) CYCLE
      
                     BFOUND = .FALSE.
                     DO KOR0 = 3, -3, -1
                        IF (TWO_D .AND. ABS(KOR0) == 2) CYCLE
                        IF (KC + AB%OFFSET(KOR0) == JC) THEN
                           BFOUND = .TRUE.
                           EXIT
                        ENDIF
                     ENDDO
                     IF (BFOUND) AB%RELAX(JC, AB%POS(JOR0)) = AB%RELAX(JC, AB%POS(JOR0)) - &
                                 AB%RELAX(IC, AB%POS(IOR0)) * AB%RELAX(KC, AB%POS(KOR0))
      
                  ENDDO
      
               ENDDO COLUMN_BANDED_LOOP
            ENDDO CELL_BANDED_LOOP
      
      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_ILU


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CLUSTER(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: L, G, MKL, ACS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I 
REAL (EB) :: TNOW
REAL (EB) :: DUMMY(1)=0.0_EB
REAL (FB) :: DUMMY_FB(1)=0.0_FB

TNOW = CURRENT_TIME()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)

      MKL => L%MKL
      ACS => G%ACS

      !>
      !> Allocate workspace for parameters and pointers needed in MKL-routine
      !>
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL_IPARM')

      IF (.NOT.ALLOCATED(MKL%CT)) THEN
         ALLOCATE(MKL%CT(64), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'CT', IERROR)
         DO I=1,64
            MKL%CT(I)%DUMMY = 0
         ENDDO
      ENDIF

      !>
      !> Define corresponding parameters
      !> Note: IPARM-vectory is allocate from 1:64, not from 0:63
      !>
      MKL%NRHS   =  1         ! one right hand side
      MKL%MAXFCT =  1         ! one matrix
      MKL%MNUM   =  1         ! number of matrix to be factorized
      MKL%ERROR  =  0         ! initialize error flag
      MKL%MSGLVL =  0         ! do not print statistical information
      IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
         MKL%MTYPE  = -2      ! Matrix type real and symmetric indefinite
      ELSE
         MKL%MTYPE  = 11      ! Matrix type real and non-symmetric
      ENDIF
      MKL%IPARM(1)  =  1      ! no solver default
      IF (N_MPI_PROCESSES > 4) THEN 
         MKL%IPARM(2) =10     ! 10 = MPI Parallel fill-in reordering from METIS. If 3 = OpenMP parallel reordering in Master
      ELSE                    ! Note IPARM(2)=10 has a bug which has been fixed from Intel MKL 2018 update 2 onwards.
         MKL%IPARM(2) = 3
      ENDIF
      MKL%IPARM(4)  = 0       ! no iterative-direct algorithm
      MKL%IPARM(5)  = 0       ! no user fill-in reducing permutation
      MKL%IPARM(6)  = 0       ! =0 solution on the first n components of x
      MKL%IPARM(8)  = 2       ! numbers of iterative refinement steps
      MKL%IPARM(10) = 13      ! perturb the pivot elements with 1E-13
      MKL%IPARM(11) = 1       ! use nonsymmetric permutation and scaling MPS  !!!!! was 1
      MKL%IPARM(13) = 1       ! maximum weighted matching algorithm is switched-off
                              !(default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
      MKL%IPARM(14) = 0       ! Output: number of perturbed pivots
      MKL%IPARM(18) = 0       ! Output: number of nonzeros in the factor LU
      MKL%IPARM(19) = 0       ! Output: Mflops for LU factorization
      MKL%IPARM(20) = 0       ! Output: Numbers of CG Iterations
      MKL%IPARM(21) = 1       ! 1x1 diagonal pivoting for symmetric indefinite matrices.
      MKL%IPARM(24) = 0
      MKL%IPARM(27) = 1       ! Check matrix
      MKL%IPARM(40) = 2       ! Matrix, solution and rhs provided in distributed assembled matrix input format.


      MKL%IPARM(41) = G%NC_OFFSET(NM) + 1                      ! first global cell number for mesh NM
      MKL%IPARM(42) = G%NC_OFFSET(NM) + G%NC_LOCAL(NM)         ! last global cell number for mesh NM
      !MKL%IPARM(39) = 2                                       ! provide matrix in distributed format
      !MKL%IPARM(40) = G%NC_OFFSET(NM)+1                       ! first global cell number for mesh NM
      !MKL%IPARM(41) = G%NC_OFFSET(NM)+G%NC_LOCAL(NM)         ! last global cell number for mesh NM

      !> 
      !> First perform only reordering and symbolic factorization
      !> Then perform only factorization
      !> 
      IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

         MKL%IPARM(28) = 1         ! single precision
         MKL%PHASE = 11
         CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      ACS%VAL_FB, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MPI_COMM_WORLD, MKL%ERROR)
         MKL%PHASE = 22
         CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      ACS%VAL_FB, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MPI_COMM_WORLD, MKL%ERROR)

      ELSE

         MKL%IPARM(28) = 0         ! double precision
         MKL%PHASE = 11
         CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      ACS%VAL, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
         MKL%PHASE = 22
         CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      ACS%VAL, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)

      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_CLUSTER


!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PARDISO(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: L, G, MKL, ACS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I, IDUMMY(1)=0
REAL (EB) :: TNOW
REAL (EB) :: DUMMY(1)=0.0_EB
REAL (FB) :: DUMMY_FB(1)=0.0_FB

TNOW = CURRENT_TIME()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID(NM, NL)
      MKL => L%MKL
      ACS => G%ACS

      !>
      !> Allocate workspace for parameters nnd pointers eeded in MKL-routine
      !>
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL_IPARM')

      IF (.NOT.ALLOCATED(MKL%PT)) THEN
         ALLOCATE(MKL%PT(64), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'PT', IERROR)
         DO I=1,64
            MKL%PT(I)%DUMMY = 0
         ENDDO
      ENDIF

      !>
      !> Define corresponding parameters
      !> Note: IPARM-vectory is allocate from 1:64, not from 0:63
      !>
      MKL%NRHS   = 1
      MKL%MAXFCT = 1
      MKL%MNUM   = 1

      MKL%IPARM(1)  =  1      ! no solver default
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
      MKL%IPARM(37) =  0      ! matrix storage in COMPACT-format
      MKL%IPARM(40) = 2       ! Matrix, solution and rhs provided in distributed assembled matrix input format.

      MKL%ERROR  =  0         ! initialize error flag
      MKL%MSGLVL =  0         ! do not print statistical information
      MKL%MTYPE  = -2         ! Matrix type real non-symmetric

      !> 
      !> First perform only reordering and symbolic factorization
      !> Then perform only factorization
      !> 
      IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

         MKL%IPARM(28) = 1         ! single precision
         MKL%PHASE = 11
         CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        ACS%VAL_FB, ACS%ROW, ACS%COL, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)
         MKL%PHASE = 22
         CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        ACS%VAL_FB, ACS%ROW, ACS%COL, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)

      ELSE

         MKL%IPARM(28) = 0         ! double precision
         MKL%PHASE = 11
         CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        ACS%VAL, ACS%ROW, ACS%COL, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
         MKL%PHASE = 22
         CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        ACS%VAL, ACS%ROW, ACS%COL, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_PARDISO
#endif

!> ------------------------------------------------------------------------------------------------
!> Define sizes for system matrix A (including extended regions related to overlaps)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_SIZES(NTYPE, NL)
USE SCARC_POINTERS, ONLY: L, G, OG, AB, AC
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM, IW

SELECT_MATRIX_SIZE: SELECT CASE (NTYPE)
   
   CASE (NSCARC_SIZE_MATRIX)
   
      MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_GRID(NM, NL)
   
         SELECT_MATRIX_TYPE: SELECT CASE (SCARC_MATRIX_LEVEL(NL))
   
            !>
            !> -------- Matrix in compact storage technique
            !>
            CASE (NSCARC_MATRIX_COMPACT)
      
               AC => G%AC
               IF (TWO_D) THEN
                  AC%N_STENCIL = 5
                  AC%POS(-3:3) = (/1,0,2,3,4,0,5/)     !> assignment of IOR settings to position in stencil
               ELSE
                  AC%N_STENCIL = 7
                  AC%POS(-3:3) = (/1,2,3,4,5,6,7/)
               ENDIF
      
               AC%N_VAL  = G%NC * AC%N_STENCIL
               AC%N_ROW  = G%NC + 1
      
               !> Determine sizes of overlapped parts for later communication with corresponding neighbors
               DO IW = 1, L%N_WALL_CELLS_EXT
                  NOM = G%WALL(IW)%NOM
                  IF (NOM /= 0) THEN
                     CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)
                     OG%AC%N_VAL = OG%AC%N_VAL + AC%N_STENCIL
                  ENDIF
               ENDDO
      
            !>
            !> -------- Matrix in banded storage technique
            !>
            CASE (NSCARC_MATRIX_BANDED)
      
               AB => G%AB
      
               IF (TWO_D) THEN
   
                  AB%N_STENCIL   = 5                      !> 5-point Laplacian
                  AB%POS(-3:3)   = (/5,0,4,3,2,0,1/)      !> assignment of IOR settings to columns in matrix array
   
                  AB%OFFSET( 3)  = -L%NX                  !> lower z
                  AB%OFFSET( 1)  = -1                     !> lower x
                  AB%OFFSET( 0)  =  0                     !> diag
                  AB%OFFSET(-1)  =  1                     !> upper x
                  AB%OFFSET(-3)  =  L%NX                  !> upper z
   
                  AB%SOURCE( 3)   =  1                    !> lower z
                  AB%SOURCE( 1)   =  1                    !> lower x
                  AB%SOURCE( 0)   =  1                    !> diag
                  AB%SOURCE(-1)   =  2                    !> upper x
                  AB%SOURCE(-3)   =  L%NX+1               !> upper z
   
                  AB%TARGET( 3)   =  L%NX+1               !> lower z
                  AB%TARGET( 1)   =  2                    !> lower x
                  AB%TARGET( 0)   =  1                    !> diag
                  AB%TARGET(-1)   =  1                    !> upper x
                  AB%TARGET(-3)   =  1                    !> upper z
   
                  AB%LENGTH( 3)  =  G%NC - L%NX           !> lower z
                  AB%LENGTH( 1)  =  G%NC - 1              !> lower x
                  AB%LENGTH( 0)  =  G%NC                  !> diag
                  AB%LENGTH(-1)  =  G%NC - 1              !> upper x
                  AB%LENGTH(-3)  =  G%NC - L%NX           !> upper z
   
               ELSE
   
                  AB%N_STENCIL   = 7                      !> 7-point Laplacian
                  AB%POS(-3:3)   = (/7,6,5,4,3,2,1/)      !> assignment of IOR settings to columns in matrix array
   
                  AB%OFFSET( 3)  = -L%NX*L%NY             !> lower z
                  AB%OFFSET( 2)  = -L%NX                  !> lower y
                  AB%OFFSET( 1)  = -1                     !> lower x
                  AB%OFFSET( 0)  =  0                     !> diag
                  AB%OFFSET(-1)  =  1                     !> upper x
                  AB%OFFSET(-2)  =  L%NX                  !> upper y
                  AB%OFFSET(-3)  =  L%NX*L%NY             !> upper z
   
                  AB%SOURCE( 3)  =  1                     !> lower z
                  AB%SOURCE( 2)  =  1                     !> lower y
                  AB%SOURCE( 1)  =  1                     !> lower x
                  AB%SOURCE( 0)  =  1                     !> diag
                  AB%SOURCE(-1)  =  2                     !> upper x
                  AB%SOURCE(-2)  =  L%NX+1                !> upper y
                  AB%SOURCE(-3)  =  L%NX*L%NY+1           !> upper z
   
                  AB%TARGET( 3)  =  L%NX*L%NY+1           !> lower z
                  AB%TARGET( 2)  =  L%NX+1                !> lower y
                  AB%TARGET( 1)  =  2                     !> lower x
                  AB%TARGET( 0)  =  1                     !> diag
                  AB%TARGET(-1)  =  1                     !> upper x
                  AB%TARGET(-2)  =  1                     !> upper y
                  AB%TARGET(-3)  =  1                     !> upper z
   
                  AB%LENGTH( 3)  =  G%NC - L%NX*L%NY      !> lower z
                  AB%LENGTH( 2)  =  G%NC - L%NX           !> lower y
                  AB%LENGTH( 1)  =  G%NC - 1              !> lower x
                  AB%LENGTH( 0)  =  G%NC                  !> diag
                  AB%LENGTH(-1)  =  G%NC - 1              !> upper x
                  AB%LENGTH(-2)  =  G%NC - L%NX           !> upper y
                  AB%LENGTH(-3)  =  G%NC - L%NX*L%NY      !> upper z
   
               ENDIF
      
               AB%N_VAL  = G%NC * AB%N_STENCIL
               AB%N_DIAG = G%NC
      
               !> Determine sizes of overlapping parts for later communication with corresponding neighbors
               DO IW = 1, L%N_WALL_CELLS_EXT
                  NOM = G%WALL(IW)%NOM
                  IF (NOM /= 0) THEN
                     CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)
                     OG%AB%N_VAL     = OG%AB%N_VAL + AB%N_STENCIL
                     OG%AB%N_DIAG    = OG%AB%N_DIAG + 1
                     OG%AB%N_STENCIL = AB%N_STENCIL
                  ENDIF
               ENDDO
      
         END SELECT SELECT_MATRIX_TYPE
   
      ENDDO MESHES_LOOP1
   
      !>
      !> -------- Exchange matrix sizes in case of a multi-mesh geometry
      !>
      IF (NMESHES > 1) THEN

         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SIZE, NL)
   
         IF (SCARC_MATRIX_LEVEL(NL) == NSCARC_MATRIX_BANDED) THEN
            DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
               CALL SCARC_POINT_TO_GRID(NM, NL)
               DO NOM = 1, NMESHES
                  IF (.NOT. ARE_NEIGHBORS(NM, NOM)) CYCLE
                  CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)
                  CALL SCARC_ALLOCATE_MATRIX_BANDED(OG%AB, NL, NSCARC_INIT_ZERO, 'OG%AB')
               ENDDO
            ENDDO
         ENDIF

      ENDIF

      !> Temporarily disabled, corresponding routines stored separatedly
      !>   CASE (NSCARC_SIZE_PROLONGATION)
      !>   CASE (NSCARC_SIZE_RESTRICTION)

END SELECT SELECT_MATRIX_SIZE

END SUBROUTINE SCARC_SETUP_MATRIX_SIZES


!> ------------------------------------------------------------------------------------------------------
!> Usual ScaRC solver routine
!> ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SOLVER(DT_CURRENT)
REAL (EB), INTENT(IN) :: DT_CURRENT
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

DT  = DT_CURRENT
DTI = 1.0_EB/DT_CURRENT

ITE_PRES = ITE_PRES + 1
ITE_GLOBAL = ICYC

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !>
   !> ---------------- Krylov method (CG/BICG) --------------------------------
   !>
   CASE (NSCARC_METHOD_KRYLOV)
   
      SELECT_KRYLOV: SELECT CASE (TYPE_KRYLOV)
         CASE (NSCARC_KRYLOV_CG)
            CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
#ifdef WITH_SCARC_BAROCLINIC
         CASE (NSCARC_KRYLOV_CGBARO)
            CALL SCARC_METHOD_BAROCLINIC (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
#endif
         CASE (NSCARC_KRYLOV_BICG)
            CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)  ! temporarily
            !CALL SCARC_SHUTDOWN(NSCARC_ERROR_BICG_DISABLED, SCARC_KRYLOV, NSCARC_NONE)
            !            CALL SCARC_METHOD_BICG(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
      END SELECT SELECT_KRYLOV
   
   !>
   !> ---------------- Multigrid method ---------------------------------------
   !>
   CASE (NSCARC_METHOD_MULTIGRID)
   
      CALL SCARC_METHOD_MULTIGRID(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   

#ifdef WITH_SCARC_MGM
   !>
   !> ---------------- McKenny-Greengard-Mayo method (MGM) --------------------
   !>
   CASE (NSCARC_METHOD_MGM)
   
      !> first solve inhomogeneous Poisson problem on structured grid with ScaRC (with Block-FFT)
      CALL SCARC_ASSIGN_GRID_TYPE(NSCARC_GRID_STRUCTURED)

      CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   
      CALL SCARC_STORE_MGM(NLEVEL_MIN, 1)
   
      CALL SCARC_MGM_INTERNAL_VELOCITY(NLEVEL_MIN)
      CALL SCARC_STORE_MGM(NLEVEL_MIN, 1)
   
      !> then solve homogeneous Poisson problem on unstructured grid with UScaRC (first with SSOR-preconditioning)
      !> later the preconditioning will be replaced by an individual LU-process
      CALL SCARC_ASSIGN_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)
      CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_HOMOGENEOUS, NLEVEL_MIN)

      CALL SCARC_STORE_MGM(NLEVEL_MIN, 2)
      CALL SCARC_STORE_MGM(NLEVEL_MIN, 3)

#endif
   
   !>
   !> ---------------- MKL method ---------------------------------------------
   !>
#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LU)
   
      SELECT_MKL: SELECT CASE (TYPE_MKL(0))
         CASE (NSCARC_MKL_GLOBAL)
            CALL SCARC_METHOD_CLUSTER(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
         CASE (NSCARC_MKL_LOCAL)
            CALL SCARC_METHOD_PARDISO(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
      END SELECT SELECT_MKL
#endif
   
END SELECT SELECT_METHOD

IF (STOP_STATUS==SETUP_STOP) RETURN


T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
CPU(MYID)%SOLVER =CPU(MYID)%SOLVER+CURRENT_TIME()-TNOW
CPU(MYID)%OVERALL=CPU(MYID)%OVERALL+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SOLVER


!> ------------------------------------------------------------------------------------------------------
!> Perform preceding FFT method to improve start solution for ScaRC
!> ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_FFT
USE MESH_POINTERS
USE POIS, ONLY: H2CZSS, H3CZSS
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
INTEGER :: NM, I, J, K
LOGICAL :: WITH_BDRY = .FALSE.

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)
   
   IF (PREDICTOR) THEN
      HP => H
   ELSE
      HP => HS
   ENDIF
   
   !>
   !> Call the Poisson solver
   !>
   IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
   IF (TWO_D .AND. .NOT.CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
   
   DO K=1,KBAR
     DO J=1,JBAR
        DO I=1,IBAR
            HP(I,J,K) = PRHS(I,J,K)
        ENDDO
      ENDDO
   ENDDO
   
   !>
   !> Apply boundary conditions to H
   !>
   IF (WITH_BDRY) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            IF (LBC==3 .OR. LBC==4)             HP(0,J,K)    = HP(1,J,K)    - DXI*BXS(J,K)
            IF (LBC==3 .OR. LBC==2 .OR. LBC==6) HP(IBP1,J,K) = HP(IBAR,J,K) + DXI*BXF(J,K)
            IF (LBC==1 .OR. LBC==2)             HP(0,J,K)    =-HP(1,J,K)    + 2._EB*BXS(J,K)
            IF (LBC==1 .OR. LBC==4 .OR. LBC==5) HP(IBP1,J,K) =-HP(IBAR,J,K) + 2._EB*BXF(J,K)
            IF (LBC==5 .OR. LBC==6)             HP(0,J,K)    = HP(1,J,K)
            IF (LBC==0) THEN
               HP(0,J,K) = HP(IBAR,J,K)
               HP(IBP1,J,K) = HP(1,J,K)
            ENDIF
         ENDDO
      ENDDO
      
      DO K=1,KBAR
         DO I=1,IBAR
            IF (MBC==3 .OR. MBC==4) HP(I,0,K)    = HP(I,1,K)    - DETA*BYS(I,K)
            IF (MBC==3 .OR. MBC==2) HP(I,JBP1,K) = HP(I,JBAR,K) + DETA*BYF(I,K)
            IF (MBC==1 .OR. MBC==2) HP(I,0,K)    =-HP(I,1,K)    + 2._EB*BYS(I,K)
            IF (MBC==1 .OR. MBC==4) HP(I,JBP1,K) =-HP(I,JBAR,K) + 2._EB*BYF(I,K)
            IF (MBC==0) THEN
               HP(I,0,K) = HP(I,JBAR,K)
               HP(I,JBP1,K) = HP(I,1,K)
            ENDIF
         ENDDO
      ENDDO
      
      DO J=1,JBAR
         DO I=1,IBAR
            IF (NBC==3 .OR. NBC==4)  HP(I,J,0)    = HP(I,J,1)    - DZETA*BZS(I,J)
            IF (NBC==3 .OR. NBC==2)  HP(I,J,KBP1) = HP(I,J,KBAR) + DZETA*BZF(I,J)
            IF (NBC==1 .OR. NBC==2)  HP(I,J,0)    =-HP(I,J,1)    + 2._EB*BZS(I,J)
            IF (NBC==1 .OR. NBC==4)  HP(I,J,KBP1) =-HP(I,J,KBAR) + 2._EB*BZF(I,J)
            IF (NBC==0) THEN
               HP(I,J,0) = HP(I,J,KBAR)
               HP(I,J,KBP1) = HP(I,J,1)
            ENDIF
         ENDDO
      ENDDO
   ENDIF

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_METHOD_FFT


!> ------------------------------------------------------------------------------------------------------
!> Check of the accuracy of the separable pressure solution, del^2 H = -del dot F - dD/dt
!> ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_POISSON
USE MESH_POINTERS
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP, RESIDUAL
REAL(EB) :: LHSS, RHSS
INTEGER :: NM, I, J, K

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)
   
   IF (PREDICTOR) THEN
      HP => H
   ELSE
      HP => HS
   ENDIF
   RESIDUAL => WORK8(1:IBAR,1:JBAR,1:KBAR)
   
   !$OMP PARALLEL DO PRIVATE(I,J,K,RHSS,LHSS) SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHSS = ( R(I-1)*FVX(I-1,J,K) - R(I)*FVX(I,J,K) )*RDX(I)*RRN(I) &
                 + (        FVY(I,J-1,K) -      FVY(I,J,K) )*RDY(J)        &    
                 + (        FVZ(I,J,K-1) -      FVZ(I,J,K) )*RDZ(K)        &    
                 - DDDT(I,J,K)
            LHSS = ((HP(I+1,J,K)-HP(I,J,K))*RDXN(I)*R(I) - (HP(I,J,K)-HP(I-1,J,K))*RDXN(I-1)*R(I-1) )*RDX(I)*RRN(I) &
                 + ((HP(I,J+1,K)-HP(I,J,K))*RDYN(J)      - (HP(I,J,K)-HP(I,J-1,K))*RDYN(J-1)        )*RDY(J)        &    
                 + ((HP(I,J,K+1)-HP(I,J,K))*RDZN(K)      - (HP(I,J,K)-HP(I,J,K-1))*RDZN(K-1)        )*RDZ(K)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   POIS_ERR = MAXVAL(RESIDUAL)

ENDDO MESHES_LOOP
   
END SUBROUTINE SCARC_CHECK_POISSON

!> ------------------------------------------------------------------------------------------------------
!> Check of how well the computed pressure satisfies the inseparable Poisson equation:
!> LHSS = del dot (1/rho) del p + del K = -del dot F - dD/dt = RHSS
!> ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ITERATE_BAROCLINIC
USE MESH_POINTERS
REAL(EB), POINTER, DIMENSION(:,:,:) :: PPP, RESIDUAL, RHOP, HP
REAL(EB) :: LHSS, RHSS
INTEGER :: NM, I, J, K

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)
   
   IF (PREDICTOR) THEN
      HP => H
      RHOP => RHO
   ELSE
      HP => HS
      RHOP => RHOS
   ENDIF

   PPP => WORK7
   PPP = RHOP*(HP-KRES)
   RESIDUAL => WORK8(1:IBAR,1:JBAR,1:KBAR)

   !$OMP PARALLEL PRIVATE(I,J,K,RHSS,LHSS,NM)
   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHSS = ( R(I-1)*(FVX(I-1,J,K)-FVX_B(I-1,J,K)) - R(I)*(FVX(I,J,K)-FVX_B(I,J,K)) )*RDX(I)*RRN(I) &
                 + (        (FVY(I,J-1,K)-FVY_B(I,J-1,K)) -      (FVY(I,J,K)-FVY_B(I,J,K)) )*RDY(J)        &
                 + (        (FVZ(I,J,K-1)-FVZ_B(I,J,K-1)) -      (FVZ(I,J,K)-FVZ_B(I,J,K)) )*RDZ(K)        &
                 - DDDT(I,J,K)
            LHSS = ((PPP(I+1,J,K)-PPP(I,J,K))*RDXN(I)*R(I)    *2._EB/(RHOP(I+1,J,K)+RHOP(I,J,K)) - &
                    (PPP(I,J,K)-PPP(I-1,J,K))*RDXN(I-1)*R(I-1)*2._EB/(RHOP(I-1,J,K)+RHOP(I,J,K)))*RDX(I)*RRN(I) &
                 + ((PPP(I,J+1,K)-PPP(I,J,K))*RDYN(J)         *2._EB/(RHOP(I,J+1,K)+RHOP(I,J,K)) - &
                    (PPP(I,J,K)-PPP(I,J-1,K))*RDYN(J-1)       *2._EB/(RHOP(I,J-1,K)+RHOP(I,J,K)))*RDY(J)        &
                 + ((PPP(I,J,K+1)-PPP(I,J,K))*RDZN(K)         *2._EB/(RHOP(I,J,K+1)+RHOP(I,J,K)) - &
                    (PPP(I,J,K)-PPP(I,J,K-1))*RDZN(K-1)       *2._EB/(RHOP(I,J,K-1)+RHOP(I,J,K)))*RDZ(K)        &
                 + ((KRES(I+1,J,K)-KRES(I,J,K))*RDXN(I)*R(I) - (KRES(I,J,K)-KRES(I-1,J,K))*RDXN(I-1)*R(I-1) )*RDX(I)*RRN(I) &
                 + ((KRES(I,J+1,K)-KRES(I,J,K))*RDYN(J)      - (KRES(I,J,K)-KRES(I,J-1,K))*RDYN(J-1)        )*RDY(J)        &
                 + ((KRES(I,J,K+1)-KRES(I,J,K))*RDZN(K)      - (KRES(I,J,K)-KRES(I,J,K-1))*RDZN(K-1)        )*RDZ(K)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
   !$OMP END PARALLEL

   PRESSURE_ERROR_MAX(NM) = MAXVAL(RESIDUAL)
   PRESSURE_ERROR_MAX_LOC(:,NM) = MAXLOC(RESIDUAL)

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_ITERATE_BAROCLINIC

!> ------------------------------------------------------------------------------------------------
!> Compute global matrix-vector product on grid level NL
!>                     Y := A*X
!> where NV1 is a reference to X and NV2 is a reference to Y
!> including data exchange along internal boundaries
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATVEC_PRODUCT(NV1, NV2, NL)
USE SCARC_POINTERS, ONLY: L, G, F, OG, GWC, AC, AB, V1, V2
INTEGER, INTENT(IN) :: NV1, NV2, NL           
REAL(EB) :: TNOW
INTEGER :: NM, NOM, IC, JC, IOR0, ICOL, INBR, ICE, ICW, ICG
INTEGER :: I, J, K, IW, IS=0, IT=0, IL=0, INUM1, INUM2
REAL(EB) :: TMP, VSAVE
#ifdef WITH_MKL
EXTERNAL :: DAXPBY
#endif

TNOW = CURRENT_TIME()

!>
!> If this call is related to a globally acting solver, exchange internal boundary values of
!> vector1 such that the ghost values contain the corresponding overlapped values of adjacent neighbor
!>
IF (TYPE_MATVEC == NSCARC_MATVEC_GLOBAL) THEN
    TYPE_VECTOR = NV1
    CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)
ENDIF

!>
!> Perform global matrix-vector product:
!> Note: - matrix already contains subdiagonal values from neighbor along internal boundaries
!>       - if vector1 contains neighboring values, then correct values of global matvec are achieved
!>
SELECT_MATRIX_TYPE: SELECT CASE (SCARC_MATRIX_LEVEL(NL))

   !>
   !> ------------- COMPACT storage technique
   !>
   CASE (NSCARC_MATRIX_COMPACT)
   
      MESHES_COMPACT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_GRID(NM, NL)
         AC => G%AC
   
         V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
         V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)
   
         DO IC = 1, G%NC
   
            ICOL = AC%ROW(IC)                                                !> diagonal entry
            JC   = AC%COL(ICOL)
            V2(IC) = AC%VAL(ICOL)* V1(JC)
            DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1                            !> subdiagonal entries
               JC = AC%COL(ICOL)
               V2(IC) =  V2(IC) + AC%VAL(ICOL) * V1(JC)
            ENDDO
         ENDDO
   
      ENDDO MESHES_COMPACT_LOOP
   
   !>
   !> ------------- BANDED storage technique
   !> matrix diagonals are supposed to be constant
   !> matrix-vector multiplication is based on daxpy-routines using the constant matrix stencil
   !> the 'wrong' entries due to boundary conditions and zero entries in subdiagonals are explicitly corrected 
   !>
   CASE (NSCARC_MATRIX_BANDED)
   
      SELECT_STENCIL_TYPE: SELECT CASE (TYPE_STENCIL)

         !>
         !> ------------- Variable entries with own implementation of daxpyv
         !> matrix-vector multiplication is based on variable matrix stencil
         !>
         CASE (NSCARC_STENCIL_VARIABLE)
         
            MESHES_BANDED_VARIABLE_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID(NM, NL)
               AB => G%AB
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               !> point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               !> point to Y-vector
               V2 = 0.0_EB
         
               !>$OMP PARALLEL default(none) PRIVATE(IOR0, IL, IS, IT) SHARED(AB, V1, V2) 
               DO IOR0 = 3, -3, -1

                  IF (AB%POS(IOR0) == 0) CYCLE

                  IL = AB%LENGTH(IOR0) 
                  IS = AB%SOURCE(IOR0)
                  IT = AB%TARGET(IOR0)

                  CALL SCARC_DAXPY_VARIABLE(IL, AB%VAL(IT:IT+IL-1,AB%POS(IOR0)), V1(IS:IS+IL-1), V2(IT:IT+IL-1))

               ENDDO
               !>$OMP END PARALLEL

               IF (NM==NMESHES) VSAVE = V2(G%NC)                        !> save value in case of pure Neumann bdry

               DO IOR0 = 3, -3, -1
                  F => L%FACE(IOR0)
                  DO INBR = 1, F%N_NEIGHBORS
                     NOM = F%NEIGHBORS(INBR)
                     CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)
                     DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                        ICW = OG%ICG_TO_ICW(ICG)
                        ICE = OG%ICG_TO_ICE(ICG)
                        V2(ICW) = V2(ICW) + F%SCAL_FACE * V1(ICE)
                     ENDDO
                  ENDDO
               ENDDO

               IF (IS_PURE_NEUMANN.AND.NM==NMESHES) V2(G%NC) = VSAVE   !> restore value in last cell

            ENDDO MESHES_BANDED_VARIABLE_LOOP

         !>
         !> ------------- Storage of variable matrix entries with local multiplications
         !> matrix-vector-multiplication according to compact storage technique
         !>
         CASE (NSCARC_STENCIL_VARIABLE_LOCAL)
         
            MESHES_BANDED_VARIABLE_LOCAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID(NM, NL)
               AB => G%AB
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               !> point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               !> point to Y-vector
         
               DO IC = 1, G%NC
                  V2(IC) = V1(IC)*AB%VAL(IC, AB%POS(0))                      !> main-diagonal contribution
                  DO IOR0 = 3, 1, -1                                         !> lower sub-diagonal contributions
                     IF (AB%POS(IOR0) == 0) CYCLE                            !> no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= 1) V2(IC) = V2(IC) + V1(JC)*AB%VAL(IC, AB%POS(IOR0))
                  ENDDO
                  DO IOR0 = -1, -3, -1                                       !> upper sub-diagonal contributions
                     IF (AB%POS(IOR0) == 0) CYCLE                            !> no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC <= G%NC) V2(IC) = V2(IC) + V1(JC)*AB%VAL(IC, AB%POS(IOR0))
                  ENDDO
               ENDDO
         
            ENDDO MESHES_BANDED_VARIABLE_LOCAL_LOOP
      
         !>
         !> ------------- Storage of constant matrix entries - with corrections at subdiagonals and diagonal (BC's)
         !>
         CASE (NSCARC_STENCIL_CONSTANT)

            MESHES_BANDED_CONSTANT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID(NM, NL)
               AB => G%AB
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               !> point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               !> point to Y-vector
               V2 = 0.0_EB
         
               DO IOR0 = 3, -3, -1

                  IF (AB%POS(IOR0) == 0) CYCLE
                  AB%AUX(1:G%NC)=0.0_EB
                  AB%AUX(1:AB%LENGTH(IOR0)) = V1(AB%SOURCE(IOR0):AB%SOURCE(IOR0)+AB%LENGTH(IOR0)-1)

                  IF (ABS(IOR0) == 1) THEN
                     DO IC = 1, G%NC
                        IF (MOD(IC,L%NX)==0) AB%AUX(IC)=0.0_EB
                     ENDDO
                  ELSE IF (ABS(IOR0) == 2) THEN
                     INUM1 = L%NX*L%NY
                     INUM2 = L%NX*L%NY - L%NX
                     DO IC = 1, G%NC
                        IF (MOD(IC,INUM1) > INUM2) AB%AUX(IC)=0.0_EB
                     ENDDO
                  ENDIF

                  IS = AB%SOURCE(IOR0)
                  IT = AB%TARGET(IOR0)
                  IL = AB%LENGTH(IOR0)

#ifdef WITH_MKL
                  CALL DAXPY(IL, AB%STENCIL(IOR0), AB%AUX(1:IL), 1, V2(IT:IT+IL-1), 1)
#else
                  WRITE(*,*) 'TODO: MATVEC: CONSTANT: NO-MKL: CHECK HERE'
                  V2(IT:IT+IL-1) = AB%STENCIL(IOR0) * AB%AUX(1:IL)
#endif

               ENDDO
               WALL_CELLS_BANDED_LOOP: DO IW = 1, G%NW 

                  IOR0 = G%WALL(IW)%IOR
                  IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE              !> cycle boundaries in y-direction for 2D-cases

                  GWC => G%WALL(IW)
                  F  => L%FACE(IOR0)
            
                  I = GWC%IXW
                  J = GWC%IYW
                  K = GWC%IZW
            
                  IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE 
            
                  IC  = G%CELL_NUMBER(I, J, K) 
            
                  !> SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE
                  IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN
                     TMP = V2(IC)
                     SELECT CASE (GWC%BTYPE)
                        CASE (DIRICHLET)
                           V2(IC) = V2(IC) - F%SCAL_BOUNDARY * V1(IC)
                        CASE (NEUMANN)
                           V2(IC) = V2(IC) + F%SCAL_BOUNDARY * V1(IC)
                     END SELECT
                  ENDIF 
            
               ENDDO WALL_CELLS_BANDED_LOOP
         
            ENDDO MESHES_BANDED_CONSTANT_LOOP
         !>
         !> ------------- Constant entries with own implementation of daxpy
         !> 'Wrong' entries due to boundary conditions and zero entries in subdiagonals are locally corrected 
         !>
         CASE (NSCARC_STENCIL_CONSTANT_LOCAL)
         
            MESHES_BANDED_CONSTANT_LOCAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID(NM, NL)
               AB => G%AB
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               !> point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               !> point to Y-vector
               V2 = 0.0_EB
         
               DO IOR0 = 3, -3, -1

                  IF (AB%POS(IOR0) == 0) CYCLE

                  IS = AB%SOURCE(IOR0)
                  IT = AB%TARGET(IOR0)
                  IL = AB%LENGTH(IOR0)

                  AB%AUX(1:G%NC)=0.0_EB
                  AB%AUX(1:IL) = V1(IS:IS+IL-1)

                  IF (ABS(IOR0) == 1) THEN
                     DO IC = 1, G%NC
                        IF (MOD(IC,L%NX)==0) AB%AUX(IC)=0.0_EB
                     ENDDO
                  ELSE IF (ABS(IOR0) == 2) THEN
                     INUM1 = L%NX*L%NY
                     INUM2 = L%NX*L%NY - L%NX
                     DO IC = 1, G%NC
                        IF (MOD(IC,INUM1) > INUM2) AB%AUX(IC)=0.0_EB
                     ENDDO
                  ENDIF

                  CALL SCARC_DAXPY_CONSTANT(IL, AB%STENCIL(IOR0), AB%AUX(1:IL), V2(IT:IT+IL-1))

               ENDDO
               WALL_CELLS_BANDED_OWN_LOOP: DO IW = 1, G%NW 

                  IOR0 = G%WALL(IW)%IOR
                  IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE              !> cycle boundaries in y-direction for 2D-cases

                  GWC => G%WALL(IW)
                  F  => L%FACE(IOR0)
            
                  I = GWC%IXW
                  J = GWC%IYW
                  K = GWC%IZW
            
                  IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE 
            
                  IC  = G%CELL_NUMBER(I, J, K) 
            
                  !> SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE
                  IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN
                  TMP = V2(IC)
                     SELECT CASE (GWC%BTYPE)
                        CASE (DIRICHLET)
                           V2(IC) = V2(IC) - F%SCAL_BOUNDARY * V1(IC)
                        CASE (NEUMANN)
                           V2(IC) = V2(IC) + F%SCAL_BOUNDARY * V1(IC)
                     END SELECT
                  ENDIF 
            
               ENDDO WALL_CELLS_BANDED_OWN_LOOP
         
            ENDDO MESHES_BANDED_CONSTANT_LOCAL_LOOP

      END SELECT SELECT_STENCIL_TYPE
END SELECT SELECT_MATRIX_TYPE

CPU(MYID)%MATVEC_PRODUCT =CPU(MYID)%MATVEC_PRODUCT+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_MATVEC_PRODUCT


!> ------------------------------------------------------------------------------------------------
!> Compute global scalar-product (including global data exchange)
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_SCALAR_PRODUCT(NV1, NV2, NL)
USE SCARC_POINTERS, ONLY: G, V1, V2
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB) :: TNOW, RANK_REAL
INTEGER :: NM
#ifdef WITH_MKL
REAL(EB) :: DDOT
EXTERNAL :: DDOT
#else
INTEGER :: IC
#endif

TNOW = CURRENT_TIME()

RANK_REAL = 0.0_EB
MESH_REAL = 0.0_EB

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
   V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

#ifdef WITH_MKL
   MESH_REAL(NM) = DDOT(G%NC, V1, 1, V2, 1)
#else
   MESH_REAL(NM) = 0.0_EB
   !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
   DO IC = 1, G%NC
      MESH_REAL(NM) = MESH_REAL(NM) + V1(IC) * V2(IC)
   ENDDO
   !$OMP END PARALLEL DO 
#endif

   RANK_REAL = RANK_REAL + MESH_REAL(NM)

ENDDO

!>
!> Compute global scalar product as sum of local scalar products
!>
IF (N_MPI_PROCESSES>1) & 
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_REAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERROR)

SCARC_SCALAR_PRODUCT = RANK_REAL

CPU(MYID)%SCALAR_PRODUCT = CPU(MYID)%SCALAR_PRODUCT + CURRENT_TIME()-TNOW
END FUNCTION SCARC_SCALAR_PRODUCT


!> ------------------------------------------------------------------------------------------------
!> Compute global L2-norm (including global data exchange)
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_L2NORM(NV1, NL)
INTEGER, INTENT(IN) :: NV1, NL
REAL(EB) :: TNOW
TNOW = CURRENT_TIME()

GLOBAL_REAL = SCARC_SCALAR_PRODUCT(NV1, NV1, NL)
GLOBAL_REAL = SQRT (GLOBAL_REAL)

SCARC_L2NORM = GLOBAL_REAL

CPU(MYID)%L2NORM =CPU(MYID)%L2NORM+CURRENT_TIME()-TNOW
END FUNCTION SCARC_L2NORM


!> ------------------------------------------------------------------------------------------------
!> Compute linear combination of two vectors for banded storage technique
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_SUM(NV1, NV2, SCAL1, SCAL2, NL)
USE SCARC_POINTERS, ONLY: G, V1, V2
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB), INTENT(IN) :: SCAL1, SCAL2
INTEGER :: NM
#ifdef WITH_MKL
EXTERNAL :: DAXPBY
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
   V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

#ifdef WITH_MKL
   CALL DAXPBY(G%NCE, SCAL1, V1, 1, SCAL2, V2, 1)
#else
   CALL SCARC_DAXPY_CONSTANT_DOUBLE(G%NCE, SCAL1, V1, SCAL2, V2)
#endif

ENDDO

END SUBROUTINE SCARC_VECTOR_SUM


!> ------------------------------------------------------------------------------------------------
!> Define vector2 to be a scaled copy of vector 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_COPY(NV1, NV2, SCAL1, NL)
USE SCARC_POINTERS, ONLY: G, V1, V2
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB), INTENT(IN) :: SCAL1
INTEGER :: NM
#ifdef WITH_MKL
EXTERNAL :: DCOPY, DSCAL
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
   V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

#ifdef WITH_MKL
   CALL DCOPY(G%NCE, V1, 1, V2, 1)
   CALL DSCAL(G%NCE, SCAL1, V2, 1)
#else
   CALL SCARC_SCALING_CONSTANT(G%NCE, SCAL1, V1, V2)
#endif

ENDDO

END SUBROUTINE SCARC_VECTOR_COPY


!> ------------------------------------------------------------------------------------------------
!> Clear vector
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_CLEAR(NV, NL)
USE SCARC_POINTERS, ONLY: VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   VC => SCARC_POINT_TO_VECTOR(NM, NL, NV)
   VC =  0.0_EB
ENDDO

END SUBROUTINE SCARC_VECTOR_CLEAR


!> ------------------------------------------------------------------------------------------------
!> Preset vector with specified value
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_RANDOM_INIT (NV, NL)
USE SCARC_POINTERS, ONLY: L, G, VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: IC, NM, I, J, K
REAL (EB) :: VAL

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   !$OMP PARALLEL DO PRIVATE(I, J, K, IC) SCHEDULE(STATIC)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
            IC = G%CELL_NUMBER(I,J,K)
            CALL RANDOM_NUMBER(VAL)
            VC(IC) = VAL
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO 

ENDDO

END SUBROUTINE SCARC_VECTOR_RANDOM_INIT


!> ------------------------------------------------------------------------------------------------
!> Preset vector with specified value
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_INIT (NV, VAL, NL)
USE SCARC_POINTERS, ONLY: L, G, VC
INTEGER, INTENT(IN) :: NV, NL
REAL (EB), INTENT(IN) :: VAL
INTEGER :: IC, NM, I, J, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   !$OMP PARALLEL DO PRIVATE(I, J, K, IC) SCHEDULE(STATIC)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
            IC = G%CELL_NUMBER(I,J,K)
            VC(IC) = VAL
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO 
   DO IC = G%NC+1, G%NCE
      VC(IC) = VAL
   ENDDO

ENDDO

END SUBROUTINE SCARC_VECTOR_INIT


!> ------------------------------------------------------------------------------------------------
!> Perform preconditioning
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELAXATION (NV1, NV2, NS, NP, NL)
USE SCARC_POINTERS, ONLY: L, G, AB, AC, FFT, V1, V2, V1_FB, V2_FB, ACS, MKL
USE POIS, ONLY: H2CZSS, H3CZSS
REAL(EB) :: AUX, OMEGA_SSOR = 1.5_EB 
REAL (EB) :: TNOW
INTEGER, INTENT(IN) :: NV1, NV2, NS, NP, NL
INTEGER :: NM, IC, JC, ICOL, ITYPE, IDIAG, IPTR, INCR, IOR0, IC0, IY, IZ

TNOW = CURRENT_TIME()
ITYPE = STACK(NS-1)%SOLVER%TYPE_RELAX

SELECT CASE (ITYPE)

   !> 
   !> --------- Preconditioning by blockwise Jacobi
   !> 
   CASE (NSCARC_RELAX_JAC)

      JACOBI_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID(NM, NL)
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         SELECT CASE (SCARC_MATRIX_LEVEL(NL))
            
            !>
            !> matrix in COMPACT storage technique
            !>
            CASE (NSCARC_MATRIX_COMPACT)

               AC => G%AC
               !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  V2(IC) = V2(IC) / AC%VAL(AC%ROW(IC))
               ENDDO
               !$OMP END PARALLEL DO

            !>
            !> matrix in BANDED storage technique
            !>
            CASE (NSCARC_MATRIX_BANDED)

               AB => G%AB
               !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  V2(IC) = V2(IC) / AB%VAL(IC, AB%POS(0))
               ENDDO
               !$OMP END PARALLEL DO

         END SELECT 

      ENDDO JACOBI_MESHES_LOOP

   !> 
   !> --------- Preconditioning by blockwise SSOR
   !> 
   CASE (NSCARC_RELAX_SSOR)

      SSOR_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID(NM, NL)
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         SELECT CASE (SCARC_MATRIX_LEVEL(NL))

            !>
            !> matrix in COMPACT storage technique
            !>
            CASE (NSCARC_MATRIX_COMPACT)

               AC => G%AC

               SSOR_FORWARD_COMPACT_LOOP: DO IC = 1, G%NC                                   !> forward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1
                     IF (AC%COL(ICOL) >= IC) CYCLE                                          !> only process lower diags
                     IF (AC%COL(ICOL) <= G%NC) THEN
                        AUX = AUX + AC%VAL(ICOL) * V2(AC%COL(ICOL))  !> ignore overlaps
                     ENDIF
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AC%VAL(AC%ROW(IC))
               ENDDO SSOR_FORWARD_COMPACT_LOOP

               SSOR_BACKWARD_COMPACT_LOOP: DO IC = G%NC-1, 1, -1                           !> backward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1
                     IF (AC%COL(ICOL) <= IC) CYCLE                                         !> only process upper diags
                     IF (AC%COL(ICOL) <= G%NC) THEN                                        !> ignore overlaps
                        AUX = AUX + AC%VAL(ICOL) * V2(AC%COL(ICOL))
                     ENDIF
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AC%VAL(AC%ROW(IC))
               ENDDO SSOR_BACKWARD_COMPACT_LOOP

         !>
         !> matrix in BANDED storage technique
         !>
         CASE (NSCARC_MATRIX_BANDED)

            AB => G%AB

            !>
            !> 2D version
            !>
            IF (TWO_D) THEN

               SSOR_FORWARD_BANDED_2D_LOOP: DO IC = 1, G%NC                 !> forward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = 1, 3, 2                                          !> only process lower x- and z-diag
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= 1 .AND. JC <= G%NC) AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_FORWARD_BANDED_2D_LOOP

               SSOR_BACKWARD_BANDED_2D_LOOP: DO IC = G%NC-1, 1, -1          !> backward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = -1, -3, -2                                       !> only process upper x- and z-diag
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC <= G%NC) THEN
                        AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                     ENDIF
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_BACKWARD_BANDED_2D_LOOP

            !>
            !> 3D version
            !>
            ELSE

               SSOR_FORWARD_BANDED_3D_LOOP: DO IC = 1, G%NC                  !> forward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = 1, 3                                             !> only process lower diags
                     IF (AB%POS(IOR0) == 0) CYCLE                            !> no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= 1 .AND. JC <= G%NC) AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_FORWARD_BANDED_3D_LOOP

               SSOR_BACKWARD_BANDED_3D_LOOP: DO IC = G%NC-1, 1, -1           !> backward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = -1, -3, -1                                       !> only process upper diags
                     IF (AB%POS(IOR0) == 0) CYCLE                            !> no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= IC .AND. JC <= G%NC) AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_BACKWARD_BANDED_3D_LOOP

            ENDIF

         END SELECT 

      ENDDO SSOR_MESHES_LOOP

   !> 
   !> --------- Preconditioning by Jacobi in matrix form
   !> 
   CASE (NSCARC_RELAX_MJAC)

      MJAC_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID(NM, NL)

         V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         SELECT CASE(TYPE_MATRIX)

            !> ------------ matrix in COMPACT storage technique
            CASE (NSCARC_MATRIX_COMPACT)
               AC => G%AC
               CALL SCARC_SCALING_VARIABLE(G%NC, AC%RELAX, V1, V2)

            !> ------------ matrix in BANDED storage technique
            CASE (NSCARC_MATRIX_BANDED)
               AB => G%AB
               CALL SCARC_SCALING_VARIABLE(G%NC, AB%RELAXD, V1, V2)

         END SELECT

      ENDDO MJAC_MESHES_LOOP

   !> 
   !> --------- Preconditioning by different matrix-form preconditioners
   !> in all cases the preconditioner is given as separate matrix which is based
   !> on the same storage technique as the matrix AC itself;
   !> two tridiagonal systems have to be solved
   !> V1 contains the RHS to be solved for, V2 will contain the solution
   !> 
   CASE (NSCARC_RELAX_MGS, NSCARC_RELAX_MSGS, NSCARC_RELAX_MSOR, NSCARC_RELAX_MSSOR, NSCARC_RELAX_ILU)

      LU_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID(NM, NL)

         V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)
      
         SELECT CASE(TYPE_MATRIX)

            !> ------------ matrix in COMPACT storage technique
            CASE (NSCARC_MATRIX_COMPACT)

               AC => G%AC
      
               !> Forward solve:   Solve V2 = L^-1 V1
               !> Compute sol(i) = rhs(i) - sum L(i,j) x sol(j)
               DO IC = 1, G%NC
                  V2(IC) = V1(IC)
                  DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
                     JC = AC%COL(IPTR)
                     IF (JC >= IC) CYCLE
                     V2(IC) = V2(IC) - AC%RELAX(IPTR) * V2(JC)
                  ENDDO
               ENDDO
      
               !> If preconditioner is not symmetric, upper matrix U is zero and nothing has to be solved
               IF (ITYPE == NSCARC_RELAX_MGS .OR. ITYPE == NSCARC_RELAX_MSOR) CYCLE
      
               !> Backward solve
               !> Compute sol: inv(U) sol
               DO IC = G%NC, 1, -1
      
                  DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
                     JC = AC%COL(IPTR)
                     IF (JC <= IC) CYCLE
                     V2(IC) = V2(IC) - AC%RELAX(IPTR) * V2(JC)
                  ENDDO
      
                  !> Compute sol(i) = sol(i)/U(i,i)
                  IDIAG = AC%ROW(IC)
                  V2(IC) = V2(IC)/AC%RELAX(IDIAG)
      
               ENDDO
      

            !>
            !> matrix in BANDED storage technique
            !>
            CASE (NSCARC_MATRIX_BANDED)

               AB => G%AB
      
               IF (TWO_D) THEN
                  INCR = -2
               ELSE 
                  INCR = -1
               ENDIF
               
               !> Forward solve:   Solve V2 = L^-1 V1
               !> Compute sol(i) = rhs(i) - sum L(i,j) x sol(j)
               !>>$OMP PARALLEL DO PRIVATE(IC, JC, IOR0) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  V2(IC) = V1(IC)
                  DO IOR0 = 3, 1, INCR
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC <= 0) CYCLE
                     V2(IC) = V2(IC) - AB%RELAX(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
               ENDDO
               !>>$OMP END PARALLEL DO 
      
               !> If preconditioner is not symmetric, upper matrix U is zero and nothing has to be solved
               IF (ITYPE == NSCARC_RELAX_MGS .OR. ITYPE == NSCARC_RELAX_MSOR) CYCLE
      
               !> Backward solve
               !> Compute sol: inv(U) sol
               !>>$OMP PARALLEL DO PRIVATE(IC, JC, IOR0) SCHEDULE(STATIC)
               DO IC = G%NC, 1, -1
      
                  DO IOR0 = -1, -3, INCR
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC > G%NC) CYCLE
                     V2(IC) = V2(IC) - AB%RELAX(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
      
                  !> Compute sol(i) = sol(i)/U(i,i)
                  V2(IC) = V2(IC)/AB%RELAX(IC, AB%POS(0))
               ENDDO
               !>>$OMP END PARALLEL DO 
      
         END SELECT
      
      ENDDO LU_MESHES_LOOP
      

   !> 
   !> --------- Preconditioning by blockwise Geometric Multigrid
   !> 
   CASE (NSCARC_RELAX_GMG)

      CALL SCARC_METHOD_MULTIGRID (NS, NP, NSCARC_RHS_DEFECT, NLEVEL_MIN)

   !> 
   !> --------- Preconditioning by blockwise FFT based on Crayfishpak
   !> 
   CASE (NSCARC_RELAX_FFT)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID(NM, NL)
         FFT => L%FFT

         V1  => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2  => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
         DO IC = 1, G%NC
            FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = V1(IC)
         ENDDO
         !$OMP END PARALLEL DO

         IF (TWO_D) THEN
            CALL H2CZSS (FFT%BXS,  FFT%BXF, FFT%BZS, FFT%BZF, FFT%ITRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ELSE
            CALL H3CZSS (FFT%BXS,  FFT%BXF, FFT%BYS, FFT%BYF, FFT%BZS, FFT%BZF, FFT%ITRN, FFT%JTRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ENDIF

         !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
         DO IC = 1, G%NC
            V2(IC) = FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) 
         ENDDO
         !$OMP END PARALLEL DO 

      ENDDO

   !> 
   !> --------- Preconditioning by blockwise overlapping FFT based on Crayfishpak 
   !>           still test-version for tunnel-shaped geometries of type Mx1
   !> 
   CASE (NSCARC_RELAX_FFTO)

      !> Exchange overlapping parts
      TYPE_VECTOR = NV1           
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID(NM, NL)
         FFT => L%FFT

         V1  => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2  => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         IF (NM == 1) THEN

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = V1(IC)
            ENDDO
            !$OMP END PARALLEL DO

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  FFT%PRHS(FFT%IBAR, IY, IZ) = V1(IC0)
               IC0 = IC0 + 1
               ENDDO
            ENDDO

         ELSE IF (NM == NMESHES) THEN

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               FFT%PRHS(G%ICX(IC)+1, G%ICY(IC), G%ICZ(IC)) = V1(IC)
            ENDDO
            !$OMP END PARALLEL DO

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  FFT%PRHS(1, IY, IZ) = V1(IC0)
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ELSE 

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               FFT%PRHS(G%ICX(IC)+1, G%ICY(IC), G%ICZ(IC)) = V1(IC)
            ENDDO
            !$OMP END PARALLEL DO

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  FFT%PRHS(1, IY, IZ) = V1(IC0)
                  IC0 = IC0 + 1
               ENDDO
            ENDDO
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  FFT%PRHS(FFT%IBAR, IY, IZ) = V1(IC0)
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ENDIF

         !> 
         !> Call corresponding FFT solver
         !> 
         IF (TWO_D) THEN
            CALL H2CZSS (FFT%BXS,  FFT%BXF, FFT%BZS, FFT%BZF, FFT%ITRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ELSE
            CALL H3CZSS (FFT%BXS,  FFT%BXF, FFT%BYS, FFT%BYF, FFT%BZS, FFT%BZF, FFT%ITRN, FFT%JTRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ENDIF

         !> 
         !> Extract computed data from FFT%PRHS
         !> 
         IF (NM == 1) THEN
            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               V2(IC) = FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) 
            ENDDO
            !$OMP END PARALLEL DO 

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  V2(IC0) = FFT%PRHS(FFT%IBAR, IY, IZ) 
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ELSE IF (NM == NMESHES) THEN

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               V2(IC) = FFT%PRHS(G%ICX(IC)+1, G%ICY(IC), G%ICZ(IC)) 
            ENDDO
            !$OMP END PARALLEL DO 

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  V2(IC0) = FFT%PRHS(1, IY, IZ) 
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ELSE

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               V2(IC) = FFT%PRHS(G%ICX(IC)+1, G%ICY(IC), G%ICZ(IC)) 
            ENDDO
            !$OMP END PARALLEL DO 

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  V2(IC0) = FFT%PRHS(1, IY, IZ) 
                  IC0 = IC0 + 1
               ENDDO
            ENDDO
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  V2(IC0) = FFT%PRHS(FFT%IBAR, IY, IZ) 
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ENDIF

      ENDDO

      !> Exchange overlapping parts
      TYPE_VECTOR = NV2           
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_MEAN, NL)


#ifdef WITH_MKL
   !> 
   !> --------- Preconditioning by LU-decomposition
   !> 
   CASE (NSCARC_RELAX_MKL)

      !> Preconditioning by Cluster Sparse Solver from MKL
      !>
      MKL_SCOPE_IF: IF (STACK(NS)%SOLVER%TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL) THEN

         MKL_SCOPE_GLOBAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

            CALL SCARC_POINT_TO_GRID(NM, NL)
            MKL => L%MKL
            ACS => G%ACS

            MKL%PHASE  = 33                            ! only solving

            V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
            V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

            IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

               V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV1)
               V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV2)

               V1_FB(1:G%NC) = REAL(V1(1:G%NC), FB)
               V2_FB(1:G%NC) = 0.0_FB

               CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                            ACS%VAL_FB, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                            MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)

               V2(1:G%NC) = REAL(V2_FB(1:G%NC), EB)

            ELSE
               CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                            ACS%VAL, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                            MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)
            ENDIF
            IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

         ENDDO MKL_SCOPE_GLOBAL_LOOP

      !>
      !> Preconditioning by Pardiso Solver from MKL
      !>
      ELSE MKL_SCOPE_IF

         MKL_SCOPE_LOCAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

            CALL SCARC_POINT_TO_GRID(NM, NL)
            MKL => L%MKL
            ACS => G%ACS

            MKL%PHASE  = 33                            ! only solving

            V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
            V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

            IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

               V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV1)
               V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV2)

               V1_FB(1:G%NC) = REAL(V1(1:G%NC), FB)
               V2_FB(1:G%NC) = 0.0_FB
   
               CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                              ACS%VAL_FB, ACS%ROW, ACS%COL, &
                              MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

               V2(1:G%NC) = REAL(V2_FB(1:G%NC), EB)
            ELSE

               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

               CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                              ACS%VAL, ACS%ROW, ACS%COL, &
                              MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1, V2, MKL%ERROR)

            ENDIF
            IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

         ENDDO MKL_SCOPE_LOCAL_LOOP

      ENDIF MKL_SCOPE_IF

#endif

END SELECT

CPU(MYID)%RELAXATION =CPU(MYID)%RELAXATION+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_RELAXATION


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CLUSTER(NSTACK, NPARENT, NLEVEL)
USE SCARC_POINTERS, ONLY: L, G, MKL, ACS, V1, V2, V1_FB, V2_FB
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NS, NP, NL
REAL (EB) :: TNOW

NS = NSTACK
NP = NPARENT
NL = NLEVEL

TNOW = CURRENT_TIME()

CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NSCARC_RHS_INHOMOGENEOUS)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)
   ACS => G%ACS

   V1 => SCARC_POINT_TO_VECTOR (NM, NL, B)
   V2 => SCARC_POINT_TO_VECTOR (NM, NL, X)

   MKL => L%MKL
   MKL%PHASE  = 33                                !> only solving

   IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

      V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, B)
      V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, X)

      V1_FB = REAL(V1, FB)
      V2_FB = 0.0_FB
      CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                   ACS%VAL_FB, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)
      V2 = REAL(V2_FB, EB)

   ELSE

      CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                   ACS%VAL, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)
   ENDIF

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

ENDDO MESHES_LOOP

TYPE_VECTOR = X
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NS, NP)

END SUBROUTINE SCARC_METHOD_CLUSTER
#endif


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_PARDISO(NSTACK, NPARENT, NLEVEL)
USE SCARC_POINTERS, ONLY: L, G, MKL, ACS, V1, V2, V1_FB, V2_FB
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NS, NP, NL
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

NS = NSTACK
NP = NPARENT
NL = NLEVEL

CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NSCARC_RHS_INHOMOGENEOUS)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)
   ACS => G%ACS

   V1 => SCARC_POINT_TO_VECTOR (NM, NL, B)
   V2 => SCARC_POINT_TO_VECTOR (NM, NL, X)

   MKL => L%MKL
   MKL%PHASE  = 33         ! only solving


   IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

      V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, B)
      V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, X)

      V1_FB = REAL(V1, FB)
      V2_FB = 0.0_FB
      CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                     ACS%VAL_FB, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

      V2 = REAL(V2_FB, EB)

   ELSE

      V2 = 0.0_EB
      CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                     ACS%VAL, ACS%ROW, ACS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, V1, V2, MKL%ERROR)
   ENDIF

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

ENDDO MESHES_LOOP

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NSTACK, NPARENT)

END SUBROUTINE SCARC_METHOD_PARDISO
#endif


!> ------------------------------------------------------------------------------------------------
!> Increase corresponding iteration count (just for visualization of convergence behavior)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INCREASE_ITERATION_COUNTS(ITE0)
INTEGER, INTENT(IN) :: ITE0

SELECT CASE (TYPE_SOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      SELECT CASE (TYPE_METHOD)
         CASE (NSCARC_METHOD_KRYLOV)
            ITE_CG = ITE0
         CASE (NSCARC_METHOD_MULTIGRID)
            ITE_MG = ITE0
         CASE (NSCARC_METHOD_LU)
            ITE_LU = ITE0
      END SELECT
   CASE (NSCARC_SOLVER_PRECON)
      ITE_MG = ITE0
   CASE (NSCARC_SOLVER_SMOOTH)
      ITE_SMOOTH = ITE0
   CASE (NSCARC_SOLVER_COARSE)
      ITE_COARSE = ITE0
END SELECT
ITE_TOTAL = ITE_TOTAL + 1

END SUBROUTINE SCARC_INCREASE_ITERATION_COUNTS


!> ------------------------------------------------------------------------------------------------
!> Perform global CG-method based on global Possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CG(NSTACK, NPARENT, NRHS, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NRHS, NLEVEL
INTEGER :: NSTATE, NS, NP, NL
REAL (EB) :: ALPHA, BETA, SIGMA, SIGMA0=0.0_EB
REAL (EB) :: TNOW, TNOWI

TNOW = CURRENT_TIME()
ITE_CG = 0

!> get current and parent stack position, and current level
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS  = NSTACK
NP  = NPARENT
NL  = NLEVEL

!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Get parameters for current scope (note: NL denotes the finest level)
!>   - Get right hand side vector and clear solution vectors
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NRHS)

!> In case of pure Neumann boundary conditions setup condensed system
IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_VECTOR_INIT (X, 0.0_EB, NL)                    
   CALL SCARC_FILTER_MEANVALUE(B, NL)                       
   CALL SCARC_SETUP_SYSTEM_CONDENSED (B, NL, 1)            
ENDIF

!> Compute initial residual 
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
CALL SCARC_MATVEC_PRODUCT (X, R, NL)                         !>  r^0 := A*x^0
CALL SCARC_VECTOR_SUM     (B, R, -1.0_EB, 1.0_EB, NL)        !>  r^0 := r^0 - b     corresponds to  A*x^0 - b

RES    = SCARC_L2NORM (R, NL)                                !>  res   := ||r^0||
RESIN  = RES                                                 !>  resin := res
NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                 !>  res < tolerance ?

!> Perform initial preconditioning
IF (NSTATE /= NSCARC_STATE_CONV0) THEN                       !>  if no convergence yet, start precon
   CALL SCARC_PRECONDITIONER(NS, NS, NL)                     !>  v^0 := Precon(r^0)
   SIGMA0 = SCARC_SCALAR_PRODUCT(R, V, NL)                   !>  SIGMA0 := (r^0,v^0)
   CALL SCARC_VECTOR_COPY (V, D, -1.0_EB, NL)                !>  d^0 := -v^0
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Perform conjugate gradient looping
!> ------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, NIT

   TNOWI = CURRENT_TIME()
   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   !TYPE_MATVEC = NSCARC_MATVEC_LOCAL
   CALL SCARC_MATVEC_PRODUCT (D, Y, NL)                      !>  y^k := A*d^k

   ALPHA = SCARC_SCALAR_PRODUCT (D, Y, NL)                   !>  alpha := (d^k,y^k)     corresponds to   (d^k,A*d^k)
   ALPHA = SIGMA0/ALPHA                                      !>  alpha := (r^k,v^k)/(d^k,A*d^k)

   CALL SCARC_VECTOR_SUM (D, X, ALPHA, 1.0_EB, NL)           !>  x^{k+1} := x^k + alpha * d^k
   CALL SCARC_VECTOR_SUM (Y, R, ALPHA, 1.0_EB, NL)           !>  r^{k+1} := r^k + alpha * y^k   ~  r^k + alpha * A * d^k

   RES = SCARC_L2NORM (R, NL)                                !>  res := ||r^{k+1}||
   NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)              !>  res < tolerance ??
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP

   CALL SCARC_PRECONDITIONER(NS, NS, NL)                     !>  v^{k+1} := Precon(r^{k+1})

   SIGMA  = SCARC_SCALAR_PRODUCT (R, V, NL)                  !>  sigma := (r^{k+1},v^{k+1})
   BETA   = SIGMA/SIGMA0                                     !>  beta  := (r^{k+1},v^{k+1})/(r^k,v^k)
   SIGMA0 = SIGMA                                            !>  save last sigma

   CALL SCARC_VECTOR_SUM (V, D, -1.0_EB, BETA, NL)           !>  d^{k+1} := -v^{k+1} + beta * d^{k+1}

   CPU(MYID)%ITERATION=MAX(CPU(MYID)%ITERATION,CURRENT_TIME()-TNOWI)

ENDDO CG_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of CG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_RESTORE_LAST_CELL(X, NL)
   CALL SCARC_FILTER_MEANVALUE(X, NL)
ENDIF

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_MAINCELLS(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
#ifdef WITH_SCARC_STANDALONE
   CALL SCARC_PRESSURE_DIFFERENCE(NLEVEL_MIN)
#endif
ENDIF

CALL SCARC_RELEASE_SOLVER(NS, NP)

END SUBROUTINE SCARC_METHOD_CG


!> -----------------------------------------------------------------------------------------------
!> Preconditioning method which is based on the following input and output convention:
!>  - the residual which has to be preconditioned is passed in via vector R
!>  - the result of preconditioning is passed out via vector V
!>  - for several variants Y and Z are used as auxiliary vectors
!>  - in the comments: call is based on current grid level l (mostly the finest one)
!>  -                  l=1 denotes the finest  grid level NLEVEL_MIN
!>  -                  l=L denotes the coarset grid level NLEVEL_MAX
!> -----------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRECONDITIONER(NS, NP, NL)
INTEGER, INTENT(IN) :: NS, NP, NL     
INTEGER :: IL

SELECT_PRECON_TYPE: SELECT CASE (TYPE_TWOLEVEL)

   !> 
   !> ---------- Classical one-level preconditioning
   !> 
   CASE (NSCARC_TWOLEVEL_NONE)

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !>  v := r
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !>  v := Relax(r)

   !> 
   !> ---------- Additive two-level preconditioning
   !> 
   CASE (NSCARC_TWOLEVEL_ADD)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !>  Use r^l as right hand side for preconditioner
      DO IL = NL, NLEVEL_MAX-1                                    !>  successively restrict to coarser levels up to coarsest
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !>  b^{l+1} := Restriction(r^l)
      ENDDO
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Z, 1.0_EB, NLEVEL_MAX)           !>  z^L := x^L

      DO IL = NLEVEL_MAX-1, NL, -1                                !>  successively interpolate to finer levels up to finest
         CALL SCARC_PROLONGATION(Z, Z, IL+1, IL)                  !>  z^l := Prolongation(z^{l+1})
      ENDDO
      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !>  v^l := r^l
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !>  v^l := Relax(r^l)
      CALL SCARC_VECTOR_SUM (Z, V, 1.0_EB, 1.0_EB, NL)            !>  v^l := z^l + v^l

   !> 
   !> ---------- Multiplicative two-level preconditioning (coarse first, fine second)
   !> 
   CASE (NSCARC_TWOLEVEL_MUL)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !>  Use r^l as right hand side for preconditioner

      DO IL = NL, NLEVEL_MAX-1
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !>  b^{l+1} := Restriction(r^l)
      ENDDO

      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !>  y^L := x^L

      DO IL = NLEVEL_MAX-1, NL, -1
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !>  y^l := Prolongation(y^{l+1})
      ENDDO
      CALL SCARC_MATVEC_PRODUCT (Y, Z, NL)                        !>  z^l := A^l * y^l

      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !>  z^l := r^l - z^l
      CALL SCARC_VECTOR_COPY (Z, V, 1.0_EB, NL)                   !>  v^l := z^l
      CALL SCARC_RELAXATION (Z, V, NS+1, NP, NL)                  !>  v^l := Relax(z^l)
      CALL SCARC_VECTOR_SUM (Y, V, 1.0_EB, 1.0_EB, NL)            !>  v^l := y^l - z^l

   !> 
   !> ---------- Multiplicative two-level preconditioning (fine first, coarse second):
   !> coarse level is one level away from finest one (one coarsening step)
   !> 
   CASE (NSCARC_TWOLEVEL_MUL2)

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !>  v^l := r^l
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !>  v^l := Relax(r^l)
      CALL SCARC_MATVEC_PRODUCT (V, Z, NL)                        !>  z^l := A^{l} * v^l

      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !>  z^l := r^l - z^l

      CALL SCARC_RESTRICTION (Z, B, NL, NL+1)                     !>  b^{l+1} := rest(R^{l})
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  x^{l+1} := A^{l+1}^{-1}(b^{l+1})
      CALL SCARC_PROLONGATION (X, Z, NL+1, NL)                    !>  v^l := Prolongation(x^{l+1})
      CALL SCARC_VECTOR_SUM (Z, V, 1.0_EB, 1.0_EB, NL)            !>  z^l := r^l - z^l

   !> 
   !> ---------- Only coarse grid preconditioner
   !> 
   CASE (NSCARC_TWOLEVEL_COARSE)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !>  Use r^l as right hand side for preconditioner
      DO IL = NL, NLEVEL_MAX-1                                    !>  successively restrict to coarser levels up to coarsest
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !>  b^{l+1} := Restriction(b^l)
      ENDDO

      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !>  y^L := x^L

      DO IL = NLEVEL_MAX-1, NL, -1                                !>  successively interpolate to finer levels up to finest
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !>  y^l := Prolongation(y^{l+1})
      ENDDO
      CALL SCARC_VECTOR_COPY (Y, V, 1.0_EB, NL)                   !>  v^l := y^l

END SELECT SELECT_PRECON_TYPE

END SUBROUTINE SCARC_PRECONDITIONER


!> ------------------------------------------------------------------------------------------------
!> Call requested coarse grid solver (iterative/direct)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_COARSE(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL

SELECT CASE (TYPE_COARSE)

   CASE (NSCARC_COARSE_ITERATIVE)
      CALL SCARC_METHOD_CG (NSTACK, NPARENT, NSCARC_RHS_DEFECT, NLEVEL)

   CASE (NSCARC_COARSE_DIRECT)
#ifdef WITH_MKL
      IF (STACK(NPARENT)%SOLVER%TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL .AND. NMESHES > 1) THEN
         CALL SCARC_METHOD_CLUSTER (NSTACK, NPARENT, NLEVEL)
      ELSE
         CALL SCARC_METHOD_PARDISO (NSTACK, NPARENT, NLEVEL)
      ENDIF
#else
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_DIRECT_NOMKL, SCARC_NONE, NLEVEL)
#endif

END SELECT

END SUBROUTINE SCARC_METHOD_COARSE


!> ------------------------------------------------------------------------------------------------
!> Perform geometric multigrid method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MULTIGRID(NSTACK, NPARENT, NRHS, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NRHS, NLEVEL
INTEGER :: NS, NP, NL
INTEGER :: NSTATE, ICYCLE
REAL (EB) :: TNOW, TNOW_COARSE

TNOW = CURRENT_TIME()
ITE_MG = 0

!> store current and parent stack position and current level
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS = NSTACK
NP = NPARENT
NL = NLEVEL

!>
!> ---------- Initialization:
!>   - Save SETTING (in case that subsequent solvers with different SETTING are called)
!>   - Define parameters for current scope (note: NL denotes the finest level)
!>   - Initialize solution, right hand side vector
!> 
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NRHS)

!> 
!> ---------- Compute initial defect:  
!>            RESIN := || B - A*X ||
!>   - Initialize cycle counts for MG-iteration
!>   - Perform initial matrix-vector product on finest level
!>   - calculate norm of initial residual on finest level
!> 
CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                  !>  V := A*X
CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                     !>  V := B - V

RES    = SCARC_L2NORM (V, NL)                                         !>  RESIN := ||V||
RESIN  = RES
NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                          !>  RES < TOL already ??

ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_SETUP, NL)

!> 
!> ---------- Perform multigrid-looping (start each iteration on finest level)
!> 
MULTIGRID_LOOP: DO ITE = 1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLING_EXIT)

      !>
      !> Presmoothing  (smoothing/restriction till coarsest level is reached)
      !> initial and final residual are passed via vector V by default
      !>
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_PRESMOOTH, NS+1, NS, NL)         !> D_fine   := Smooth(defect)
         CALL SCARC_RESTRICTION (V, B, NL, NL+1)                              !> B_coarse := Rest(D_fine)
         CALL SCARC_VECTOR_CLEAR (X, NL+1)                                    !> use zero initial guess on coarse level
         NL = NL + 1                                                          !> set coarser level
      ENDDO PRESMOOTHING_LOOP

      !>
      !> Coarse grid solver
      !>
      TNOW_COARSE = CURRENT_TIME()
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)                          !> X_coarse := exact_sol(.)
      CPU(MYID)%COARSE =CPU(MYID)%COARSE+CURRENT_TIME()-TNOW_COARSE

      !>
      !> Postsmoothing (smoothing/restriction till finest level is reached again)
      !>
      POSTSMOOTHING_LOOP: DO WHILE (NL > NLEVEL_MIN)
         NL=NL-1
         CALL SCARC_PROLONGATION (X, V, NL+1, NL)                             !> V_fine := Prol(X_coarse)
         CALL SCARC_VECTOR_SUM (V, X, 1.0_EB, 1.0_EB, NL)                     !> X_fine := V_fine + X_fine
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_POSTSMOOTH, NS+1, NS, NL)        !> V_fine := Smooth(defect)
         ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_PROCEED, NL)           !> perform requested cycle
         IF (ICYCLE /= NSCARC_CYCLING_POSTSMOOTH) CYCLE CYCLE_LOOP
      ENDDO POSTSMOOTHING_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MULTIGRID_LEVEL, SCARC_NONE, NL)

   !> 
   !> Compute norm of new residual on finest level and  leave loop correspondingly
   !> 
   CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                       !> V := A*X
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                          !> V := F - V

   RES = SCARC_L2NORM (V, NL)                                                 !> RES := ||V||
   NSTATE = SCARC_CONVERGENCE_STATE(0, NS, NL)                                !> convergence ?
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT MULTIGRID_LOOP

ENDDO MULTIGRID_LOOP

!> 
!> ---------- Determine convergence rate and print corresponding information:
!> In case of MG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries (consistency!)
!> 
CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)

SELECT CASE (TYPE_SOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      CALL SCARC_UPDATE_MAINCELLS(NLEVEL_MIN)
      CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
#ifdef WITH_SCARC_STANDALONE
      CALL SCARC_PRESSURE_DIFFERENCE(NLEVEL_MIN)
#endif
   CASE (NSCARC_SOLVER_PRECON)
      CALL SCARC_UPDATE_PRECONDITIONER(NLEVEL_MIN)
END SELECT

CALL SCARC_RELEASE_SOLVER(NS, NP)

END SUBROUTINE SCARC_METHOD_MULTIGRID


!> ------------------------------------------------------------------------------------------------
!> Control multigrid cycling (F/V/W)
!> Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CYCLING_CONTROL(NTYPE, NL)
USE SCARC_POINTERS, ONLY: MG
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, IL, ICYCLE

SELECT CASE (NTYPE)

   !> 
   !> initialize cycle counts at beginning of multigrid method
   !> 
   CASE (NSCARC_CYCLING_SETUP)
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         MG => SCARC(NM)%LEVEL(NL)%MG
         MG%CYCLING(2)=1
   
         DO IL = NLEVEL_MIN+1, NLEVEL_MAX - 1
            MG => SCARC(NM)%LEVEL(IL)%MG
            IF (TYPE_CYCLING==NSCARC_CYCLING_F) THEN
               MG%CYCLING(2)=2
            ELSE
               MG%CYCLING(2)=TYPE_CYCLING
            ENDIF
         ENDDO
      ENDDO
   
      ICYCLE = NSCARC_CYCLING_NEXT
   
   !> 
   !> reset cycle counts at beginning of each new multigrid iteration
   !> 
   CASE (NSCARC_CYCLING_RESET)
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         DO IL = NLEVEL_MIN, NLEVEL_MAX
            MG => SCARC(NM)%LEVEL(IL)%MG
            MG%CYCLING(1)=MG%CYCLING(2)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLING_NEXT
   
   !> 
   !> determine where to proceed with cycling
   !> 
   CASE (NSCARC_CYCLING_PROCEED)
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         MG => SCARC(NM)%LEVEL(NL)%MG
         MG%CYCLING(1)=MG%CYCLING(1)-1
   
         IF (MG%CYCLING(1)==0) THEN
            IF (TYPE_CYCLING==NSCARC_CYCLING_F) THEN
               MG%CYCLING(1)=1
            ELSE
               MG%CYCLING(1)=MG%CYCLING(2)
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
INTEGER, INTENT(IN) :: NTYPE, NSTACK, NPARENT, NLEVEL
INTEGER :: NSTATE=0, NS, NP, NL
REAL(EB) :: TNOW
LOGICAL :: BMATVEC, BL2NORM, BVERBOSE

!>
!> ---------- Initialization
!>
TNOW = CURRENT_TIME()
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS = NSTACK
NP = NPARENT
NL = NLEVEL

CALL SCARC_SETUP_SOLVER(NS, NP)

!>
!> Calculate initial defect on l2-norm on level NL (only if BMATVEC and Bl2NORM are set to .TRUE.)
!> Because initial vector in MG is set to zero, this defect corresponds to F
!>
ITE = 0
BVERBOSE = .TRUE.
IF (BVERBOSE) THEN
   BL2NORM  = .TRUE.
   BMATVEC  = .TRUE.
ELSE
   BL2NORM  = .FALSE.
   BMATVEC  = .FALSE.
ENDIF

IF (BMATVEC) THEN
   CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                  !>  v := A*x
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                     !>  v := b - v    corresponds to   b - A*x
ENDIF

IF (BL2NORM) THEN
   RESIN = SCARC_L2NORM (V, NL)                                          !>  resin := ||v||
ELSE
   RESIN = SCARC_RESIDUAL
ENDIF
IF (BVERBOSE) NSTATE = SCARC_CONVERGENCE_STATE(NTYPE, NS, NL)

!>
!> ---------- Smoothing loop - only temporarily
!>
IF (NTYPE == NSCARC_CYCLING_PRESMOOTH) THEN
   NIT = SCARC_MULTIGRID_PRESMOOTH
ELSE
   NIT = SCARC_MULTIGRID_POSTSMOOTH
ENDIF

SMOOTH_LOOP: DO ITE=1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

#ifdef WITH_MKL
   IF (TYPE_SMOOTH == NSCARC_RELAX_MKL) THEN
      CALL SCARC_VECTOR_COPY(V, Z, 1.0_EB, NL)                          !>  use additional auxiliary vector Z
      CALL SCARC_RELAXATION (Z, V, NS, NP, NL)                          !>  v := Relax(z)
   ELSE
      CALL SCARC_RELAXATION (V, V, NS, NP, NL)                          !>  v := Relax(v)
   ENDIF
#else
   CALL SCARC_RELAXATION (V, V, NS, NP, NL)                             !>  v := Relax(v)
#endif

   CALL SCARC_VECTOR_SUM      (V, X, OMEGA, 1.0_EB, NL)                 !>  x := omega * v + x

   CALL SCARC_MATVEC_PRODUCT  (X, V, NL)                                !>  v := A*x
   CALL SCARC_VECTOR_SUM      (B, V, 1.0_EB, -1.0_EB, NL)               !>  v := b - v

   IF (BL2NORM) THEN
      RES = SCARC_L2NORM (V, NL)                                        !>  res := ||v||
      IF (BVERBOSE) THEN
         NSTATE = SCARC_CONVERGENCE_STATE(NTYPE, NS, NL)
         IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP
      ENDIF
   ENDIF

ENDDO SMOOTH_LOOP

CALL SCARC_RELEASE_SOLVER(NS, NP)

CPU(MYID)%SMOOTHER = CPU(MYID)%SMOOTHER + CURRENT_TIME() - TNOW
END SUBROUTINE SCARC_SMOOTHER


!> ------------------------------------------------------------------------------------------------
!> Setup environement in every solver CALL (i.e. set pointers to used vectors) related to NSTACK
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SOLVER(NS, NP)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NS, NP                          !> references to current stack and parent

!>
!> if not first solver in stack, store last iteration parameters of parent solver NP
!>
IF (NP > 0) THEN
   SVP => STACK(NP)%SOLVER
   SVP%ITE   = ITE
   SVP%RES   = RES
   SVP%RESIN = RESIN
   SVP%ERR   = ERR
   SVP%CAPPA = CAPPA
ENDIF

!>
!> set new environment for solver on stack position NS
!>
SV => STACK(NS)%SOLVER

CNAME = SV%CNAME
ITE   = 0
NIT   = SV%NIT
EPS   = SV%EPS
OMEGA = SV%OMEGA
RESIN = SV%RESIN
CAPPA = -1.0

TYPE_PARENT = NP

TYPE_ACCURACY   = SV%TYPE_ACCURACY
TYPE_COARSE     = SV%TYPE_COARSE
TYPE_COARSENING = SV%TYPE_COARSENING
TYPE_CYCLING    = SV%TYPE_CYCLING
TYPE_GRID       = SV%TYPE_GRID
TYPE_EXCHANGE   = SV%TYPE_EXCHANGE
TYPE_INTERPOL   = SV%TYPE_INTERPOL
TYPE_KRYLOV     = SV%TYPE_KRYLOV
TYPE_LEVEL      = SV%TYPE_LEVEL
TYPE_MATRIX     = SV%TYPE_MATRIX
TYPE_METHOD     = SV%TYPE_METHOD
TYPE_MKL        = SV%TYPE_MKL
TYPE_MULTIGRID  = SV%TYPE_MULTIGRID
TYPE_PARENT     = SV%TYPE_PARENT
TYPE_MKL_PRECISION  = SV%TYPE_MKL_PRECISION
TYPE_PRECON     = SV%TYPE_PRECON
TYPE_RELAX      = SV%TYPE_RELAX
TYPE_SCOPE      = SV%TYPE_SCOPE
TYPE_SMOOTH     = SV%TYPE_SMOOTH
TYPE_SOLVER     = SV%TYPE_SOLVER
TYPE_STAGE      = SV%TYPE_STAGE
TYPE_STENCIL    = SV%TYPE_STENCIL
TYPE_TWOLEVEL   = SV%TYPE_TWOLEVEL
TYPE_VECTOR     = SV%TYPE_VECTOR

X = SV%X
B = SV%B
D = SV%D
R = SV%R
V = SV%V
Y = SV%Y
Z = SV%Z

#ifdef WITH_SCARC_DEBUG
E = SV%E
#endif

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) ITE_TOTAL = 0

#ifdef SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*)
#endif
END SUBROUTINE SCARC_SETUP_SOLVER


!> ------------------------------------------------------------------------------------------------
!> Reset SETTING of calling CURRENT-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELEASE_SOLVER(NS, NP)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN)  :: NS, NP                            !> references to current stack and parent

SV  => STACK(NS)%SOLVER

!>
!> store convergence information of preceding solver for FDS dump routine
!>
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   SCARC_CAPPA      = CAPPA
   SCARC_RESIDUAL   = RES
   SCARC_ITERATIONS = ITE
   !CALL SCARC_CLOSE_CSV_FILE()
ENDIF

SV%RESIN = RESIN
SV%RES   = RES
SV%ITE   = ITE
SV%CAPPA = CAPPA
SV%ERR   = ERR

!>
!> if not first solver in stack, reset environment of parent (calling) routine
!>
IF (NP > 0) THEN

   SVP => STACK(NP)%SOLVER

   ITE   = SVP%ITE
   NIT   = SVP%NIT
   EPS   = SVP%EPS
   RESIN = SVP%RESIN
   RES   = SVP%RES
   OMEGA = SVP%OMEGA
   CAPPA = SVP%CAPPA

   TYPE_ACCURACY   = SVP%TYPE_ACCURACY
   TYPE_COARSE     = SVP%TYPE_COARSE
   TYPE_COARSENING = SVP%TYPE_COARSENING
   TYPE_CYCLING    = SVP%TYPE_CYCLING
   TYPE_GRID       = SVP%TYPE_GRID
   TYPE_EXCHANGE   = SVP%TYPE_EXCHANGE
   TYPE_INTERPOL   = SVP%TYPE_INTERPOL
   TYPE_KRYLOV     = SVP%TYPE_KRYLOV
   TYPE_LEVEL      = SVP%TYPE_LEVEL
   TYPE_MATRIX     = SVP%TYPE_MATRIX
   TYPE_METHOD     = SVP%TYPE_METHOD
   TYPE_MKL        = SVP%TYPE_MKL
   TYPE_MULTIGRID  = SVP%TYPE_MULTIGRID
   TYPE_PARENT     = SVP%TYPE_PARENT
   TYPE_MKL_PRECISION  = SVP%TYPE_MKL_PRECISION
   TYPE_PRECON     = SVP%TYPE_PRECON
   TYPE_RELAX      = SVP%TYPE_RELAX
   TYPE_SCOPE      = SVP%TYPE_SCOPE
   TYPE_SMOOTH     = SVP%TYPE_SMOOTH
   TYPE_SOLVER     = SVP%TYPE_SOLVER
   TYPE_STAGE      = SVP%TYPE_STAGE
   TYPE_STENCIL    = SVP%TYPE_STENCIL
   TYPE_TWOLEVEL   = SVP%TYPE_TWOLEVEL
   TYPE_VECTOR     = SVP%TYPE_VECTOR

   X = SVP%X
   B = SVP%B
   D = SVP%D
   R = SVP%R
   V = SVP%V
   Y = SVP%Y
   Z = SVP%Z

#ifdef WITH_SCARC_DEBUG
   E = SVP%E
#endif

ENDIF

END SUBROUTINE SCARC_RELEASE_SOLVER


!> ----------------------------------------------------------------------------------------------------
!> Set initial solution corresponding to boundary data in BXS, BXF, ...
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WORKSPACE(NS, NL, NRHS)
USE SCARC_POINTERS, ONLY: M, L, F, G, SV, ST, STP, GWC, PRHS, HP
#ifdef WITH_SCARC_MGM
USE SCARC_POINTERS, ONLY: MGM
#endif
#ifdef WITH_SCARC_STANDALONE
USE SCARC_POINTERS, ONLY: P
#endif
INTEGER, INTENT(IN) :: NS, NL, NRHS
INTEGER :: NM, IW, IW1, IW2, IOR0, I, J, K, IC
REAL(EB) :: VAL

SV  => STACK(NS)%SOLVER

SELECT_SOLVER_TYPE: SELECT CASE (SV%TYPE_SOLVER)

   !>
   !> --------------- If used as main solver use values from pressure-routine as initialization
   !>
   CASE (NSCARC_SOLVER_MAIN)
   
      MAIN_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         !CALL SCARC_POINT_TO_GRID(NM, NL)
         CALL SCARC_POINT_TO_GRID(NM, NL)
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
   
         PRHS => M%PRHS
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF

#ifdef WITH_SCARC_STANDALONE
         IF (PREDICTOR) THEN
            P%H_OLD = P%H_NEW
         ELSE
            P%H_OLD = P%HS_NEW
         ENDIF
         P%B_OLD = ST%B
#endif
   
         !> get right hand side (PRHS from pres.f90) and initial vector (H or HS from last time step)
         SELECT_RHS_TYPE: SELECT CASE (NRHS)
   
            !> Solve original problem with inhomegeneous boundary conditions
            CASE (NSCARC_RHS_INHOMOGENEOUS)
      
               !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  ST%X(IC) = HP(G%ICX(IC), G%ICY(IC), G%ICZ(IC))        !> use last iterate as initial solution
                  ST%B(IC) = PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC))      !> get new RHS from surrounding code
               ENDDO                         
               !$OMP END PARALLEL DO
      
               !!$OMP PARALLEL 
               MAIN_INHOMOGENEOUS_LOOP: DO IOR0 = -3, 3, 1 
      
                  IF (IOR0 == 0) CYCLE
      
                  F => SCARC(NM)%LEVEL(NL)%FACE(IOR0)
                  
                  IW1 = F%NCW0
                  IW2 = F%NCW0 + F%NCW - 1
      
                  !!$OMP DO PRIVATE(IW, GWC, I, J, K, IC, VAL) SCHEDULE(STATIC)
                  FACE_INHOMOGENEOUS_LOOP: DO IW = IW1, IW2
      
                     GWC => G%WALL(IW)
         
                     I = GWC%IXW
                     J = GWC%IYW
                     K = GWC%IZW
         
                     IF (TWO_D .AND. J /= 1) CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_INDEX, SCARC_NONE, J)
         
                     IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
         
                     IC = G%CELL_NUMBER(I,J,K)
         
                     !> Dirichlet BC's:
                     !> these are based on the SETTING in BTYPE
                     !> in the structured case this corresponds to the face-wise SETTING according to the FFT
                     !> (this allows to use local FFT's as preconditioners)
                     !> in the unstructured case only open boundary cells lead to Dirichlet BC's
                     IF_DIRICHLET: IF (GWC%BTYPE == DIRICHLET) THEN
         
                        SELECT CASE (IOR0)
                        CASE (1)
                           VAL =  M%BXS(J,K)
                        CASE (-1)
                           VAL =  M%BXF(J,K)
                        CASE (2)
                           VAL =  M%BYS(I,K)
                        CASE (-2)
                           VAL =  M%BYF(I,K)
                        CASE (3)
                           VAL =  M%BZS(I,J)
                        CASE (-3)
                           VAL =  M%BZF(I,J)
                        END SELECT
         
                        ST%B(IC) = ST%B(IC) + F%SCAL_DIRICHLET * VAL
         
                     ENDIF IF_DIRICHLET
         
                     !> Neumann BC's:
                     !> Note for the unstructured case only:
                     !> Here, the matrix also contains Neumann BC's for those cells which have a
                     !> PRESSURE_BC_INDEX == DIRICHLET but are NOT open; these cells must be excluded below,
                     !> because BXS, BXF, ... contain the Dirichlet information from pres.f90 there;
                     !> excluding them corresponds to a homogeneous Neumann condition for these cells
                     IF_NEUMANN: IF (GWC%BTYPE == NEUMANN) THEN
         
                        IF (IS_UNSTRUCTURED .AND. M%WALL(IW)%PRESSURE_BC_INDEX /= NEUMANN) CYCLE
         
                        SELECT CASE (IOR0)
                        CASE (1)
                           VAL =  M%BXS(J,K)
                        CASE (-1)
                           VAL =  M%BXF(J,K)
                        CASE (2)
                           VAL =  M%BYS(I,K)
                        CASE (-2)
                           VAL =  M%BYF(I,K)
                        CASE (3)
                           VAL =  M%BZS(I,J)
                        CASE (-3)
                           VAL =  M%BZF(I,J)
                        END SELECT
         
                        ST%B(IC) = ST%B(IC) + F%SCAL_NEUMANN * VAL
         
                     ENDIF IF_NEUMANN
      
                  ENDDO FACE_INHOMOGENEOUS_LOOP
                  !!$OMP END DO

               ENDDO MAIN_INHOMOGENEOUS_LOOP
               !!$OMP END PARALLEL 
   
#ifdef WITH_SCARC_MGM

            !> Solve problem with homegeneous boundary conditions (MGM only)
            CASE (NSCARC_RHS_HOMOGENEOUS)
   
               ST%B = 0.0_EB                                    !> set RHS to zero
               ST%X = 0.0_EB                                    !> use zero as initial vector
   
               MGM => L%MGM

               !$OMP PARALLEL DO PRIVATE(IW, GWC, I, J, K, IOR0, IC) SCHEDULE(STATIC)
               DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
   
                  GWC => G%WALL(IW)
   
                  I = GWC%IXW
                  J = GWC%IYW
                  K = GWC%IZW
   
                  IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
   
                  IOR0 = GWC%IOR
                  IC   = G%CELL_NUMBER(I,J,K)
   
                  SELECT CASE (ABS(IOR0))
                     CASE(1)
                        ST%B(IC) = L%DXI * DTI * MGM%US(IW)
                     CASE(2)
                        ST%B(IC) = L%DYI * DTI * MGM%VS(IW)
                     CASE(3)
                        ST%B(IC) = L%DZI * DTI * MGM%WS(IW)
                  END SELECT
   
               ENDDO
               !$OMP END PARALLEL DO
#endif
   
         END SELECT SELECT_RHS_TYPE
 
#ifdef WITH_SCARC_STANDALONE
         P%B_NEW = ST%B
#endif
   
      ENDDO MAIN_MESHES_LOOP
      
      !> In case of a Krylov method clear overlapping parts of auxiliary vectors
      IF (IS_CG.OR.HAS_TWO_LEVELS) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            ST%D = 0.0_EB
            ST%R = 0.0_EB
            ST%V = 0.0_EB
            ST%Y = 0.0_EB
            ST%Z = 0.0_EB
         ENDDO
      ENDIF
   
      !> In case of a multigrid method as main solver clear
      !> overlapping parts of auxiliary vectors and coarse grid solver vectors
      IF (IS_GMG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            ST%V = 0.0_EB
            ST%Z = 0.0_EB
         ENDDO
      ENDIF
   
      !> In case of pure Neumann or periodic BCs, broadcast RHS(end) from last mesh
      !> to all and store it on all meshes
      IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
         IF (UPPER_MESH_INDEX == NMESHES) THEN
            L  => SCARC(NMESHES)%LEVEL(NL)
            ST => SCARC(NMESHES)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            MESH_REAL = ST%B(G%NC)
         ELSE
            MESH_REAL = 0.0_EB
         ENDIF
         IF (N_MPI_PROCESSES > 1) &
            CALL MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MESH_REAL, 1, MPI_DOUBLE_PRECISION,&
                               MPI_COMM_WORLD, IERROR)
         DO NM = 1, NMESHES
            SCARC(NM)%RHS_END = MESH_REAL(NMESHES)
         ENDDO
      ENDIF
   
   !>
   !> --------------- If MG is used as Krylov-preconditioner, vector G of main Krylov is the RHS for MG
   !>
   CASE (NSCARC_SOLVER_PRECON)
   
      IF (IS_CG_GMG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L   => SCARC(NM)%LEVEL(NL)
            ST  => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            STP => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)
            ST%X = 0.0_EB
            ST%B = STP%R
            ST%V = 0.0_EB
            ST%Z = 0.0_EB
         ENDDO
      ENDIF
   
   !>
   !> --------------- If used as coarse grid solver start with zero initialization
   !>
   CASE (NSCARC_SOLVER_COARSE)
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
         ST%X = 0.0_EB
         ST%D = 0.0_EB
         ST%R = 0.0_EB
         ST%V = 0.0_EB
         ST%Y = 0.0_EB
         ST%Z = 0.0_EB
      ENDDO
   
END SELECT SELECT_SOLVER_TYPE

END SUBROUTINE SCARC_SETUP_WORKSPACE


!> ------------------------------------------------------------------------------------------------
!> Check if solver converges or diverges and print out residual information
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CONVERGENCE_STATE(ISM, NS, NL)
INTEGER, INTENT(IN) :: NL, NS, ISM
INTEGER :: NSTATE

NSTATE = NSCARC_STATE_PROCEED

SELECT CASE (TYPE_ACCURACY)
   CASE (NSCARC_ACCURACY_RELATIVE)
      IF (RES <= RESIN*EPS .OR. RES <= NSCARC_THRESHOLD_CONVERGENCE) NSTATE = NSCARC_STATE_CONV
   CASE (NSCARC_ACCURACY_ABSOLUTE)
      IF (RES <= EPS .AND. RES <= RESIN) THEN
         IF (ITE == 0) THEN
            NSTATE = NSCARC_STATE_CONV0
         ELSE
            NSTATE = NSCARC_STATE_CONV
         ENDIF
         NIT = 0
      ENDIF
END SELECT
IF (RES > NSCARC_THRESHOLD_DIVGERGENCE) NSTATE = NSCARC_STATE_DIVG

SCARC_CONVERGENCE_STATE = NSTATE

IF (HAS_CSV_DUMP) CALL SCARC_DUMP_CSV(ISM, NS, NL)

#ifdef WITH_SCARC_VERBOSE
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) &
   WRITE(MSG%LU_VERBOSE,1100) STACK(NS)%SOLVER%CNAME, NL, ITE, RES
1100 FORMAT (A30,': Level=',i4,': Iteration = ',i4,': Residual =',e12.4)
#endif

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG, 1000) STACK(NS)%SOLVER%CNAME, NL, ITE, RES
1000 FORMAT (A30,': Level=',i4,': Iteration = ',i4,': Residual =',e25.16)
#endif

END FUNCTION SCARC_CONVERGENCE_STATE


!> ------------------------------------------------------------------------------------------------
!> Compute convergence rate and print out residual information for final loop
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)
INTEGER, INTENT(IN) :: NSTATE, NS, NL

IF (NSTATE == NSCARC_STATE_DIVG) THEN
   ITE   = - 1
   CAPPA = 1.0_EB
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

CALL SCARC_DUMP_CSV(0, NS, NL)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,2000) STACK(NS)%SOLVER%CNAME, ITE, CAPPA
#endif

#ifdef WITH_SCARC_VERBOSE
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) &
   WRITE(MSG%LU_VERBOSE,2000) STACK(NS)%SOLVER%CNAME, ITE, CAPPA
#endif

#if defined (WITH_SCARC_DEBUG) || defined (WITH_SCARC_VERBOSE)
2000 FORMAT (A30,': Iterations: ',i6,':   Convergence Rate =',e14.6,/)
#endif

END SUBROUTINE SCARC_CONVERGENCE_RATE


!> ------------------------------------------------------------------------------------------------
!> Perform restriction from finer to coarser grid in multigrid method
!>    - 'VF' corresponds to vector on fine   grid
!>    - 'VC' corresponds to vector on coarse grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTRICTION (NVB, NVC, NLF, NLC)
USE SCARC_POINTERS, ONLY: LF, LC, GF, GC, VF, VC
INTEGER, INTENT(IN) :: NVB, NVC, NLF, NLC
INTEGER :: NM
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-2:2,-2:2)=0
INTEGER :: IXC, IYC, IZC, ICC

!>
!> ------------------ Twolevel-CG or Geometric multigrid (as main solver or preconditioner) --------------
!>
IF (HAS_MULTI_LEVELS) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LF => SCARC(NM)%LEVEL(NLF)
      LC => SCARC(NM)%LEVEL(NLC)

      SELECT CASE(TYPE_GRID)
         CASE(NSCARC_GRID_STRUCTURED)
            GC => LC%STRUCTURED
            GF => LF%STRUCTURED
         CASE(NSCARC_GRID_UNSTRUCTURED)
            GC => LC%UNSTRUCTURED
            GF => LF%UNSTRUCTURED
      END SELECT

      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)
      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)

      IF (TWO_D) THEN

         SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)

            !> ---------- Constant Interpolation
            CASE (NSCARC_INTERPOL_CONSTANT)
   
               !$OMP PARALLEL DO PRIVATE(IXC, IZC, IXF, IZF, ICC, ICF) SCHEDULE(STATIC)
               DO IZC = 1, LC%NZ
                  DO IXC = 1, LC%NX
   
                     IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, 1, IZC)) CYCLE
   
                     IXF = 2*IXC
                     IZF = 2*IZC
   
                     ICC = GC%CELL_NUMBER(IXC, 1, IZC)
   
                     ICF(1) = GF%CELL_NUMBER(IXF-1, 1, IZF-1)
                     ICF(2) = GF%CELL_NUMBER(IXF-1, 1, IZF  )
                     ICF(3) = GF%CELL_NUMBER(IXF  , 1, IZF-1)
                     ICF(4) = GF%CELL_NUMBER(IXF  , 1, IZF  )
   
                     VC(ICC) = 0.25_EB * (  VF(ICF(1)) &
                                          + VF(ICF(2)) &
                                          + VF(ICF(3)) &
                                          + VF(ICF(4)) )
                  ENDDO
               ENDDO
               !$OMP END PARALLEL DO
   
            !> ---------- Bilinear Interpolation
            CASE (NSCARC_INTERPOL_BILINEAR)
   
               VC=0.0_EB
   
               !$OMP PARALLEL DO PRIVATE(IXC, IYC, IZC, IXF, IYF, IZF, ICC, ICFB) SCHEDULE(STATIC)
               DO IZC = 1, LC%NZ
                  DO IXC = 1, LC%NX
   
                     IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, 1, IZC)) CYCLE
   
                     IXF = 2*IXC
                     IZF = 2*IZC
   
                     ICC = GC%CELL_NUMBER(IXC, 1, IZC)
   
                     ICFB(-2,-2) = GF%CELL_NUMBER(IXF-2, 1, IZF-2)
                     ICFB(-1,-2) = GF%CELL_NUMBER(IXF-1, 1, IZF-2)
                     ICFB( 1,-2) = GF%CELL_NUMBER(IXF  , 1, IZF-2)
                     ICFB( 2,-2) = GF%CELL_NUMBER(IXF+1, 1, IZF-2)
   
                     ICFB(-2,-1) = GF%CELL_NUMBER(IXF-2, 1, IZF-1)
                     ICFB(-1,-1) = GF%CELL_NUMBER(IXF-1, 1, IZF-1)
                     ICFB( 1,-1) = GF%CELL_NUMBER(IXF  , 1, IZF-1)
                     ICFB( 2,-1) = GF%CELL_NUMBER(IXF+1, 1, IZF-1)
   
                     ICFB(-2, 1) = GF%CELL_NUMBER(IXF-2, 1, IZF)
                     ICFB(-1, 1) = GF%CELL_NUMBER(IXF-1, 1, IZF)
                     ICFB( 1, 1) = GF%CELL_NUMBER(IXF  , 1, IZF)
                     ICFB( 2, 1) = GF%CELL_NUMBER(IXF+1, 1, IZF)
   
                     ICFB(-2, 2) = GF%CELL_NUMBER(IXF-2, 1, IZF+1)
                     ICFB(-1, 2) = GF%CELL_NUMBER(IXF-1, 1, IZF+1)
                     ICFB( 1, 2) = GF%CELL_NUMBER(IXF  , 1, IZF+1)
                     ICFB( 2, 2) = GF%CELL_NUMBER(IXF+1, 1, IZF+1)
   
                     IF (IXC==1.AND.IZC==1) THEN
                        VC(ICC) = SCALR*( &
                           W4 *VF(ICFB(-1, 2)) + W3 *VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) + &
                           W12*VF(ICFB(-1, 1)) + W9 *VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) + &
                           W16*VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) + W4*VF(ICFB(2,-1)) )
                     ELSE IF (IXC==LC%NX.AND.IZC==  1) THEN
                        VC(ICC) = SCALR*( &
                           W1 *VF(ICFB(-2, 2)) + W3 *VF(ICFB(-1, 2)) + W4 *VF(ICFB(1, 2)) + &
                           W3 *VF(ICFB(-2, 1)) + W9 *VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) + &
                           W4 *VF(ICFB(-2,-1)) + W12*VF(ICFB(-1,-1)) + W16*VF(ICFB(1,-1)) )
                     ELSE IF (IXC==  1.AND.IZC==LC%NZ) THEN
                        VC(ICC) = SCALR*( &
                           W16*VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) + W4*VF(ICFB(2, 1)) + &
                           W12*VF(ICFB(-1,-1)) + W9 *VF(ICFB(1,-1)) + W3*VF(ICFB(2,-1)) + &
                           W4 *VF(ICFB(-1,-2)) + W3 *VF(ICFB(1,-2)) + W1*VF(ICFB(2,-2)) )
                     ELSE IF (IXC==LC%NX.AND.IZC==LC%NZ) THEN
                        VC(ICC) = SCALR*( &
                           W4 *VF(ICFB(-2, 1)) + W12*VF(ICFB(-1, 1)) + W16*VF(ICFB(1, 1)) + &
                           W3 *VF(ICFB(-2,-1)) + W9 *VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) + &
                           W1 *VF(ICFB(-2,-2)) + W3 *VF(ICFB(-1,-2)) + W4 *VF(ICFB(1,-2)) )
                     ELSE IF (IZC==  1) THEN
                        VC(ICC) = SCALR*( &
                           W1*VF(ICFB(-2, 2)) + W3 *VF(ICFB(-1, 2)) + W3 *VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) + &
                           W3*VF(ICFB(-2, 1)) + W9 *VF(ICFB(-1, 1)) + W9 *VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) + &
                           W4*VF(ICFB(-2,-1)) + W12*VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) + W4*VF(ICFB(2,-1)) )
                     ELSE IF (IZC==LC%NZ) THEN
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
                     ELSE IF (IXC==LC%NX) THEN
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
               !$OMP END PARALLEL DO
   
         END SELECT SELECT_INTERPOL

      !> ---------- Constant Interpolation (Note: 3D-bilinear case is still missing)
      ELSE

         !$OMP PARALLEL DO PRIVATE(IXC, IYC, IZC, IXF, IYF, IZF, ICC, ICF) SCHEDULE(STATIC)
         DO IZC = 1, LC%NZ
            DO IYC = 1, LC%NY
               DO IXC = 1, LC%NX

                  IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, IYC, IZC)) CYCLE

                  IXF = 2*IXC
                  IYF = 2*IYC
                  IZF = 2*IZC

                  ICC = GC%CELL_NUMBER(IXC, IYC, IZC)

                  ICF(1) = GF%CELL_NUMBER(IXF-1, IYF-1, IZF-1)
                  ICF(2) = GF%CELL_NUMBER(IXF-1, IYF-1, IZF  )
                  ICF(3) = GF%CELL_NUMBER(IXF-1, IYF  , IZF-1)
                  ICF(4) = GF%CELL_NUMBER(IXF-1, IYF  , IZF  )
                  ICF(5) = GF%CELL_NUMBER(IXF  , IYF-1, IZF-1)
                  ICF(6) = GF%CELL_NUMBER(IXF  , IYF-1, IZF  )
                  ICF(7) = GF%CELL_NUMBER(IXF  , IYF  , IZF-1)
                  ICF(8) = GF%CELL_NUMBER(IXF  , IYF  , IZF  )

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
         !$OMP END PARALLEL DO

      ENDIF
   ENDDO
ENDIF

END SUBROUTINE SCARC_RESTRICTION


!> ------------------------------------------------------------------------------------------------
!> Perform prolongation from coarser to finer grid
!>    - 'VC' corresponds to coarser grid
!>    - 'VF' corresponds to finer   grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PROLONGATION (NVC, NVB, NLC, NLF)
USE SCARC_POINTERS, ONLY: LF, LC, GF, GC, VF, VC
INTEGER, INTENT(IN) :: NVC, NVB, NLC, NLF
INTEGER :: NM, I
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-1:1,-1:1)=0
INTEGER :: IXC, IYC, IZC, ICC

!>
!> ------------------ Twolevel CG or Geometric Multigrid 
!>
IF (HAS_MULTI_LEVELS) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LC => SCARC(NM)%LEVEL(NLC)
      LF => SCARC(NM)%LEVEL(NLF)

      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)

      SELECT CASE(TYPE_GRID)
         CASE(NSCARC_GRID_STRUCTURED)
            GC => LC%STRUCTURED
            GF => LF%STRUCTURED
         CASE(NSCARC_GRID_UNSTRUCTURED)
            GC => LC%UNSTRUCTURED
            GF => LF%UNSTRUCTURED
      END SELECT

      IF (TWO_D) THEN

         SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)

            CASE (NSCARC_INTERPOL_CONSTANT)
   
               !>$OMP PARALLEL DO PRIVATE(IXC, IZC, IXF, IZF, ICC, ICF) SCHEDULE(STATIC)
               DO IZC = 1, LC%NZ
                  DO IXC = 1, LC%NX
   
                     IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, 1, IZC)) CYCLE
   
                     IXF = 2*IXC
                     IZF = 2*IZC
   
                     ICC = GC%CELL_NUMBER(IXC, 1, IZC)
   
                     ICF(1) = GF%CELL_NUMBER(IXF-1, 1, IZF-1)
                     ICF(2) = GF%CELL_NUMBER(IXF-1, 1, IZF  )
                     ICF(3) = GF%CELL_NUMBER(IXF  , 1, IZF-1)
                     ICF(4) = GF%CELL_NUMBER(IXF  , 1, IZF  )
   
                     DO I = 1, 4
                        VF(ICF(I)) = VC(ICC)
                     ENDDO
                  ENDDO
               ENDDO
               !>$OMP END PARALLEL DO 
   
            CASE (NSCARC_INTERPOL_BILINEAR)
   
               !>$OMP PARALLEL DO PRIVATE(IXC, IZC, IXF, IZF, ICC, ICFB) SCHEDULE(STATIC)
               DO IZC = 1, LC%NZ
                  DO IXC = 1, LC%NX
   
                     IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, 1, IZC)) CYCLE
   
                     IXF = 2*IXC
                     IZF = 2*IZC
   
                     ICC = GC%CELL_NUMBER(IXC, 1, IZC)
   
                     ICFB(-1,-1) = GF%CELL_NUMBER(IXF-1, 1, IZF-1)
                     ICFB(-1, 1) = GF%CELL_NUMBER(IXF-1, 1, IZF  )
                     ICFB( 1,-1) = GF%CELL_NUMBER(IXF  , 1, IZF-1)
                     ICFB( 1, 1) = GF%CELL_NUMBER(IXF  , 1, IZF  )
   
                     IF (IXC==1.AND.IZC==1) THEN
                        VF(ICFB(-1,-1)) = VC(ICC)
                        VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+LC%NX))
                        VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                        VF(ICFB( 1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX+1))
                     ELSE IF (IXC==1 .AND. IZC==LC%NZ) THEN
                        VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-LC%NX))
                        VF(ICFB(-1, 1)) = VC(ICC)
                        VF(ICFB( 1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX+1))
                        VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                     ELSE IF (IXC==LC%NX .AND. IZC==1) THEN
                        VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                        VF(ICFB(-1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX-1))
                        VF(ICFB( 1,-1)) = VC(ICC)
                        VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+LC%NX))
                     ELSE IF (IXC==LC%NX .AND. IZC==LC%NZ) THEN
                        VF(ICFB(-1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX-1))
                        VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                        VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-LC%NX))
                        VF(ICFB( 1, 1)) = VC(ICC)
                     ELSE IF (IZC==1) THEN
                        VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                        VF(ICFB(-1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX-1))
                        VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                        VF(ICFB( 1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX+1))
                     ELSE IF (IZC==LC%NZ) THEN
                        VF(ICFB(-1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX-1))
                        VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                        VF(ICFB( 1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX+1))
                        VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                     ELSE IF (IXC==1) THEN
                        VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-LC%NX))
                        VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+LC%NX))
                        VF(ICFB( 1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX+1))
                        VF(ICFB( 1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX+1))
                     ELSE IF (IXC==LC%NX) THEN
                        VF(ICFB(-1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX-1))
                        VF(ICFB(-1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX-1))
                        VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-LC%NX))
                        VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+LC%NX))
                     ELSE
                        VF(ICFB(-1,-1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX-1))
                        VF(ICFB(-1, 1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX-1))
                        VF(ICFB( 1,-1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX+1))
                        VF(ICFB( 1, 1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX+1))
                     ENDIF
                  ENDDO
               ENDDO
               !>$OMP END PARALLEL DO 
   
         END SELECT SELECT_INTERPOL

      ELSE

         !> ---------- Constant Interpolation (Note: 3D-bilinear case is still missing)
         !>$OMP PARALLEL DO PRIVATE(IXC, IYC, IZC, IXF, IYF, IZF, ICC, ICF) SCHEDULE(STATIC)
         DO IZC = 1, LC%NZ
            DO IYC = 1, LC%NY
               DO IXC = 1, LC%NX

                  IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, IYC, IZC)) CYCLE

                  IXF = 2*IXC
                  IYF = 2*IYC
                  IZF = 2*IZC

                  ICC = GC%CELL_NUMBER(IXC, IYC, IZC)

                  ICF(1) = GF%CELL_NUMBER(IXF-1, IYF-1, IZF-1)
                  ICF(2) = GF%CELL_NUMBER(IXF-1, IYF-1, IZF  )
                  ICF(3) = GF%CELL_NUMBER(IXF-1, IYF  , IZF-1)
                  ICF(4) = GF%CELL_NUMBER(IXF-1, IYF  , IZF  )
                  ICF(5) = GF%CELL_NUMBER(IXF  , IYF-1, IZF-1)
                  ICF(6) = GF%CELL_NUMBER(IXF  , IYF-1, IZF  )
                  ICF(7) = GF%CELL_NUMBER(IXF  , IYF  , IZF-1)
                  ICF(8) = GF%CELL_NUMBER(IXF  , IYF  , IZF  )

                  DO I = 1, 8
                     VF(ICF(I)) = VC(ICC)
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
         !>$OMP END PARALLEL DO 

      ENDIF
   ENDDO
ENDIF

END SUBROUTINE SCARC_PROLONGATION


!> ------------------------------------------------------------------------------------------------
!> Copy final solution from GMG (as preconditioner) to corresponding vector of CG (as main solver)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRECONDITIONER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%V = SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%X
ENDDO
END SUBROUTINE SCARC_UPDATE_PRECONDITIONER


!> ------------------------------------------------------------------------------------------------
!> Finalize data
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_MAINCELLS(NL)
USE SCARC_POINTERS, ONLY: M, G, ST, HP
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC !, I, J, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)
   ST => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

   !DO K = 1, M%KBAR
   !   DO J = 1, M%JBAR
   !      DO I = 1, M%IBAR
   !         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
   !         IC = G%CELL_NUMBER(I,J,K)
   !         HP(I, J, K) = ST%X(IC)
   !      ENDDO
   !   ENDDO
   !ENDDO

   !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
   DO IC = 1, G%NC
      HP (G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = ST%X(IC)
   ENDDO
   !>$OMP END PARALLEL DO 

ENDDO

!IF (MYID == 0) THEN
!WRITE(*,*) 'UPDATE_MAIN:'
!WRITE(*,'(5E25.16)') (((HP(I, J, K), I=5,9), J=1,4), K=1,4)
!ENDIF

END SUBROUTINE SCARC_UPDATE_MAINCELLS


!> ------------------------------------------------------------------------------------------------
!> Set correct boundary values at external and internal boundaries
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_GHOSTCELLS(NL)
USE SCARC_POINTERS, ONLY: M, L, G, GWC, HP
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IXG, IYG, IZG, IXW, IYW, IZW !, I, J, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

   !>
   !> Compute ghost cell values
   !>
   !$OMP PARALLEL DO SHARED(HP, M, L, G) PRIVATE(IW, IXG, IYG, IZG, IXW, IYW, IZW, IOR0, GWC) SCHEDULE(STATIC)
   WALL_CELLS_LOOP: DO IW = 1, L%N_WALL_CELLS_EXT

      GWC => G%WALL(IW)

      IXG = GWC%IXG
      IYG = GWC%IYG
      IZG = GWC%IZG

      IXW = GWC%IXW
      IYW = GWC%IYW
      IZW = GWC%IZW

      IOR0 = GWC%IOR

      SELECT CASE (IOR0)
         CASE ( 1)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXG,IYW,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXS(IYW,IZW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) - L%DX *M%BXS(IYW,IZW)
            ENDIF
         CASE (-1)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXG,IYW,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXF(IYW,IZW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) + L%DX *M%BXF(IYW,IZW)
            ENDIF
         CASE ( 2)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYS(IXW,IZW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) - L%DY *M%BYS(IXW,IZW)
            ENDIF
         CASE (-2)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYF(IXW,IZW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) + L%DY *M%BYF(IXW,IZW)
            ENDIF
         CASE ( 3)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZS(IXW,IYW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) - L%DZ *M%BZS(IXW,IYW)
            ENDIF
         CASE (-3)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZF(IXW,IYW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) + L%DZ *M%BZF(IXW,IYW)
            ENDIF
      END SELECT

   ENDDO WALL_CELLS_LOOP
   !$OMP END PARALLEL DO

#ifdef WITH_SCARC_DEBUG2
   WRITE(MSG%LU_DEBUG,*) 'Updating ghostcells '
   WRITE(MSG%LU_DEBUG,*) 'HP'
   DO IZG=0,M%KBP1
      DO IYG=0,M%JBP1
         WRITE(MSG%LU_DEBUG,'(8E12.4)') (HP(IXG,IYG,IZG), IXG=0,M%IBP1)
      ENDDO
   ENDDO
#endif

ENDDO

!> -----------------------------------------------------------------------------------------------
!> Perform data exchange to achieve consistency of ghost values along internal boundaries
!> Note: this is no longer necessary because MESH_EXCHANGE(5) is used after the call of ScaRC
!> -----------------------------------------------------------------------------------------------
CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_PRESSURE, NL)

!IF (MYID == 0) THEN
!WRITE(*,*) 'UPDATE_MAIN2:'
!WRITE(*,'(5E25.16)') (((HP(I, J, K), I=5,9), J=1,4), K=1,4)
!ENDIF

END SUBROUTINE SCARC_UPDATE_GHOSTCELLS


!> ------------------------------------------------------------------------------------------------
!>  Perform data exchange corresponding to requested exchange type (CALL receive and send-routines)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE (NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL

N_REQ = 0
TYPE_EXCHANGE = NTYPE

CALL SCARC_EXCHANGE_RECV(NL)
CALL SCARC_EXCHANGE_SEND(NL)

END SUBROUTINE SCARC_EXCHANGE


!> ------------------------------------------------------------------------------------------------
!>  Receive data from neighbors (corresponds to POST_RECEIVES)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_RECV (NL)
USE SCARC_POINTERS, ONLY: OS
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM

RECEIVE_MESH_INDEX: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   RECEIVE_NBR_INDEX: DO NOM = 1, NMESHES

      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_NBR_INDEX

      IF (.NOT. ARE_NEIGHBORS(NM, NOM)) CYCLE RECEIVE_NBR_INDEX

      CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)

      SELECT_EXCHANGE_TYPE: SELECT CASE (TYPE_EXCHANGE)

         !> 
         !> Exchange information about neighboring step size along internal boundary
         !> 
         CASE (NSCARC_EXCHANGE_SETUP)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT0(1), 1, MPI_INTEGER, SNODE, &
                           TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)

         !> 
         !> Exchange information about neighboring step size along internal boundary
         !> 
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT0(1), 3, MPI_INTEGER, SNODE, &
                           TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)

         !> 
         !> Exchange information about neighboring step size along internal boundary
         !> 
         CASE (NSCARC_EXCHANGE_CELL_WIDTH)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_REAL0(1), 6, MPI_DOUBLE_PRECISION, SNODE, &
                           TAG, MPI_COMM_WORLD,REQ(N_REQ), IERROR)

         !> 
         !> Exchange information about neighboring cell numbers 
         !> 
         CASE (NSCARC_EXCHANGE_CELL_NUMBER)

            N_REQ = N_REQ+1
            OS%RECV_INT1 = 0
            CALL MPI_IRECV(OS%RECV_INT1(1), SIZE(OS%RECV_INT1), MPI_INTEGER, SNODE, &
                           TAG, MPI_COMM_WORLD,REQ(N_REQ), IERROR)

         !> 
         !> Exchange information about overlapping stencils of neighboring matrices
         !> 
         CASE (NSCARC_EXCHANGE_MATRIX_STENCIL)

            N_REQ = N_REQ+1
            OS%RECV_REAL7 = 0
            CALL MPI_IRECV(OS%RECV_REAL7(1), SIZE(OS%RECV_REAL7), MPI_DOUBLE_PRECISION, SNODE, &
                           TAG, MPI_COMM_WORLD,REQ(N_REQ), IERROR)

         !> 
         !> Exchange information about mean values of overlapping parts
         !> 
         CASE (NSCARC_EXCHANGE_VECTOR_MEAN)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_REAL2(1), SIZE(OS%RECV_REAL2), MPI_DOUBLE_PRECISION, SNODE, &
                           TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)

         !> 
         !> Exchange information about pressure and general vectors as well as  matrix subdiagonals
         !> 
         CASE (NSCARC_EXCHANGE_PRESSURE, NSCARC_EXCHANGE_VECTOR, NSCARC_EXCHANGE_MATRIX_SUBDIAG)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_REAL1(1), SIZE(OS%RECV_REAL1), MPI_DOUBLE_PRECISION, SNODE, &
                           TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)

         CASE DEFAULT

            CALL SCARC_SHUTDOWN(NSCARC_ERROR_EXCHANGE_RECV, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SELECT_EXCHANGE_TYPE
   ENDDO RECEIVE_NBR_INDEX
ENDDO RECEIVE_MESH_INDEX

END SUBROUTINE SCARC_EXCHANGE_RECV


!> ------------------------------------------------------------------------------------------------
!> Different data exchanges between mesh neighbors corresponding to requested TYPE_EXCHANGE
!>
!> NSCARC_EXCHANGE_SETUP:           Setup of initial communication data
!> NSCARC_EXCHANGE_PRESSURE:        Pressure vector (H or HS) from predictor/corrector 
!> NSCARC_EXCHANGE_VECTOR:          Specified working vector for mat-vec product
!> NSCARC_EXCHANGE_VECTOR_MEAN:     Specified working vector for mat-vec product with summing up
!> NSCARC_EXCHANGE_CELL_WIDTH:      Cell width information of overlapping parts 
!> NSCARC_EXCHANGE_CELL_NUMBER:     Cell number information of overlapped parts 
!> NSCARC_EXCHANGE_MATRIX_SIZE:     Size of overlapped matrix
!> NSCARC_EXCHANGE_MATRIX_SUBDIAG:  Overlapping diagonal 'leg' of neighboring matrix
!> NSCARC_EXCHANGE_MATRIX_STENCIL:  Overlapping stencils of neighboring matrix
!>
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_SEND (NL)
USE SCARC_POINTERS, ONLY: OS, L, OL, G, OG, AC, VC, &
                          RECV_INT0, RECV_INT1, RECV_REAL0, RECV_REAL1, RECV_REAL7
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, NOM
INTEGER  :: IWG, ICOL, ICE, ICW, ICG, IOR0, ID, LL
INTEGER  :: IXW, IYW, IZW, IXG, IYG, IZG
REAL(EB) :: TNOW

!> 
!> ------------------------------ First part ------------------------------------------------------
!> Collect data for requested exchange type in corresponding SEND-buffer
!> 
TNOW = CURRENT_TIME()

MESHES_PACK_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   OMESHES_PACK_LOOP: DO NOM = 1, NMESHES

      IF (.NOT. ARE_NEIGHBORS(NM, NOM)) CYCLE OMESHES_PACK_LOOP
      CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)

      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      OMESHES_PACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

         !> 
         !> ---------- Pack size of exchange vector
         !> 
         CASE (NSCARC_EXCHANGE_SETUP)

            OS%SEND_INT0(1)=OG%NCG

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT0(1), 1, MPI_INTEGER, SNODE,&
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> 
         !> ---------- Pack face data pressure vector (H or HS) from predictor/corrector 
         !> 
         CASE (NSCARC_EXCHANGE_PRESSURE)

            IF (PREDICTOR) THEN
               HVC => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_H)
            ELSE
               HVC => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_HS)
            ENDIF

            DO IOR0 = -3, 3
               IF (OL%NCG(IOR0) == 0) CYCLE

               !>$OMP PARALLEL DO PRIVATE(ICG, IWG) SCHEDULE(STATIC)
               DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                  IWG = OG%ICG_TO_IWG(ICG)
                  OS%SEND_REAL1(ICG) = HVC(G%WALL(IWG)%IXW, G%WALL(IWG)%IYW, G%WALL(IWG)%IZW)
               ENDDO
               !>$OMP END PARALLEL DO 
            ENDDO

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL1(1), SIZE(OS%SEND_REAL1), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> 
         !> ---------- Pack face data of specified vector
         !> 
         CASE (NSCARC_EXCHANGE_VECTOR)

            VC => SCARC_POINT_TO_VECTOR(NM, NL, TYPE_VECTOR)

            IF (IS_STRUCTURED) THEN

               DO IOR0 = -3, 3
                  IF (OL%NCG(IOR0) == 0) CYCLE
                  !$OMP PARALLEL DO PRIVATE(ICG, ICW) SCHEDULE(STATIC)
                  DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                     ICW = OG%ICG_TO_ICW(ICG)
                     OS%SEND_REAL1(ICG) = VC(ICW)
                  ENDDO
                  !$OMP END PARALLEL DO 
               ENDDO

            ELSE

               DO IOR0 = -3, 3
                  IF (OL%NCG(IOR0) == 0) CYCLE
                  !$OMP PARALLEL DO PRIVATE(ICG, ICW) SCHEDULE(STATIC)
                  DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                     ICW = OG%ICG_TO_ICW(ICG)
                     IF (ICW > 0) THEN
                        OS%SEND_REAL1(ICG) = VC(ICW)
                     ELSE
                        OS%SEND_REAL1(ICG) = NSCARC_HUGE_REAL
                     ENDIF
                  ENDDO
                  !$OMP END PARALLEL DO 
               ENDDO

            ENDIF

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL1(1), SIZE(OS%SEND_REAL1), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> 
         !> ---------- Pack face data of specified vector
         !> 
         CASE (NSCARC_EXCHANGE_VECTOR_MEAN)

            VC => SCARC_POINT_TO_VECTOR(NM, NL, TYPE_VECTOR)

            IF (IS_STRUCTURED) THEN

               LL = 1
               DO IOR0 = -3, 3
                  IF (OL%NCG(IOR0) == 0) CYCLE
                  DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                     ICW = OG%ICG_TO_ICW(ICG)
                     OS%SEND_REAL2(LL) = VC(ICW)
                     LL = LL + 1
                  ENDDO
                  DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                     ICE = OG%ICG_TO_ICE(ICG)
                     OS%SEND_REAL2(LL) = VC(ICE)
                     LL = LL + 1
                  ENDDO
               ENDDO

            ELSE

               DO IOR0 = -3, 3
                  IF (OL%NCG(IOR0) == 0) CYCLE

                  !$OMP PARALLEL DO PRIVATE(ICG, ICW) SCHEDULE(STATIC)
                  DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                     ICE = OG%ICG_TO_ICE(ICG)
                     IF (ICW > 0) THEN
                        OS%SEND_REAL2(ICG) = VC(ICE)
                     ELSE
                        OS%SEND_REAL2(ICG) = NSCARC_HUGE_REAL
                     ENDIF
                  ENDDO
                  !$OMP END PARALLEL DO 
               ENDDO

            ENDIF

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL2(1), SIZE(OS%SEND_REAL2), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> 
         !> ---------- Pack cell width information of overlapping parts 
         !> 
         CASE (NSCARC_EXCHANGE_CELL_WIDTH)

            OS%SEND_REAL0(1) = L%DXL(0)
            OS%SEND_REAL0(2) = L%DXL(L%NX)
            OS%SEND_REAL0(3) = L%DYL(0)
            OS%SEND_REAL0(4) = L%DYL(L%NY)
            OS%SEND_REAL0(5) = L%DZL(0)
            OS%SEND_REAL0(6) = L%DZL(L%NZ)

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL0(1), 6, MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> 
         !> ---------- Pack cell number information of overlapped parts 
         !> 
         CASE (NSCARC_EXCHANGE_CELL_NUMBER)

            DO IOR0 = -3, 3

               IF (OL%NCG(IOR0) == 0) CYCLE

               !$OMP PARALLEL DO PRIVATE(ICG, IWG, IXW, IYW, IZW) SCHEDULE(STATIC)
               DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                  IWG = OG%ICG_TO_IWG(ICG)
                  IXW = G%WALL(IWG)%IXW
                  IYW = G%WALL(IWG)%IYW
                  IZW = G%WALL(IWG)%IZW
                  IF (L%IS_SOLID(IXW, IYW, IZW)) THEN
                     OS%SEND_INT1(ICG) = -G%CELL_NUMBER(IXW, IYW, IZW)
                  ELSE
                     OS%SEND_INT1(ICG) =  G%CELL_NUMBER(IXW, IYW, IZW)
                  ENDIF
               ENDDO
               !$OMP END PARALLEL DO 

            ENDDO

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT1(1),SIZE(OS%SEND_INT1),MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF

         !>
         !> ---------- Pack size of overlapping matrix parts
         !> 
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)
   
            SELECT CASE (SCARC_MATRIX_LEVEL(NL))
               CASE (NSCARC_MATRIX_BANDED)
                  OS%SEND_INT0(1) = OG%AB%N_VAL
                  OS%SEND_INT0(2) = OG%AB%N_DIAG
                  OS%SEND_INT0(3) = OG%AB%N_STENCIL
               CASE (NSCARC_MATRIX_COMPACT)
                  OS%SEND_INT0(1) = OG%AC%N_VAL
                  OS%SEND_INT0(2) = OG%AC%N_ROW
                  OS%SEND_INT0(3) = OG%AC%N_STENCIL
            END SELECT
   
            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT0(1), 3, MPI_INTEGER, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> 
         !> ---------- Pack overlapping subdiagonal 'leg' for neighboring banded matrix
         !> 
         CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG)

            SELECT CASE (SCARC_MATRIX_LEVEL(NL))

               !> Banded matrix case
               CASE (NSCARC_MATRIX_BANDED)                  !> Only structured discretization possible
                  AB => G%AB
                  DO IOR0 = -3, 3
                     IF (OL%NCG(IOR0) == 0) CYCLE
                     !$OMP PARALLEL DO PRIVATE(ICG, ICW, ICE, ICOL) SCHEDULE(STATIC)
                     DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                        ICW = OG%ICG_TO_ICW(ICG)
                        OS%SEND_REAL1(ICG) = L%FACE(IOR0)%SCAL_FACE
                     ENDDO
                     !$OMP END PARALLEL DO 
                  ENDDO

               !> Compact matrix case
               CASE (NSCARC_MATRIX_COMPACT)
                  AC => G%AC
                  IF (IS_STRUCTURED) THEN
      
                     DO IOR0 = -3, 3
                        IF (OL%NCG(IOR0) == 0) CYCLE
                        !$OMP PARALLEL DO PRIVATE(ICG, ICW, ICE, ICOL) SCHEDULE(STATIC)
                        DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                           ICW = OG%ICG_TO_ICW(ICG)
                           ICE = OG%ICG_TO_ICE(ICG)
                           DO ICOL = AC%ROW(ICW)+1, AC%ROW(ICW+1)-1
                              IF (AC%COL(ICOL) == ICE) THEN
                                 OS%SEND_REAL1(ICG) = AC%VAL(ICOL)
                                 EXIT
                              ENDIF
                           ENDDO
                        ENDDO
                        !$OMP END PARALLEL DO 
                     ENDDO
      
                  ELSE
      
                     DO IOR0 = -3, 3
                        IF (OL%NCG(IOR0) == 0) CYCLE
                        !$OMP PARALLEL DO PRIVATE(ICG, ICW, ICE, ICOL) SCHEDULE(STATIC)
                        DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                           ICW = OG%ICG_TO_ICW(ICG)
                           IF (ICW > 0) THEN
                              ICE = OG%ICG_TO_ICE(ICG)
                              DO ICOL = AC%ROW(ICW)+1, AC%ROW(ICW+1)-1
                                 IF (AC%COL(ICOL) == ICE) THEN
                                    OS%SEND_REAL1(ICG) = AC%VAL(ICOL)
                                    EXIT
                                 ENDIF
                              ENDDO
                           ELSE
                              OS%SEND_REAL1(ICG) = NSCARC_HUGE_REAL
                           ENDIF
                        ENDDO
                        !$OMP END PARALLEL DO 
                     ENDDO
      
                  ENDIF
            END SELECT
      
            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL1(1), SIZE(OS%SEND_REAL1), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> 
         !> ---------- Pack overlapping subdiagonal 'leg' for neighboring compact matrix
         !> 
         CASE (NSCARC_EXCHANGE_MATRIX_STENCIL)

            LL = 1
            SELECT CASE (SCARC_MATRIX_LEVEL(NL))

               !> Banded matrix case
               CASE (NSCARC_MATRIX_BANDED)                 !> Only structured discretization possible

                  AB => G%AB
                  DO IOR0 = -3, 3
                     IF (OL%NCG(IOR0) == 0) CYCLE
                     !$OMP PARALLEL DO PRIVATE(ICG, ICW, ICE, ICOL) SCHEDULE(STATIC)
                     DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                        ICW = OG%ICG_TO_ICW(ICG)
                        DO ID = 1, AB%N_STENCIL
                           OS%SEND_REAL7(LL) = AB%VAL(ICW, ID)
                           LL = LL + 1
                        ENDDO
                     ENDDO
                     !$OMP END PARALLEL DO 
                  ENDDO
      
 
               !> Compact matrix case
               CASE (NSCARC_MATRIX_COMPACT)

                  AC => G%AC
                  IF (IS_STRUCTURED) THEN
      
                     DO IOR0 = -3, 3
                        IF (OL%NCG(IOR0) == 0) CYCLE
                        !$OMP PARALLEL DO PRIVATE(ICG, ICW, ICE, ICOL) SCHEDULE(STATIC)
                        DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                           ICW = OG%ICG_TO_ICW(ICG)
                           DO ICOL = AC%ROW(ICW), AC%ROW(ICW+1)-1
                              OS%SEND_REAL7(LL) = AC%VAL(ICOL)
                              LL = LL + 1
                           ENDDO
                        ENDDO
                        !$OMP END PARALLEL DO 
                     ENDDO
      
                  ELSE
      
                     DO IOR0 = -3, 3
                        IF (OL%NCG(IOR0) == 0) CYCLE

                        !$OMP PARALLEL DO PRIVATE(ICG, ICW, ICE, ICOL) SCHEDULE(STATIC)
                        DO ICG= OL%NCG0(IOR0), OL%NCG(IOR0)
                           ICW = OG%ICG_TO_ICW(ICG)
                           IF (ICW > 0) THEN
                              DO ICOL = AC%ROW(ICW), AC%ROW(ICW+1)-1
                                 OS%SEND_REAL7(LL) = AC%VAL(ICOL)
                                 LL = LL + 1
                              ENDDO
                           ELSE
                              OS%SEND_REAL7(LL) = NSCARC_HUGE_REAL
                              LL = LL + 1
                           ENDIF
                        ENDDO
                        !$OMP END PARALLEL DO 
                     ENDDO
      
                  ENDIF
            END SELECT

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL7(1), SIZE(OS%SEND_REAL7), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         CASE DEFAULT

            CALL SCARC_SHUTDOWN(NSCARC_ERROR_EXCHANGE_SEND, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT OMESHES_PACK_SELECT
   ENDDO OMESHES_PACK_LOOP
ENDDO MESHES_PACK_LOOP

CPU(MYID)%BUFFER_PACKING = CPU(MYID)%BUFFER_PACKING   + CURRENT_TIME() - TNOW

!> 
!> ------------------------------ Second part -----------------------------------------------------
!> Wait for all meshes to have sent and received their data
!> 
IF (N_MPI_PROCESSES>1.AND.N_REQ/=0) CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),MPI_STATUSES_IGNORE,IERROR)

!> 
!> ------------------------------ Third part ------------------------------------------------------
!> Extract exchanged data from corresponding RECEIVE-buffers
!> 
TNOW = CURRENT_TIME()
MESHES_UNPACK_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   OMESHES_UNPACK_LOOP: DO NOM=1,NMESHES

      SNODE  = PROCESS(NM)
      RNODE  = PROCESS(NOM)

      IF (.NOT. ARE_NEIGHBORS(NM, NOM)) CYCLE 
      CALL SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)

      OMESHES_UNPACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

         !> 
         !> ---------- Unpack information about neighboring exchange size
         !> 
         CASE (NSCARC_EXCHANGE_SETUP)

            CALL SCARC_POINT_TO_BUFFER (NM, NOM, RNODE, SNODE, NSCARC_BUFFER_INT0)
            OG%NCG = RECV_INT0(1)

         !> 
         !> ---------- Unpack information about neighboring cell widths
         !> 
         CASE (NSCARC_EXCHANGE_CELL_WIDTH)

            CALL SCARC_POINT_TO_BUFFER (NM, NOM, RNODE, SNODE, NSCARC_BUFFER_REAL0)

            IF (OL%NCG(-1) /= 0) L%DXL(0)    = 0.5_EB*(RECV_REAL0(1) + L%DXL(0))
            IF (OL%NCG( 1) /= 0) L%DXL(L%NX) = 0.5_EB*(RECV_REAL0(2) + L%DXL(L%NX))
            IF (OL%NCG(-2) /= 0) L%DYL(0)    = 0.5_EB*(RECV_REAL0(3) + L%DYL(0))
            IF (OL%NCG( 2) /= 0) L%DYL(L%NY) = 0.5_EB*(RECV_REAL0(4) + L%DYL(L%NY))
            IF (OL%NCG(-3) /= 0) L%DZL(0)    = 0.5_EB*(RECV_REAL0(5) + L%DZL(0))
            IF (OL%NCG( 3) /= 0) L%DZL(L%NZ) = 0.5_EB*(RECV_REAL0(6) + L%DZL(L%NZ))
             
         !> 
         !> ---------- Unpack information about neighboring cell numbers
         !> 
         CASE (NSCARC_EXCHANGE_CELL_NUMBER)

            CALL SCARC_POINT_TO_BUFFER (NM, NOM, RNODE, SNODE, NSCARC_BUFFER_INT1)

            LL = 1
            DO IOR0 = -3, 3
               IF (OL%NCG(IOR0) == 0) CYCLE

               !$OMP PARALLEL DO PRIVATE(ICG, IWG, ICE, IXG, IYG, IZG) SCHEDULE(STATIC)
               DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
   
                  ICE = OG%ICG_TO_ICE(ICG)
                  G%ICE_TO_ICN(ICE) = RECV_INT1(LL)

                  IWG = OG%ICG_TO_IWG(ICG)
                  IXG = G%WALL(IWG)%IXG
                  IYG = G%WALL(IWG)%IYG
                  IZG = G%WALL(IWG)%IZG
   
                  IF (RECV_INT1(ICG) < 0) THEN
                     L%IS_SOLID(IXG, IYG, IZG) = .TRUE.
                  ELSE
                     L%IS_SOLID(IXG, IYG, IZG) = .FALSE.
                  ENDIF

                  LL = LL + 1
               ENDDO
               !$OMP END PARALLEL DO 
            ENDDO


         !> 
         !> ---------- Unpack overlapping parts of pressure vectors H or HS from predictor/corrector
         !> 
         CASE (NSCARC_EXCHANGE_PRESSURE)

            CALL SCARC_POINT_TO_BUFFER (NM, NOM, RNODE, SNODE, NSCARC_BUFFER_REAL1)

            IF (PREDICTOR) THEN
               HVC => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_H)
            ELSE
               HVC => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_HS)
            ENDIF

            LL = 1
            DO IOR0 = -3, 3
               IF (OL%NCG(IOR0) == 0) CYCLE
               UNPACK_PRESSURE: DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                  IWG = OG%ICG_TO_IWG(ICG)
                  IXG=G%WALL(IWG)%IXG
                  IYG=G%WALL(IWG)%IYG
                  IZG=G%WALL(IWG)%IZG
                  HVC(IXG, IYG, IZG) = RECV_REAL1(LL)
                  LL = LL + 1
               ENDDO UNPACK_PRESSURE
            ENDDO

         !>
         !> ---------- Unpack overlapping parts of a specified vector and overwrite own positions
         !> 
         CASE (NSCARC_EXCHANGE_VECTOR)

            CALL SCARC_POINT_TO_BUFFER (NM, NOM, RNODE, SNODE, NSCARC_BUFFER_REAL1)

            VC => SCARC_POINT_TO_VECTOR (NM, NL, TYPE_VECTOR)

            LL = 1
            DO IOR0 = -3, 3
               IF (OL%NCG(IOR0) == 0) CYCLE
               DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                  VC(OG%ICG_TO_ICE(ICG)) = RECV_REAL1(LL)
                  LL = LL + 1
               ENDDO
            ENDDO

         !>
         !> ---------- Unpack overlapping parts of a specified vector and build mean values with own positions
         !> 
         CASE (NSCARC_EXCHANGE_VECTOR_MEAN)

            CALL SCARC_POINT_TO_BUFFER (NM, NOM, RNODE, SNODE, NSCARC_BUFFER_REAL2)

            VC => SCARC_POINT_TO_VECTOR (NM, NL, TYPE_VECTOR)

            LL = 1
            DO IOR0 = -3, 3
               IF (OL%NCG(IOR0) == 0) CYCLE

               DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                  !VC(OG%ICG_TO_ICE(ICG)) = 0.5_EB * (VC(OG%ICG_TO_ICE(ICG)) + RECV_REAL2(LL))
                  VC(OG%ICG_TO_ICE(ICG)) = 1.0_EB/3.0_EB * VC(OG%ICG_TO_ICE(ICG)) + 2.0_EB/3.0_EB * RECV_REAL2(LL)
                  LL = LL + 1
               ENDDO
               DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                  !VC(OG%ICG_TO_ICW(ICG)) = 0.5_EB * (VC(OG%ICG_TO_ICW(ICG)) + RECV_REAL2(LL))
                  VC(OG%ICG_TO_ICW(ICG)) = 2.0_EB/3.0_EB * VC(OG%ICG_TO_ICW(ICG)) + 1.0_EB/3.0_EB * RECV_REAL2(LL)
                  LL = LL + 1
               ENDDO
            ENDDO

         !> 
         !> ---------- Unpack information about neighboring matrix sizes
         !> 
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            CALL SCARC_POINT_TO_BUFFER (NM, NOM, RNODE, SNODE, NSCARC_BUFFER_INT0)

            SELECT CASE (SCARC_MATRIX_LEVEL(NL))
               CASE (NSCARC_MATRIX_BANDED)
                  OG%AB%N_VAL     = RECV_INT0(1)
                  OG%AB%N_DIAG    = RECV_INT0(2)
                  OG%AB%N_STENCIL = RECV_INT0(3)
               CASE (NSCARC_MATRIX_COMPACT)
                  OG%AC%N_VAL     = RECV_INT0(1)
                  OG%AC%N_ROW     = RECV_INT0(2)
                  OG%AC%N_STENCIL = RECV_INT0(3)
            END SELECT

         !> 
         !> ---------- Unpack information about neighboring subdiagonal overlapping leg
         !> 
         CASE (NSCARC_EXCHANGE_MATRIX_SUBDIAG)

            CALL SCARC_POINT_TO_BUFFER (NM, NOM, RNODE, SNODE, NSCARC_BUFFER_REAL1)
            SELECT CASE (SCARC_MATRIX_LEVEL(NL))

               !> Banded matrix case
               CASE (NSCARC_MATRIX_BANDED)                 !> Only structured discretization possible

                  LL = 1
                  DO IOR0 = -3, 3
                     IF (OL%NCG(IOR0) == 0) CYCLE
                     DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                        ICE = OG%ICG_TO_ICE(ICG)
                        G%ICE_TO_VAL(ICE) = RECV_REAL1(LL)
                        LL = LL + 1
                     ENDDO
                  ENDDO

               CASE (NSCARC_MATRIX_COMPACT)

                  LL = 1
                  DO IOR0 = -3, 3
                     IF (OL%NCG(IOR0) == 0) CYCLE
                     DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                        ICE = OG%ICG_TO_ICE(ICG)
                        G%ICE_TO_VAL(ICE) = RECV_REAL1(LL)
                        LL = LL + 1
                     ENDDO
                  ENDDO

            END SELECT

         !> 
         !> ---------- Unpack information about neighboring matrix stencils
         !> 
         CASE (NSCARC_EXCHANGE_MATRIX_STENCIL)

            CALL SCARC_POINT_TO_BUFFER (NM, NOM, RNODE, SNODE, NSCARC_BUFFER_REAL7)
            SELECT CASE (SCARC_MATRIX_LEVEL(NL))

               !> Banded matrix case
               CASE (NSCARC_MATRIX_BANDED)                 !> Only structured discretization possible

                  LL = 1
                  DO IOR0 = -3, 3
                     IF (OL%NCG(IOR0) == 0) CYCLE
                     DO ICG = OL%NCG0(IOR0), OL%NCG(IOR0)
                        DO ID = 1, OG%AB%N_STENCIL
                           OG%AB%VAL(ICG, ID) = RECV_REAL7(LL)
                           LL = LL + 1
                        ENDDO
                     ENDDO
                  ENDDO

               CASE (NSCARC_MATRIX_COMPACT)

                  LL = 1
                  !WRITE(*,*) 'UNPACK EXCHANGE_MATRIX_STENCIL, COMPACT, CURRENTLY NOT NEEDED !!!'
                  !IF (IS_STRUCTURED) THEN
                  !  ...
                  !ELSE
                  !   ...
                  !ENDIF

            END SELECT

         CASE DEFAULT

            CALL SCARC_SHUTDOWN(NSCARC_ERROR_EXCHANGE_SEND, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT OMESHES_UNPACK_SELECT
   ENDDO OMESHES_UNPACK_LOOP
ENDDO MESHES_UNPACK_LOOP
CPU(MYID)%BUFFER_UNPACKING = CPU(MYID)%BUFFER_UNPACKING   + CURRENT_TIME() - TNOW

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


!> ----------------------------------------------------------------------------------------------------
!> Filter out mean value
!> --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_FILTER_MEANVALUE(NV, NL)
USE SCARC_POINTERS, ONLY: L, G, VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM, IC, I, J, K

MESH_REAL = 0.0_EB
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)
   VC => SCARC_POINT_TO_VECTOR(NM, NL, NV)

   DO IC = 1, G%NC
      MESH_REAL(NM) = MESH_REAL(NM) + VC(IC)
   ENDDO

ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,MESH_REAL,COUNTS,DISPLS,&
                       MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

GLOBAL_REAL = SUM(MESH_REAL(1:NMESHES))/REAL(NC_GLOBAL(NL))

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)
   VC  => SCARC_POINT_TO_VECTOR(NM, NL, NV)

   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
            IC = G%CELL_NUMBER(I,J,K)
            VC(IC) = VC(IC) - GLOBAL_REAL
         ENDDO
      ENDDO
   ENDDO

ENDDO

END SUBROUTINE SCARC_FILTER_MEANVALUE


!> ------------------------------------------------------------------------------------------------
!> Restore last cell of last mesh
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_LAST_CELL (XX, NL)
USE SCARC_POINTERS, ONLY: S, VC
INTEGER, INTENT(IN) :: XX, NL

IF (UPPER_MESH_INDEX /= NMESHES .OR. TYPE_RELAX == NSCARC_RELAX_FFT) RETURN
S => SCARC(UPPER_MESH_INDEX)

VC => SCARC_POINT_TO_VECTOR (UPPER_MESH_INDEX, NL, XX)
VC(S%NC) = S%RHS_END

END SUBROUTINE SCARC_RESTORE_LAST_CELL


!> ================================================================================================
!>  Set of different pointer functions which return pointers to specified structures
!> ================================================================================================

!> -----------------------------------------------------------------------------
!> Point to mesh
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_MESH(NM)
USE SCARC_POINTERS, ONLY: M, S
INTEGER, INTENT(IN) :: NM
M => MESHES(NM)
S => SCARC(NM)
END SUBROUTINE SCARC_POINT_TO_MESH

!> -----------------------------------------------------------------------------
!> This routine sets the requested combination of mesh and level 
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_LEVEL(NM, NL)
USE SCARC_POINTERS, ONLY: M, S, L
INTEGER, INTENT(IN) :: NM, NL
M => MESHES(NM)
S => SCARC(NM)
L => S%LEVEL(NL)
END SUBROUTINE SCARC_POINT_TO_LEVEL


!> -----------------------------------------------------------------------------
!> This routine sets the requested combination of mesh and level 
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_NEIGHBOR(NM, NOM, NL)
USE SCARC_POINTERS, ONLY : OS, OL, OG
INTEGER, INTENT(IN) :: NM, NOM, NL

OS => SCARC(NM)%OSCARC(NOM)
OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      OG => OL%STRUCTURED
   CASE (NSCARC_GRID_UNSTRUCTURED)
      OG => OL%UNSTRUCTURED
END SELECT

END SUBROUTINE SCARC_POINT_TO_NEIGHBOR


!> ----------------------------------------------------------------------------------------------------
!> This routine sets the requested combination of mesh, level and discretization
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_GRID(NM, NL)
USE SCARC_POINTERS, ONLY: M, S, L, G, W, P
INTEGER, INTENT(IN) ::  NM, NL

M => MESHES(NM)
S => SCARC(NM)
L => S%LEVEL(NL)
SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      G => L%STRUCTURED
      G%NW = L%N_WALL_CELLS_EXT 
   CASE (NSCARC_GRID_UNSTRUCTURED)
      G => L%UNSTRUCTURED
      G%NW = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
END SELECT
W => G%WALL
P => SCARC(NM)%LEVEL(NL)%PRES

END SUBROUTINE SCARC_POINT_TO_GRID


!> -----------------------------------------------------------------------------
!> Point to requested receive buffer
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_BUFFER(NM, NOM, RNODE, SNODE, ITYPE)
USE SCARC_POINTERS, ONLY: RECV_INT0, RECV_INT1, RECV_REAL0, RECV_REAL1, RECV_REAL2, RECV_REAL7
INTEGER, INTENT(IN) ::  NM, NOM, RNODE, SNODE, ITYPE

SELECT CASE (ITYPE)
   CASE (NSCARC_BUFFER_INT0)
      IF (RNODE/=SNODE) THEN
         RECV_INT0  => SCARC(NM)%OSCARC(NOM)%RECV_INT0
      ELSE
         RECV_INT0  => SCARC(NOM)%OSCARC(NM)%SEND_INT0
      ENDIF
   CASE (NSCARC_BUFFER_INT1)
      IF (RNODE/=SNODE) THEN
         RECV_INT1  => SCARC(NM)%OSCARC(NOM)%RECV_INT1
      ELSE
         RECV_INT1  => SCARC(NOM)%OSCARC(NM)%SEND_INT1
      ENDIF
   CASE (NSCARC_BUFFER_REAL0)
      IF (RNODE/=SNODE) THEN
         RECV_REAL0 => SCARC(NM)%OSCARC(NOM)%RECV_REAL0
      ELSE
         RECV_REAL0 => SCARC(NOM)%OSCARC(NM)%SEND_REAL0
      ENDIF
   CASE (NSCARC_BUFFER_REAL1)
      IF (RNODE/=SNODE) THEN
         RECV_REAL1 => SCARC(NM)%OSCARC(NOM)%RECV_REAL1
      ELSE
         RECV_REAL1 => SCARC(NOM)%OSCARC(NM)%SEND_REAL1
      ENDIF
   CASE (NSCARC_BUFFER_REAL2)
      IF (RNODE/=SNODE) THEN
         RECV_REAL2 => SCARC(NM)%OSCARC(NOM)%RECV_REAL2
      ELSE
         RECV_REAL2 => SCARC(NOM)%OSCARC(NM)%SEND_REAL2
      ENDIF
   CASE (NSCARC_BUFFER_REAL7)
      IF (RNODE/=SNODE) THEN
         RECV_REAL7 => SCARC(NM)%OSCARC(NOM)%RECV_REAL7
      ELSE
         RECV_REAL7 => SCARC(NOM)%OSCARC(NM)%SEND_REAL7
      ENDIF
END SELECT

END SUBROUTINE SCARC_POINT_TO_BUFFER


!> -----------------------------------------------------------------------------
!> Point to exchange structure
!> -----------------------------------------------------------------------------
LOGICAL FUNCTION ARE_NEIGHBORS(NM, NOM)
USE SCARC_POINTERS, ONLY: OM
INTEGER, INTENT(IN) :: NM, NOM

ARE_NEIGHBORS = .TRUE.

OM => MESHES(NM)%OMESH(NOM)
IF (OM%NIC_R == 0 .AND. OM%NIC_S == 0) ARE_NEIGHBORS = .FALSE.

END FUNCTION ARE_NEIGHBORS

!> ------------------------------------------------------------------------------------------------
!> Set vector pointer on given combination of mesh and level
!> ------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_VECTOR (NM, NL, NV)
USE SCARC_POINTERS, ONLY: L
REAL(EB), POINTER, DIMENSION(:) :: SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NM, NL, NV

L => SCARC(NM)%LEVEL(NL)
SELECT CASE (NV)

   !>
   !> Stage one vectors
   !>
   CASE (NSCARC_VECTOR_ONE_X)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%X
   CASE (NSCARC_VECTOR_ONE_B)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%B
   CASE (NSCARC_VECTOR_ONE_D)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%D
   CASE (NSCARC_VECTOR_ONE_R)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%R
   CASE (NSCARC_VECTOR_ONE_V)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%V
   CASE (NSCARC_VECTOR_ONE_Y)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%Y
   CASE (NSCARC_VECTOR_ONE_Z)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%Z
#ifdef WITH_SCARC_DEBUG
   CASE (NSCARC_VECTOR_ONE_E)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%E
#endif

   !>
   !> Stage two vectors
   !>
   CASE (NSCARC_VECTOR_TWO_X)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%X
   CASE (NSCARC_VECTOR_TWO_B)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%B
   CASE (NSCARC_VECTOR_TWO_D)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%D
   CASE (NSCARC_VECTOR_TWO_R)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%R
   CASE (NSCARC_VECTOR_TWO_V)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%V
   CASE (NSCARC_VECTOR_TWO_Y)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%Y
   CASE (NSCARC_VECTOR_TWO_Z)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%Z
#ifdef WITH_SCARC_DEBUG
   CASE (NSCARC_VECTOR_TWO_E)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%E
#endif
END SELECT

END FUNCTION SCARC_POINT_TO_VECTOR


!> ------------------------------------------------------------------------------------------------
!> Single precision: Set vector pointer on given combination of mesh and level
!> ------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_VECTOR_FB(NM, NL, NV)
USE SCARC_POINTERS, ONLY: L
REAL(FB), POINTER, DIMENSION(:) :: SCARC_POINT_TO_VECTOR_FB
INTEGER, INTENT(IN) :: NM, NL, NV

L => SCARC(NM)%LEVEL(NL)
SELECT CASE (NV)
   CASE (NSCARC_VECTOR_ONE_X)
      SCARC_POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%X_FB
   CASE (NSCARC_VECTOR_ONE_B)
      SCARC_POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%B_FB
   CASE (NSCARC_VECTOR_ONE_R)
      SCARC_POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%R_FB
   CASE (NSCARC_VECTOR_ONE_V)
      SCARC_POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%V_FB
END SELECT

END FUNCTION SCARC_POINT_TO_VECTOR_FB


!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for banded storage technique
!> ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_HVECTOR(NM, NV)
REAL(EB), POINTER, DIMENSION(:,:,:) :: POINT_TO_HVECTOR
INTEGER, INTENT(IN) :: NM, NV
SELECT CASE (NV)
   CASE (NSCARC_VECTOR_H)
      POINT_TO_HVECTOR => MESHES(NM)%H
   CASE (NSCARC_VECTOR_HS)
      POINT_TO_HVECTOR => MESHES(NM)%HS
END SELECT
END FUNCTION POINT_TO_HVECTOR



!> ================================================================================================
!> Routines for different vector-operations, based on use of OpenMP
!> ================================================================================================
!> ------------------------------------------------------------------------------------------------
!> Vector multiplied with a constant scalar is added to another vector 
!>     DY(I) = DA * DX(I) + DY(I) 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DAXPY_CONSTANT(N, DA, DX, DY)
REAL(EB), INTENT(IN):: DA
REAL(EB), INTENT(IN), DIMENSION(:):: DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) = DY(I) + DA * DX(I)
ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE SCARC_DAXPY_CONSTANT


!> ------------------------------------------------------------------------------------------------
!> Vector multiplied with a constant scalar is added to vector multiplied with another scalar
!>     DY(I) = DA1 * DX(I) + DA2 * DY(I) 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DAXPY_CONSTANT_DOUBLE(N, DA1, DX, DA2, DY)
REAL(EB), INTENT(IN):: DA1, DA2
REAL(EB), INTENT(IN), DIMENSION(:):: DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) = DA1 * DX(I) + DA2 * DY(I)
ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE SCARC_DAXPY_CONSTANT_DOUBLE


!> ------------------------------------------------------------------------------------------------
!> Vector multiplied with variable scalars (componentwise) is added to another vector 
!>     DY(I) = DA(I)*DX(I) + DY(I)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DAXPY_VARIABLE(N, DA, DX, DY)
REAL(EB), INTENT(IN), DIMENSION(:):: DA, DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) = DY(I) + DA(I) * DX(I)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A, I6, 3E14.6)') 'I, DX, DA, DY:', I, DX(I), DA(I), DY(I)
#endif
ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE SCARC_DAXPY_VARIABLE


!> ------------------------------------------------------------------------------------------------
!> Vector is multiplied with a constant scalar 
!>     DY(I) = DA(I)*DX(I) 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SCALING_CONSTANT(N, SCAL, DX, DY)
REAL(EB), INTENT(IN), DIMENSION(:):: DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
REAL(EB), INTENT(IN) :: SCAL
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) =  SCAL * DX(I)
ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE SCARC_SCALING_CONSTANT


!> ------------------------------------------------------------------------------------------------
!> Vector is multiplied with variable scalars (componentwise)
!>     DY(I) = DA(I)*DX(I) 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SCALING_VARIABLE(N, DA, DX, DY)
REAL(EB), INTENT(IN), DIMENSION(:):: DA, DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) = DA(I) * DX(I)
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE SCARC_SCALING_VARIABLE



!> ================================================================================================
!> Collection of routines for the (de)allocation of different working spaces
!> ================================================================================================
!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT1(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT1', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
   END SELECT
#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
ENDIF
#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating INT1     array  ',A20,' in length (',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_INT1


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 2
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT2(WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT2', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
   END SELECT
#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating INT2     array  ',A20,' in length (',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_INT2


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 3
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT3', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
   END SELECT
#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2, NL3, NR3
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
      SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating INT3     array  ',A20,' in length (',I8,':',I8,' , ',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_INT3


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize Logical array of dimension 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_LOG1(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT1', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_TRUE)
         WORKSPACE = .TRUE.
      CASE (NSCARC_INIT_FALSE)
         WORKSPACE = .FALSE.
   END SELECT
#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating LOG3     array  ',A20,' in length (',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_LOG1


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize Logical array of dimension 2
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_LOG2(WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT3', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_TRUE)
         WORKSPACE = .TRUE.
      CASE (NSCARC_INIT_FALSE)
         WORKSPACE = .FALSE.
   END SELECT
#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating LOG3     array  ',A20,' in length (',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_LOG2


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize Logical array of dimension 3
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_LOG3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT3', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_TRUE)
         WORKSPACE = .TRUE.
      CASE (NSCARC_INIT_FALSE)
         WORKSPACE = .FALSE.
   END SELECT
#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2, NL3, NR3
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
      SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating LOG3     array  ',A20,' in length (',I8,':',I8,' , ',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_LOG3


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL1(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT

IF (.NOT.ALLOCATED(WORKSPACE)) THEN

   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL1', CTEXT, IERROR)

   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
   END SELECT

#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
ENDIF

#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating REAL1    array  ',A20,' in length (',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_REAL1


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL1_FB(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(FB), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL1_FB', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_FB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_FB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_FB
   END SELECT
#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
ENDIF
#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating REAL1_FB array  ',A20,' in length (',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_REAL1_FB


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 2
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL2(WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL2', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
   END SELECT
#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating REAL2    array    ',A20,' in length (',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_REAL2


!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 3
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL3', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
   END SELECT
#ifdef WITH_SCARC_VERBOSE2
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2, NL3, NR3
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
      SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating REAL3    array    ',A20,' in length (',I8,':',I8,' , ',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_REAL3


!> ------------------------------------------------------------------------------------------------
!> COMPACT storage format:
!> Allocate matrix with corresponding pointer and length structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_MATRIX_COMPACT(A, NL, NPREC, NINIT, CNAME)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: A      
INTEGER, INTENT(IN) :: NPREC, NINIT, NL
CHARACTER(*), INTENT(IN) :: CNAME
CHARACTER(40) :: CINFO

A%CNAME = TRIM(CNAME)
#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000) CNAME
#endif

IF (NPREC == NSCARC_PRECISION_SINGLE) THEN

   WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CNAME),'_LEV',NL,'.VAL'
   CALL SCARC_ALLOCATE_REAL1_FB(A%VAL_FB, 1, A%N_VAL, NINIT, CINFO)

ELSE

   WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CNAME),'_LEV',NL,'.VAL'
   CALL SCARC_ALLOCATE_REAL1(A%VAL, 1, A%N_VAL, NINIT, CINFO)

ENDIF

WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CNAME),'_LEV',NL,'.ROW'
CALL SCARC_ALLOCATE_INT1(A%ROW, 1, A%N_ROW, NINIT, CINFO)

WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CNAME),'_LEV',NL,'.COL'
CALL SCARC_ALLOCATE_INT1(A%COL, 1, A%N_VAL, NINIT, CINFO)

#ifdef WITH_MKL
IF (IS_MKL_LEVEL(NL)) THEN
   WRITE(CINFO,'(A,A)') CNAME,'.COL_GLOBAL'
   CALL SCARC_ALLOCATE_INT1 (A%COL_GLOBAL, 1, A%N_VAL, NINIT, CINFO)
ENDIF
#endif

#ifdef WITH_SCARC_VERBOSE2
1000 FORMAT('Allocating COMPACT  matrix ',A20)
#endif
END SUBROUTINE SCARC_ALLOCATE_MATRIX_COMPACT

!> ------------------------------------------------------------------------------------------------
!> COMPACT storage format:
!> Deallocate matrix with corresponding pointer and length structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_MATRIX_COMPACT(A)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: A      

A%N_STENCIL   = 0
A%N_CONDENSED = 0
A%N_ROW       = 0
A%N_VAL       = 0

A%STENCIL = 0
A%POS     = 0

IF (ALLOCATED(A%VAL)) DEALLOCATE(A%VAL)
IF (ALLOCATED(A%RELAX)) DEALLOCATE(A%RELAX)
A%STENCIL = 0

IF (ALLOCATED(A%ROW)) DEALLOCATE(A%ROW)
IF (ALLOCATED(A%COL)) DEALLOCATE(A%COL)
IF (ALLOCATED(A%COL_GLOBAL)) DEALLOCATE(A%COL_GLOBAL)

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000)
1000 FORMAT('Removing   COMPACT matrix ')
#endif
END SUBROUTINE SCARC_DEALLOCATE_MATRIX_COMPACT


!> ------------------------------------------------------------------------------------------------
!> COMPACT storage format:
!> Reduce size of matrix to specified size
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESIZE_MATRIX_COMPACT(A, NSIZE, CNAME)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: A      
INTEGER, INTENT(IN) :: NSIZE
CHARACTER(*), INTENT(IN) :: CNAME
CHARACTER(40) :: CINFO
REAL(EB), ALLOCATABLE, DIMENSION(:) :: VAL
INTEGER , ALLOCATABLE, DIMENSION(:) :: COL

IF (A%N_VAL == NSIZE) THEN
   RETURN                                  !> matrix has already desired size
ELSE IF (A%N_VAL > NSIZE) THEN

   A%N_VAL = NSIZE

   WRITE(CINFO,'(A,A)') TRIM(CNAME),'.VAL'

   CALL SCARC_ALLOCATE_REAL1(VAL, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   VAL(1:NSIZE) = A%VAL(1:NSIZE)
   DEALLOCATE(A%VAL)

   CALL SCARC_ALLOCATE_REAL1(A%VAL, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   A%VAL(1:NSIZE) = VAL(1:NSIZE)
   DEALLOCATE(VAL)

   WRITE(CINFO,'(A,A)') TRIM(CNAME),'.COL'

   CALL SCARC_ALLOCATE_INT1(COL, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   COL(1:NSIZE) = A%COL(1:NSIZE)
   DEALLOCATE(A%COL)

   CALL SCARC_ALLOCATE_INT1(A%COL, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   A%COL(1:NSIZE) = COL(1:NSIZE)
   DEALLOCATE(COL)

ELSE
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SIZE, SCARC_NONE, NSIZE)
ENDIF

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000) A%CNAME, NSIZE
1000 FORMAT('Resizing   COMPACT matrix ',A20,' to size ',I8)
#endif
END SUBROUTINE SCARC_RESIZE_MATRIX_COMPACT


!> ------------------------------------------------------------------------------------------------
!> COMPACT storage format:
!> Copy matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_COPY_MATRIX_COMPACT(A, B, CB)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: A, B
CHARACTER(*), INTENT(IN) :: CB

B%CNAME       = CB

B%N_VAL       = A%N_VAL
B%N_ROW       = A%N_ROW
B%N_STENCIL   = A%N_STENCIL
B%N_CONDENSED = A%N_CONDENSED

B%STENCIL = A%STENCIL
B%POS     = A%POS

IF (.NOT.ALLOCATED(B%VAL)) CALL SCARC_ALLOCATE_REAL1(B%VAL, 1, B%N_VAL, NSCARC_INIT_ZERO, 'B%VAL')
B%VAL = A%VAL

IF (.NOT.ALLOCATED(B%COL)) CALL SCARC_ALLOCATE_INT1 (B%COL, 1, B%N_VAL, NSCARC_INIT_ZERO, 'B%COL')
B%COL = A%COL

IF (.NOT.ALLOCATED(B%ROW)) CALL SCARC_ALLOCATE_INT1 (B%ROW, 1, B%N_ROW, NSCARC_INIT_ZERO, 'B%ROW')
B%ROW = A%ROW

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000) CA, CB
1000 FORMAT('Copying   COMPACT matrix ',A20,' to ',A20)
#endif
END SUBROUTINE SCARC_COPY_MATRIX_COMPACT


!> ------------------------------------------------------------------------------------------------
!> BANDED storage format:
!> Allocate matrix with corresponding pointer and length structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_MATRIX_BANDED(A, NL, NINIT, CNAME)
TYPE (SCARC_MATRIX_BANDED_TYPE), INTENT(INOUT) :: A      
CHARACTER(*), INTENT(IN) :: CNAME
INTEGER, INTENT(IN) :: NINIT, NL
CHARACTER(40) :: CINFO

A%CNAME = CNAME

WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CNAME),'_LEV',NL,'.AUX'
CALL SCARC_ALLOCATE_REAL1(A%AUX, 1, A%N_DIAG, NINIT, CINFO)
WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CNAME),'_LEV',NL,'.VAL'
CALL SCARC_ALLOCATE_REAL2(A%VAL, 1, A%N_DIAG, 1, A%N_STENCIL, NINIT, CINFO)

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000) CNAME, A%N_VAL, A%N_DIAG
1000 FORMAT('Allocating BANDED  matrix ',A20,' with NA=',I8,' and NDIAG=',I8)
#endif
END SUBROUTINE SCARC_ALLOCATE_MATRIX_BANDED


!> ------------------------------------------------------------------------------------------------
!> BANDED storage format:
!> Deallocate matrix with corresponding pointer and length structures
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_MATRIX_BANDED(A)
TYPE (SCARC_MATRIX_BANDED_TYPE), INTENT(INOUT) :: A      

A%N_STENCIL   = 0
A%N_CONDENSED = 0
A%N_VAL       = 0
A%N_DIAG      = 0

IF (ALLOCATED(A%AUX))    DEALLOCATE(A%AUX)
IF (ALLOCATED(A%VAL))    DEALLOCATE(A%VAL)
IF (ALLOCATED(A%RELAX))  DEALLOCATE(A%RELAX)
IF (ALLOCATED(A%RELAXD)) DEALLOCATE(A%RELAXD)
A%STENCIL = 0

A%OFFSET = 0
A%LENGTH = 0
A%SOURCE = 0
A%TARGET = 0

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000) 
1000 FORMAT('Removing   COMPACT matrix ')
#endif
END SUBROUTINE SCARC_DEALLOCATE_MATRIX_BANDED


!> ====================================================================================================
!> Routines for the output of different information
!> ====================================================================================================
!> ----------------------------------------------------------------------------------------------------
!> Dump residual information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_CSV(ISM, NS, NL)
INTEGER, INTENT(IN) :: ISM, NS, NL

IF (.NOT.HAS_CSV_DUMP .OR. MYID /= 0) RETURN
IF (ITE_TOTAL == 0 .AND. TYPE_SOLVER /= NSCARC_SOLVER_MAIN) RETURN
IF (TYPE_SOLVER == NSCARC_SOLVER_COARSE) RETURN
WRITE(MSG%LU_STAT,1000) ITE_PRES, NS, ITE_TOTAL, ITE_CG, ITE_MG, NL, ITE_SMOOTH, ISM, ITE_COARSE, ITE_LU, RES, CAPPA

1000 FORMAT(10(I8,','), E12.4,',',E12.4)
END SUBROUTINE SCARC_DUMP_CSV


!> ----------------------------------------------------------------------------------------------------
!> Dump CPU times
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_TIMERS
INTEGER, PARAMETER :: LINE_LENGTH = 5 + 12*11
INTEGER :: N, STATUS(MPI_STATUS_SIZE)
CHARACTER(LEN=LINE_LENGTH) :: LINE
CHARACTER(LEN=LINE_LENGTH), DIMENSION(0:N_MPI_PROCESSES-1) :: LINE_ARRAY

! All MPI processes except root send their timings to the root process. The root process then writes them out to a file.
WRITE(LINE,'(I5,12(",",ES10.3))') MYID,                       &
                                  CPU(MYID)%OVERALL,          &
                                  CPU(MYID)%SETUP,            &
                                  CPU(MYID)%SOLVER,           &
                                  CPU(MYID)%ITERATION,        &
                                  CPU(MYID)%MATVEC_PRODUCT,   &
                                  CPU(MYID)%SCALAR_PRODUCT,   &
                                  CPU(MYID)%RELAXATION,       &
                                  CPU(MYID)%SMOOTHER,         &
                                  CPU(MYID)%COARSE,           &
                                  CPU(MYID)%EXCHANGE,         &
                                  CPU(MYID)%BUFFER_PACKING,   &
                                  CPU(MYID)%BUFFER_UNPACKING

IF (MYID>0) THEN
   CALL MPI_SEND(LINE,LINE_LENGTH,MPI_CHARACTER,0,MYID,MPI_COMM_WORLD,IERROR)
ELSE
   LINE_ARRAY(0) = LINE
   DO N=1,N_MPI_PROCESSES-1
      CALL MPI_RECV(LINE_ARRAY(N),LINE_LENGTH,MPI_CHARACTER,N,N,MPI_COMM_WORLD,STATUS,IERROR)
   ENDDO
   MSG%FILE_CPU = TRIM(CHID)//'_scarc_cpu.csv'
   OPEN (MSG%LU_CPU, FILE=MSG%FILE_CPU, STATUS='REPLACE',FORM='FORMATTED')
   WRITE(MSG%LU_CPU,'(A,A)') 'Rank,OVERALL,SETUP,SOLVER,ITERATION,MATVEC_PRODUCT,SCALAR_PRODUCT,',&
                             'RELAXATION,SMOOTHER,COARSE,EXCHANGE,PACKING,UNPACKING'
   DO N=0,N_MPI_PROCESSES-1
      WRITE(MSG%LU_CPU,'(A)') LINE_ARRAY(N)
   ENDDO
   CLOSE(MSG%LU_CPU)
ENDIF

END SUBROUTINE SCARC_DUMP_TIMERS



#ifdef WITH_SCARC_MGM
!> ====================================================================================================
!> Begin  SCARC_WITH_MGM - Part
!> Collection of routines to experiment with McKeeney-Greengard-Mayo method 
!> ====================================================================================================
!> ----------------------------------------------------------------------------------------------------
!> Store preliminary solution vector in MGM structure
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_STORE_MGM(NL, NTYPE)
USE SCARC_POINTERS, ONLY: ST, MGM
INTEGER, INTENT(IN) :: NL, NTYPE
INTEGER :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   ST=> SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)
   MGM => SCARC(NM)%LEVEL(NL)%MGM

   SELECT CASE(NTYPE)
      CASE (1)
         MGM%H1 = ST%X
      CASE (2)
         MGM%H2 = ST%X
      CASE (3)
         ST%X = MGM%H1 + MGM%H2
   END SELECT

ENDDO

END SUBROUTINE SCARC_STORE_MGM


!> ------------------------------------------------------------------------------------------------
!> Set internal boundary conditions for unstructured, homogeneous part of MGM method
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_INTERNAL_VELOCITY(NL)
USE SCARC_POINTERS, ONLY: M, L, MGM, UU, VV, WW, HP
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, I, J, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M   => MESHES(NM)
   L   => SCARC(NM)%LEVEL(NL)
   MGM => L%MGM

   IF (PREDICTOR) THEN
      UU => M%U
      VV => M%V
      WW => M%W
      HP => M%H
   ELSE
      UU => M%US
      VV => M%VS
      WW => M%WS
      HP => M%HS
   ENDIF

   DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT

      GWC => L%STRUCTURED%WALL(IW)

      I = GWC%IXW
      J = GWC%IYW
      K = GWC%IZW

      MGM%US(IW) = UU(I,J,K) - DT*( M%FVX(I,J,K) + M%RDXN(I)*(HP(I+1,J  ,K  ) - HP(I,J,K)) )
      MGM%VS(IW) = VV(I,J,K) - DT*( M%FVY(I,J,K) + M%RDYN(J)*(HP(I  ,J+1,K  ) - HP(I,J,K)) )
      MGM%WS(IW) = WW(I,J,K) - DT*( M%FVZ(I,J,K) + M%RDZN(K)*(HP(I  ,J  ,K+1) - HP(I,J,K)) )

#ifdef WITH_SCARC_DEBUG
      WRITE(MSG%LU_DEBUG,'(a,4i5,5F12.6)') 'MGM: I,J,K,IW,US,WS,DT,RDXN, RDZN:', &
         I, J, K, IW, MGM%US(IW), MGM%WS(IW), DT, M%RDXN(I), M%RDZN(I)
#endif

   ENDDO

ENDDO

END SUBROUTINE SCARC_MGM_INTERNAL_VELOCITY


!> ----------------------------------------------------------------------------------------------------
!> Allocate velocity vectors for the setting of internal BC's in MGM-method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: L, G, MGM, AC
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, IW, IP
LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_PROCESSED

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO NL = NLMIN, NLMAX

      L => SCARC(NM)%LEVEL(NL)
      G => L%UNSTRUCTURED
      MGM => L%MGM

      MGM%NWE = L%N_WALL_CELLS_EXT
      MGM%NWI = L%N_WALL_CELLS_INT

      MGM%NW1 = L%N_WALL_CELLS_EXT + 1
      MGM%NW2 = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT

      CALL SCARC_ALLOCATE_REAL1(MGM%H1, 1, G%NC, NSCARC_INIT_ZERO, 'H1')
      CALL SCARC_ALLOCATE_REAL1(MGM%H2, 1, G%NC, NSCARC_INIT_ZERO, 'H2')

      CALL SCARC_ALLOCATE_REAL1(MGM%US, MGM%NW1, MGM%NW2, NSCARC_INIT_ZERO, 'UMGM')
      CALL SCARC_ALLOCATE_REAL1(MGM%VS, MGM%NW1, MGM%NW2, NSCARC_INIT_ZERO, 'VMGM')
      CALL SCARC_ALLOCATE_REAL1(MGM%WS, MGM%NW1, MGM%NW2, NSCARC_INIT_ZERO, 'WMGM')

      CALL SCARC_ALLOCATE_INT1 (MGM%PERM , 1, G%NC, NSCARC_INIT_ZERO, 'PERMUTATION')

      CALL SCARC_ALLOCATE_REAL2(MGM%A , 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'A')
      CALL SCARC_ALLOCATE_REAL2(MGM%L , 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'L')
      CALL SCARC_ALLOCATE_REAL2(MGM%U , 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'U')

      CALL SCARC_ALLOCATE_REAL2(MGM%IA, 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'IA')
      CALL SCARC_ALLOCATE_REAL2(MGM%IL, 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'IL')
      CALL SCARC_ALLOCATE_REAL2(MGM%IU, 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'IU')

      CALL SCARC_ALLOCATE_REAL1(MGM%B, 1, G%NC, NSCARC_INIT_ZERO, 'B')
      CALL SCARC_ALLOCATE_REAL1(MGM%Y, 1, G%NC, NSCARC_INIT_ZERO, 'Y')
      CALL SCARC_ALLOCATE_REAL1(MGM%X, 1, G%NC, NSCARC_INIT_ZERO, 'X')

      CALL SCARC_ALLOCATE_LOG1(IS_PROCESSED, 1, G%NC, NSCARC_INIT_FALSE, 'IS_PROCESSED')
      AC => G%AC

      IP = G%NC - MGM%NWI + 1
      DO IW = MGM%NW1, MGM%NW2
         IC = G%WALL(IW)%ICW
         IS_PROCESSED(IC) = .TRUE.
         MGM%PERM(IP) = IC
         IP = IP + 1
      ENDDO

      IP = 1
      DO IC = 1, G%NC
         IF (IS_PROCESSED(IC)) CYCLE
         IS_PROCESSED(IC) = .TRUE.
         MGM%PERM(IP) = IC
         IP = IP + 1
      ENDDO

      DEALLOCATE (IS_PROCESSED)

      CALL SCARC_SETUP_MGM_LU(NM, NL)
      CALL SCARC_METHOD_MGM_LU(NM, NL)

      CALL SCARC_SETUP_MGM_ILU(NM, NL)
      CALL SCARC_METHOD_MGM_ILU(NM, NL)

   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_MGM


!> ---------------------------------------------------------------------------------------------
!> Setup LU-decomposition for MGM-method
!> ---------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_LU(NM, NL)
USE SCARC_POINTERS, ONLY: G, MGM, AC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, JC, ICOL, IP, IVERSION
REAL (EB) :: SCAL
REAL (EB), DIMENSION(:,:), POINTER :: A, L, U

G   => SCARC(NM)%LEVEL(NL)%UNSTRUCTURED
MGM => SCARC(NM)%LEVEL(NL)%MGM
AC  => G%AC

A => MGM%A
L => MGM%L
U => MGM%U

!>
!> Temporarily extract full matrix from compact storage technique - just for proof of concept
!> Consider permutation in MGM%PERM
!>
DO JC = 1, G%NC
   IC = MGM%PERM(JC)
   IC = JC

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'PERMUTATION  --- JC =', JC, '-----------> IC=',IC
#endif

   DO IP = AC%ROW(IC), AC%ROW(IC+1)-1
      ICOL = AC%COL(IP)

#ifdef WITH_SCARC_DEBUG
      WRITE(MSG%LU_DEBUG,*) 'PERMUTATION ICOL ', ICOL, AC%VAL(IP), ' TO ', JC, ICOL
#endif

      MGM%A(JC,ICOL) = AC%VAL(IP)
   ENDDO
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '------- MGM%A - Copy (1:24)'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%A(IC, JC), JC=1, 24)
ENDDO
#endif

IVERSION = 1

!> Version 1
!IF (IVERSION == 1) THEN
DO I = 1, G%NC
   U(I,I) = 1.0_EB
ENDDO

DO J = 1, G%NC
   DO I = J, G%NC
      SCAL = 0.0_EB
      DO K = 1, J-1
         SCAL = SCAL + L(I,K) * U(K,J)
      ENDDO
      L(I,J) = A(I,J) - SCAL
   ENDDO
   DO I = J, G%NC
      SCAL = 0.0_EB
      DO K = 1, J-1
         SCAL = SCAL + L(J,K) * U(K,I)
      ENDDO
      U(J,I) = (A(J,I) - SCAL) / L(J,J)
   ENDDO
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== MGM%A (1:24) '
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%A(IC, JC), JC=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) 'MGM%L'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%L(IC, JC), JC=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) 'MGM%U'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%U(IC, JC), JC=1, 24)
ENDDO
#endif

!> Version 2
!ELSE IF (IVERSION == 2) THEN

DO I = 1, G%NC

   L(I,I) = 1.0_EB

   DO J = I, G%NC

      SCAL = 0.0_EB
      DO K = 1, I-1
         SCAL = SCAL + L(I,K) * U(K,J)
      ENDDO
      U(I,J) = A(I,J) - SCAL

      SCAL = 0.0_EB
      DO K = 1, I-1
         SCAL = SCAL + L(J,K) * U(K,I)
      ENDDO
      L(J,I) = (A(J,I) - SCAL)/U(I,I)

   ENDDO
ENDDO

#ifdef WITH_SCARC_MGM
WRITE(MSG%LU_DEBUG,*) '=============================== MGM%A (1:24) '
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%A(IC, JC), JC=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) 'MGM%L'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%L(IC, JC), JC=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) 'MGM%U'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%U(IC, JC), JC=1, 24)
ENDDO
#endif

!> Version 3
!ELSE IF (IVERSION == 3) THEN

!   DO I = 1, G%NC
!      L(I,1) = A(I,1)
!      U(I,I) = 1.0_EB
!   ENDDO
!
!   DO J = 2, G%NC
!      U(1,J) = A(1,J)/L(1,1)
!   ENDDO
!
!   DO I = 2, G%NC
!
!      DO J = 2, I
!         SCAL = 0.0_EB
!         DO K = 1, J-1
!            SCAL = SCAL + L(I,K)*U(K,J)
!         ENDDO
!         L (I,J) = A(I,J) - SCAL
!      ENDDO
!
!      DO J = I+1, G%NC
!         SCAL = 0.0_EB
!         DO K = 1, I-1
!            SCAL = SCAL + L(I,K)*U(K,J)/L(I,I)
!         ENDDO
!         U(I,J) = A(I,J) - SCAL
!      ENDDO
!   ENDDO

!ENDIF

DO I = 1, G%NC
   DO J = 1, I-1
      A(I,J) = L(I,J)
   ENDDO
   DO J = I, G%NC
      A(I,J) = U(I,J)
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_MGM_LU


!> ---------------------------------------------------------------------------------------------
!> Setup ILU-decomposition for MGM-method
!> ---------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_ILU(NM, NL)
USE SCARC_POINTERS, ONLY: G, MGM, AC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, JC, ICOL, IP
REAL (EB), DIMENSION(:,:), POINTER :: IA, IL, IU

G   => SCARC(NM)%LEVEL(NL)%UNSTRUCTURED
MGM => SCARC(NM)%LEVEL(NL)%MGM
AC  => G%AC

IA => MGM%IA
IL => MGM%IL
IU => MGM%IU

!>
!> Temporarily extract full matrix from compact storage technique - just for proof of concept
!> Consider permutation in MGM%PERM
!>
DO JC = 1, G%NC
   !   IC = MGM%PERM(JC)
   IC = JC
   DO IP = AC%ROW(IC), AC%ROW(IC+1)-1
      ICOL = AC%COL(IP)
      IA(JC,ICOL) = AC%VAL(IP)
   ENDDO
ENDDO

DO I = 2, G%NC
   DO K = 1, I-1
      IF (IA(I,K) /= 0 .AND. IA(K,K) /= 0) THEN
         IA(I,K) = IA(I,K)/IA(K,K)
         DO J = K+1, G%NC
            IF (IA(I,J) /= 0) THEN
               IA(I,J) = IA(I,J) - IA(I,K)*IA(K,J)
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_MGM_ILU


!> ---------------------------------------------------------------------------------------------
!> Setup LU-decomposition for MGM-method
!> ---------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MGM_LU(NM, NL)
USE SCARC_POINTERS, ONLY: G, MGM, AC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, IW, N
REAL (EB), DIMENSION(:,:), POINTER :: A, L, U
REAL (EB), DIMENSION(:),   POINTER :: B, X, Y

G   => SCARC(NM)%LEVEL(NL)%UNSTRUCTURED
MGM => SCARC(NM)%LEVEL(NL)%MGM
AC  => G%AC

A => MGM%A
L => MGM%L
U => MGM%U

B => MGM%B
X => MGM%X
Y => MGM%Y

N = G%NC

DO IW = MGM%NW1, MGM%NW2
   IC = G%WALL(IW)%ICW
   B(IC) = 1.0_EB
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== A'
DO I = 1, N
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (A(I, J), J=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) '=============================== L'
DO I = 1, N
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (L(I, J), J=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) '=============================== U'
DO I = 1, N
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (U(I, J), J=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) '=============================== B'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (B(I), I=1, G%NC)
#endif

DO J = 1, N
   Y(J) = B(J)
   DO K = 1, J-1
      Y(J) = Y(J) - A(J,K)*Y(K)
   ENDDO
ENDDO

DO J = N, 1, -1
   X(J) = Y(J)
   DO K = J+1, N
      X(J) = X(J) - A(J,K)*X(K)
   ENDDO
   X(J) = X(J)/A(J,J)
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== Y'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (Y(I), I=1, G%NC)
WRITE(MSG%LU_DEBUG,*) '=============================== X'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (X(I), I=1, G%NC)
#endif

END SUBROUTINE SCARC_METHOD_MGM_LU


!> ---------------------------------------------------------------------------------------------
!> Setup ILU-decomposition for MGM-method
!> ---------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MGM_ILU(NM, NL)
USE SCARC_POINTERS, ONLY: G, MGM, AC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, IW, N
REAL (EB), DIMENSION(:,:), POINTER :: IA, IL, IU
REAL (EB), DIMENSION(:),   POINTER :: B, X, Y

G   => SCARC(NM)%LEVEL(NL)%UNSTRUCTURED
MGM => SCARC(NM)%LEVEL(NL)%MGM
AC  => G%AC

IA => MGM%IA
IL => MGM%IL
IU => MGM%IU

B => MGM%B
X => MGM%X
Y => MGM%Y

N = G%NC

DO IW = MGM%NW1, MGM%NW2
   IC = G%WALL(IW)%ICW
   B(IC) = 1.0_EB
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== A'
DO I = 1, N
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (IA(I, J), J=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) '=============================== B'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (B(I), I=1, G%NC)
#endif

DO J = 1, N
   Y(J) = B(J)
   DO K = 1, J-1
      Y(J) = Y(J) - IA(J,K)*Y(K)
   ENDDO
ENDDO

DO J = N, 1, -1
   X(J) = Y(J)
   DO K = J+1, N
      X(J) = X(J) - IA(J,K)*X(K)
   ENDDO
   X(J) = X(J)/IA(J,J)
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== Y'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (Y(I), I=1, G%NC)
WRITE(MSG%LU_DEBUG,*) '=============================== X'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (X(I), I=1, G%NC)
#endif

END SUBROUTINE SCARC_METHOD_MGM_ILU

!> ====================================================================================================
!> End SCARC_WITH_MGM - Part
!> ====================================================================================================
#endif



#ifdef WITH_SCARC_DEBUG
!> ================================================================================================
!> Start  WITH_SCARC_DEBUG  - Part
!> Collection of routines which print out different quantities or allow to preset them
!> ================================================================================================
!> ------------------------------------------------------------------------------------------------
!> Preset right hand side in such a way that exact solution is known
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_EXACT (NE, NL)
USE SCARC_POINTERS, ONLY: M, L, G, VC, XMID, ZMID
INTEGER, INTENT(IN) :: NE, NL
INTEGER :: IC, NM, I, K

IF (ITE_TOTAL == 0) WRITE(*,*) 'TODO: PRESET_EXACT is active !!!'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NE)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I,1,K)) CYCLE
         IC = G%CELL_NUMBER(I,1,K)
         IF (NL == NLEVEL_MIN) THEN
            XMID => M%XC
            ZMID => M%ZC
         ELSE
            XMID => L%XMID
            ZMID => L%ZMID
         ENDIF
         VC(IC) = EXACT(XMID(I),ZMID(K))
         WRITE(MSG%LU_DEBUG,'(A,i3,a,e10.2,a,e10.2,a,e12.4)') 'IC=',IC,':X=',XMID(i),':Z=',ZMID(k),': RHS=',VC(IC)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_EXACT


!> ------------------------------------------------------------------------------------------------
!> Preset vector
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_VECTOR (NV, NL)
USE SCARC_POINTERS, ONLY: M, G, VC, XMID, ZMID
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: IC, NM, I, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I,1,K)) CYCLE
         IC = G%CELL_NUMBER(I,1,K)
         IF (NL == NLEVEL_MIN) THEN
            XMID => M%XC
            ZMID => M%ZC
         ELSE
            XMID => L%XMID
            ZMID => L%ZMID
         ENDIF
         VC(IC) = XMID(I)
         WRITE(MSG%LU_DEBUG,'(A,i3,a,e10.2,a,e10.2,a,e12.4)') 'IC=',IC,':X=',XMID,':Z=',ZMID,': RHS=',VC(IC)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_VECTOR


!> ------------------------------------------------------------------------------------------------
!> Preset right hand side in such a way that exact solution is known
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_RHS (NV, NL)
USE SCARC_POINTERS, ONLY: M, L, G, VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: IC, NM, I, K
REAL (EB) :: X, Z

IF (NL > NLEVEL_MIN) WRITE(*,*) 'Wrong level for presetting RHS '

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   M%BXS = 0.0_EB
   M%BXF = 0.0_EB
   M%BZS = 0.0_EB
   M%BZF = 0.0_EB

   CALL SCARC_POINT_TO_GRID(NM, NL)
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I,1,K)) CYCLE
         IC = G%CELL_NUMBER(I,1,K)
         X  = M%XC(I)
         Z  = M%ZC(K)
         WRITE(MSG%LU_DEBUG,'(A,i3,a,e10.2,a,e10.2,a,e12.4)') 'IC=',IC,':X=',X,':Z=',Z,': RHS=',VC(IC)
         VC(IC) = RHS(X,Z)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_RHS


!> ------------------------------------------------------------------------------------------------
!> Set exact solution
!> ------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION EXACT(X,Z)
REAL (EB), INTENT(IN) :: X, Z
!EXACT = (X**2 - X**4) * (Z**4 - Z**2)                                    !> FUNCTION 1
!EXACT = (X**2 - 1) * (Z**2 - 1)                                         !> FUNCTION 2
!EXACT =  625.0_EB/16.0_EB * X * (0.8_EB - X) * Z * (0.8_EB - Z)        !> FUNCTION 3
EXACT = - X * (0.8_EB - X) * Z * (0.8_EB - Z)        !> FUNCTION 3
END FUNCTION EXACT


!> ------------------------------------------------------------------------------------------------
!> Set right hand side
!> ------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION RHS(X,Z)
REAL (EB), INTENT(IN) :: X, Z
!RHS = 2.0_EB*((1.0_EB - 6.0_EB*X**2)*Z**2*(1.0_EB-Z**2)+(1.0_EB-6.0_EB*Z**2)*X**2*(1.0_EB-X**2))
!RHS = -X**2 - Z**2 +2
!RHS = 625.0_EB/8.0_EB * (X * (0.8_EB - X) + Z * (0.8_EB - Z))
RHS = 2.0_EB * (X * (0.8_EB - X) + Z * (0.8_EB - Z))
END FUNCTION RHS


!> ------------------------------------------------------------------------------------------------
!> Save dump of vector in dump-directory
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_QUANTITY (NV, CNAME, ISM, NS, NL)
USE SCARC_POINTERS, ONLY: G, VC
INTEGER, INTENT(IN) :: NV, NS, NL, ISM
CHARACTER(*), INTENT(IN) :: CNAME
CHARACTER(80) :: FN_DUMP, CDIR
INTEGER :: LU_DUMP, IC, NM

IF (TRIM(CNAME) == 'RESIDUAL') THEN
   CDIR = 'res'
ELSE IF (TRIM(CNAME) == 'ERROR') THEN
   CDIR = 'err'
ELSE IF (TRIM(CNAME) == 'RHS') THEN
   CDIR = 'rhs'
ELSE IF (TRIM(CNAME) == 'EXACT') THEN
   CDIR = 'exa'
ELSE IF (TRIM(CNAME) == 'DISCRET') THEN
   CDIR = 'dis'
ENDIF

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)
   CALL SCARC_POINT_TO_GRID(NM, NL)
   if (ISM == 0) THEN
      WRITE (FN_DUMP, '(A,A3,A,i3.3,A,I3.3,A,I3.3,A,I3.3,A,I1,A,A)') &
         'dump/',CDIR,'_t',ITE_TOTAL,'_ALL',ITE_SMOOTH,'_mg',ITE_MG,'_cg',ITE_CG,'_level',NL,'_',&
         TRIM(STACK(NS)%SOLVER%CNAME)
   else if (ISM == NSCARC_CYCLING_PRESMOOTH) THEN
      WRITE (FN_DUMP, '(A,A3,A,i3.3,A,I3.3,A,I3.3,A,I3.3,A,I1,A,A)') &
         'dump/',CDIR,'_t',ITE_TOTAL,'_PRE',ITE_SMOOTH,'_mg',ITE_MG,'_cg',ITE_CG,'_level',NL,'_',&
         TRIM(STACK(NS)%SOLVER%CNAME)
   ELSE IF (ISM == NSCARC_CYCLING_POSTSMOOTH) THEN
      WRITE (FN_DUMP, '(A,A3,A,i3.3,A,I3.3,A,I3.3,A,I3.3,A,I1,A,A)') &
         'dump/',CDIR,'_t',ITE_TOTAL,'_POST',ITE_SMOOTH,'_mg',ITE_MG,'_cg',ITE_CG,'_level',NL,'_',&
         TRIM(STACK(NS)%SOLVER%CNAME)
   ENDIF
   LU_DUMP = GET_FILE_NUMBER()
   OPEN (LU_DUMP, FILE=FN_DUMP)
   DO IC = 1, G%NC
      !WRITE(LU_DUMP,'(F25.16)') VC(IC)
      WRITE(LU_DUMP,*) VC(IC)
   ENDDO
   CLOSE(LU_DUMP)
ENDDO

FN_DUMP = CNAME
END SUBROUTINE SCARC_DUMP_QUANTITY


!> ------------------------------------------------------------------------------------------------
!> Print out vector information on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_LEVEL (NV, CVEC, NL)
USE SCARC_POINTERS, ONLY: L, G, VC
INTEGER, INTENT(IN) :: NV, NL
REAL (EB) :: VALUES(0:100)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CVEC

!IF (TYPE_SOLVER /= NSCARC_SOLVER_MAIN) RETURN
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   NX8=MIN(10,L%NX)
   NY8=MIN(10,L%NY)
   NZ8=MIN(8,L%NZ)

   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   WRITE(MSG%LU_DEBUG,2001) CVEC, NM, NL
   WRITE(MSG%LU_DEBUG,2002) G%NC, NX8, NY8, NZ8, NV, SIZE(VC)
   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   !IF (NL == NLEVEL_MIN) THEN
   DO KK = NZ8, 1, - 1
      DO JJ = NY8, 1, - 1
         DO II=1,NX8
            IF (IS_UNSTRUCTURED.AND.L%IS_SOLID(II,JJ,KK)) THEN
               VALUES(II)=-999999.0_EB
            ELSE
               IC=G%CELL_NUMBER(II,JJ,KK)
               IF (ABS(VC(IC))<1.0E-14_EB) THEN
                  VALUES(II)=0.0_EB
               ELSE
                  VALUES(II)=VC(IC)
               ENDIF
            ENDIF
         ENDDO
         !WRITE(MSG%LU_DEBUG, '(5E23.14)') (VALUES(II), II=4, NX8)
         WRITE(MSG%LU_DEBUG, '(12E14.6)') (VALUES(II), II=1, NX8)
         !WRITE(MSG%LU_DEBUG, '(12E11.3)') (VALUES(II), II=1, NX8)
      ENDDO
      IF (.NOT. TWO_D) WRITE(MSG%LU_DEBUG, *) '----------------'
   ENDDO
   !ENDIF
   WRITE(MSG%LU_DEBUG, *) '---------------- Overlap ----------------'
   WRITE(MSG%LU_DEBUG, '(4E14.6)') (VC(IC), IC = G%NC+1, G%NCE)
ENDDO

!CALL SCARC_MATLAB_VECTOR(NV, CVEC, NL)

!2000 FORMAT('=== ',A,' : ',A,' on mesh ',I8,' on level ',I8, ': NX, NY, NZ=',3I8,': NV=',I8)
2001 FORMAT('=== ',A,' on mesh ',I8,' on level ',I8)
2002 FORMAT('=== NC = ',I6, ': NX, NY, NZ=',3I6,': NV=',I6,': Size=',I8)
END SUBROUTINE SCARC_DEBUG_LEVEL


!> ------------------------------------------------------------------------------------------------
!> Print out vector information on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_LEVEL2 (NV, CVEC, NL)
USE SCARC_POINTERS, ONLY:  G, VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM, IC
CHARACTER (*), INTENT(IN) :: CVEC

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   WRITE(MSG%LU_DEBUG,*) CVEC, NM, NL
   WRITE(MSG%LU_DEBUG,*) G%NC, NV, SIZE(VC)
   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   DO IC=1, G%NC
      WRITE(MSG%LU_DEBUG, '(I4,E25.16)') IC, VC(IC)
   ENDDO
ENDDO

END SUBROUTINE SCARC_DEBUG_LEVEL2


!> ------------------------------------------------------------------------------------------------
!> Debug specified quantity
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_QUANTITY(NTYPE, NL, CQUANTITY)
#ifdef WITH_MKL
USE SCARC_POINTERS, ONLY: M, L, OL, G, OG, SV, AB, AC, ACS
#else
USE SCARC_POINTERS, ONLY: M, L, OL, G, OG, SV, AB, AC
#endif
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM, IP, IC, ID, IW, I, J, IOR0, INBR, III, JJJ, KKK, IWG, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CQUANTITY

SELECT CASE (NTYPE)

   !> ------------------------------------------------------------------------------------------------
   !> Debug system matrix A (corresponding to system type)
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_PRESSURE)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID(NM, NL)

         NX8=MIN(8,L%NX)
         NY8=MIN(8,L%NY)
         NZ8=MIN(8,L%NZ)

         WRITE(MSG%LU_DEBUG,*) '========= PRESSURE ========= ', NM, NL, CQUANTITY
         IF (PREDICTOR) THEN
            WRITE(MSG%LU_DEBUG,*) 'RHO'
            DO KKK = NZ8+1, 0, -1
               WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%RHO(III, 1, KKK), III=0,NX8+1)
            ENDDO
            WRITE(MSG%LU_DEBUG,*) 'H'
            DO KKK = NZ8+1, 0, -1
               WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%H(III, 1, KKK), III=0,NX8+1)
            ENDDO
         ELSE
            WRITE(MSG%LU_DEBUG,*) 'RHOS'
            DO KKK = NZ8+1, 0, -1
               WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%RHOS(III, 1, KKK), III=0,NX8+1)
            ENDDO
            WRITE(MSG%LU_DEBUG,*) 'HS'
            DO KKK = NZ8+1, 0, -1
               WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%HS(III, 1, KKK), III=0,NX8+1)
            ENDDO
         ENDIF
         WRITE(MSG%LU_DEBUG,*) 'FVX'
         DO KKK = NZ8+1, 0, -1
            WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%FVX(III, 1, KKK), III=0,NX8+1)
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'FVY'
         DO KKK = NZ8+1, 0, -1
            WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%FVY(III, 1, KKK), III=0,NX8+1)
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'FVZ'
         DO KKK = NZ8+1, 0, -1
            WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%FVZ(III, 1, KKK), III=0,NX8+1)
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'KRES'
         DO KKK = NZ8+1, 0, -1
            WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%KRES(III, 1, KKK), III=0,NX8+1)
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'M%PRHS'
         DO KKK = NZ8+1, 1, -1
            WRITE(MSG%LU_DEBUG,'(I4,9E14.6)') KKK, (M%PRHS(III, 1, KKK), III=1,NX8+1)
         ENDDO
         !WRITE(MSG%LU_DEBUG,*) 'P%PRHS'
         !DO KKK = NZ8+1, 1, -1
         !   WRITE(MSG%LU_DEBUG,'(I4,9E14.6)') KKK, (P%PRHS(III, 1, KKK), III=1,NX8+1)
         !ENDDO

      ENDDO

   !> ------------------------------------------------------------------------------------------------
   !> Debug system matrix A (corresponding to system type)
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_STACK)

      WRITE(MSG%LU_DEBUG,*) 'N_STACK_TOTAL=',N_STACK_TOTAL
      DO I = 1, N_STACK_TOTAL
         SV  => STACK(I)%SOLVER
         WRITE(MSG%LU_DEBUG,*) '===================== STACK ', I,' ======================'
         WRITE(MSG%LU_DEBUG,*) '-------------------- SOLVER:'
         WRITE(MSG%LU_DEBUG,*) 'NAME=',SV%CNAME
         WRITE(MSG%LU_DEBUG,*) '-- SETTING:'
         WRITE(MSG%LU_DEBUG,*) 'EPS   = ',SV%EPS
         WRITE(MSG%LU_DEBUG,*) 'RES   = ',SV%RES
         WRITE(MSG%LU_DEBUG,*) 'RESIN = ',SV%RESIN
         WRITE(MSG%LU_DEBUG,*) 'OMEGA = ',SV%OMEGA
         WRITE(MSG%LU_DEBUG,*) 'ITE   = ',SV%ITE
         WRITE(MSG%LU_DEBUG,*) 'NIT   = ',SV%NIT
         WRITE(MSG%LU_DEBUG,*) '-- TYPES:'
         WRITE(MSG%LU_DEBUG,*) 'TYPE_PARENT   = ',SV%TYPE_PARENT
         WRITE(MSG%LU_DEBUG,*) 'TYPE_SOLVER   = ',SV%TYPE_SOLVER
         WRITE(MSG%LU_DEBUG,*) 'TYPE_STAGE    = ',SV%TYPE_STAGE
         WRITE(MSG%LU_DEBUG,*) 'TYPE_SCOPE(0) = ',SV%TYPE_SCOPE(0)
         WRITE(MSG%LU_DEBUG,*) 'TYPE_PRECON   = ',SV%TYPE_RELAX
         WRITE(MSG%LU_DEBUG,*) 'TYPE_ACCURACY = ',SV%TYPE_ACCURACY
         WRITE(MSG%LU_DEBUG,*) 'TYPE_INTERPOL = ',SV%TYPE_INTERPOL
         WRITE(MSG%LU_DEBUG,*) 'TYPE_CYCLING  = ',SV%TYPE_CYCLING
         WRITE(MSG%LU_DEBUG,*) 'TYPE_TWOLEVEL = ',SV%TYPE_TWOLEVEL
         WRITE(MSG%LU_DEBUG,*) '-- POINTERS:'
         WRITE(MSG%LU_DEBUG,*) 'X   = ',SV%X
         WRITE(MSG%LU_DEBUG,*) 'B   = ',SV%B
         WRITE(MSG%LU_DEBUG,*) 'D   = ',SV%D
         WRITE(MSG%LU_DEBUG,*) 'E   = ',SV%E
         WRITE(MSG%LU_DEBUG,*) 'R   = ',SV%R
         WRITE(MSG%LU_DEBUG,*) 'V   = ',SV%V
         WRITE(MSG%LU_DEBUG,*) 'Y   = ',SV%Y
         WRITE(MSG%LU_DEBUG,*) 'Z   = ',SV%Z
      ENDDO

   !> ------------------------------------------------------------------------------------------------
   !> Debug system matrix A (corresponding to system type)
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIX)

      SELECT CASE (SCARC_MATRIX_LEVEL(NL))

         CASE (NSCARC_MATRIX_COMPACT)

            DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
               CALL SCARC_POINT_TO_GRID(NM, NL)
               AC => G%AC
               WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
               WRITE(MSG%LU_DEBUG,*) '----------- SHOWING FULL COMPACT MATRIX ENTRIES'
               WRITE(MSG%LU_DEBUG,*) 'G%NC =',G%NC
               WRITE(MSG%LU_DEBUG,*) '=========================================='
               WRITE(MSG%LU_DEBUG,*) 'AC%N_VAL =',AC%N_VAL
               WRITE(MSG%LU_DEBUG,*) 'AC%N_ROW =',AC%N_ROW
               WRITE(MSG%LU_DEBUG,*) 'SIZE(AC%VAL) =',SIZE(AC%VAL)
               WRITE(MSG%LU_DEBUG,*) 'SIZE(AC%COL) =',SIZE(AC%COL)
               WRITE(MSG%LU_DEBUG,*) 'SIZE(AC%ROW) =',SIZE(AC%ROW)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC%POS:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AC%POS(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC%STENCIL:'
               WRITE(MSG%LU_DEBUG,'(7F10.3)') (AC%STENCIL(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC%ROW:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AC%ROW(IC), IC=1,AC%N_ROW)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC%COL:'
               DO IC = 1, AC%N_ROW-1
                  WRITE(MSG%LU_DEBUG,'(i5,a,20i9)') IC,':',(AC%COL(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
               ENDDO
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC(',NM,',',NL,'):'
               DO IC = 1, AC%N_ROW-1
                  WRITE(MSG%LU_DEBUG,'(i5,a,20f8.1)') IC,':',(AC%VAL(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
               ENDDO
               !IF (ALLOCATED(AC%RELAX)) THEN
               !   WRITE(MSG%LU_DEBUG,*) '---------------------- RELAX:'
               !   DO IC = 1, AC%N_ROW-1
               !      WRITE(MSG%LU_DEBUG,'(i5,a,20f15.6)') IC,':',(AC%RELAX(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
               !   ENDDO
               !ENDIF
#ifdef WITH_MKL
               !IF ((TYPE_METHOD == NSCARC_METHOD_LU) .AND. (TYPE_MKL(0) == NSCARC_MKL_GLOBAL)) THEN
               !   WRITE(MSG%LU_DEBUG,*) '---------------------- AG_COL:'
               !   DO IC = 1, G%NC
               !      WRITE(MSG%LU_DEBUG,'(i5,a,20i9)') IC,':',(AC%COLG(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
               !      WRITE(MSG%LU_DEBUG,*)  IC,':',(AC%COLG(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
               !   ENDDO
               !ENDIF
               !DO IC=1,AC%N_ROW-1
               !   DO IP=AC%ROW(IC),AC%ROW(IC+1)-1
               !       WRITE(MSG%LU_DEBUG,'(2i10,F24.12)') IC,AC%COL(IP),AC%VAL(IP)
               !       WRITE(MSG%LU_DEBUG,*) IC,AC%COL(IP),AC%VAL(IP)
               !   ENDDO
               !   WRITE(MSG%LU_DEBUG,*)
               !ENDDO
               !DO IC=1,AC%N_ROW-1
               !   DO IP=AC%ROW(IC),AC%ROW(IC+1)-1
               !       IF (IC == AC%COL(IP)) WRITE(MSG%LU_DEBUG,'(2i10,F24.12)') IC,AC%COL(IP),AC%VAL(IP)
               !   ENDDO
               !ENDDO
#endif
               IF (IS_STRUCTURED) THEN
                  CALL SCARC_MATLAB_MATRIX(AC%VAL, AC%ROW, AC%COL, G%NC, G%NC, NM, NL, 'A_struct')
               ELSE
                  CALL SCARC_MATLAB_MATRIX(AC%VAL, AC%ROW, AC%COL, G%NC, G%NC, NM, NL, 'A_unstruct')
               ENDIF

            ENDDO

         CASE (NSCARC_MATRIX_BANDED)

            DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
               CALL SCARC_POINT_TO_GRID(NM, NL)
               AB => G%AB
               WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
               WRITE(MSG%LU_DEBUG,*) '----------- SHOWING FULL BANDED MATRIX ENTRIES'
               WRITE(MSG%LU_DEBUG,*) 'NC =',G%NC
               WRITE(MSG%LU_DEBUG,*) 'N_VAL  =',AB%N_VAL
               WRITE(MSG%LU_DEBUG,*) 'N_DIAG =',AB%N_DIAG
               WRITE(MSG%LU_DEBUG,*) 'N_STENCIL =',AB%N_STENCIL
               WRITE(MSG%LU_DEBUG,*) 'SIZE(AB%VAL) =',SIZE(AB%VAL)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%POS:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%POS(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%STENCIL:'
               WRITE(MSG%LU_DEBUG,'(7F10.3)') (AB%STENCIL(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%OFFSET:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%OFFSET(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%LENGTH:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%LENGTH(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%SOURCE:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%SOURCE(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%TARGET:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%TARGET(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%VAL:'
               DO IC = 1, AB%N_DIAG
                 WRITE(MSG%LU_DEBUG,'(I5,A,7(F12.3))') IC,':',(AB%VAL(IC, ID),ID=1,AB%N_STENCIL)
               ENDDO
            ENDDO

      END SELECT

   !> ------------------------------------------------------------------------------------------------
   !> Debug symmetric system matrix AS
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIX_SYM)

      IF (NL /= NLEVEL_MAX) RETURN

#ifdef WITH_MKL
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         CALL SCARC_POINT_TO_GRID(NM, NL)
         ACS => G%ACS
         WRITE(MSG%LU_DEBUG,*) '----------- SHOWING SYMMETRIC MATRIX ENTRIES'
         WRITE(MSG%LU_DEBUG,*) 'G%NC=', G%NC
         WRITE(MSG%LU_DEBUG,*) 'ACS%ROW(G%NC)=', ACS%ROW(G%NC)
         WRITE(MSG%LU_DEBUG,*) 'ACS%ROW(G%NC+1)=', ACS%ROW(G%NC+1)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(ASC) =',    SIZE(ACS%VAL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(ACS%COL) =',SIZE(ACS%COL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(ACS%ROW) =',SIZE(ACS%ROW)
         WRITE(MSG%LU_DEBUG,*) '---------------------- ACS%ROW:', G%NC
         WRITE(MSG%LU_DEBUG,*) (ACS%ROW(IC), IC=1,G%NC+1)
         WRITE(MSG%LU_DEBUG,*) '---------------------- ACS%COL:'
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(ACS%COL(IP),IP=ACS%ROW(IC),ACS%ROW(IC+1)-1)
         ENDDO
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(ACS%VAL(IP),IP=ACS%ROW(IC),ACS%ROW(IC+1)-1)
         ENDDO
         DO IC=1,G%NC
            DO IP=ACS%ROW(IC),ACS%ROW(IC+1)-1
               WRITE(MSG%LU_DEBUG,*) IC,ACS%COL(IP),ACS%VAL(IP)
            ENDDO
            WRITE(MSG%LU_DEBUG,*)
         ENDDO

      ENDDO
#endif

   !> ------------------------------------------------------------------------------------------------
   !> Debug FACEINFO
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_FACE)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         DO IOR0 = -3, 3
            IF (IOR0 == 0) CYCLE
            WRITE(MSG%LU_DEBUG,*) '========================================='
            WRITE(MSG%LU_DEBUG,*) '============= DEBUGGING FACE(',IOR0,'): FOR LEVEL ', NL
            WRITE(MSG%LU_DEBUG,*) '========================================='
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%DH: '
            WRITE(MSG%LU_DEBUG,'(12F8.2)')  L%FACE(IOR0)%DH
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NOP:', L%FACE(IOR0)%NOP
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NCW:', L%FACE(IOR0)%NCW
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NX: ', L%FACE(IOR0)%NX
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NY: ', L%FACE(IOR0)%NY
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NZ: ', L%FACE(IOR0)%NZ
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NW0:', L%FACE(IOR0)%NCW0
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%N_NEIGHBORS: ', L%FACE(IOR0)%N_NEIGHBORS
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%SCAL_BOUNDARY: ', L%FACE(IOR0)%SCAL_BOUNDARY
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%SCAL_DIRICHLET: ', L%FACE(IOR0)%SCAL_DIRICHLET
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%SCAL_NEUMANN: ', L%FACE(IOR0)%SCAL_NEUMANN
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%SCAL_FACE: ', L%FACE(IOR0)%SCAL_FACE
            WRITE(MSG%LU_DEBUG,*) '----------------------------------------------'
            WRITE(MSG%LU_DEBUG,*)
            IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
               DO INBR=1,L%FACE(IOR0)%N_NEIGHBORS
                  NOM = L%FACE(IOR0)%NEIGHBORS(INBR)
                  OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  SELECT CASE(TYPE_GRID)
                  CASE (NSCARC_GRID_STRUCTURED)
                     OG => OL%STRUCTURED
                  CASE (NSCARC_GRID_UNSTRUCTURED)
                     OG => OL%UNSTRUCTURED
                  END SELECT
                  WRITE(MSG%LU_DEBUG,*) 'N_NEIGHBORS:', L%FACE(IOR0)%N_NEIGHBORS
                  WRITE(MSG%LU_DEBUG,*) 'NOM:', NOM
                  WRITE(MSG%LU_DEBUG,*) 'SIZE(OG%ICG_TO_IWG)=',SIZE(OG%ICG_TO_IWG)
                  WRITE(MSG%LU_DEBUG,*) 'SIZE(OG%ICG_TO_ICE)=',SIZE(OG%ICG_TO_ICE)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a,2i8)') '---OG(',NOM,')%NCG(.):',OG%NCG
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OG(',NOM,')%ICG_TO_IWG:'
                  DO IW = 1, OG%NCG
                     WRITE(MSG%LU_DEBUG,'(16i8)') OG%ICG_TO_IWG(IW)
                  ENDDO
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OG(',NOM,')%ICG_TO_ICE:'
                  DO IW = 1, OG%NCG
                     WRITE(MSG%LU_DEBUG,'(16i8)') OG%ICG_TO_ICE(IW)
                  ENDDO
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OG(',NOM,')%ICG_TO_ICW:'
                  DO IW = 1, OG%NCG
                     WRITE(MSG%LU_DEBUG,'(16i8)') OG%ICG_TO_ICW(IW)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

   !> ------------------------------------------------------------------------------------------------
   !> Debug WALLINFO
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_WALL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_GRID(NM, NL)
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         WRITE(MSG%LU_DEBUG,*) 'SIZE(G%ICE_TO_IWG)=',SIZE(G%ICE_TO_IWG)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(G%ICE_TO_ICG)=',SIZE(G%ICE_TO_ICG)
         WRITE(MSG%LU_DEBUG,*) 'NM  =',NM
         WRITE(MSG%LU_DEBUG,*) 'NL  =',NL
         WRITE(MSG%LU_DEBUG,*) 'NC  =',G%NC
         WRITE(MSG%LU_DEBUG,*) 'NCE =',G%NCE
         WRITE(MSG%LU_DEBUG,*) 'N_WALL_CELLS  =',L%N_WALL_CELLS
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------G%ICE_TO_IWG:'
         WRITE(MSG%LU_DEBUG,'(8I8)') (G%ICE_TO_IWG(IW), IW = G%NC+1, G%NCE)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------G%ICE_TO_ICG:'
         WRITE(MSG%LU_DEBUG,'(8I8)') (G%ICE_TO_ICG(IW), IW = G%NC+1, G%NCE)
         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
            WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------G%ICE_TO_ICN:'
            WRITE(MSG%LU_DEBUG,'(8I8)') (G%ICE_TO_ICN(IW), IW = G%NC+1, G%NCE)
            WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------G%ICE_TO_VAL:'
            WRITE(MSG%LU_DEBUG,'(8e12.4)') (G%ICE_TO_VAL(IW), IW = G%NC+1, G%NCE)
         ENDIF
         WRITE(MSG%LU_DEBUG,*)
         IF (NL == 1) THEN
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXG:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IXG, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYG:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IYG, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZG:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IZG, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXW:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IXW, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYW:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IYW, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZW:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IZW, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXN(1):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IXN(1), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYN(1):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IYN(1), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZN(1):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IZN(1), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXN(2):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IXN(2), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYN(2):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IYN(2), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZN(2):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IZN(2), IW=1,L%N_WALL_CELLS)
         ENDIF
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%BTYPE:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%BTYPE, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IOR:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IOR, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%NOM:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%NOM, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICW:', NM
         DO IW=1,L%N_WALL_CELLS
            WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',G%WALL(IW)%ICW
         ENDDO
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICE:', NM
         DO IW=1,L%N_WALL_CELLS
            IF (G%WALL(IW)%NOM /=0) &
               WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',G%WALL(IW)%ICE
         ENDDO
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICG:', NM
         DO IW=1,L%N_WALL_CELLS
            IF (G%WALL(IW)%NOM /=0) &
               WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',G%WALL(IW)%ICG
         ENDDO
         WRITE(MSG%LU_DEBUG,*) '====================================================='
         WRITE(MSG%LU_DEBUG,*) ' Plotting out M%WALL-structure'
         WRITE(MSG%LU_DEBUG,*) '====================================================='
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%WALL_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%WALL_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%SURF_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%SURF_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%BACK_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%BACK_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%BOUNDARY_TYPE'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%BOUNDARY_TYPE, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%OBST_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%OBST_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%PRESSURE_BC_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%PRESSURE_BC_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%PRESSURE_ZONE'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%PRESSURE_ZONE, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%VENT_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%VENT_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%NOM'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%EXTERNAL_WALL(IW)%NOM, IW=1,L%N_WALL_CELLS_EXT)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%II'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%II, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%JJ'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%JJ, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%KK'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%KK, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%IIG'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%IIG, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%JJG'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%JJG, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%KKG'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%KKG, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%N_LAYER_CELLS'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%IOR, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%IOR'
         !WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%N_LAYER_CELLS, IW=1,L%N_WALL_CELLS)

      ENDDO
      !ENDIF

   !> ------------------------------------------------------------------------------------------------
   !> Debug complete grid information
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_GRID)
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_GRID(NM, NL)
         WRITE(MSG%LU_DEBUG,*) 'M%N_OBST=',M%N_OBST
         WRITE(MSG%LU_DEBUG,*) 'M%OBST ... I1, I2, J1, J2, K1, K2'
         DO IWG = 1, M%N_OBST
            WRITE(MSG%LU_DEBUG,'(6I8)') M%OBSTRUCTION(IWG)%I1,M%OBSTRUCTION(IWG)%I2,&
               M%OBSTRUCTION(IWG)%J1,M%OBSTRUCTION(IWG)%J2,&
               M%OBSTRUCTION(IWG)%K1,M%OBSTRUCTION(IWG)%K2
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'M%N_OBST=',M%N_OBST
         WRITE(MSG%LU_DEBUG,*) 'M%N_WALL_CELLS=',M%N_WALL_CELLS
         WRITE(MSG%LU_DEBUG,*) 'M%CELL_INDEX:'
         DO JJJ = L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (M%CELL_INDEX(III,JJJ,KKK), III=0,L%NX+1)
            ENDDO
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'M%WALL(.)%BOUNDARY_TYPE:'
         WRITE(MSG%LU_DEBUG,'(16i6)') (M%WALL(IWG)%BOUNDARY_TYPE, IWG=1, L%N_WALL_CELLS)

         WRITE(MSG%LU_DEBUG,*) 'CELL_COUNT:', CELL_COUNT(NM)

         WRITE(MSG%LU_DEBUG,*) 'M%WALL_INDEX:'
         WRITE(MSG%LU_DEBUG,'(i8,a,6i8)') ( IWG, ' : ', &
            M%WALL_INDEX(IWG, 1),  &
            M%WALL_INDEX(IWG,-1),  &
            M%WALL_INDEX(IWG, 2),  &
            M%WALL_INDEX(IWG,-2),  &
            M%WALL_INDEX(IWG, 3),  &
            M%WALL_INDEX(IWG,-3),  &
            IWG=1, CELL_COUNT(NM))

         WRITE(MSG%LU_DEBUG,*) 'M%WALL(.)%ONE_D% IOR,II,JJ,KK, BOUNDARY_TYPE, BTYPE, PRESSURE_BC_INDEX:'
         DO IWG = 1, L%N_WALL_CELLS
            WRITE(MSG%LU_DEBUG,'(9I8)') &
               IWG,M%WALL(IWG)%ONE_D%IOR,M%WALL(IWG)%ONE_D%II,M%WALL(IWG)%ONE_D%JJ,M%WALL(IWG)%ONE_D%KK,&
               M%WALL(IWG)%BOUNDARY_TYPE, G%WALL(IWG)%BTYPE, M%WALL(IWG)%PRESSURE_BC_INDEX
         ENDDO

         WRITE(MSG%LU_DEBUG,*) 'GRIG%CELL_NUMBER(...)'
         DO JJJ=L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (G%CELL_NUMBER(III,JJJ,KKK),III=0,L%NX+1)
            ENDDO
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'GRIL%IS_SOLID(...)'
         DO JJJ=L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (L%IS_SOLID(III,JJJ,KKK),III=0,L%NX+1)
            ENDDO
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'M%CELL_INDEX:'
         DO JJJ=L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (M%CELL_INDEX(III,JJJ,KKK),III=0,L%NX+1)
            ENDDO
         ENDDO
      ENDDO

   !> ------------------------------------------------------------------------------------------------
   !> Debug PRESSURE_BC_INDEX
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_BDRY)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         WRITE(MSG%LU_DEBUG, '(16I6)') (G%WALL(J)%BTYPE, J=1,L%N_WALL_CELLS)
      ENDDO

END SELECT

1000 FORMAT('======================================================================================',/, &
   '=== ', A30,' for mesh ',i3,' on level ', i3, /, &
   '======================================================================================')

END SUBROUTINE SCARC_DEBUG_QUANTITY


!> ------------------------------------------------------------------------------------------------
!> Print out vector information on level NL for matlab
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATLAB_VECTOR (NV, CVEC, NL)
USE SCARC_POINTERS, ONLY: G, VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM 
CHARACTER (*), INTENT(IN) :: CVEC
INTEGER :: JC, MVEC
CHARACTER(60) :: CNAME, CFORM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)
   WRITE (CNAME, '(A,A1,A,i2.2,A,i2.2,A)') 'matlab/',CVEC,'_mesh',NM,'_level',NL,'_vec.txt'
   WRITE (CFORM, '(I3, A)' ) G%NC-1, "(F7.2,;),F7.2"
   MVEC=GET_FILE_NUMBER()

   OPEN(MVEC,FILE=CNAME)
   WRITE(MVEC, *) CVEC, ' = ['
   WRITE(MVEC,'(8F12.2)') (VC(JC),JC=1,G%NC)
   WRITE(MVEC, *) ' ]'
   CLOSE(MVEC)

ENDDO

END SUBROUTINE SCARC_MATLAB_VECTOR


!> ------------------------------------------------------------------------------------------------
!> Print out matrix information on level NL for matlab
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATLAB_MATRIX(VAL, ROW, COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: VAL
INTEGER, DIMENSION(:), INTENT(IN) :: ROW
INTEGER, DIMENSION(:), INTENT(IN) :: COL
INTEGER, INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(*), INTENT(IN) :: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60) :: CFILE, CFORM
REAL(EB) :: MATRIX_LINE(1000)

WRITE (CFILE, '(A,A,A,i2.2,A,i2.2,A)') 'matlab/',TRIM(CNAME),'_mesh',NM,'_level',NL,'_mat.txt'
!WRITE (CFORM, '(A,I3, 2A)' ) "(", NC2-1, "(F9.3,','),F9.3,';')"
!WRITE (CFORM, '(A,I3, 2A)' ) "(", NC2-1, "(F9.3,' '),F9.3,' ')"
WRITE (CFORM, '(A,I3, 2A)' ) "(", NC2-1, "(F9.3,' '),F9.3,' ')"
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=CFILE)
!WRITE(MMATRIX, *) CNAME, ' = ['
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= ROW(IC), ROW(IC+1)-1
         !IF (COL(ICOL)==JC) MATRIX_LINE(JC)=VAL(ICOL)/25.0_EB
         IF (COL(ICOL)==JC) MATRIX_LINE(JC)=VAL(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX, CFORM) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO
!WRITE(MMATRIX, *) ' ]'
CLOSE(MMATRIX)

END SUBROUTINE SCARC_MATLAB_MATRIX

!> ================================================================================================
!> End  WITH_SCARC_DEBUG  - Part
!> ================================================================================================
#endif


#ifdef WITH_SCARC_STANDALONE
!> ================================================================================================
!> Begin  WITH_SCARC_STANDALONE  - PART
!> ================================================================================================
!> ------------------------------------------------------------------------------------------------
!> Dump matrix and vectors belonging to pressure system 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_SYSTEM (ITYPE)
USE SCARC_POINTERS, ONLY: G, P, AB
INTEGER, INTENT(IN) :: ITYPE
INTEGER :: NM, IC, ID, IX, IY, IZ

!IF (MOD(ITE_TOTAL,10) /= 0) RETURN

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)

   SELECT CASE(ITYPE)

      CASE(NSCARC_DUMP_AB)            ! print out banded matrix

         AB => G%AB
         MSG%LU_CNN = GET_FILE_NUMBER()
         WRITE (MSG%FILE_CNN, '(A,I5.5)') 'cnn/A_t',ITE_GLOBAL
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN)
         DO ID = 1, AB%N_STENCIL
            DO IC = 1, AB%N_DIAG
               WRITE(MSG%LU_CNN,*) AB%VAL(IC, ID)
            ENDDO
         ENDDO
         CLOSE(MSG%LU_CNN)

      CASE(NSCARC_DUMP_B_OLD)

         MSG%LU_CNN = GET_FILE_NUMBER()
         WRITE (MSG%FILE_CNN, '(A,I5.5)') 'cnn/b_old_t',ITE_GLOBAL
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN)
         DO IC = 1, G%NC
            WRITE(MSG%LU_CNN,*) P%B_OLD(IC)
         ENDDO
         CLOSE(MSG%LU_CNN)

      CASE(NSCARC_DUMP_B_NEW)

         MSG%LU_CNN = GET_FILE_NUMBER()
         WRITE (MSG%FILE_CNN, '(A,I5.5)') 'cnn/b_new_t',ITE_GLOBAL
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN)
         DO IC = 1, G%NC
            WRITE(MSG%LU_CNN,*) P%B_NEW(IC)
         ENDDO
         CLOSE(MSG%LU_CNN)

      CASE(NSCARC_DUMP_H_OLD)

         MSG%LU_CNN = GET_FILE_NUMBER()
         WRITE (MSG%FILE_CNN, '(A,I5.5)') 'cnn/h_old_t',ITE_GLOBAL
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN)
         DO IZ = 1, L%NZ
            DO IY = 1, L%NY
               DO IX = 1, L%NX
                  WRITE(MSG%LU_CNN,*) P%H_OLD(IX, IY, IZ)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(MSG%LU_CNN)

      CASE(NSCARC_DUMP_H_NEW)

         MSG%LU_CNN = GET_FILE_NUMBER()
         WRITE (MSG%FILE_CNN, '(A,I5.5)') 'cnn/h_new_t',ITE_GLOBAL
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN)
         DO IZ = 1, L%NZ
            DO IY = 1, L%NY
               DO IX = 1, L%NX
                  WRITE(MSG%LU_CNN,*) P%H_NEW(IX, IY, IZ)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(MSG%LU_CNN)

      CASE(NSCARC_DUMP_HS_OLD)

         MSG%LU_CNN = GET_FILE_NUMBER()
         WRITE (MSG%FILE_CNN, '(A,I5.5)') 'cnn/hs_old_t',ITE_GLOBAL
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN)
         DO IZ = 1, L%NZ
            DO IY = 1, L%NY
               DO IX = 1, L%NX
                  WRITE(MSG%LU_CNN,*) P%HS_OLD(IX, IY, IZ)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(MSG%LU_CNN)

      CASE(NSCARC_DUMP_HS_NEW)

         MSG%LU_CNN = GET_FILE_NUMBER()
         WRITE (MSG%FILE_CNN, '(A,I5.5)') 'cnn/hs_new_t',ITE_GLOBAL
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN)
         DO IZ = 1, L%NZ
            DO IY = 1, L%NY
               DO IX = 1, L%NX
                  WRITE(MSG%LU_CNN,*) P%HS_NEW(IX, IY, IZ)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(MSG%LU_CNN)

   END SELECT

ENDDO

END SUBROUTINE SCARC_DUMP_SYSTEM

!> ------------------------------------------------------------------------------------------------
!> Dump several arrays and structures needed for the standalone version of ScaRC 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_ENVIRONMENT
USE SCARC_POINTERS, ONLY: M
INTEGER :: NM, NOM, I, J, K, IW

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)

   ! Open SCARC-file for mesh NM
   WRITE (MSG%FILE_SCARC, '(A,A,A,i3.3)') 'dump/',TRIM(CHID),'_elk.dump',NM
   MSG%LU_SCARC = GET_FILE_NUMBER()
   OPEN (MSG%LU_SCARC, FILE=MSG%FILE_SCARC, ACTION = 'readwrite')

   !> Current process number
   WRITE(MSG%LU_SCARC,*) 'PROCESS(NM)'
   WRITE(MSG%LU_SCARC,*) PROCESS(NM)

   !> Dump information about all necessary MESHES-quantities
   WRITE(MSG%LU_SCARC,*) 'M%IBAR, M%JBAR, M%KBAR'
   WRITE(MSG%LU_SCARC,*) M%IBAR, M%JBAR, M%KBAR
   WRITE(MSG%LU_SCARC,*) 'M%IBP1, M%JBP1, M%KBP1'
   WRITE(MSG%LU_SCARC,*) M%IBP1, M%JBP1, M%KBP1
   WRITE(MSG%LU_SCARC,*) 'M%ITRN, M%JTRN, M%KTRN'
   WRITE(MSG%LU_SCARC,*) M%ITRN, M%JTRN, M%KTRN
   WRITE(MSG%LU_SCARC,*) 'M%LBC, M%MBC, M%NBC'
   WRITE(MSG%LU_SCARC,*) M%LBC, M%MBC, M%NBC
   WRITE(MSG%LU_SCARC,*) 'M%N_WALL_CELLS'
   WRITE(MSG%LU_SCARC,*) M%N_WALL_CELLS
   WRITE(MSG%LU_SCARC,*) 'M%N_EXTERNAL_WALL_CELLS'
   WRITE(MSG%LU_SCARC,*) M%N_EXTERNAL_WALL_CELLS
   WRITE(MSG%LU_SCARC,*) 'M%N_INTERNAL_WALL_CELLS'
   WRITE(MSG%LU_SCARC,*) M%N_INTERNAL_WALL_CELLS
   WRITE(MSG%LU_SCARC,*) 'M%IPS, M%LSAVE, M%LWORK'
   WRITE(MSG%LU_SCARC,*) M%IPS, M%LSAVE, M%LWORK
   WRITE(MSG%LU_SCARC,*) 'M%POIS_PTB'
   WRITE(MSG%LU_SCARC,*) M%POIS_PTB
   WRITE(MSG%LU_SCARC,*) 'M%XS, M%XF, M%YS, M%YF, M%ZS, M%ZF'
   WRITE(MSG%LU_SCARC,*) M%XS, M%XF, M%YS, M%YF, M%ZS, M%ZF
   WRITE(MSG%LU_SCARC,*) 'M%DXI, M%DETA, M%DZETA'
   WRITE(MSG%LU_SCARC,*) M%DXI, M%DETA, M%DZETA
   WRITE(MSG%LU_SCARC,*) 'M%RDXI, M%RDETA, M%RDZETA'
   WRITE(MSG%LU_SCARC,*) M%RDXI, M%RDETA, M%RDZETA
   WRITE(MSG%LU_SCARC,*) 'M%X'
   WRITE(MSG%LU_SCARC,*) (M%X(I), I=0, M%IBAR)
   WRITE(MSG%LU_SCARC,*) 'M%Y'
   WRITE(MSG%LU_SCARC,*) (M%Y(I), I=0, M%JBAR)
   WRITE(MSG%LU_SCARC,*) 'M%Z'
   WRITE(MSG%LU_SCARC,*) (M%Z(I), I=0, M%KBAR)
   WRITE(MSG%LU_SCARC,*) 'M%XC'
   WRITE(MSG%LU_SCARC,*) (M%XC(I), I=0, M%IBP1)
   WRITE(MSG%LU_SCARC,*) 'M%YC'
   WRITE(MSG%LU_SCARC,*) (M%YC(I), I=0, M%JBP1)
   WRITE(MSG%LU_SCARC,*) 'M%ZC'
   WRITE(MSG%LU_SCARC,*) (M%ZC(I), I=0, M%KBP1)
   WRITE(MSG%LU_SCARC,*) 'M%HX'
   WRITE(MSG%LU_SCARC,*) (M%HX(I), I=0, M%IBP1)
   WRITE(MSG%LU_SCARC,*) 'M%HY'
   WRITE(MSG%LU_SCARC,*) (M%HY(I), I=0, M%JBP1)
   WRITE(MSG%LU_SCARC,*) 'M%HZ'
   WRITE(MSG%LU_SCARC,*) (M%HZ(I), I=0, M%KBP1)

   DO NOM = 1, NMESHES
      WRITE(MSG%LU_SCARC,*) 'M%OMESH(',NOM,')%NIC_R, NICS'
      WRITE(MSG%LU_SCARC,*) M%OMESH(NOM)%NIC_R, M%OMESH(NOM)%NIC_S
   ENDDO

   IF (PREDICTOR) THEN
      WRITE(MSG%LU_SCARC,*) 'M%H'
      WRITE(MSG%LU_SCARC,*) (((M%H(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)
   ELSE
      WRITE(MSG%LU_SCARC,*) 'M%HS'
      WRITE(MSG%LU_SCARC,*) (((M%HS(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)
   ENDIF

   WRITE(MSG%LU_SCARC,*) 'M%PRHS'
   IF (TWO_D) THEN
      WRITE(MSG%LU_SCARC,*) ((M%PRHS(I,1,K), I=1, M%IBP1), K=1, M%KBP1)
   ELSE
      WRITE(MSG%LU_SCARC,*) (((M%PRHS(I,J,K), I=1, M%IBP1), J=1, M%JBP1), K=1, M%KBP1)
   ENDIF

   IF (TWO_D) THEN
      WRITE(MSG%LU_SCARC,*) 'M%BXS'
      WRITE(MSG%LU_SCARC,*) (M%BXS(1,K), K=1, M%KBP1)
      WRITE(MSG%LU_SCARC,*) 'M%BXF'
      WRITE(MSG%LU_SCARC,*) (M%BXF(1,K), K=1, M%KBP1)
   ELSE
      WRITE(MSG%LU_SCARC,*) 'M%BXS'
      WRITE(MSG%LU_SCARC,*) ((M%BXS(J,K), J=1, M%JBP1), K=1, M%KBP1)
      WRITE(MSG%LU_SCARC,*) 'M%BXF'
      WRITE(MSG%LU_SCARC,*) ((M%BXF(J,K), J=1, M%JBP1), K=1, M%KBP1)
   ENDIF

   WRITE(MSG%LU_SCARC,*) 'M%BYS'
   WRITE(MSG%LU_SCARC,*) ((M%BYS(I,K), I=1, M%IBP1), K=1, M%KBP1)
   WRITE(MSG%LU_SCARC,*) 'M%BYF'
   WRITE(MSG%LU_SCARC,*) ((M%BYF(I,K), I=1, M%IBP1), K=1, M%KBP1)

   IF (TWO_D) THEN
      WRITE(MSG%LU_SCARC,*) 'M%BZS'
      WRITE(MSG%LU_SCARC,*) (M%BZS(I,1), I=1, M%IBP1)
      WRITE(MSG%LU_SCARC,*) 'M%BZF'
      WRITE(MSG%LU_SCARC,*) (M%BZF(I,1), I=1, M%IBP1)
   ELSE
      WRITE(MSG%LU_SCARC,*) 'M%BZS'
      WRITE(MSG%LU_SCARC,*) ((M%BZS(I,J), I=1, M%IBP1), J=1, M%JBP1)
      WRITE(MSG%LU_SCARC,*) 'M%BZF'
      WRITE(MSG%LU_SCARC,*) ((M%BZF(I,J), I=1, M%IBP1), J=1, M%JBP1)
   ENDIF

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%PRESSURE_BC_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%PRESSURE_BC_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%BOUNDARY_TYPE'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%BOUNDARY_TYPE, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%WALL_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%WALL_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%SURF_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%SURF_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%BACK_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%BACK_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%OBST_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%OBST_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%VENT_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%VENT_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%IOR'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%IOR, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%II'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%II, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%JJ'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%JJ, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%KK'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%KK, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%IIG'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%IIG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%JJG'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%JJG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%KKG'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%KKG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%NOM'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%NOM, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%IIO_MIN'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%IIO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%IIO_MAX'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%IIO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%JJO_MIN'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%JJO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%JJO_MAX'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%JJO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%KKO_MIN'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%KKO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%KKO_MAX'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%KKO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%N_OBST'
   WRITE(MSG%LU_SCARC,*) M%N_OBST
   IF (M%N_OBST > 0) THEN
      WRITE(MSG%LU_SCARC,*) 'M%OBSTRUCTION(.)%I1, I2, J1, J2, K1, K2'
      DO I = 1, M%N_OBST
         WRITE(MSG%LU_SCARC,*) M%OBSTRUCTION(I)%I1, M%OBSTRUCTION(I)%I2, &
                               M%OBSTRUCTION(I)%J1, M%OBSTRUCTION(I)%J2, &
                               M%OBSTRUCTION(I)%K1, M%OBSTRUCTION(I)%K2
      ENDDO
   ENDIF

   WRITE(MSG%LU_SCARC,*) 'CELL_COUNT(',NM,')'
   WRITE(MSG%LU_SCARC,*) CELL_COUNT(NM)

   WRITE(MSG%LU_SCARC,*) 'M%SOLID'
   WRITE(MSG%LU_SCARC,*) (M%SOLID(I), I=1, CELL_COUNT(NM))

   WRITE(MSG%LU_SCARC,*) 'M%CELL_INDEX'
   WRITE(MSG%LU_SCARC,*) (((M%CELL_INDEX(I,J,K), I=0, M%IBAR+1), J=0, M%JBAR+1), K=0, M%KBAR+1)

   WRITE(MSG%LU_SCARC,*) 'M%WALL_INDEX'
   DO I = 0, CELL_COUNT(NM)
      WRITE(MSG%LU_SCARC,'(7I8)') (M%WALL_INDEX(I,J), J=-3,3)
   ENDDO

   WRITE(MSG%LU_SCARC,*) 'M%SAVE1'
   WRITE(MSG%LU_SCARC,*) (M%SAVE1(I), I=-3, M%LSAVE)

   WRITE(MSG%LU_SCARC,*) 'M%WORK'
   WRITE(MSG%LU_SCARC,*) (M%WORK(I), I=1, M%LWORK)

   CLOSE (MSG%LU_SCARC)

ENDDO

SCARC_DUMP = .FALSE.
END SUBROUTINE SCARC_DUMP_ENVIRONMENT


!> ------------------------------------------------------------------------------------------------
!> Dump several arrays and structures needed for ScaRC - standalone version
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_ENVIRONMENT
USE SCARC_POINTERS, ONLY: M
INTEGER :: NM, NOM, I, J, K, IW, LU, IERR=0

!CPATH = 'D:\GIT\HHP\A_ScaRC\VisualStudio\Cases\'
!CPATH = '../VisualStudio/Cases/'

ALLOCATE(CELL_COUNT(NMESHES), STAT= IERR)
CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'CELL_COUNT', IERR)

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)

   MSG%LU_SCARC = 99
   WRITE (MSG%FILE_SCARC, '(A,A,A,i3.3)') 'dump/',TRIM(CHID),'.dump',NM
   write(*,*) MSG%FILE_SCARC
   MSG%FILE_SCARC = TRIM(MSG%FILE_SCARC)

   LU = MSG%LU_SCARC
   OPEN (LU, FILE=MSG%FILE_SCARC)

   !> Current process number
   READ(LU,*)
   READ(LU,*) PROCESS(NM)

   !> Dump information about all necessary MESHES-quantities
   READ(LU,*)
   READ(LU,*) M%IBAR, M%JBAR, M%KBAR
   READ(LU,*) 
   READ(LU,*) M%IBP1, M%JBP1, M%KBP1
   READ(LU,*) 
   READ(LU,*) M%ITRN, M%JTRN, M%KTRN
   READ(LU,*) 
   READ(LU,*) M%LBC, M%MBC, M%NBC
   READ(LU,*) 
   READ(LU,*) M%N_WALL_CELLS
   READ(LU,*)
   READ(LU,*) M%N_EXTERNAL_WALL_CELLS
   READ(LU,*)
   READ(LU,*) M%N_INTERNAL_WALL_CELLS
   READ(LU,*)
   READ(LU,*) M%IPS, M%LSAVE, M%LWORK
   READ(LU,*) 
   READ(LU,*) M%POIS_PTB
   READ(LU,*)
   READ(LU,*) M%XS, M%XF, M%YS, M%YF, M%ZS, M%ZF
   READ(LU,*) 
   READ(LU,*) M%DXI, M%DETA, M%DZETA
   READ(LU,*) 
   READ(LU,*) M%RDXI, M%RDETA, M%RDZETA

   ALLOCATE(M%X(0:M%IBAR), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%X', IERR)

   ALLOCATE(M%Y(0:M%JBAR), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%Y', IERR)

   ALLOCATE(M%Z(0:M%KBAR), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%Z', IERR)

   ALLOCATE(M%XC(0:M%IBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%XC', IERR)

   ALLOCATE(M%YC(0:M%JBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%YC', IERR)

   ALLOCATE(M%ZC(0:M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%ZC', IERR)

   ALLOCATE(M%HX(0:M%IBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%HX', IERR)

   ALLOCATE(M%HY(0:M%JBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%HY', IERR)

   ALLOCATE(M%HZ(0:M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%HZ', IERR)

   READ(LU,*)
   READ(LU,*) (M%X(I), I=0, M%IBAR)
   READ(LU,*)
   READ(LU,*) (M%Y(I), I=0, M%JBAR)
   READ(LU,*)
   READ(LU,*) (M%Z(I), I=0, M%KBAR)
   READ(LU,*)
   READ(LU,*) (M%XC(I), I=0, M%IBP1)
   READ(LU,*)
   READ(LU,*) (M%YC(I), I=0, M%JBP1)
   READ(LU,*)
   READ(LU,*) (M%ZC(I), I=0, M%KBP1)

   READ(LU,*) 
   READ(LU,*) (M%HX(I), I=0, M%IBP1)
   READ(LU,*) 
   READ(LU,*) (M%HY(I), I=0, M%JBP1)
   READ(LU,*) 
   READ(LU,*) (M%HZ(I), I=0, M%KBP1)

   ALLOCATE(M%OMESH(NMESHES), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%OMESH', IERR)

   DO NOM = 1, NMESHES
      READ(LU,*)
      READ(LU,*) M%OMESH(NOM)%NIC_R, M%OMESH(NOM)%NIC_S
   ENDDO

   IF (PREDICTOR) THEN

      ALLOCATE(M%H(0:M%IBP1, 0:M%JBP1, 0:M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%H', IERR)
      M%H = 0.0_EB

      READ(LU,*)
      READ(LU,*) (((M%H(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)

   ELSE

      ALLOCATE(M%HS(0:M%IBP1, 0:M%JBP1, 0:M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%H', IERR)
      M%HS = 0.0_EB

      READ(LU,*)
      READ(LU,*) (((M%HS(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)

   ENDIF
   IF (TWO_D) THEN

      ALLOCATE(M%PRHS(M%IBP1, 1, M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%H', IERR)

      READ(LU,*)
      READ(LU,*) ((M%PRHS(I,1,K), I=1, M%IBP1), K=1, M%KBP1)

   ELSE

      ALLOCATE(M%PRHS(M%IBP1, M%JBP1, M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%H', IERR)

      READ(LU,*)
      READ(LU,*) (((M%PRHS(I,J,K), I=1, M%IBP1), J=1, M%JBP1), K=1, M%KBP1)

   ENDIF

   IF (TWO_D) THEN

      ALLOCATE(M%BXS(1,M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BXS', IERR)

      ALLOCATE(M%BXF(1,M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BXF', IERR)

      ALLOCATE(M%BZS(M%IBP1,1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BZS', IERR)

      ALLOCATE(M%BZF(M%IBP1,1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BZF', IERR)

   ELSE

      ALLOCATE(M%BXS(M%JBP1,M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BXS', IERR)

      ALLOCATE(M%BXF(M%JBP1,M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BXF', IERR)

      ALLOCATE(M%BZS(M%IBP1,M%JBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BZS', IERR)

      ALLOCATE(M%BZF(M%IBP1,M%JBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BZF', IERR)

   ENDIF
   
   ALLOCATE(M%BYS(M%IBP1,M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%BYS', IERR)

   ALLOCATE(M%BYF(M%IBP1,M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%BYF', IERR)

   IF (TWO_D) THEN
      READ(LU,*)
      READ(LU,*) (M%BXS(1,K), K=1, M%KBP1)
      READ(LU,*)
      READ(LU,*) (M%BXF(1,K), K=1, M%KBP1)
   ELSE
      READ(LU,*)
      READ(LU,*) ((M%BXS(J,K), J=1, M%JBP1), K=1, M%KBP1)
      READ(LU,*)
      READ(LU,*) ((M%BXF(J,K), J=1, M%JBP1), K=1, M%KBP1)
   ENDIF
   
   READ(LU,*)
   READ(LU,*) ((M%BYS(I,K), I=1, M%IBP1), K=1, M%KBP1)
   READ(LU,*)
   READ(LU,*) ((M%BYF(I,K), I=1, M%IBP1), K=1, M%KBP1)

   IF (TWO_D) THEN
      READ(LU,*)
      READ(LU,*) (M%BZS(I,1), I=1, M%IBP1)
      READ(LU,*)
      READ(LU,*) (M%BZF(I,1), I=1, M%IBP1)
   ELSE
      READ(LU,*)
      READ(LU,*) ((M%BZS(I,J), I=1, M%IBP1), J=1, M%JBP1)
      READ(LU,*)
      READ(LU,*) ((M%BZF(I,J), I=1, M%IBP1), J=1, M%JBP1)
   ENDIF

   ALLOCATE(M%WALL(M%N_WALL_CELLS), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%WALL', IERR)

   ALLOCATE(M%EXTERNAL_WALL(M%N_EXTERNAL_WALL_CELLS), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%EXTERNAL_WALL', IERR)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%PRESSURE_BC_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%BOUNDARY_TYPE, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%WALL_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%SURF_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%BACK_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%OBST_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%VENT_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%IOR, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%II, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%JJ, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%KK, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%IIG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%JJG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%KKG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*) 
   READ(LU,*) (M%EXTERNAL_WALL(IW)%NOM, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%EXTERNAL_WALL(IW)%IIO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%EXTERNAL_WALL(IW)%IIO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*) 
   READ(LU,*) (M%EXTERNAL_WALL(IW)%JJO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*) 
   READ(LU,*) (M%EXTERNAL_WALL(IW)%JJO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*) 
   READ(LU,*) (M%EXTERNAL_WALL(IW)%KKO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%EXTERNAL_WALL(IW)%KKO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) M%N_OBST

   IF (M%N_OBST > 0) THEN

      ALLOCATE(M%OBSTRUCTION(M%N_OBST), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%OBSTRUCTION', IERR)

      READ(LU,*)
      DO I = 1, M%N_OBST
         READ(LU,*) M%OBSTRUCTION(I)%I1, M%OBSTRUCTION(I)%I2, &
                    M%OBSTRUCTION(I)%J1, M%OBSTRUCTION(I)%J2, &
                    M%OBSTRUCTION(I)%K1, M%OBSTRUCTION(I)%K2
      ENDDO
   ENDIF

   READ(LU,*)
   READ(LU,*) CELL_COUNT(NM)

   ALLOCATE(M%SOLID(0:CELL_COUNT(NM)), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%SOLID', IERR)
   M%SOLID=.FALSE.

   ALLOCATE(M%CELL_INDEX(0:M%IBP1,0:M%JBP1,0:M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%CELL_INDEX', IERR)
   M%CELL_INDEX = 0

   ALLOCATE(M%WALL_INDEX(0:CELL_COUNT(NM),-3:3), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%WALL_INDEX', IERR)

   READ(LU,*)
   READ(LU,*) (M%SOLID(I), I=1, CELL_COUNT(NM))

   READ(LU,*)
   READ(LU,*) (((M%CELL_INDEX(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)

   READ(LU,*)
   DO I = 0, CELL_COUNT(NM)
      READ(LU,'(7I8)') (M%WALL_INDEX(I,J), J=-3,3)
   ENDDO

   ALLOCATE(M%SAVE1(-3:M%LSAVE), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%SAVE1', IERR)

   ALLOCATE(M%WORK(M%LWORK), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%WORK', IERR)

   READ(LU,*) 
   READ(LU,*) (M%SAVE1(I), I=-3, M%LSAVE)

   READ(LU,*) 
   READ(LU,*) (M%WORK(I), I=1, M%LWORK)

   CLOSE (LU)

ENDDO

END SUBROUTINE SCARC_RESTORE_ENVIRONMENT


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors pressure diagnostics - only for developping purposes
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PRESSURE()
USE SCARC_POINTERS, ONLY: L, P
INTEGER :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)

   CALL SCARC_ALLOCATE_REAL1(P%B_OLD, 1, G%NC, NSCARC_INIT_ZERO, 'P%B_OLD')
   CALL SCARC_ALLOCATE_REAL1(P%B_NEW, 1, G%NC, NSCARC_INIT_ZERO, 'P%B_NEW')

   CALL SCARC_ALLOCATE_REAL3(P%H_OLD, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'P%H_OLD')
   CALL SCARC_ALLOCATE_REAL3(P%H_NEW, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'P%H_NEW')

   CALL SCARC_ALLOCATE_REAL3(P%HS_OLD, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'P%HS_OLD')
   CALL SCARC_ALLOCATE_REAL3(P%HS_NEW, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'P%HS_NEW')

   !CALL SCARC_ALLOCATE_REAL3(P%PRHS, 1, L%NX+1, 1, L%NY+1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'P%PRHS')

ENDDO

END SUBROUTINE SCARC_SETUP_PRESSURE


!> ------------------------------------------------------------------------------------------------
!> Compute Differences between old and new pressure solutions - only for developping purposes
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESSURE_DIFFERENCE(NL)
USE SCARC_POINTERS, ONLY: L, P
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IX, IY, IZ

!> -----------------------------------------------------------------------------------------------
!> Store new pressure vector for comparison in next time step
!> -----------------------------------------------------------------------------------------------
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NL)

   IF (PREDICTOR) THEN
   
      P%H_NEW  = M%H
      P%DIFF_H = 0.0_EB
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
               P%DIFF_H = P%DIFF_H + (P%H_NEW(IX, IY, IZ) - P%H_OLD(IX, IY, IZ))**2         
            ENDDO
         ENDDO
      ENDDO
      P%DIFF_H = P%DIFF_H / REAL(L%N_CELLS, EB)

   ELSE

      P%HS_NEW = M%HS
      P%DIFF_HS = 0.0_EB
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
               P%DIFF_HS = P%DIFF_HS + (P%HS_NEW(IX, IY, IZ) - P%HS_OLD(IX, IY, IZ))**2         
            ENDDO
         ENDDO
      ENDDO
      P%DIFF_HS = P%DIFF_HS / REAL(L%N_CELLS, EB)

   ENDIF

ENDDO

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'Differences of pressure vectors on mesh ', NM,' : ', P%DIFF_H, P%DIFF_HS
#endif

END SUBROUTINE SCARC_PRESSURE_DIFFERENCE

!> ================================================================================================
!> End  WITH_SCARC_STANDALONE  - PART
!> ================================================================================================
#endif

END MODULE SCRC




