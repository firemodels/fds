!===========================================================================================================================
!
!> \brief Scalable Recursive Clustering (ScaRC): Collection of alternative pressure solvers for the FDS pressure equation
!
!===========================================================================================================================
#define WITH_SCARC_AMG
#undef WITH_SCARC_VERBOSE


!=======================================================================================================================
!
! MODULE SCARC_CONSTANTS
!
!> \brief Define all constants needed in for the different solution strategies in ScaRC/UScaRC
!
!=======================================================================================================================
MODULE SCARC_CONSTANTS

USE PRECISION_PARAMETERS

IMPLICIT NONE (TYPE,EXTERNAL)

INTEGER, PARAMETER :: NSCARC_ACCURACY_ABSOLUTE       =  1        !< Type of requested accuracy of method: absolute
INTEGER, PARAMETER :: NSCARC_ACCURACY_RELATIVE       =  2        !< Type of requested accuracy of method: relative
                   
INTEGER, PARAMETER :: NSCARC_BUFFER_BASIC            =  1        !< Length of exchange buffer: basic initialization
INTEGER, PARAMETER :: NSCARC_BUFFER_FULL             =  2        !< Length of exchange buffer: full length
INTEGER, PARAMETER :: NSCARC_BUFFER_LAYER1           =  3        !< Length of exchange buffer: one ghost cell layer
INTEGER, PARAMETER :: NSCARC_BUFFER_LAYER2           =  4        !< Length of exchange buffer: two ghost cell layers
INTEGER, PARAMETER :: NSCARC_BUFFER_LAYER4           =  5        !< Length of exchange buffer: four ghost cell layers
INTEGER, PARAMETER :: NSCARC_BUFFER_STENCIL          =  6        !< Length of exchange buffer: stencil size

INTEGER, PARAMETER :: NSCARC_COARSE_ITERATIVE        =  1        !< Type of coarse grid solver: iterative solver
INTEGER, PARAMETER :: NSCARC_COARSE_DIRECT           =  2        !< Type of coarse grid solver: direct solver
INTEGER, PARAMETER :: NSCARC_COARSE_LEVEL            =  3        !< Type of coarse grid solver: only on specified level
INTEGER, PARAMETER :: NSCARC_COARSE_MACRO            =  4        !< Type of coarse grid solver: on a separate macro

INTEGER, PARAMETER :: NSCARC_COARSENING_AGGREGATED   =  1        !< Type of grid coarsening: aggregation-based 
INTEGER, PARAMETER :: NSCARC_COARSENING_AGGREGATED_S =  2        !< Type of grid coarsening: aggregation-based, staggered 
INTEGER, PARAMETER :: NSCARC_COARSENING_CUBIC        =  3        !< Type of grid coarsening: cubic zones 
INTEGER, PARAMETER :: NSCARC_COARSENING_GMG          =  4        !< Type of grid coarsening: GMG-like 

INTEGER, PARAMETER :: NSCARC_CYCLING_F               =  0        !< Type of MG-cycling: F-cycle 
INTEGER, PARAMETER :: NSCARC_CYCLING_V               =  1        !< Type of MG-cycling: V-cycle 
INTEGER, PARAMETER :: NSCARC_CYCLING_W               =  2        !< Type of MG-cycling: W-cycle 
INTEGER, PARAMETER :: NSCARC_CYCLING_FMG             =  3        !< Type of MG-cycling: Full multigrid grid cycle 
INTEGER, PARAMETER :: NSCARC_CYCLING_SETUP           =  4        !< State of MG-cycling: initialize cycle counts
INTEGER, PARAMETER :: NSCARC_CYCLING_RESET           =  5        !< State of MG-cycling: reset cycle counts
INTEGER, PARAMETER :: NSCARC_CYCLING_PROCEED         =  6        !< State of MG-cycling: proceed cycle counts
INTEGER, PARAMETER :: NSCARC_CYCLING_NEXT            =  7        !< State of MG-cycling: perform next cycling loop
INTEGER, PARAMETER :: NSCARC_CYCLING_EXIT            =  8        !< State of MG-cycling: exit cycling loop
INTEGER, PARAMETER :: NSCARC_CYCLING_PRESMOOTH       = -1        !< State of MG-cycling: presmoothing cycle
INTEGER, PARAMETER :: NSCARC_CYCLING_POSTSMOOTH      =  1        !< State of MG-cycling: postsmoothing cycle

INTEGER, PARAMETER :: NSCARC_DATA_BMATRIX            =  1        !< Type of allocated structure: bandwise stored matrix
INTEGER, PARAMETER :: NSCARC_DATA_CMATRIX            =  2        !< Type of allocated structure: compactly stored matrix
INTEGER, PARAMETER :: NSCARC_DATA_INTEGER            =  3        !< Type of allocated structure: integer array
INTEGER, PARAMETER :: NSCARC_DATA_LOGICAL            =  4        !< Type of allocated structure: integer array
INTEGER, PARAMETER :: NSCARC_DATA_REAL_EB            =  5        !< Type of allocated structure: double precision array
INTEGER, PARAMETER :: NSCARC_DATA_REAL_FB            =  6        !< Type of allocated structure: single precision array

INTEGER, PARAMETER :: NSCARC_ERROR_BOUNDARY_SUM      =  3        !< Type of error message: wrong sum of boundary elements
INTEGER, PARAMETER :: NSCARC_ERROR_BOUNDARY_TYPE     =  4        !< Type of error message: wrong boundary type
INTEGER, PARAMETER :: NSCARC_ERROR_DIRECT_NOMKL      =  5        !< Type of error message: MKL for direct solver missing
INTEGER, PARAMETER :: NSCARC_ERROR_EXCHANGE_RECV     =  6        !< Type of error message: wrong receive exchange structure
INTEGER, PARAMETER :: NSCARC_ERROR_EXCHANGE_DIAG     =  7        !< Type of error message: wrong receive exchange structure
INTEGER, PARAMETER :: NSCARC_ERROR_EXCHANGE_SEND     =  8        !< Type of error message: wrong send exchange structure
INTEGER, PARAMETER :: NSCARC_ERROR_FFT_GRID          =  9        !< Type of error message: no unstructured FFT possible
INTEGER, PARAMETER :: NSCARC_ERROR_GRID_INDEX        = 10        !< Type of error message: error with grid index
INTEGER, PARAMETER :: NSCARC_ERROR_GRID_NUMBER       = 11        !< Type of error message: uneven cell number
INTEGER, PARAMETER :: NSCARC_ERROR_GRID_NUMBERX      = 12        !< Type of error message: uneven cell number in x 
INTEGER, PARAMETER :: NSCARC_ERROR_GRID_NUMBERY      = 13        !< Type of error message: uneven cell number in y 
INTEGER, PARAMETER :: NSCARC_ERROR_GRID_NUMBERZ      = 14        !< Type of error message: uneven cell number in z
INTEGER, PARAMETER :: NSCARC_ERROR_GRID_RESOLUTION   = 15        !< Type of error message: error with grid resolution
INTEGER, PARAMETER :: NSCARC_ERROR_NEIGHBOR_NUMBER   = 16        !< Type of error message: wrong neighbor number
INTEGER, PARAMETER :: NSCARC_ERROR_NEIGHBOR_TYPE     = 17        !< Type of error message: wrong neighbor type
INTEGER, PARAMETER :: NSCARC_ERROR_MATRIX_ALLOCATION = 18        !< Type of error message: error in matrix allocation
INTEGER, PARAMETER :: NSCARC_ERROR_MATRIX_COPY       = 19        !< Type of error message: subdiagonal missing
INTEGER, PARAMETER :: NSCARC_ERROR_MATRIX_SETUP      = 20        !< Type of error message: error in matrix setup
INTEGER, PARAMETER :: NSCARC_ERROR_MATRIX_SIZE       = 21        !< Type of error message: error in matrix size
INTEGER, PARAMETER :: NSCARC_ERROR_MATRIX_SUBDIAG    = 22        !< Type of error message: subdiagonal missing
INTEGER, PARAMETER :: NSCARC_ERROR_MATRIX_SYMMETRY   = 23        !< Type of error message: matrix not symmetric
INTEGER, PARAMETER :: NSCARC_ERROR_MGM_PARDISO       = 24        !< Type of error message: no PARDISO solver for MGM-Laplace 
INTEGER, PARAMETER :: NSCARC_ERROR_MGM_PERMUTATION   = 25        !< Type of error message: error in MGM LU permutation
INTEGER, PARAMETER :: NSCARC_ERROR_MKL_CLUSTER       = 26        !< Type of error message: CLUSTER_SPARSE_SOLVER missing
INTEGER, PARAMETER :: NSCARC_ERROR_MKL_INTERNAL      = 27        !< Type of error message: internal error in MKL routine
INTEGER, PARAMETER :: NSCARC_ERROR_MKL_PARDISO       = 28        !< Type of error message: PARDISO solver missing
INTEGER, PARAMETER :: NSCARC_ERROR_MKL_STACK         = 29        !< Type of error message: MKL method not available on this stack
INTEGER, PARAMETER :: NSCARC_ERROR_MKL_STORAGE       = 30        !< Type of error message: wrong storage scheme in MKL
INTEGER, PARAMETER :: NSCARC_ERROR_MULTIGRID_LEVEL   = 31        !< Type of error message: wrong multigrid level
INTEGER, PARAMETER :: NSCARC_ERROR_PARSE_INPUT       = 32        !< Type of error message: wrong input parameter
INTEGER, PARAMETER :: NSCARC_ERROR_STACK_MESSAGE     = 33        !< Type of error message: error with stack message
INTEGER, PARAMETER :: NSCARC_ERROR_STACK_SOLVER      = 34        !< Type of error message: error in solver stack
INTEGER, PARAMETER :: NSCARC_ERROR_STENCIL           = 35        !< Type of error message: error in matrix stencil
INTEGER, PARAMETER :: NSCARC_ERROR_VECTOR_LENGTH     = 36        !< Type of error message: error in vector length

INTEGER, PARAMETER :: NSCARC_EXCHANGE_AUXILIARY      =  1        !< Type of data exchange: various auxiliary data 
INTEGER, PARAMETER :: NSCARC_EXCHANGE_BASIC_SIZES    =  2        !< Type of data exchange: basic sizes during setup
INTEGER, PARAMETER :: NSCARC_EXCHANGE_CELL_NEIGHBORS =  3        !< Type of data exchange: neighboring cells
INTEGER, PARAMETER :: NSCARC_EXCHANGE_CELL_NUMBERS   =  4        !< Type of data exchange: neighboring cell numbers
INTEGER, PARAMETER :: NSCARC_EXCHANGE_CELL_SIZES     =  5        !< Type of data exchange: neighboring cell sizes
INTEGER, PARAMETER :: NSCARC_EXCHANGE_LAYER2_NUMS    =  6        !< Type of data exchange: numbers of second layer cells
INTEGER, PARAMETER :: NSCARC_EXCHANGE_LAYER2_VALS    =  7        !< Type of data exchange: values of second layer cells
INTEGER, PARAMETER :: NSCARC_EXCHANGE_MATRIX_COLS    =  8        !< Type of data exchange: (local) columns of Poisson matrix
INTEGER, PARAMETER :: NSCARC_EXCHANGE_MATRIX_COLSG   =  9        !< Type of data exchange: global columns of Poisson matrix
INTEGER, PARAMETER :: NSCARC_EXCHANGE_MATRIX_DIAGS   = 10        !< Type of data exchange: diagonal of Poisson matrix
INTEGER, PARAMETER :: NSCARC_EXCHANGE_MATRIX_SIZES   = 11        !< Type of data exchange: size of Poisson matrix 
INTEGER, PARAMETER :: NSCARC_EXCHANGE_MATRIX_VALS    = 12        !< Type of data exchange: values of Poisson matrix
INTEGER, PARAMETER :: NSCARC_EXCHANGE_MGM_SINGLE     = 13        !< Type of data exchange: MGM - Mean value interface
INTEGER, PARAMETER :: NSCARC_EXCHANGE_MGM_DOUBLE     = 14        !< Type of data exchange: MGM - True approximate interface
INTEGER, PARAMETER :: NSCARC_EXCHANGE_MGM_VELO       = 15        !< Type of data exchange: MGM - Velocity interface
INTEGER, PARAMETER :: NSCARC_EXCHANGE_MGM_TRUE       = 16        !< Type of data exchange: MGM - Velocity interface
INTEGER, PARAMETER :: NSCARC_EXCHANGE_NULLSPACE      = 17        !< Type of data exchange: nullspace entries (AMG only)
INTEGER, PARAMETER :: NSCARC_EXCHANGE_PRESSURE       = 18        !< Type of data exchange: pressure values
INTEGER, PARAMETER :: NSCARC_EXCHANGE_SOLIDS         = 19        !< Type of data exchange: pressure values
INTEGER, PARAMETER :: NSCARC_EXCHANGE_VECTOR_MEAN    = 20        !< Type of data exchange: mean values of a vector
INTEGER, PARAMETER :: NSCARC_EXCHANGE_VECTOR_PLAIN   = 21        !< Type of data exchange: plain values of a vector
INTEGER, PARAMETER :: NSCARC_EXCHANGE_ZONE_NEIGHBORS = 22        !< Type of data exchange: aggregation zones (AMG only)
INTEGER, PARAMETER :: NSCARC_EXCHANGE_ZONE_TYPES     = 23        !< Type of data exchange: aggregation zones types (AMG only)

INTEGER, PARAMETER :: NSCARC_GRID_STRUCTURED         =  1        !< Type of discretization: structured 
INTEGER, PARAMETER :: NSCARC_GRID_UNSTRUCTURED       =  2        !< Type of discretization: unstructured 

INTEGER, PARAMETER :: NSCARC_INIT_UNDEF              =-999       !< Type of data allocation: initialize as undefined
INTEGER, PARAMETER :: NSCARC_INIT_NONE               =  -2       !< Type of data allocation: do not initialize
INTEGER, PARAMETER :: NSCARC_INIT_MINUS              =  -1       !< Type of data allocation: initialize with minus one
INTEGER, PARAMETER :: NSCARC_INIT_ZERO               =   0       !< Type of data allocation: initialize with zero
INTEGER, PARAMETER :: NSCARC_INIT_ONE                =   1       !< Type of data allocation: initialize with one
INTEGER, PARAMETER :: NSCARC_INIT_TRUE               =   2       !< Type of data allocation: initialize with .TRUE.
INTEGER, PARAMETER :: NSCARC_INIT_FALSE              =   3       !< Type of data allocation: initialize with .FALSE.
INTEGER, PARAMETER :: NSCARC_INIT_HUGE               =   4       !< Type of data allocation: initialize with .FALSE.

INTEGER, PARAMETER :: NSCARC_INTERPOL_BILINEAR       =  1        !< Type of grid interpolation: bilinear 
INTEGER, PARAMETER :: NSCARC_INTERPOL_CONSTANT       =  2        !< Type of grid interpolation: constant 
INTEGER, PARAMETER :: NSCARC_INTERPOL_CLASSICAL      =  3        !< Type of grid interpolation: classical
INTEGER, PARAMETER :: NSCARC_INTERPOL_DIRECT         =  4        !< Type of grid interpolation: direct 
INTEGER, PARAMETER :: NSCARC_INTERPOL_STANDARD       =  5        !< Type of grid interpolation: standard 
INTEGER, PARAMETER :: NSCARC_INTERPOL_AMG            =  6        !< Type of grid interpolation: AMG-defined
INTEGER, PARAMETER :: NSCARC_INTERPOL_BILINEAR2      =  7        !< Type of grid interpolation: bilinear 

INTEGER, PARAMETER :: NSCARC_KRYLOV_MAIN             =  1        !< Type of Krylov solver: use it as main solver
INTEGER, PARAMETER :: NSCARC_KRYLOV_COARSE           =  2        !< Type of Krylov solver: use it as coarse grid solver

INTEGER, PARAMETER :: NSCARC_LEVEL_MIN               =  0        !< Range of multigrid levels: minimum level
INTEGER, PARAMETER :: NSCARC_LEVEL_MAX               = 10        !< Range of multigrid levels: maximum level
INTEGER, PARAMETER :: NSCARC_LEVEL_SINGLE            =  1        !< Type of multigrid levels: only one level needed
INTEGER, PARAMETER :: NSCARC_LEVEL_MULTI             =  2        !< Type of multigrid levels: multiple levels needed

INTEGER, PARAMETER :: NSCARC_MATRIX_BANDWISE         =  1        !< Type of matrix storage technique: bandwise
INTEGER, PARAMETER :: NSCARC_MATRIX_COMPACT          =  2        !< Type of matrix storage technique: compact
INTEGER, PARAMETER :: NSCARC_MATRIX_CONDENSED        =  3        !< Flag for matrix treatment: condensing applied
INTEGER, PARAMETER :: NSCARC_MATRIX_CONNECTION       =  4        !< Flag for matrix selection: strength of connection matrix
INTEGER, PARAMETER :: NSCARC_MATRIX_FULL             =  5        !< Flag for matrix allocation: maximum possible size
INTEGER, PARAMETER :: NSCARC_MATRIX_LIGHT            =  7        !< Flag for matrix allocation: reduced size (no COLG info)
INTEGER, PARAMETER :: NSCARC_MATRIX_MINIMAL          =  8        !< Flag for matrix allocation: reduced size (only ROW and COL)
INTEGER, PARAMETER :: NSCARC_MATRIX_POISSON          =  9        !< Flag for matrix selection: Poisson 
INTEGER, PARAMETER :: NSCARC_MATRIX_POISSON_PROL     = 10        !< Flag for matrix selection: Poisson times Prolongation 
INTEGER, PARAMETER :: NSCARC_MATRIX_POISSON_SYM      = 11        !< Flag for matrix selection: symmetric Poisson 
INTEGER, PARAMETER :: NSCARC_MATRIX_LAPLACE          = 12        !< Flag for matrix selection: Poisson 
INTEGER, PARAMETER :: NSCARC_MATRIX_LAPLACE_SYM      = 13        !< Flag for matrix selection: Poisson 
INTEGER, PARAMETER :: NSCARC_MATRIX_PROLONGATION     = 14        !< Flag for matrix selection: Prolongation (AMG only)
INTEGER, PARAMETER :: NSCARC_MATRIX_RESTRICTION      = 15        !< Flag for matrix selection: Restriction (AMG only)
INTEGER, PARAMETER :: NSCARC_MATRIX_ZONES            = 16        !< Flag for matrix selection: Aggregation zones (AMG only)
INTEGER, PARAMETER :: NSCARC_MATRIX_LOWER            = 17        !< Flag for matrix selection: Lower matrix for MGM-LU decomposition
INTEGER, PARAMETER :: NSCARC_MATRIX_UPPER            = 18        !< Flag for matrix selection: Upper matrix for MGM-LU decomposition

INTEGER, PARAMETER :: NSCARC_MATVEC_GLOBAL           =  1        !< Scope of matrix-vector product: globally 
INTEGER, PARAMETER :: NSCARC_MATVEC_LOCAL            =  2        !< Scope of matrix-vector product: locally

INTEGER, PARAMETER :: NSCARC_MAX_FACE_NEIGHBORS      = 10        !< Maximum settings: Number of administrable face neighbors
INTEGER, PARAMETER :: NSCARC_MAX_STENCIL             =  7        !< Maximum settings: Number of legs in Poisson stencil
INTEGER, PARAMETER :: NSCARC_MAX_BUFFER0             = 10        !< Maximum settings: Buffer size for initial data exchanges

INTEGER, PARAMETER :: NSCARC_STORAGE_CREATE          =  1        !< Type of memory operation: create array
INTEGER, PARAMETER :: NSCARC_STORAGE_REMOVE          =  2        !< Type of memory operation: remove array
INTEGER, PARAMETER :: NSCARC_STORAGE_RESIZE          =  3        !< Type of memory operation: resize array
INTEGER, PARAMETER :: NSCARC_STORAGE_MAX             = 10000     !< Current maximum of allocatable arrays (may be increased)
                   
INTEGER, PARAMETER :: NSCARC_METHOD_KRYLOV           =  1        !< Global ScaRC method: Krylov solver
INTEGER, PARAMETER :: NSCARC_METHOD_MULTIGRID        =  2        !< Global ScaRC method: Multigrid solver
INTEGER, PARAMETER :: NSCARC_METHOD_LU               =  3        !< Global ScaRC method: LU-decomposition based on MKL
INTEGER, PARAMETER :: NSCARC_METHOD_MGM              =  4        !< Global ScaRC method: McKeeney-Greengard-Mayo solver

INTEGER, PARAMETER :: NSCARC_MKL_NONE                =  0        !< Type of MKL method: no use of MKL 
INTEGER, PARAMETER :: NSCARC_MKL_LOCAL               =  1        !< Type of MKL method: local LU-decompositions 
INTEGER, PARAMETER :: NSCARC_MKL_GLOBAL              =  2        !< Type of MKL method: global LU-decomposition
INTEGER, PARAMETER :: NSCARC_MKL_COARSE              =  3        !< Type of MKL method: only coarse grid level

INTEGER, PARAMETER :: NSCARC_MGM_POISSON             =  1        !< Type of MGM pass: First (inhomogeneous Poisson)
INTEGER, PARAMETER :: NSCARC_MGM_LAPLACE             =  2        !< Type of MGM pass: Second (homogeneous Laplace)
INTEGER, PARAMETER :: NSCARC_MGM_SCARC               =  3        !< Type of MGM pass: Process structured ScaRC solution
INTEGER, PARAMETER :: NSCARC_MGM_USCARC              =  4        !< Type of MGM pass: Process unstructured UScaRC solution
INTEGER, PARAMETER :: NSCARC_MGM_MERGE               =  5        !< Type of MGM pass: Merge first and second pass
INTEGER, PARAMETER :: NSCARC_MGM_TERMINATE           =  6        !< Type of MGM pass: Terminate current MGM  
INTEGER, PARAMETER :: NSCARC_MGM_FAILURE             = 11        !< Type of MGM convergence: Failure
INTEGER, PARAMETER :: NSCARC_MGM_SUCCESS             = 12        !< Type of MGM convergence: Success
INTEGER, PARAMETER :: NSCARC_MGM_BC_EXPOL            = 21        !< Type of internal MGM boundary: Linear extrapolatioln
INTEGER, PARAMETER :: NSCARC_MGM_BC_MEAN             = 22        !< Type of internal MGM boundary: Simple mean value 
INTEGER, PARAMETER :: NSCARC_MGM_BC_TAYLOR           = 23        !< Type of internal MGM boundary: Taylor expansion
INTEGER, PARAMETER :: NSCARC_MGM_BC_TRUE             = 24        !< Type of internal MGM boundary: True approximate
INTEGER, PARAMETER :: NSCARC_MGM_DSCARC_TO_UHL       = 31        !< Type of MGM copy: diff ScaRC-UScaRC to unstruct hom Laplace
INTEGER, PARAMETER :: NSCARC_MGM_DSCARC_TO_UHL2      = 32        !< Type of MGM copy: diff ScaRC-UScaRC to prev unstruct hom Laplace
INTEGER, PARAMETER :: NSCARC_MGM_SCARC_TO_SIP        = 33        !< Type of MGM copy: ScaRC to struct inhom Poisson
INTEGER, PARAMETER :: NSCARC_MGM_USCARC_TO_UIP       = 34        !< Type of MGM copy: UScarc to unstructured inhom Poisson
INTEGER, PARAMETER :: NSCARC_MGM_SIP_TO_UIP          = 35        !< Type of MGM copy: struct Poisson to unstruct Poisson (setup)
INTEGER, PARAMETER :: NSCARC_MGM_UHL_TO_UHL2         = 36        !< Type of MGM copy: unstruct hom Laplace to prev 
INTEGER, PARAMETER :: NSCARC_MGM_OUHL_TO_OUHL2       = 37        !< Type of MGM copy: other unstruct hom Laplace to prev
INTEGER, PARAMETER :: NSCARC_MGM_USCARC_VS_SCARC     = 41        !< Type of MGM difference: UScaRC vs ScaRC
INTEGER, PARAMETER :: NSCARC_MGM_UHL_VS_DSCARC       = 42        !< Type of MGM difference: UScaRC-ScaRC vs unstruct hom Laplace
INTEGER, PARAMETER :: NSCARC_MGM_UIP_VS_USCARC       = 43        !< Type of MGM difference: UScaRC vs unstruct inhom Poisson
INTEGER, PARAMETER :: NSCARC_MGM_LAPLACE_AMG         = 51        !< Type of MGM Laplace solver: AMG method
INTEGER, PARAMETER :: NSCARC_MGM_LAPLACE_CG          = 52        !< Type of MGM Laplace solver: CG method
INTEGER, PARAMETER :: NSCARC_MGM_LAPLACE_PARDISO     = 53        !< Type of MGM Laplace solver: IntelMKL-Pardiso solver
INTEGER, PARAMETER :: NSCARC_MGM_LAPLACE_FFT         = 54        !< Type of MGM Laplace solver: Crayfishpak FFT solver
INTEGER, PARAMETER :: NSCARC_MGM_LAPLACE_LU          = 55        !< Type of MGM Laplace solver: Own LU method
INTEGER, PARAMETER :: NSCARC_MGM_LAPLACE_LUPERM      = 56        !< Type of MGM Laplace solver: Own permuted LU method
INTEGER, PARAMETER :: NSCARC_MGM_LAPLACE_OPTIMIZED   = 57        !< Type of MGM Laplace solver: Best optimized method (FFT/PARDISO)
INTEGER, PARAMETER :: NSCARC_MGM_INTERPOL_LINEAR     = 61        !< Type of MGM interpolation: Use linear interpolation for BC's
INTEGER, PARAMETER :: NSCARC_MGM_INTERPOL_SQUARE     = 62        !< Type of MGM interpolation: Use quadratic interpolation for BC's

INTEGER, PARAMETER :: NSCARC_MULTIGRID_GEOMETRIC     =  1        !< Type of multigrid method: geometric multigrid
INTEGER, PARAMETER :: NSCARC_MULTIGRID_ALGEBRAIC     =  2        !< Type of multigrid method: algebraic multigrid
INTEGER, PARAMETER :: NSCARC_MULTIGRID_MAIN          =  1        !< Type of multigrid method: used as main solver
INTEGER, PARAMETER :: NSCARC_MULTIGRID_PRECON        =  2        !< Type of multigrid method: used as preconditioner

INTEGER, PARAMETER :: NSCARC_ORDER_ACTIVE            =  1        !< Order of aggregation: mesh is active
INTEGER, PARAMETER :: NSCARC_ORDER_LOCKED            = -1        !< Order of aggregation: mesh is locked
INTEGER, PARAMETER :: NSCARC_ORDER_UNASSIGNED        =  0        !< Order of aggregation: mesh is unassigned 

INTEGER, PARAMETER :: NSCARC_PRECISION_SINGLE        =  1        !< Type of data precision: single
INTEGER, PARAMETER :: NSCARC_PRECISION_DOUBLE        =  2        !< Type of data precision: double 

INTEGER, PARAMETER :: NSCARC_RELAX_FFT               =  1        !< Type of preconditioner: FFT-methods
INTEGER, PARAMETER :: NSCARC_RELAX_FFTO              =  2        !< Type of preconditioner: FFTO-methods (including overlap)
INTEGER, PARAMETER :: NSCARC_RELAX_GMG               =  3        !< Type of preconditioner: multigrid methods
INTEGER, PARAMETER :: NSCARC_RELAX_ILU               =  4        !< Type of preconditioner: ILU-decompositions (own)
INTEGER, PARAMETER :: NSCARC_RELAX_JAC               =  5        !< Type of preconditioner: JACOBI-methods
INTEGER, PARAMETER :: NSCARC_RELAX_LU                =  6        !< Type of preconditioner: LU-decompositions (own)
INTEGER, PARAMETER :: NSCARC_RELAX_MGS               =  7        !< Type of preconditioner: MGS-methods (matrix form)
INTEGER, PARAMETER :: NSCARC_RELAX_MJAC              =  8        !< Type of preconditioner: MJAC-methods (matrix form)
INTEGER, PARAMETER :: NSCARC_RELAX_MKL               =  9        !< Type of preconditioner: LU-decompositions (MKL)
INTEGER, PARAMETER :: NSCARC_RELAX_MSGS              = 10        !< Type of preconditioner: MSGS-methods (matrix form)
INTEGER, PARAMETER :: NSCARC_RELAX_MSOR              = 11        !< Type of preconditioner: MSOR-methods (matrix form)
INTEGER, PARAMETER :: NSCARC_RELAX_MSSOR             = 12        !< Type of preconditioner: MSSOR-methods (matrix form)
INTEGER, PARAMETER :: NSCARC_RELAX_OPTIMIZED         = 13        !< Type of preconditioner: Optimized mixture of MKL and FFT
INTEGER, PARAMETER :: NSCARC_RELAX_SSOR              = 14        !< Type of preconditioner: SSOR-methods

INTEGER, PARAMETER :: NSCARC_SCOPE_GLOBAL            =  0        !< Scope of current solver: global
INTEGER, PARAMETER :: NSCARC_SCOPE_LOCAL             =  1        !< Scope of current solver: local

INTEGER, PARAMETER :: NSCARC_SOLVER_MAIN             =  1        !< Type of solver: main solver
INTEGER, PARAMETER :: NSCARC_SOLVER_PRECON           =  2        !< Type of solver: preconditioner
INTEGER, PARAMETER :: NSCARC_SOLVER_SMOOTH           =  3        !< Type of solver: smoother
INTEGER, PARAMETER :: NSCARC_SOLVER_COARSE           =  4        !< Type of solver: coarse grid solver
INTEGER, PARAMETER :: NSCARC_SOLVER_MGM              =  5        !< Type of solver: second pass in MGM solver

INTEGER, PARAMETER :: NSCARC_STACK_ZERO              =   0       !< Order in solver stack: zero position
INTEGER, PARAMETER :: NSCARC_STACK_ROOT              =   1       !< Order in solver stack: root position
INTEGER, PARAMETER :: NSCARC_STACK_MAX               =  10       !< Order in solver stack: maximum position
INTEGER, PARAMETER :: NSCARC_STACK_NOPARENT          = -99       !< Order in solver stack: no parent available

INTEGER, PARAMETER :: NSCARC_STAGE_ONE               =  1        !< Stage of administration for current method: primary stage 
INTEGER, PARAMETER :: NSCARC_STAGE_TWO               =  2        !< Stage of administration for current method: secondary stage 

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED           =  0        !< State of multigrid: proceed loop
INTEGER, PARAMETER :: NSCARC_STATE_CONV_INITIAL      =  1        !< State of multigrid: check initial residual
INTEGER, PARAMETER :: NSCARC_STATE_CONV              =  2        !< State of multigrid: check residual
INTEGER, PARAMETER :: NSCARC_STATE_DIVG              =  3        !< State of multigrid: check divergence

INTEGER, PARAMETER :: NSCARC_STENCIL_CONSTANT        =  1        !< Type of matrix stencil: constant matrix entries
INTEGER, PARAMETER :: NSCARC_STENCIL_VARIABLE        =  2        !< Type of matrix stencil: variable matrix entries

INTEGER, PARAMETER :: NSCARC_TWOLEVEL_NONE           =  0        !< Type of two-level method: only one level
INTEGER, PARAMETER :: NSCARC_TWOLEVEL_COARSE_ADD     =  1        !< Type of two-level method: additive 2-level 
INTEGER, PARAMETER :: NSCARC_TWOLEVEL_COARSE_MUL1    =  2        !< Type of two-level method: multiplicative1 2-level (fine last)
INTEGER, PARAMETER :: NSCARC_TWOLEVEL_COARSE_MUL2    =  3        !< Type of two-level method: multiplicative2 2-level (fine first)
INTEGER, PARAMETER :: NSCARC_TWOLEVEL_COARSE_ONLY    =  4        !< Type of two-level method: only coarse grid
INTEGER, PARAMETER :: NSCARC_TWOLEVEL_MACRO          =  5        !< Type of two-level method: use macro solver 
INTEGER, PARAMETER :: NSCARC_TWOLEVEL_XMEAN_ADD      =  6        !< Type of two-level method: additive xmean 
INTEGER, PARAMETER :: NSCARC_TWOLEVEL_XMEAN_MUL1     =  7        !< Type of two-level method: multiplicative1 xmean (fine last)
INTEGER, PARAMETER :: NSCARC_TWOLEVEL_XMEAN_MUL2     =  8        !< Type of two-level method: multiplicative2 xmean (fine first)

INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_X            =  1        !< Flag for 1D-vector on stage 1: X
INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_B            =  2        !< Flag for 1D-vector on stage 1: B
INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_D            =  3        !< Flag for 1D-vector on stage 1: D
INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_R            =  4        !< Flag for 1D-vector on stage 1: R
INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_V            =  5        !< Flag for 1D-vector on stage 1: V
INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_Y            =  6        !< Flag for 1D-vector on stage 1: Y
INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_Z            =  7        !< Flag for 1D-vector on stage 1: Z
INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_E            =  8        !< Flag for 1D-vector on stage 1: E
INTEGER, PARAMETER :: NSCARC_VECTOR_TWO_X            =  9        !< Flag for 1D-vector on stage 2: X
INTEGER, PARAMETER :: NSCARC_VECTOR_TWO_B            = 10        !< Flag for 1D-vector on stage 2: B
INTEGER, PARAMETER :: NSCARC_VECTOR_TWO_D            = 11        !< Flag for 1D-vector on stage 2: D                
INTEGER, PARAMETER :: NSCARC_VECTOR_TWO_R            = 12        !< Flag for 1D-vector on stage 2: R
INTEGER, PARAMETER :: NSCARC_VECTOR_TWO_V            = 13        !< Flag for 1D-vector on stage 2: V
INTEGER, PARAMETER :: NSCARC_VECTOR_TWO_Y            = 14        !< Flag for 1D-vector on stage 2: Y                
INTEGER, PARAMETER :: NSCARC_VECTOR_TWO_Z            = 15        !< Flag for 1D-vector on stage 2: Z
INTEGER, PARAMETER :: NSCARC_VECTOR_TWO_E            = 16        !< Flag for 1D-vector on stage 2: E
INTEGER, PARAMETER :: NSCARC_VECTOR_H                = 17        !< Flag for 3D-vector H 
INTEGER, PARAMETER :: NSCARC_VECTOR_HS               = 18        !< Flag for 3D-vector HS

INTEGER, PARAMETER :: NSCARC_WARNING_NO_MKL_PRECON   =  1        !< Type of warning message: No MKL preconditioner available
INTEGER, PARAMETER :: NSCARC_WARNING_NO_MKL_SMOOTH   =  2        !< Type of warning message: No MKL smoother available
INTEGER, PARAMETER :: NSCARC_WARNING_NO_MKL_LU       =  3        !< Type of warning message: Using LU instead of MKL
INTEGER, PARAMETER :: NSCARC_WARNING_NO_GLOBAL_SCOPE =  4        !< Type of warning message: No global scope solver available
INTEGER, PARAMETER :: NSCARC_WARNING_ONLY_SSOR_PRECON=  5        !< Type of warning message: Only SSOR preconditioner available
INTEGER, PARAMETER :: NSCARC_WARNING_ONLY_SSOR_SMOOTH=  6        !< Type of warning message: Only SSOR smoother available

INTEGER, PARAMETER :: NSCARC_UNDEF_INT               = -1        !< Flag for undefined integer value
INTEGER, PARAMETER :: NSCARC_ZERO_INT                =  0        !< Flag for zero integer value
INTEGER, PARAMETER :: NSCARC_ONE_INT                 =  1        !< Flag for one integer value

REAL(EB), PARAMETER:: NSCARC_UNDEF_REAL_EB           = -1.0_EB   !< Flag for undefined double precision value
REAL(EB), PARAMETER:: NSCARC_ZERO_REAL_EB            =  0.0_EB   !< Flag for zero double precision value
REAL(EB), PARAMETER:: NSCARC_ONE_REAL_EB             =  1.0_EB   !< Flag for one double precision value

REAL(FB), PARAMETER:: NSCARC_UNDEF_REAL_FB           = -1.0_FB   !< Flag for undefined single precision value
REAL(FB), PARAMETER:: NSCARC_ZERO_REAL_FB            =  0.0_FB   !< Flag for zero single precision value
REAL(FB), PARAMETER:: NSCARC_ONE_REAL_FB             =  1.0_FB   !< Flag for one single precision value

CHARACTER(40), PARAMETER :: SCARC_NONE = 'NONE'                  !< Flag for a dummy character value 
INTEGER, PARAMETER  :: NSCARC_NONE = -123456789                  !< Flag for a dummy integer value 

REAL(FB), PARAMETER :: NSCARC_HUGE_REAL_FB = -999999999.0_FB     !< Flag for an undefined double precision value
REAL(EB), PARAMETER :: NSCARC_HUGE_REAL_EB = -999999999.0_EB     !< Flag for an undefined double precision value
INTEGER, PARAMETER  :: NSCARC_HUGE_INT     = -999999999          !< Flag for an undefined integer value

REAL(EB), PARAMETER :: NSCARC_THRESHOLD_CONVERGENCE = 1.0E-12_EB !< Threshold for convergence
REAL(EB), PARAMETER :: NSCARC_THRESHOLD_DIVGERGENCE = 1.0E+15_EB !< Threshold for divergence

REAL(EB), PARAMETER :: SCALR  = 0.015625_EB                      !< Scaling parameter for geometric multigrid method
REAL(EB), PARAMETER :: SCALP  = 0.0625_EB                        !< Scaling parameter for geometric multigrid method
REAL(EB), PARAMETER :: W1     =  1.0_EB                          !< Weighting parameter for geometric multigrid method
REAL(EB), PARAMETER :: W3     =  3.0_EB                          !< Weighting parameter for geometric multigrid method
REAL(EB), PARAMETER :: W4     =  4.0_EB                          !< Weighting parameter for geometric multigrid method
REAL(EB), PARAMETER :: W9     =  9.0_EB                          !< Weighting parameter for geometric multigrid method
REAL(EB), PARAMETER :: W12    = 12.0_EB                          !< Weighting parameter for geometric multigrid method
REAL(EB), PARAMETER :: W16    = 16.0_EB                          !< Weighting parameter for geometric multigrid method

END MODULE SCARC_CONSTANTS


!=======================================================================================================================
!
! MODULE SCARC_TYPES
!
!> \brief Collection of data types used for the different ScaRC/UScaRC solvers
!
!=======================================================================================================================
MODULE SCARC_TYPES

USE PRECISION_PARAMETERS
USE SCARC_CONSTANTS

#ifdef WITH_MKL
USE MKL_PARDISO
USE MKL_CLUSTER_SPARSE_SOLVER
#endif

IMPLICIT NONE (TYPE,EXTERNAL)

!> \brief Detailed information about arrays created within the ScaRC memory manager

TYPE SCARC_ALLOCATION_TYPE

   INTEGER :: NTYPE  = NSCARC_INIT_NONE               !< Data type of array
   INTEGER :: NDIM   = NSCARC_INIT_NONE               !< Dimension of array
   INTEGER :: NINIT  = NSCARC_INIT_NONE               !< Initialization type of array
   INTEGER :: NRANK  = NSCARC_INIT_NONE               !< Rank of array (order of allocation)
   INTEGER :: NSTATE = NSCARC_INIT_NONE               !< State of array (allocated/removed)

   INTEGER :: LBND(3) = NSCARC_INIT_NONE              !< Left bounds of array for x-,y- and z-direction
   INTEGER :: RBND(3) = NSCARC_INIT_NONE              !< Right bounds of array for x-,y- and z-direction

   CHARACTER(60) :: CNAME                             !< Name of array
   CHARACTER(60) :: CSCOPE                            !< Name of allocating routine 

END TYPE SCARC_ALLOCATION_TYPE

!> \brief ScaRC memory manager type
 
TYPE SCARC_STORAGE_TYPE

   TYPE (SCARC_ALLOCATION_TYPE), ALLOCATABLE, DIMENSION(:) :: ALLOCATION_LIST  !< Administrative list of allocated structures

   INTEGER :: IP                                           !< Pointer to current array entry
   INTEGER :: IRANK                                        !< Rank of memory accesses

   INTEGER :: NSUM_INT     = NSCARC_INIT_ZERO              !< Total sum of all allocated integer data
   INTEGER :: NSUM_REAL_EB = NSCARC_INIT_ZERO              !< Total sum of all allocated double precision data
   INTEGER :: NSUM_REAL_FB = NSCARC_INIT_ZERO              !< Total sum of all allocated single precision data

   INTEGER :: N_ARRAYS  = NSCARC_INIT_ZERO                 !< Number of allocated arrays
   INTEGER :: N_BMATRIX = NSCARC_INIT_ZERO                 !< Number of allocated bandwise stored matrices
   INTEGER :: N_CMATRIX = NSCARC_INIT_ZERO                 !< Number of allocated compactly stored matrices
   INTEGER :: N_INT     = NSCARC_INIT_ZERO                 !< Number of allocated integer data
   INTEGER :: N_LOG     = NSCARC_INIT_ZERO                 !< Number of allocated logical data
   INTEGER :: N_REAL_EB = NSCARC_INIT_ZERO                 !< Number of allocated double precision data
   INTEGER :: N_REAL_FB = NSCARC_INIT_ZERO                 !< Number of allocated single precision data

   INTEGER :: NWORK_INT     = NSCARC_INIT_ZERO             !< Workspace occupied by integer arrays
   INTEGER :: NWORK_LOG     = NSCARC_INIT_ZERO             !< Workspace occupied by logical arrays
   INTEGER :: NWORK_REAL_EB = NSCARC_INIT_ZERO             !< Workspace occupied by double precision arrays
   INTEGER :: NWORK_REAL_FB = NSCARC_INIT_ZERO             !< Workspace occupied by single precision arrays

END TYPE SCARC_STORAGE_TYPE

!> \brief Messaging and debugging mechanisms
 
TYPE SCARC_MESSAGE_TYPE

   CHARACTER(60) :: FILE_CPU                               !< Output file name for CPU measurements
   CHARACTER(60) :: FILE_MEM                               !< Output file name for memory management information
   CHARACTER(60) :: FILE_STAT                              !< Output file name for convergence statistcis
   INTEGER :: LU_CPU                                       !< Logical unit for CPU measurements
   INTEGER :: LU_MEM                                       !< Logical unit for memory management information
   INTEGER :: LU_STAT                                      !< Logical unit for convergence statistics

END TYPE SCARC_MESSAGE_TYPE

!> \brief Information about the neighborship structur within the mesh decomposition 
 
TYPE SCARC_SUBDIVISION_TYPE

   INTEGER, ALLOCATABLE, DIMENSION (:)   :: N_NEIGHBORS    !< Number of meshes for in complete subdivision
   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: NEIGHBORS      !< Global neighborship structure between meshes
   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ORDER          !< Order of meshes during aggregation process
   INTEGER :: N_CYCLES                                     !< Number of needed cycles during aggregation process
   INTEGER :: N_NEIGHBORS_TOTAL                            !< Sum of all neighbors

END TYPE SCARC_SUBDIVISION_TYPE

!> \brief Face information related to wall cells and neighbors
 
TYPE SCARC_FACE_TYPE

   REAL(EB), POINTER, DIMENSION(:) :: DH              !< Step size vector between adjacent faces

   REAL(EB) :: INCR_BOUNDARY                          !< Increment for boundary conditions
   REAL(EB) :: INCR_FACE                              !< Increment for matrix subdiagonals on face
   REAL(EB) :: INCR_INSIDE                            !< Increment for matrix subdiagonals between opposite faces

   REAL(EB) :: SCAL_DIRICHLET                         !< Scaling factor for Dirichlet BC's
   REAL(EB) :: SCAL_NEUMANN                           !< Scaling factor for Neumann BC's
   REAL(EB) :: SCAL_MGM                               !< Scaling factor for MGM method

   INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS    !< Adjacent neighbors at that face
   INTEGER  :: N_NEIGHBORS = 0                        !< Number of adjacent neighbors 

   INTEGER  :: NOP = 0                                !< Number of cells between opposite faces
   INTEGER  :: NX = 0, NY = 0, NZ = 0                 !< Cells in different directions on that face
   INTEGER  :: NCW0 = 0, NCW = 0                      !< Number of first wall cell and total number of wall cells
   INTEGER  :: INCRX = 0, INCRY = 0, INCRZ = 0        !< Increments to next internal cell in that face direction

   INTEGER  :: INCR_STENCIL (-3:3)                    !< Increments within stencil for HB-vector on that face

END TYPE SCARC_FACE_TYPE
 
!> \brief Wall information related to neighbors and BC's
 
TYPE SCARC_WALL_TYPE

   ! different properties of wall cell
   INTEGER :: BTYPE = 0                               !< Type of wall cell (Dirichlet/Neumann/Internal)
   INTEGER :: BOUNDARY_TYPE = 0                       !< Type of boundary for wall cell (Solid/Interpolated/Open))
   INTEGER :: IOR = 0                                 !< Orientation of wall cell
   INTEGER :: NOM = 0                                 !< Adjacent neighbor at wall cell

   INTEGER :: ICW = NSCARC_UNDEF_INT                  !< Internal wall cell for IW
   INTEGER :: IXG, IYG, IZG                           !< Coordinate indices of ghost cells
   INTEGER :: IXW, IYW, IZW                           !< Coordinate indices of (internal) wall cells
   INTEGER :: IXN(2), IYN(2), IZN(2)                  !< Coordinate indices of neighboring cells
   INTEGER :: ICE, ICG                                !< Flag for externa and ghost cell

END TYPE SCARC_WALL_TYPE

!> \brief Obstruction information
 
TYPE SCARC_OBST_TYPE
   INTEGER :: I1, I2, J1, J2, K1, K2                  !< Cell indices of obstructions
END TYPE SCARC_OBST_TYPE

!> \brief Preconditioning information for unidirection meanvalue preconditioning
 
TYPE SCARC_PRECON_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:) :: DIAG        !< Main diagonal workspace for Twolevel-Xmean preconditioning
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: AUX         !< Auxiliary main diagonal workspace for Twolevel-Xmean preconditioning
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: LOW         !< Lower diagonal workspace for Twolevel-Xmean preconditioning
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: UP          !< Upper workspace for Twolevel-Xmean preconditioning
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: RHS         !< RHS workspace for Twolevel-Xmean preconditioning

END TYPE SCARC_PRECON_TYPE

!> \brief Compact matrix entries which will be exchanged during generation of condensed system
 
TYPE SCARC_MATRIX_COMPACT_CONDENSED_TYPE

   REAL(EB) :: VAL1(NSCARC_MAX_STENCIL) = 0.0_EB      !< Original values (double precision)
   REAL(EB) :: VAL2(NSCARC_MAX_STENCIL) = 0.0_EB      !< Condensed values (double precision)

   INTEGER :: COL(NSCARC_MAX_STENCIL) = 0             !< Column pointers
   INTEGER :: PTR(NSCARC_MAX_STENCIL) = 0             !< Storage pointer
   INTEGER :: N_COL

END TYPE SCARC_MATRIX_COMPACT_CONDENSED_TYPE

!> \brief Bandwise matrix entries which will exchanged during generation of condensed system
 
TYPE SCARC_MATRIX_BANDWISE_CONDENSED_TYPE

   REAL(EB) :: VAL1(NSCARC_MAX_STENCIL) = 0.0_EB      !< Original values (double precision)
   REAL(EB) :: VAL2(NSCARC_MAX_STENCIL) = 0.0_EB      !< Condensed values (double precision)
   INTEGER  :: IOR0 = 0                               !< Position pointer
   INTEGER  :: ICO = 0                                !< Cell pointer

END TYPE SCARC_MATRIX_BANDWISE_CONDENSED_TYPE

!> \brief Compact sparse row (COMPACT) storage technique for matrices
! Is based on three arrays:
!    - non-zero matrix values
!    - corresponding columns pointers
!    - row pointers
 
TYPE SCARC_CMATRIX_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:) :: VAL         !< Values of matrix (real precision)
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: ILU         !< ILU-decomposition
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: DIAG         !< ILU-decomposition
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: RELAX       !< Workspace for relaxation
   REAL(EB), DIMENSION (-3:3)           :: STENCIL     !< Store basic stencil information in single precision

   REAL(FB), ALLOCATABLE, DIMENSION (:) :: VAL_FB      !< Values of matrix (single precision)
   REAL(FB), ALLOCATABLE, DIMENSION (:) :: RELAX_FB    !< Workspace for relaxation
   REAL(FB), DIMENSION (-3:3)           :: STENCIL_FB  !< Store basic stencil information in single precision

   TYPE (SCARC_MATRIX_COMPACT_CONDENSED_TYPE) :: CONDENSED(NSCARC_MAX_STENCIL)

   INTEGER,  ALLOCATABLE, DIMENSION (:) :: ROW         !< Row index vector
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: COL         !< Local column index vectors
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: COLG        !< Global column index vectors

   INTEGER :: POS(-3:3) = 0                            !< Position of IOR's in STENCIL
   INTEGER :: N_CONDENSED = 0                          !< Number of condensed entries in matrix
   INTEGER :: N_VAL = 0                                !< Number of matrix values
   INTEGER :: N_ROW = 0                                !< Number of matrix rows
   INTEGER :: N_STENCIL = 0                            !< Number of points in matrix stencil
   INTEGER :: N_STENCIL_MAX = 0                        !< Max stencil size (AMG only)
   INTEGER :: NTYPE = 0                                !< Matrix type
   INTEGER :: NPREC = 0                                !< Precision type

   CHARACTER(40) :: CNAME                              !< Name of matrix

END TYPE SCARC_CMATRIX_TYPE

!> \brief Bandwise storage technique for matrices
! The entries are stored one diagonal after the other
! Missing entries of subdiagonals are filled with zero
! Is based on two arrays: non-zero matrix entries diagonal-wise and  the offsets from the main diagonal
 
TYPE SCARC_BMATRIX_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: AUX           !< Auxiliary vector (double precision)
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: VAL           !< Values of matrix (double precision)
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: RELAX         !< Workspace for relaxation
   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: RELAXD        !< Workspace for relaxation - only for diagonal scaling
   REAL(EB), DIMENSION (-3:3)             :: STENCIL       !< Store basic stencil information (double precision)

   REAL(FB), ALLOCATABLE, DIMENSION (:)   :: AUX_FB        !< Auxiliary vector (single precision)
   REAL(FB), ALLOCATABLE, DIMENSION (:,:) :: VAL_FB        !< Values of matrix (single precision)
   REAL(FB), ALLOCATABLE, DIMENSION (:,:) :: RELAX_FB      !< Workspace for relaxation
   REAL(FB), DIMENSION (-3:3)             :: STENCIL_FB    !< Store basic stencil information (single precision)

   TYPE (SCARC_MATRIX_BANDWISE_CONDENSED_TYPE) :: CONDENSED(NSCARC_MAX_STENCIL)

   CHARACTER(40) :: CNAME                                  !< Name of matrix

   INTEGER,  DIMENSION (-3:3) :: OFFSET                    !< Offset pointers
   INTEGER,  DIMENSION (-3:3) :: LENGTH                    !< Relevant diagonal length 
   INTEGER,  DIMENSION (-3:3) :: SOURCE                    !< Source address in corresponding diagonal
   INTEGER,  DIMENSION (-3:3) :: TARGET                    !< Target address in corresponding diagonal

   INTEGER :: POS(-3:3) = 0                                !< position of IOR's in STENCIL and in matrix storage array
   INTEGER :: N_STENCIL = 0                                !< Number of points in matrix stencil
   INTEGER :: N_CONDENSED = 0                              !< Number of condensed entries in matrix
   INTEGER :: N_VAL = 0                                    !< Number of matrix values in general and symmetric cass
   INTEGER :: N_DIAG = 0                                   !< Length of main diagonal

END TYPE SCARC_BMATRIX_TYPE


!> \brief Information for ScaRC-internal instances of Crayfishpak-based FFT methods, used as preconditioners and smoothers
  
TYPE SCARC_FFT_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)       :: SAVE1     !< Saving area 
   REAL(EB), ALLOCATABLE, DIMENSION (:)       :: WORK      !< Workspace 
   REAL(EB), ALLOCATABLE, DIMENSION (:)       :: HX        !< Grid stretching vector 
   REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BXS       !< Boundary conditions along XS 
   REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BXF       !< Boundary conditions along XF 
   REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BYS       !< Boundary conditions along YS 
   REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BYF       !< Boundary conditions along YF 
   REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BZS       !< Boundary conditions along ZS 
   REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BZF       !< Boundary conditions along ZF 
   REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: PRHS      !< Current right hand side

   REAL(EB) :: XS                                          !< Minimum x-coordinate of mesh (may include overlap)
   REAL(EB) :: XF                                          !< Maximum x-coordinate of mesh (may include overlap)
   REAL(EB) :: YS                                          !< Minimum y-coordinate of mesh (may include overlap)
   REAL(EB) :: YF                                          !< Maximum y-coordinate of mesh (may include overlap)
   REAL(EB) :: ZS                                          !< Minimum z-coordinate of mesh (may include overlap)
   REAL(EB) :: ZF                                          !< Maximum z-coordinate of mesh (may include overlap)

   REAL(EB) :: POIS_PTB = 0.0_EB                           !< Perturbation parameter
   REAL(EB) :: XLM = 0.0_EB                                !< No Helmholtz equation used

   INTEGER :: LSAVE                                        !< Length of saving area 
   INTEGER :: LWORK                                        !< Length of workspace areas

   INTEGER :: LBC                                          !< Boundary type in x-direction
   INTEGER :: MBC                                          !< Boundary type in y-direction
   INTEGER :: NBC                                          !< Boundary type in z-direction

   INTEGER :: ITRN                                         !< Number of nodes in x-direction
   INTEGER :: JTRN                                         !< Number of nodes in y-direction
   INTEGER :: KTRN                                         !< Number of nodes in z-direction

   INTEGER :: IBAR                                         !< Number of cells in x-direction
   INTEGER :: JBAR                                         !< Number of cells in y-direction
   INTEGER :: KBAR                                         !< Number of cells in z-direction

END TYPE SCARC_FFT_TYPE

#ifdef WITH_MKL
!> \brief MKL information needed for IntelMKL PARDISO and CLUSTER_SPARSE_SOLVER solvers
  
TYPE SCARC_MKL_TYPE

   CHARACTER(40) :: CNAME                                                     !< Name of matrix

   INTEGER, ALLOCATABLE :: IPARM(:)                                           !< Parameter vector
   INTEGER :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR, MSGLVL                 !< Various MKL specific settings 
   INTEGER :: PERM(1)                                                         !< Permutation parameter 

   TYPE (MKL_PARDISO_HANDLE),               ALLOCATABLE :: PT_H(:), PT(:)     !< Handles for PARDISO 
   TYPE (MKL_CLUSTER_SPARSE_SOLVER_HANDLE), ALLOCATABLE :: CT_H(:), CT(:)     !< Handles for CLUSTER_SPARSE_SOLVER 

END TYPE SCARC_MKL_TYPE
#endif

!> \brief Different scopes for solution, rhs and auxiliary vectors of different solvers
  
TYPE SCARC_STAGE_TYPE

   REAL (EB), ALLOCATABLE, DIMENSION (:) :: X               !< Solution vector in double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: B               !< Right hand side vector in double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: R               !< Residual vector in double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: D               !< Auxiliary vector in double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: V               !< Auxiliary vector in double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: Y               !< Auxiliary vector in double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: Z               !< Auxiliary vector in double precision

   REAL (FB), ALLOCATABLE, DIMENSION (:) :: X_FB            !< Solution vector vector in single precision
   REAL (FB), ALLOCATABLE, DIMENSION (:) :: B_FB            !< Right hand side vector in single precision
   REAL (FB), ALLOCATABLE, DIMENSION (:) :: R_FB            !< Residual vector in single precision
   REAL (FB), ALLOCATABLE, DIMENSION (:) :: V_FB            !< Auxiliary vector in single precision


END TYPE SCARC_STAGE_TYPE

!> \brief Multigrid type - to be extended for algebraic multigrid
  
TYPE SCARC_MULTIGRID_TYPE

   REAL(EB) :: OMEGA = 1.0_EB                              !< Relaxation parameter

   INTEGER :: CYCLING(2) = 0                               !< Counter for multigrid cycling
   INTEGER :: N_PRESMOOTH, N_POSTSMOOTH                    !< Number of pre- and post-processing steps

#ifdef WITH_SCARC_AMG
   REAL(EB) :: APPROX_SPECTRAL_RADIUS = 2.0_EB             !< Relaxation parameter (AMG only)
   REAL(EB) :: AMG_TOL = 0.25_EB                           !< Tolerance for coarsening (AMG only)
   REAL(EB) :: THETA = 0.10_EB                             !< Threshold for aggregation process (AMG only)
#endif
END TYPE SCARC_MULTIGRID_TYPE

!> \brief Information related to discretization type (structured/unstructured)
  
TYPE SCARC_GRID_TYPE

   ! Basic boundary and interface information
   TYPE (SCARC_WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: WALL   !< wall information

   ! Matrices in different storage types
   TYPE (SCARC_BMATRIX_TYPE) :: POISSONB                       !< Poisson matrix in bandwise storage technique
   TYPE (SCARC_CMATRIX_TYPE) :: POISSON                        !< Poisson matrix in compact storage technique (default)
   TYPE (SCARC_CMATRIX_TYPE) :: LAPLACE                        !< Laplace matrix in compact storage technique

#ifdef WITH_MKL
   TYPE (SCARC_CMATRIX_TYPE) :: POISSON_SYM                    !< Symmetric part of compact Poisson matrix 
   TYPE (SCARC_CMATRIX_TYPE) :: LAPLACE_SYM                    !< Symmetric part of compact Laplace matrix (MGM only)
   TYPE (SCARC_CMATRIX_TYPE) :: GALERKIN_SYM                   !< Galerkin matrix symmetric version
#endif

   TYPE (SCARC_CMATRIX_TYPE) :: PROLONGATION                   !< Prolongation matrix
   TYPE (SCARC_CMATRIX_TYPE) :: RESTRICTION                    !< Restriction matrix

   TYPE (SCARC_CMATRIX_TYPE) :: LOWER                          !< Lower matrix for MGM-LU decompostion
   TYPE (SCARC_CMATRIX_TYPE) :: UPPER                          !< Upper matrix for MGM-LU decompostion

   INTEGER,  ALLOCATABLE, DIMENSION (:) :: LOCAL_TO_GLOBAL     !< Mapping from local to global numbering 

   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_IWG         !< Mapping from ICE to IWG
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_ICN         !< Mapping from ICE to ICN

   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ICG_TO_ICW         !< Mapping from ICG to ICW
   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ICG_TO_ICE         !< Mapping from ICG to ICE 
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_IWG         !< Mapping from ICG to IWG

   ! Forward and backward permutation vectors for grid renumbering strategies (MGM only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: PERM_FW            !< Permutation vector for LU - forward direction
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: PERM_BW            !< Permutation vector for LU - backward direction
   INTEGER :: NONZERO = 1                                      !< Index of first nonzero entry in RHS vector for permuted LU

   ! Assignment of cell coordinates
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICX                !< I-coordinate of cell IC
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICY                !< J-coordinate of cell IC
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICZ                !< J-coordinate of cell IC

   ! Number and state of single cells
   INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: CELL_NUMBER      !< Numbering of single cells

   ! Cell numbers of all meshes and offsets between meshes
   INTEGER, ALLOCATABLE, DIMENSION (:) :: NC_LOCAL             !< Number of cells in local meshes
   INTEGER, ALLOCATABLE, DIMENSION (:) :: NC_OFFSET            !< Offset in cell numbering between meshes
   INTEGER :: NC_GLOBAL   = NSCARC_ZERO_INT                    !< Global number of cells in all meshes
   INTEGER :: NC_GALERKIN = NSCARC_ZERO_INT                    !< Number of cells with influence on Galerkin matrix

   ! Local numbers of internal, extended and ghost cells
   INTEGER :: NC   = NSCARC_ZERO_INT                           !< Number of cells needed for matrix
   INTEGER :: NW   = NSCARC_ZERO_INT                           !< Number of wall cells
   INTEGER :: NCG  = NSCARC_ZERO_INT                           !< Number of ghost cells 
   INTEGER :: NCGI = NSCARC_ZERO_INT                           !< Number of ghost zones in internal direction
   INTEGER :: NCGE = NSCARC_ZERO_INT                           !< Number of ghost zones in external direction
   INTEGER :: NZG  = NSCARC_ZERO_INT                           !< Number of ghost zones 
   INTEGER :: NCE  = NSCARC_ZERO_INT                           !< Number of extended cells
   INTEGER :: NCE2 = NSCARC_ZERO_INT                           !< Number of extended cells, second layer
   INTEGER :: ICE2 = NSCARC_ZERO_INT                           !< Counter for extended cells, second layer

   ! Number of Dirichlet and Neumann boundary cells
   INTEGER :: N_DIRIC   = NSCARC_ZERO_INT                      !< Number of Dirichlet BCs
   INTEGER :: N_NEUMANN = NSCARC_ZERO_INT                      !< Number of Neumann BCs

   ! Pointer variables and arrays for data exchange with neighbors
   INTEGER :: ICG  = NSCARC_ZERO_INT                           !< Ghost cell pointer for first layer
   INTEGER :: ICG2 = NSCARC_ZERO_INT                           !< Ghost cell pointer for second layer
   INTEGER :: ICE  = NSCARC_ZERO_INT                           !< Ghost cell pointer for extended cells

   INTEGER :: NLEN_BUFFER_LAYER1  = NSCARC_ZERO_INT            !< Length for single layer length exchange on that level
   INTEGER :: NLEN_BUFFER_LAYER2  = NSCARC_ZERO_INT            !< Length for double layer length exchange on that level
   INTEGER :: NLEN_BUFFER_LAYER4  = NSCARC_ZERO_INT            !< Length for fourfold length exchange on that level
   INTEGER :: NLEN_BUFFER_STENCIL = NSCARC_ZERO_INT            !< Length for stencil layer length exchange on that level
   INTEGER :: NLEN_BUFFER_FULL    = NSCARC_ZERO_INT            !< Length for full length exchange on that level

#ifdef WITH_SCARC_AMG
   TYPE (SCARC_CMATRIX_TYPE) :: GALERKIN                       !< Galerkin matrix (AMG only)
   TYPE (SCARC_CMATRIX_TYPE) :: CONNECTION                     !< Strength of connection matrix (AMG only)
   TYPE (SCARC_CMATRIX_TYPE) :: ZONES                          !< Aggregation Zones matrix (AMG only)
   TYPE (SCARC_CMATRIX_TYPE) :: POISSON_PROL                   !< Poisson times Prolongation matrix (AMG only)

   REAL(EB), ALLOCATABLE, DIMENSION (:)  :: MEASURES           !< Measure for grid coarsening (AMG only)
   REAL(EB), ALLOCATABLE, DIMENSION (:)  :: NULLSPACE          !< Nullspace vector (AMG only)
   REAL(EB), ALLOCATABLE, DIMENSION (:)  :: AUX1, AUX2         !< Auxiliary vectors (AMG only)
   REAL(EB), ALLOCATABLE, DIMENSION (:)  :: DIAG               !< Matrix diagonal, possible inverted (AMG only)
   REAL(EB), ALLOCATABLE, DIMENSION (:)  :: QQ, RR             !< workspace for QR-decompostion (AMG only)

   INTEGER,  ALLOCATABLE, DIMENSION (:)  :: CELLS_LOCAL        !< Local coarse cells which influence Galerkin matrix (AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:)  :: CELLS_GLOBAL       !< Global coarse cells which influence Galerkin matrix (AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:)  :: ORDER              !< Search order for aggregation (AMG only)

   INTEGER,  ALLOCATABLE, DIMENSION (:)  :: ZONES_GLOBAL       !< Global zone numbers (AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:)  :: ZONES_LOCAL        !< Local  zone numbers (AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:)  :: ZONE_CENTERS       !< Zone centers for grid coarsening (AMG only)

   REAL(EB), ALLOCATABLE, DIMENSION (:)  :: ELAYER2_VALS       !< Values of cells in external second layer(AMG only)
   REAL(EB), ALLOCATABLE, DIMENSION (:)  :: ILAYER2_VALS       !< Values of cells in internal second layer(AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:)  :: ELAYER2_NUMS       !< Number of cells in external second layer (AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:)  :: ILAYER2_NUMS       !< Number of cells in internal second layer (AMG only)

   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_ECELL       !< Mapping from ICG to global neighboring cell (AMG only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_ICELL       !< Mapping from ICG to local neighboring cell (AMG only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_GCELL       !< Mapping from ICG to global neighboring cell (AMG only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_OCELL       !< Mapping from ICG to local neighboring cell (AMG only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_IZONE       !< Mapping from ICG to internal own zone (AMG only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_EZONE       !< Mapping from ICG to external own zone (AMG only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_GZONE       !< Mapping from ICG to global neighboring zone (AMG only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_OZONE       !< Mapping from ICG to local neighboring zone (AMG only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_ELAYER2     !< Mapping from ICG to external second layer (AMG only)
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_ILAYER2     !< Mapping from ICG to internal second layer (AMG only)

   INTEGER :: N_FINE    = NSCARC_ZERO_INT                      !< Number of fine cells (AMG only)
   INTEGER :: N_COARSE  = NSCARC_ZERO_INT                      !< Number of coarse cells (AMG only)
   INTEGER :: N_ZONES   = NSCARC_ZERO_INT                      !< Number of zones (AMG only)
   INTEGER :: N_STENCIL_MAX  = 25                              !< Max stencil size (AMG only)
#endif

END TYPE SCARC_GRID_TYPE

!> \brief McKenney-Greengard-Mayo method - still experimental
  
TYPE SCARC_MGM_TYPE

   TYPE (SCARC_FFT_TYPE) :: FFT                                    !< FFT preconditioner for local Laplace problems (if structured)

   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: SIP                 !< structured inhomogeneous Poisson MGM solution 
   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: UIP                 !< unstructured Poisson MGM solution 
   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: UHL                 !< unstructured Laplace MGM solution 
   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: UHL2                !< unstructured Laplace MGM solution II

   REAL(EB), ALLOCATABLE, DIMENSION (:)     :: OUIP                !< other unstructured Poisson MGM solution 
   REAL(EB), ALLOCATABLE, DIMENSION (:)     :: OUHL                !< other unstructured Laplace MGM solution 
   REAL(EB), ALLOCATABLE, DIMENSION (:)     :: OUHL2               !< other unstructured Laplace MGM solution II

   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: SCARC               !< structured ScaRC solution 
   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: USCARC              !< unstructured ScaRC solution 
   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: DSCARC              !< difference of structured and unstructured ScaRC solutions

   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: UIP_VS_USCARC       !< difference of unstructured Poisson MGM and UScaRC solutions
   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: UHL_VS_DSCARC       !< difference of unstructured Laplace MGM and diff ScaRC-UScarC

   REAL(EB), ALLOCATABLE, DIMENSION (:)     :: BC                  !< Boundary conditions along internal mesh interfaces
   REAL(EB), ALLOCATABLE, DIMENSION (:)     :: X, Y, B             !< solution, auxiliary and RHS vectors for local Laplace solvers
   REAL(FB), ALLOCATABLE, DIMENSION (:)     :: X_FB, B_FB          !< solution and RHS vectors for local Laplace solvers (SP)
   REAL(EB), ALLOCATABLE, DIMENSION (:,:)   :: AAA                 !< Matrix for LU-decomposition (experimental)
   REAL(EB), ALLOCATABLE, DIMENSION (:,:)   :: LLL                 !< Lower part of LU-decomposition (experimental)
   REAL(EB), ALLOCATABLE, DIMENSION (:,:)   :: UUU                 !< Upper part of LU-decomposition (experimental)
   REAL(EB), ALLOCATABLE, DIMENSION (:)     :: WEIGHT              !< Scaling weights for true boundary setting

   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: UVEL                !< u-velocity vector component
   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: VVEL                !< v-velocity vector component
   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: WVEL                !< w-velocity vector component

   REAL(EB), ALLOCATABLE, DIMENSION (:)     :: OUVEL               !< other u-velocity component
   REAL(EB), ALLOCATABLE, DIMENSION (:)     :: OVVEL               !< other v-velocity component
   REAL(EB), ALLOCATABLE, DIMENSION (:)     :: OWVEL               !< other w-velocity component 

   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: U1, V1, W1          !< Velocity vectors predictor
   REAL(EB), ALLOCATABLE, DIMENSION (:,:,:) :: U2, V2, W2          !< Velocity vectors predictor

   REAL(EB), ALLOCATABLE, DIMENSION (:,:)   :: BXS                 !< Boundary values at x-min face
   REAL(EB), ALLOCATABLE, DIMENSION (:,:)   :: BXF                 !< Boundary values at x-max face
   REAL(EB), ALLOCATABLE, DIMENSION (:,:)   :: BYS                 !< Boundary values at y-min face
   REAL(EB), ALLOCATABLE, DIMENSION (:,:)   :: BYF                 !< Boundary values at y-max face
   REAL(EB), ALLOCATABLE, DIMENSION (:,:)   :: BZS                 !< Boundary values at z-min face
   REAL(EB), ALLOCATABLE, DIMENSION (:,:)   :: BZF                 !< Boundary values at z-max face

   REAL(EB):: CAPPA_POISSON = 0.0_EB                               !< Convergence rate of Poisson solution
   REAL(EB):: CAPPA_LAPLACE = 0.0_EB                               !< Max convergence rate of Laplace solutions

   REAL(EB):: VELOCITY_ERROR = 0.0_EB                              !< Velocity error in single Laplace iteration
   REAL(EB):: VELOCITY_ERROR_MAX = 0.0_EB                          !< Maximum achieved velocity error over all Laplace iterations
   REAL(EB):: VELOCITY_TOLERANCE = 1.E-4_EB                        !< Requested velocity tolerance for MGM Laplace problems

   REAL(EB):: THRESHOLD = 1.E-4_EB                                 !< Default requested velocity tolerance Laplace problems
   INTEGER :: NIT                                                  !< Default max number of iterations for Laplace problems

   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: BTYPE                  !< Boundary type of an interface cell

   INTEGER :: NCS, NCU                                             !< Number of cells for structured/unstructured grid
   INTEGER :: NW1, NW2, NWI, NWE                                   !< Range of IW's with non-zero B-values

   INTEGER :: ITE = 0                            
   INTEGER :: ITE_POISSON = 0
   INTEGER :: ITE_LAPLACE = 1

   LOGICAL :: HAS_OBSTRUCTIONS = .FALSE.

END TYPE SCARC_MGM_TYPE

!> \brief Collection of grid level related information on single mesh
  
TYPE SCARC_LEVEL_TYPE

   ! Administrative structures for different components based on given grid level
   TYPE (SCARC_FACE_TYPE), ALLOCATABLE, DIMENSION(:) :: FACE   !< Face information
   TYPE (SCARC_OBST_TYPE), ALLOCATABLE, DIMENSION(:) :: OBST   !< Obstruction information
   TYPE (SCARC_STAGE_TYPE), DIMENSION(2) :: STAGE              !< Hierarchy of solvers and related working vectors
   TYPE (SCARC_GRID_TYPE)      :: STRUCTURED, UNSTRUCTURED     !< Structured and unstructured grid information
   TYPE (SCARC_FFT_TYPE)       :: FFT                          !< FFT preconditioner based on CRAYFISHPAK
   TYPE (SCARC_MGM_TYPE)       :: MGM                          !< McKenney-Greengard-Mayo method 
   TYPE (SCARC_MULTIGRID_TYPE) :: MG                           !< Multigrid method information
   TYPE (SCARC_PRECON_TYPE)    :: PRECON                       !< Preconditioning information

#ifdef WITH_MKL
   TYPE (SCARC_MKL_TYPE)       :: MKL                          !< MKL preconditioner based on Intel MKL
#endif

   ! Coordinate information
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: XCOR, YCOR, ZCOR    !< Coordinate vectors in x-, y- and z-direction
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: XMID, YMID, ZMID    !< Midpoint vectors in x-, y- and z-direction
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: DXL, DYL, DZL       !< Step size vectors in x-, y- and z-direction
   REAL(EB) :: DX , DY , DZ                                    !< Step sizes in x-, y- and z-direction
   REAL(EB) :: DX2, DY2, DZ2                                   !< Half step sizes in x-, y- and z-direction
   REAL(EB) :: DXH, DYH, DZH                                   !< Half step sizes in x-, y- and z-direction
   REAL(EB) :: DXI, DYI, DZI                                   !< Inversed of step sizes in x-, y- and z-direction
   REAL(EB) :: DXI2, DYI2, DZI2                                !< Squared and inversed step sizes in x-, y- and z-direction

   ! Cell and wall index information:
   !    - on the finest level, the original arrays from FDS are used
   !    - separate arrays will only be allocated for coarser levels
   !    - to address them on all levels, corresponding pointers are used
   INTEGER, POINTER, DIMENSION (:,:,:)     :: CELL_INDEX_PTR   !< Pointer to cell index array
   INTEGER, POINTER, DIMENSION (:,:)       :: WALL_INDEX_PTR   !< Pointer to wall index array
   INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: CELL_INDEX       !< Cell index list (only allocated for coarser levels)
   INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: WALL_INDEX       !< Wall index list (only allocated for coarser levels)

   LOGICAL, ALLOCATABLE, DIMENSION (:,:,:) :: IS_SOLID         !< State of single cells (.TRUE. if solid/.FALSE. otherwise)

   ! Orientation and different cell related lengths
   INTEGER :: GHOST_FIRSTW(-3:3) = 0                           !< First internal ghost cell numbers for all faces
   INTEGER :: GHOST_LASTW(-3:3)  = 0                           !< Last internal ghost cell numbers for all faces
   INTEGER :: GHOST_FIRSTE(-3:3) = 0                           !< First external ghost cell numbers for all faces
   INTEGER :: GHOST_LASTE(-3:3)  = 0                           !< Last external ghost cell numbers for all faces
   INTEGER :: NX = 0                                           !< Number of grid cells in x-direction
   INTEGER :: NY = 0                                           !< Number of grid cells in y-direction
   INTEGER :: NZ = 0                                           !< Number of grid cells in z-direction

   ! Number of discretizations and obstructions
   INTEGER :: N_DISCRET = 0                                    !< Number of discretization types used
   INTEGER :: N_OBST = 0                                       !< Number of obstructions
   INTEGER :: N_CELL_INDEX = 0                                 !< Number of entries in CELL_INDEX array
   INTEGER :: N_CELLS = 0                                      !< Number of cells in structured discretization
   INTEGER :: N_GHOST_ZONES = 0                                !< Number of adjacent ghost zones
   INTEGER :: N_LAYER2 = 0                                     !< Number of cells in second layer to neighbor (AMG only)
   INTEGER :: N_LAYER2_TOTAL = 0                               !< Total number of cells in second layer (AMG only) 

   ! Different wall related lengths
   INTEGER :: N_WALL_CELLS = 0                                 !< Number of wall cells
   INTEGER :: N_WALL_CELLS_EXT = 0                             !< Number of external wall cells
   INTEGER :: N_WALL_CELLS_INT = 0                             !< Number of internal wall cells
   INTEGER :: N_WALL_CELLS_LOCAL = 0                           !< Number of local wall cells

   INTEGER :: L2PTR = 0                                        !< Pointer to current second layer cell

   LOGICAL :: HAS_OBSTRUCTIONS = .FALSE.                       !< Information whether current mesh has obstructions or not

END TYPE SCARC_LEVEL_TYPE

!> \brief Sample sequence of used solvers in stack
  
TYPE SCARC_SOLVER_TYPE

   CHARACTER(30) :: CNAME = 'NONE'                             !< Name of current solver

   ! Types of different solver components
   INTEGER :: TYPE_ACCURACY      = NSCARC_ACCURACY_ABSOLUTE    !< Type of requested accuracy
   INTEGER :: TYPE_COARSE        = NSCARC_COARSE_DIRECT        !< Type of coarse grid solver for multilevel methods
   INTEGER :: TYPE_COARSENING    = NSCARC_COARSENING_CUBIC     !< Type of grid coarsening 
   INTEGER :: TYPE_CYCLING       = NSCARC_CYCLING_V            !< Type of cycling for multigrid method
   INTEGER :: TYPE_GRID          = NSCARC_GRID_STRUCTURED      !< Type of discretization
   INTEGER :: TYPE_EXCHANGE      = NSCARC_UNDEF_INT            !< Type of data exchange
   INTEGER :: TYPE_INTERPOL      = NSCARC_INTERPOL_CONSTANT    !< Type of interpolation method
   INTEGER :: TYPE_LEVEL(0:2)    = NSCARC_UNDEF_INT            !< Type of levels
   INTEGER :: TYPE_MATRIX        = NSCARC_MATRIX_COMPACT       !< Type of storage for matrix
   INTEGER :: TYPE_METHOD        = NSCARC_METHOD_KRYLOV        !< Type of ScaRC method
   INTEGER :: TYPE_MKL(0:10)     = NSCARC_UNDEF_INT            !< Type of MKL for single levels
   INTEGER :: TYPE_MKL_PRECISION = NSCARC_PRECISION_DOUBLE     !< Type of precision for MKL solver
   INTEGER :: TYPE_MULTIGRID     = NSCARC_MULTIGRID_GEOMETRIC  !< Type of multigrid method
   INTEGER :: TYPE_PARENT        = NSCARC_UNDEF_INT            !< Type of parent (calling) solver
   INTEGER :: TYPE_PRECON        = NSCARC_UNDEF_INT            !< Type of preconditioner for iterative solver
   INTEGER :: TYPE_RELAX         = NSCARC_UNDEF_INT            !< Type of preconditioner for iterative solver
   INTEGER :: TYPE_SCOPE(0:2)    = NSCARC_SCOPE_GLOBAL         !< Type of solver scopes
   INTEGER :: TYPE_SMOOTH        = NSCARC_UNDEF_INT            !< Type of smoother for multigrid method
   INTEGER :: TYPE_SOLVER        = NSCARC_SOLVER_MAIN          !< Type of surrounding solver stage
   INTEGER :: TYPE_STAGE         = NSCARC_STAGE_ONE            !< Type of surrounding solver stage
   INTEGER :: TYPE_STENCIL       = NSCARC_STENCIL_CONSTANT     !< Type of storage for matrix
   INTEGER :: TYPE_TWOLEVEL      = NSCARC_TWOLEVEL_NONE        !< Type of two-level method
   INTEGER :: TYPE_VECTOR        = NSCARC_UNDEF_INT            !< Type of vector to point to


   ! References to different vectors which are needed for the current solver
   INTEGER :: X = NSCARC_UNDEF_INT                             !< Reference to local X-vector, double precision
   INTEGER :: B = NSCARC_UNDEF_INT                             !< Reference to local B-vector, double precision
   INTEGER :: D = NSCARC_UNDEF_INT                             !< Reference to local D-vector, double precision
   INTEGER :: R = NSCARC_UNDEF_INT                             !< Reference to local R-vector, double precision
   INTEGER :: V = NSCARC_UNDEF_INT                             !< Reference to local V-vector, double precision
   INTEGER :: Y = NSCARC_UNDEF_INT                             !< Reference to local Y-vector, double precision
   INTEGER :: Z = NSCARC_UNDEF_INT                             !< Reference to local Z-vector, double precision

   INTEGER :: X_FB = NSCARC_UNDEF_INT                          !< Reference to local X-vector, single precision
   INTEGER :: B_FB = NSCARC_UNDEF_INT                          !< Reference to local B-vector, single precision
   INTEGER :: R_FB = NSCARC_UNDEF_INT                          !< Reference to local R-vector, single precision
   INTEGER :: V_FB = NSCARC_UNDEF_INT                          !< Reference to local V-vector, single precision


   ! Convergence requirements for current solver

   INTEGER  :: NIT   = NSCARC_UNDEF_INT                        !< Maximum iteration number
   INTEGER  :: ITE   = NSCARC_UNDEF_INT                        !< Current iteration number

   REAL(EB) :: EPS   = NSCARC_UNDEF_REAL_EB                    !< Required accuracy
   REAL(EB) :: RES   = NSCARC_UNDEF_REAL_EB                    !< Current residual
   REAL(EB) :: RESIN = NSCARC_UNDEF_REAL_EB                    !< Initial residual
   REAL(EB) :: ERR   = NSCARC_UNDEF_REAL_EB                    !< Initial residual
   REAL(EB) :: OMEGA = NSCARC_UNDEF_REAL_EB                    !< Relaxation parameter
   REAL(EB) :: CAPPA = NSCARC_UNDEF_REAL_EB                    !< Convergence rate

   REAL(EB) :: VELOCITY_ERROR_GLOBAL = NSCARC_UNDEF_REAL_EB    !< Global velocity error
   REAL(EB) :: VELOCITY_ERROR_MGM = NSCARC_UNDEF_REAL_EB       !< MGM velocity error

END TYPE SCARC_SOLVER_TYPE

!> \brief Stack type
  
TYPE SCARC_STACK_TYPE
   TYPE (SCARC_SOLVER_TYPE), POINTER :: SOLVER                 !< Parameters for different solvers
   INTEGER :: NSTAGE                                           !< Stage of current solver
END TYPE SCARC_STACK_TYPE

!> \brief Administration other mesh data needed for the coupling of adjacent neighbors
  
TYPE SCARC_OSCARC_TYPE

   TYPE (SCARC_LEVEL_TYPE), ALLOCATABLE, DIMENSION(:) :: LEVEL     !< Level related information

   REAL(EB) :: SEND_BUFFER_REAL0(1:NSCARC_MAX_BUFFER0)             !< Constant length send buffer for setup
   REAL(EB) :: RECV_BUFFER_REAL0(1:NSCARC_MAX_BUFFER0)             !< Constant length receive buffer for setup
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: SEND_BUFFER_REAL        !< Real send buffer 
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: RECV_BUFFER_REAL        !< Real receive buffer 

   INTEGER :: SEND_BUFFER_INT0(NSCARC_MAX_BUFFER0)                 !< Constant length send buffer for setup
   INTEGER :: RECV_BUFFER_INT0(NSCARC_MAX_BUFFER0)                 !< Constant length receive buffer for setup
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: SEND_BUFFER_INT         !< Integer send buffer 
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: RECV_BUFFER_INT         !< Integer receive buffer 

   INTEGER :: NLEN_MAX_BUFFER_LAYER1  = NSCARC_HUGE_INT            !< Length for single layer length exchange
   INTEGER :: NLEN_MAX_BUFFER_LAYER2  = NSCARC_HUGE_INT            !< Length for double layer length exchange
   INTEGER :: NLEN_MAX_BUFFER_LAYER4  = NSCARC_HUGE_INT            !< Length for double layer length exchange
   INTEGER :: NLEN_MAX_BUFFER_STENCIL = NSCARC_HUGE_INT            !< Length for stencil layer length exchange
   INTEGER :: NLEN_MAX_BUFFER_FULL    = NSCARC_HUGE_INT            !< Length for max length exchange

END TYPE SCARC_OSCARC_TYPE

!> \brief Measurement of CPU times
  
TYPE SCARC_CPU_TYPE
   REAL(EB) :: BUFFER_PACKING   = 0.0_EB                           !< Time for data exchange
   REAL(EB) :: BUFFER_UNPACKING = 0.0_EB                           !< Time for data exchange
   REAL(EB) :: AMG              = 0.0_EB                           !< Time for algebraic multigrid solver
   REAL(EB) :: COARSE           = 0.0_EB                           !< Time for coarse grid solver
   REAL(EB) :: EXCHANGE         = 0.0_EB                           !< Time for data exchange
   REAL(EB) :: ITERATION        = 0.0_EB                           !< Time for Krylov solver
   REAL(EB) :: L2NORM           = 0.0_EB                           !< Time for l2-norm
   REAL(EB) :: MATVEC_PRODUCT   = 0.0_EB                           !< Time for matrix vector multiplication
   REAL(EB) :: OVERALL          = 0.0_EB                           !< Complete time for ScaRC
   REAL(EB) :: RELAXATION       = 0.0_EB                           !< Time for relaxation
   REAL(EB) :: SCALAR_PRODUCT   = 0.0_EB                           !< Time for scalar product
   REAL(EB) :: SETUP            = 0.0_EB                           !< Time for setup of requested ScaRC solver
   REAL(EB) :: SMOOTHER         = 0.0_EB                           !< Time for smoothing
   REAL(EB) :: SOLVER           = 0.0_EB                           !< Time for solver 
   INTEGER  :: N_TIMER          = 13                               !< Total number of timers
END TYPE SCARC_CPU_TYPE

!> \brief Basic administration type for ScaRC-method
 
TYPE SCARC_TYPE

   TYPE (SCARC_OSCARC_TYPE), ALLOCATABLE, DIMENSION(:) :: OSCARC   !< ScaRC type on other mesh
   TYPE (SCARC_LEVEL_TYPE) , ALLOCATABLE, DIMENSION(:) :: LEVEL    !< Level related information

   REAL(EB) :: XS, XF, YS, YF, ZS, ZF                              !< x-, y- and z-bounds of grid
   REAL(EB) :: RHS_END = 0.0_EB                                    !< Very last RHS entry, needed for matrix condensing

   INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS                 !< List of adjacent neighbors of whole mesh
   INTEGER :: N_NEIGHBORS = 0                                      !< Number of adjacent neighbors of whole mesh
   INTEGER :: NC = 0                                               !< Total number of cells on that mesh
   INTEGER :: IBAR, JBAR, KBAR                                     !< Number of cells (corresponding to main prg)

END TYPE SCARC_TYPE

END MODULE SCARC_TYPES


!=======================================================================================================================
!
! MODULE SCARC_VARIABLES
!
!> \brief Define the variables used in the different routines of ScaRC/UScaRC
!
!=======================================================================================================================
MODULE SCARC_VARIABLES

USE PRECISION_PARAMETERS
USE SCARC_CONSTANTS
USE SCARC_TYPES
USE MPI_F08

IMPLICIT NONE (TYPE,EXTERNAL)

! ---------- Basic definitions
 
CHARACTER(40) :: SCARC_GRID               = 'STRUCTURED'         !< Type of discretization (STRUCTURED/UNSTRUCTURED)
CHARACTER(40) :: SCARC_METHOD             = 'NONE'               !< Type of global ScaRC solver (Krylov/MULTIGRID)
CHARACTER(40) :: SCARC_MATRIX             = 'NONE'               !< Type of matrix storage (COMPACT/BANDWISE)
CHARACTER(40) :: SCARC_STENCIL            = 'VARIABLE'           !< Type of matrix stencil (CONSTANT/VARIABLE)
CHARACTER(40) :: SCARC_TWOLEVEL           = 'NONE'               !< Type of two-level method (NONE/ADDITIVE/MULTIPLICATIVE)

! ---------- General iteration parameters
 
CHARACTER(40) :: SCARC_ACCURACY           = 'ABSOLUTE'           !< Type of accuracy type (ABSOLUTE/RELATIVE)
REAL (EB)     :: SCARC_CAPPA              =  0.0_EB              !< Convergence rate of selected ScarC solver
INTEGER       :: SCARC_ITERATIONS         =  0                   !< Number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_RESIDUAL           =  0.0_EB              !< Residual of globally selected ScaRC solver

! ---------- Parameters for coarse grid method
 
CHARACTER(40) :: SCARC_COARSE             = 'DIRECT'             !< Type of coarse grid solver (ITERATIVE/DIRECT)
REAL (EB)     :: SCARC_COARSE_ACCURACY    = 1.E-14_EB            !< Requested accuracy for iterative solver
INTEGER       :: SCARC_COARSE_ITERATIONS  = 100                  !< Max number of iterations for iterative solver
INTEGER       :: SCARC_COARSE_LEVEL       =  1                   !< Coarse grid level for twolevel-Krylov method
REAL (EB)     :: SCARC_COARSE_OMEGA       = 0.80E+0_EB           !< Relaxation parameter

CHARACTER(40) :: SCARC_COARSENING  = 'GMG'                       !< Coarsening strategy (CUBIC/AGGREGATED/GMG)

! ---------- Parameters for Krylov type methods
 
REAL (EB)     :: SCARC_KRYLOV_ACCURACY    = 1.E-8_EB             !< Requested accuracy for convergence
CHARACTER(40) :: SCARC_KRYLOV_INTERPOL    = 'CONSTANT'           !< Twolevel-interpolation (CONSTANT/BILINEAR)
INTEGER       :: SCARC_KRYLOV_ITERATIONS  = 1000                 !< Max number of iterations

! ---------- Parameters for multigrid-type methods
 
CHARACTER(40) :: SCARC_MULTIGRID             = 'GEOMETRIC'       !< Type of MG method (GEOMETRIC/ALGEBRAIC)
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY    = 1.E-8_EB          !< Requested accuracy for convergence
CHARACTER(3)  :: SCARC_MULTIGRID_CYCLE       = 'V'               !< Cycling type  (F/V/W/FULL)
CHARACTER(40) :: SCARC_MULTIGRID_INTERPOL    = 'CONSTANT'        !< Interpolation strategy (CONSTANT/BILINEAR)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS  = 100               !< Max number of iterations
INTEGER       :: SCARC_MULTIGRID_LEVEL       = -1                !< User defined number of MG-levels (optionally, maximum else)
INTEGER       :: SCARC_MULTIGRID_PRESMOOTH   = 4                 !< Number of presmoothing iterations
INTEGER       :: SCARC_MULTIGRID_POSTSMOOTH  = 4                 !< Number of postsmoothing iterations
LOGICAL       :: SCARC_MULTIGRID_RELAXING    = .TRUE.            !< Relaxing of nullspace (AMG only)
REAL (EB)     :: SCARC_MULTIGRID_THETA       = 0.10E+0_EB        !< Threshold for strength of connection matrix (AMG only)

! ---------- Parameters for MGM method

REAL(EB)      :: SCARC_MGM_ACCURACY         = 1.E-4_EB           !< Requested accuracy for interface velocity error 
CHARACTER(40) :: SCARC_MGM_BC               = 'MEAN'             !< Type of interface boundary condition for local Laplace problems
CHARACTER(40) :: SCARC_MGM_INTERPOLATION    = 'LINEAR'           !< Type of interpolation for Laplace BC settings
CHARACTER(40) :: SCARC_MGM_LAPLACE_SOLVER   = 'LUPERM'           !< Type of solver for local Laplace problems
INTEGER       :: SCARC_MGM_ITERATIONS       = 20                 !< Maximum allowed number of Laplace iterations 
LOGICAL       :: SCARC_MGM_CHECK_LAPLACE    = .FALSE.            !< Requested check of Laplace solutions against ScaRC-UScaRC diff
LOGICAL       :: SCARC_MGM_EXACT_INITIAL    = .FALSE.            !< Exact solutions for initialization of interface BCs are used
LOGICAL       :: SCARC_MGM_USE_FFT          = .FALSE.            !< If local mesh is structured, use FFT as Laplace solver

! ---------- Parameters for smoothing method (used in multigrids-methods)
 
CHARACTER(40) :: SCARC_SMOOTH             = 'SSOR'               !< Smoother for MG (JACOBI/SSOR)
REAL (EB)     :: SCARC_SMOOTH_ACCURACY    = 1.E-8_EB             !< Requested accuracy for convergence
INTEGER       :: SCARC_SMOOTH_ITERATIONS  = 4                    !< Max number of iterations
REAL (EB)     :: SCARC_SMOOTH_OMEGA       = 0.80E+0_EB           !< Relaxation parameter
CHARACTER(40) :: SCARC_SMOOTH_SCOPE       = 'GLOBAL'             !< Scope of action (LOCAL/GLOBAL)

! ---------- Parameters for preconditioning method (used in Krylov methods)
 
CHARACTER(40) :: SCARC_PRECON             = 'NONE'               !< Preconditioner for CG (JACOBI/SSOR/FFT/PARDISO/MG)
REAL (EB)     :: SCARC_PRECON_ACCURACY    = 1.E-10_EB            !< Requested accuracy for convergence
INTEGER       :: SCARC_PRECON_ITERATIONS  = 100                  !< Max number of iterations
REAL (EB)     :: SCARC_PRECON_OMEGA       = 1.50E+0_EB           !< Relaxation parameter
CHARACTER(40) :: SCARC_PRECON_SCOPE       = 'LOCAL'              !< Scope of action (LOCAL/GLOBAL)

! ---------- Parameter for MKL solver
 
CHARACTER(40) :: SCARC_MKL_SCOPE          = 'GLOBAL'             !< Scope of MKL solver (LOCAL/GLOBAL)
CHARACTER(40) :: SCARC_MKL_MTYPE          = 'SYMMETRIC'          !< Type of MKL matrix (SYMMETRIC/UNSYMMETRIC)
CHARACTER(6)  :: SCARC_MKL_PRECISION      = 'DOUBLE'             !< Single/double precision for MKL solver

! ---------- Dump out of error information and error handling
 
LOGICAL :: SCARC_ERROR_FILE  = .FALSE.                           !< Print ScaRC statistics into chid_scarc.csv (TRUE/FALSE)
INTEGER :: IERROR  = 0                                           !< General error flag - used at different positions


! ---------- Logical indicators for different methods and mechanisms
  
LOGICAL :: IS_STRUCTURED         = .FALSE.                       !< Flag for structured discretization
LOGICAL :: IS_UNSTRUCTURED       = .FALSE.                       !< Flag for unstructured discretization
LOGICAL :: IS_PURE_NEUMANN       = .FALSE.                       !< Flag for pure Neumann system
LOGICAL :: IS_MG                 = .FALSE.                       !< Flag for Multigrid-method
LOGICAL :: IS_AMG                = .FALSE.                       !< Flag for Algebraic Multigrid-method
LOGICAL :: IS_GMG                = .FALSE.                       !< Flag for Geometric Multigrid-method
LOGICAL :: IS_CG                 = .FALSE.                       !< Flag for Krylov method
LOGICAL :: IS_CG_ADD             = .FALSE.                       !< Flag for additive twolevel-Krylov method
LOGICAL :: IS_CG_AMG             = .FALSE.                       !< Flag for Krylov method with AMG-preconditioning
LOGICAL :: IS_CG_GMG             = .FALSE.                       !< Flag for Krylov method with GMG-preconditioning
LOGICAL :: IS_CG_COARSE          = .FALSE.                       !< Flag for only coarse grid solver
LOGICAL :: IS_CG_MACRO           = .FALSE.                       !< Flag for macro coarse grid solver
LOGICAL :: IS_CG_MG              = .FALSE.                       !< Flag for Krylov method with MG-preconditioning
LOGICAL :: IS_CG_MUL1            = .FALSE.                       !< Flag for multiplicative1 Twolevel-Krylov method (fine last)
LOGICAL :: IS_CG_MUL2            = .FALSE.                       !< Flag for multiplicative2 Twolevel-Krylov method (fine first) 
LOGICAL :: IS_FFT                = .FALSE.                       !< Flag for FFT-method
LOGICAL :: IS_FFTO               = .FALSE.                       !< Flag for FFTO-method
LOGICAL :: IS_LAPLACE            = .FALSE.                       !< Flag for use of Laplace matrix (MGM only)
LOGICAL :: IS_POISSON            = .TRUE.                        !< Flag for use of Poisson matrix (MGM only)
LOGICAL :: IS_MKL                = .FALSE.                       !< Flag for MKL-method
LOGICAL :: IS_MKL_LEVEL(10)      = .FALSE.                       !< Flag for level-dependent MKL method
LOGICAL :: IS_MGM                = .FALSE.                       !< Flag for McKeeney-Greengard-Mayo method

LOGICAL :: HAS_CSV_DUMP          = .FALSE.                       !< Flag for CSV-file to be dumped out
LOGICAL :: HAS_MULTIPLE_GRIDS    = .FALSE.                       !< Flag for multiple discretization types
LOGICAL :: HAS_TWO_LEVELS        = .FALSE.                       !< Flag for two grid levels
LOGICAL :: HAS_XMEAN_LEVELS      = .FALSE.                       !< Flag for two-level xmean values preconditioning
LOGICAL :: HAS_MULTIPLE_LEVELS   = .FALSE.                       !< Flag for multiple grid levels
LOGICAL :: HAS_AMG_LEVELS        = .FALSE.                       !< Flag for AMG-based grid levels
LOGICAL :: HAS_GMG_LEVELS        = .FALSE.                       !< Flag for GMG-based grid levels

! ---------- Globally used types for different solvers parameters
  
INTEGER :: TYPE_ACCURACY            = NSCARC_ACCURACY_ABSOLUTE   !< Type of requested accuracy
INTEGER :: TYPE_COARSE              = NSCARC_COARSE_DIRECT       !< Type of coarse grid solver 
INTEGER :: TYPE_COARSENING          = NSCARC_COARSENING_CUBIC    !< Type of grid coarsening 
INTEGER :: TYPE_CYCLING             = NSCARC_CYCLING_V           !< Type of cycling for multigrid method
INTEGER :: TYPE_GRID                = NSCARC_GRID_STRUCTURED     !< Type of discretization 
INTEGER :: TYPE_EXCHANGE            = NSCARC_UNDEF_INT           !< Type of data exchange
INTEGER :: TYPE_EXCHANGE_MATRIX     = NSCARC_MATRIX_POISSON      !< Type of matrix for exchange
INTEGER :: TYPE_INTERPOL            = NSCARC_INTERPOL_CONSTANT   !< Type of interpolation method
INTEGER :: TYPE_LEVEL(0:2)          = NSCARC_UNDEF_INT           !< Type of levels
INTEGER :: TYPE_MATRIX              = NSCARC_MATRIX_COMPACT      !< Type of storage for matrix
INTEGER :: TYPE_MATVEC              = NSCARC_MATVEC_GLOBAL       !< Type of matrix-vector multiplication
INTEGER :: TYPE_METHOD              = NSCARC_METHOD_KRYLOV       !< Type of ScaRC method
INTEGER :: TYPE_MGM_BC              = NSCARC_MGM_BC_MEAN         !< Type of internal MGM boundary conditions
INTEGER :: TYPE_MGM_LAPLACE         = NSCARC_MGM_LAPLACE_LUPERM  !< Type of solver for local Laplace problems
INTEGER :: TYPE_MGM_INTERPOL        = NSCARC_MGM_INTERPOL_LINEAR !< Type of internal MGM boundary conditions
INTEGER :: TYPE_MKL(0:10)           = NSCARC_MKL_NONE            !< Type of use of MKL solvers
INTEGER :: TYPE_MKL_PRECISION       = NSCARC_PRECISION_DOUBLE    !< Type of double precision MKL solver
INTEGER :: TYPE_MULTIGRID           = NSCARC_MULTIGRID_GEOMETRIC !< Type of multigrid method 
INTEGER :: TYPE_PARENT              = NSCARC_UNDEF_INT           !< Type of parent (calling) solver
INTEGER :: TYPE_PRECON              = NSCARC_UNDEF_INT           !< Type of preconditioner for iterative solver
INTEGER :: TYPE_RELAX               = NSCARC_UNDEF_INT           !< Type of preconditioner for iterative solver
INTEGER :: TYPE_SCOPE(0:2)          = NSCARC_SCOPE_GLOBAL        !< Type of method scopes
INTEGER :: TYPE_SMOOTH              = NSCARC_UNDEF_INT           !< Type of smoother for multigrid method
INTEGER :: TYPE_SOLVER              = NSCARC_SOLVER_MAIN         !< Type of surrounding solver stage
INTEGER :: TYPE_STAGE               = NSCARC_STAGE_ONE           !< Type of surrounding solver stage
INTEGER :: TYPE_STENCIL             = NSCARC_STENCIL_VARIABLE    !< Type of storage for matrix
INTEGER :: TYPE_TWOLEVEL            = NSCARC_TWOLEVEL_NONE       !< Type of twolevel method
INTEGER :: TYPE_VECTOR              = NSCARC_UNDEF_INT           !< Type of vector to point to

! ---------- Globally used parameters
 
INTEGER :: NLEVEL_MIN, NLEVEL_MAX                           !< Minimum and maximum number of multigrid levels
INTEGER :: NC_GLOBAL(20)  = 0                               !< Number of global cells
INTEGER :: NX_GLOBAL      = 0                               !< Number of cells in x-direction on finest level
INTEGER :: N_DIRIC_GLOBAL(20)  = 0                          !< Global number of Dirichlet BCs
INTEGER :: N_STACK_TOTAL  = 0                               !< Maximum number of used solvers in stack
INTEGER :: N_STACK_LAPLACE = 0                              !< Stack position of local MGM Laplace problems

INTEGER :: N_REQ, N_EXCHANGES, TAG                          !< Information for data exchange
INTEGER :: SNODE, RNODE                                     !< Process Indicator for data exchange

TYPE (MPI_REQUEST),  ALLOCATABLE, DIMENSION (:)  :: REQ     !< Request array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: COUNTS             !< Counter array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: DISPLS             !< Displacement array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: MESH_INT           !< Local integer data array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: NX_OFFSET          !< Offset array of offsets in x-direction
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: MESH_REAL          !< Local real data array for data exchange

INTEGER  :: GLOBAL_INT,  RANK_INT
REAL(EB) :: GLOBAL_REAL, RANK_REAL

INTEGER :: FACE_ORIENTATION(6)  = (/1,-1,2,-2,3,-3/)        !< Coordinate direction related order of mesh faces

CHARACTER(60) :: CNAME, CROUTINE

! ---------- Public variables
  
PUBLIC :: SCARC_ACCURACY                  !< Requested accuracy for ScaRC solver
PUBLIC :: SCARC_CAPPA                     !< Resulting convergence rate of ScaRC solver
PUBLIC :: SCARC_COARSENING                !< Selection parameter for AMG coarsening strategy (aggregation/cubic)
PUBLIC :: SCARC_ERROR_FILE                !< Flag to print additional convergence information about current ScaRC call
PUBLIC :: SCARC_GRID                      !< Selection parameter for requested grid variant (structured/unstructured)
PUBLIC :: SCARC_ITERATIONS                !< Final number of needed iterations for ScaRC solver
PUBLIC :: SCARC_MATRIX                    !< Selection parameter for requested matrix storage technique (compact/bandwise)
PUBLIC :: SCARC_METHOD                    !< Selection parameter for requested ScaRC variant (Krylov/Multigrid/LU)
PUBLIC :: SCARC_MKL_PRECISION             !< Selection parameter for requested MKL precision (double/single)
PUBLIC :: SCARC_RESIDUAL                  !< Final residual after call of ScaRC solver
PUBLIC :: SCARC_TWOLEVEL                  !< Selection parameter for possible twolevel variant (additive/multiplicative)

PUBLIC :: SCARC_COARSE                    !< Selection parameter for type of coarse grid solver (iterative/direct)
PUBLIC :: SCARC_COARSE_ACCURACY           !< Requested accuracy for coarse grid solver
PUBLIC :: SCARC_COARSE_ITERATIONS         !< Maximum number of allowed coarse grid iterations (direct variant only)
PUBLIC :: SCARC_COARSE_OMEGA              !< Relaxation parameter for coarse grid solver (direct variant only)
PUBLIC :: SCARC_COARSE_LEVEL              !< Grid refinement level for coarse grid solver

PUBLIC :: SCARC_KRYLOV_ACCURACY           !< Requested accuracy for Krylov solver
PUBLIC :: SCARC_KRYLOV_ITERATIONS         !< Maximum number of allowed Krylov iterations
PUBLIC :: SCARC_KRYLOV_INTERPOL           !< Selection parameter for interpolation type in case of a twolevel Krylov variant

PUBLIC :: SCARC_MULTIGRID                 !< Selection parameter for multigrid method (geometric/algebraic)
PUBLIC :: SCARC_MULTIGRID_ACCURACY        !< Requested accuracy for multigrid method
PUBLIC :: SCARC_MULTIGRID_ITERATIONS      !< Maximum number of allowed iterations for multigrid method
PUBLIC :: SCARC_MULTIGRID_INTERPOL        !< Selection parameter for interpolation type between multigrid levels
PUBLIC :: SCARC_MULTIGRID_CYCLE           !< Selection parameter of multigrid cycling type (V/W/F)
PUBLIC :: SCARC_MULTIGRID_LEVEL           !< Coarse grid level for multigrid method
PUBLIC :: SCARC_MULTIGRID_PRESMOOTH       !< Number of presmoothing iterations in multigrid method
PUBLIC :: SCARC_MULTIGRID_POSTSMOOTH      !< Number of postesmoothing iterations in multigrid method
PUBLIC :: SCARC_MULTIGRID_RELAXING        !< Relaxing of nullspace
PUBLIC :: SCARC_MULTIGRID_THETA           !< Optional relaxation parameter for multigrid

PUBLIC :: SCARC_MGM_BC                    !< Interface boundary conditions for Laplace problems of MGM method
PUBLIC :: SCARC_MGM_ACCURACY              !< Requested accuracy for velocity error of MGM method
PUBLIC :: SCARC_MGM_INTERPOLATION         !< Interpolation type for BC definition
PUBLIC :: SCARC_MGM_ITERATIONS            !< Maximum number of allowed Laplace iterations for MGM method
PUBLIC :: SCARC_MGM_CHECK_LAPLACE         !< Requested check of Laplace solutions against ScaRC-UScaRC difference 
PUBLIC :: SCARC_MGM_EXACT_INITIAL         !< Exact Laplace solution for initialization of interface BC's is used
PUBLIC :: SCARC_MGM_USE_FFT               !< If local mesh is structured, use FFT as Laplace solver

PUBLIC :: SCARC_PRECON                    !< Selection parameter for preconditioner
PUBLIC :: SCARC_PRECON_ACCURACY           !< Requested accuracy for preconditioner 
PUBLIC :: SCARC_PRECON_ITERATIONS         !< Maximum number of allowed iterations for preconditioner
PUBLIC :: SCARC_PRECON_OMEGA              !< Relaxation parameter for preconditioner
PUBLIC :: SCARC_PRECON_SCOPE              !< Scope of activity for preconditioner (global/local)

PUBLIC :: SCARC_SMOOTH                    !< Selection parameter for smoother
PUBLIC :: SCARC_SMOOTH_ACCURACY           !< Requested accuracy for smoother
PUBLIC :: SCARC_SMOOTH_ITERATIONS         !< Maximum number of allowed iterations for smoother
PUBLIC :: SCARC_SMOOTH_OMEGA              !< Relaxation parameter for smoother
PUBLIC :: SCARC_SMOOTH_SCOPE              !< Scope of activity for smoother (global/local)

! ---------- Type declarations
  
TYPE (SCARC_TYPE)      , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: SCARC       !< Main ScaRC data structure
TYPE (SCARC_STACK_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: STACK       !< Stack structure for consecutive solvers
TYPE (SCARC_CPU_TYPE),   SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: CPU         !< Structure to handle CPU time measurements

TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: POISSON_SOLVER                       !< Poisson solver 
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: POISSON_SOLVER_STRUCTURED            !< Structured Poisson solver
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: POISSON_SOLVER_UNSTRUCTURED          !< Unstructured Poisson solver
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: POISSON_PRECON                       !< Preconditioner for Poisson solver
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: POISSON_PRECON_STRUCTURED            !< Preconditioner for structured Poisson solver
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: POISSON_PRECON_UNSTRUCTURED          !< Preconditioner for unstructured Poisson solver
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: POISSON_SMOOTH                       !< Smoother for Poisson solver
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: LAPLACE_SOLVER_UNSTRUCTURED          !< Unstructured Laplace solver
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: LAPLACE_PRECON_UNSTRUCTURED          !< Preconditioner for unstructured Laplace solver
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: COARSE_SOLVER                        !< Coarse grid solver
TYPE (SCARC_SOLVER_TYPE),  SAVE, TARGET :: COARSE_PRECON                        !< Preconditioner for coarse grid solver

TYPE (SCARC_STORAGE_TYPE), SAVE, TARGET :: STORAGE                              !< Structure for ScaRC storage administration
TYPE (SCARC_MESSAGE_TYPE), SAVE, TARGET :: MSG                                  !< Structure for verbose and debug messages
TYPE (SCARC_SUBDIVISION_TYPE), SAVE, TARGET :: SUBDIVISION                      !< Structure for information about subdivision

END MODULE SCARC_VARIABLES


!=======================================================================================================================
!
! MODULE SCARC_POINTERS
!
!> \brief Define and organize a series of pointers to specify the different meshes, grid levels, 
!   discretizations and matrices, etc. in combination with corresponding methods to set them
!
!=======================================================================================================================
MODULE SCARC_POINTERS

USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MESH_VARIABLES, ONLY: MESHES, MESH_TYPE, OMESH_TYPE, WALL_TYPE, EXTERNAL_WALL_TYPE
USE SCARC_CONSTANTS
USE SCARC_TYPES
USE SCARC_VARIABLES

IMPLICIT NONE (TYPE,EXTERNAL)

TYPE (MESH_TYPE), POINTER:: M=>NULL()                  !< Pointer to specified mesh (based on MESHES from base code)
TYPE (OMESH_TYPE), POINTER:: OM=>NULL()                !< Pointer to specified neighboring mesh (based on OMESH from base code)
TYPE (WALL_TYPE), POINTER:: MWC=>NULL()                !< Pointer to specified wall cell (based on WALL from base code)

TYPE (EXTERNAL_WALL_TYPE), POINTER:: EWC=>NULL()       !< Pointer to specified external wall cell

TYPE (SCARC_ALLOCATION_TYPE), POINTER:: AL=>NULL()     !< Pointer to allocated structure within ScaRC memory management

TYPE (SCARC_TYPE), POINTER:: S=>NULL()                 !< Pointer to ScaRC-structure on a specified mesh
TYPE (SCARC_OSCARC_TYPE), POINTER:: OS=>NULL()         !< Pointer to ScaRC-structure on a specified neighboring

TYPE (SCARC_LEVEL_TYPE), POINTER:: L=>NULL()           !< Pointer to specified level
TYPE (SCARC_LEVEL_TYPE), POINTER:: LF=>NULL()          !< Pointer to specified fine level
TYPE (SCARC_LEVEL_TYPE), POINTER:: LC=>NULL()          !< Pointer to specified coarse level
TYPE (SCARC_LEVEL_TYPE), POINTER:: OL=>NULL()          !< Pointer to specified level on neighboring mesh
TYPE (SCARC_LEVEL_TYPE), POINTER:: OLF=>NULL()         !< Pointer to specified fine level on neighboring mesh
TYPE (SCARC_LEVEL_TYPE), POINTER:: OLC=>NULL()         !< Pointer to specified coarse level on neighboring mesh

TYPE (SCARC_GRID_TYPE), POINTER:: G=>NULL()            !< Pointer to specified grid discretization
TYPE (SCARC_GRID_TYPE), POINTER:: GF=>NULL()           !< Pointer to specified fine grid discretization
TYPE (SCARC_GRID_TYPE), POINTER:: GC=>NULL()           !< Pointer to specified coarse grid discretization
TYPE (SCARC_GRID_TYPE), POINTER:: OG=>NULL()           !< Pointer to specified grid discretization
TYPE (SCARC_GRID_TYPE), POINTER:: OGF=>NULL()          !< Pointer to specified fine grid discretization on neighboring mesh
TYPE (SCARC_GRID_TYPE), POINTER:: OGC=>NULL()          !< Pointer to specified coarse grid discretization on neighboring mesh

TYPE (SCARC_FACE_TYPE), POINTER:: F=>NULL()            !< Pointer to specified face of grid
TYPE (SCARC_FACE_TYPE), POINTER:: FF=>NULL()           !< Pointer to specified face of fine grid level
TYPE (SCARC_FACE_TYPE), POINTER:: FC=>NULL()           !< Pointer to specified face of coarse grid level

TYPE (SCARC_OBST_TYPE), POINTER:: OB=>NULL()           !< Pointer to specified obstruction
TYPE (SCARC_OBST_TYPE), POINTER:: OBF=>NULL()          !< Pointer to specified obstruction of fine grid level
TYPE (SCARC_OBST_TYPE), POINTER:: OBC=>NULL()          !< Pointer to specified obstruction of coarse grid level

TYPE (SCARC_WALL_TYPE), POINTER:: GWC=>NULL()          !< Pointer to specified wall cell

TYPE (SCARC_WALL_TYPE), DIMENSION(:), POINTER:: W=>NULL()    !< Pointer to complete wall structure
TYPE (SCARC_WALL_TYPE), DIMENSION(:), POINTER:: WF=>NULL()   !< Pointer to wall structure on fine grid level
TYPE (SCARC_WALL_TYPE), DIMENSION(:), POINTER:: WC=>NULL()   !< Pointer to wall structure on coarse grid level

TYPE (SCARC_SUBDIVISION_TYPE), POINTER:: SUB=>NULL()         !< Pointer to subdivision structure (only shortcut)

TYPE (SCARC_SOLVER_TYPE), POINTER:: SV=>NULL()               !< Pointer to ScaRC solver structure
TYPE (SCARC_SOLVER_TYPE), POINTER:: SVP=>NULL()              !< Pointer to parent ScaRC solver structure

TYPE (SCARC_STAGE_TYPE),  POINTER:: ST=>NULL()               !< Pointer to solver stage structure
TYPE (SCARC_STAGE_TYPE),  POINTER:: STP=>NULL()              !< Pointer to parent solver stage structure

TYPE (SCARC_FFT_TYPE),    POINTER:: FFT=>NULL()              !< Pointer to FFT structure
TYPE (SCARC_MGM_TYPE),    POINTER:: MGM=>NULL()              !< Pointer to McKeeney-Greengard-Mayo structure
TYPE (SCARC_PRECON_TYPE), POINTER:: PRE=>NULL()              !< Pointer to preconditioning structure

TYPE (SCARC_BMATRIX_TYPE), POINTER:: AB=>NULL()       !< Pointer to bandwise matrix structure
TYPE (SCARC_BMATRIX_TYPE), POINTER:: OAB=>NULL()      !< Pointer to neighboring bandwise matrix structure

TYPE (SCARC_CMATRIX_TYPE), POINTER:: A=>NULL()        !< Pointer to compactly stored matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: A1=>NULL()       !< Pointer to compactly stored matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: A2=>NULL()       !< Pointer to compactly stored matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: AC=>NULL()       !< Pointer to compactly stored coarse matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: AF=>NULL()       !< Pointer to compactly stored fine matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OA=>NULL()       !< Pointer to compactly stored neighboring matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OAC=>NULL()      !< Pointer to compactly stored coarse neighboring matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OAF=>NULL()      !< Pointer to compactly stored fine neighboring matrix

TYPE (SCARC_CMATRIX_TYPE), POINTER:: LO=>NULL()       !< Pointer to compactly stored lower matrix from MGM-LU-decomposition
TYPE (SCARC_CMATRIX_TYPE), POINTER:: UP=>NULL()       !< Pointer to compactly stored upper matrix from MGM-LU-decomposition

TYPE (SCARC_CMATRIX_TYPE), POINTER:: P=>NULL()        !< Pointer to compactly stored Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: PC=>NULL()       !< Pointer to compactly stored coarse Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: PF=>NULL()       !< Pointer to compactly stored fine matrix on coarse grid
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OP=>NULL()       !< Pointer to compactly stored neighboring Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OPC=>NULL()      !< Pointer to compactly stored coarse neighboring Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OPF=>NULL()      !< Pointer to compactly stored fine neighboring Prolongation matrix

TYPE (SCARC_CMATRIX_TYPE), POINTER:: R=>NULL()        !< Pointer to compactly stored Restriction matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: RC=>NULL()       !< Pointer to compactly stored coarse Restriction matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: RF=>NULL()       !< Pointer to compactly stored fine Restriction matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OR=>NULL()       !< Pointer to compactly stored neighboring Restriction matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: ORC=>NULL()      !< Pointer to compactly stored coarse neighboring Restriction matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: ORF=>NULL()      !< Pointer to compactly stored fine neighboring Restriction matrix

TYPE (SCARC_CMATRIX_TYPE), POINTER:: C=>NULL()        !< Pointer to compactly stored connection matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: CC=>NULL()       !< Pointer to compactly stored coarse connection matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: CF=>NULL()       !< Pointer to compactly stored fine connection matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OC=>NULL()       !< Pointer to compactly stored neighboring connection matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OCC=>NULL()      !< Pointer to compactly stored coarse neighboring connection matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OCF=>NULL()      !< Pointer to compactly stored fine neighboring connection matrix

TYPE (SCARC_CMATRIX_TYPE), POINTER:: Z=>NULL()        !< Pointer to compactly stored Zones matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: ZC=>NULL()       !< Pointer to compactly stored coarse Zones matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: ZF=>NULL()       !< Pointer to compactly stored fine Zones matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OZ=>NULL()       !< Pointer to compactly stored neighboring Zones matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OZC=>NULL()      !< Pointer to compactly stored neighboring coarse Zones matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OZF=>NULL()      !< Pointer to compactly stored neighboring fine Zones matrix

TYPE (SCARC_CMATRIX_TYPE), POINTER:: PP=>NULL()       !< Pointer to compactly stored Poisson-times-Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: PPF=>NULL()      !< Pointer to compactly stored fine Poisson-times-Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: PPC=>NULL()      !< Pointer to compactly stored coarse Poisson-times-Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OPP=>NULL()      !< Pointer to compactly stored neighboring Poisson-times-Prolongation matrix

TYPE (SCARC_MATRIX_COMPACT_CONDENSED_TYPE),  POINTER:: ACO =>NULL()    !< Pointer to compactly stored condensed Poisson matrix
TYPE (SCARC_MATRIX_BANDWISE_CONDENSED_TYPE), POINTER:: ABCO=>NULL()    !< Pointer to bandwise stored condensed Poisson matrix

TYPE (SCARC_MULTIGRID_TYPE), POINTER:: MG =>NULL()                     !< Pointer to multigrid type

REAL(EB), POINTER, DIMENSION(:):: XCOR=>NULL()              !< Pointer to vector of node coordinates in x-direction
REAL(EB), POINTER, DIMENSION(:):: YCOR=>NULL()              !< Pointer to vector of node coordinates in x-direction
REAL(EB), POINTER, DIMENSION(:):: ZCOR=>NULL()              !< Pointer to vector of node coordinates in x-direction

REAL(EB), POINTER, DIMENSION(:):: XMID=>NULL()              !< Pointer to vector of cell midpoints in x-direction
REAL(EB), POINTER, DIMENSION(:):: YMID=>NULL()              !< Pointer to vector of cell midpoints in y-direction
REAL(EB), POINTER, DIMENSION(:):: ZMID=>NULL()              !< Pointer to vector of cell midpoints in z-direction

REAL(EB), POINTER, DIMENSION(:,:):: BXS=>NULL()             !< Pointer to boundary values at x-min face
REAL(EB), POINTER, DIMENSION(:,:):: BXF=>NULL()             !< Pointer to boundary values at x-max face
REAL(EB), POINTER, DIMENSION(:,:):: BYS=>NULL()             !< Pointer to boundary values at y-min face
REAL(EB), POINTER, DIMENSION(:,:):: BYF=>NULL()             !< Pointer to boundary values at y-max face
REAL(EB), POINTER, DIMENSION(:,:):: BZS=>NULL()             !< Pointer to boundary values at z-min face
REAL(EB), POINTER, DIMENSION(:,:):: BZF=>NULL()             !< Pointer to boundary values at z-max face

INTEGER,  POINTER, DIMENSION(:,:):: BTYPE=>NULL()           !< Boundary type of face cell

REAL(EB), POINTER, DIMENSION(:):: VX=>NULL()                !< Pointer to solution vector 
REAL(EB), POINTER, DIMENSION(:):: VB=>NULL()                !< Pointer to right hand side vector
REAL(EB), POINTER, DIMENSION(:):: VC=>NULL()                !< Pointer to vector on coarse grid level
REAL(EB), POINTER, DIMENSION(:):: VF=>NULL()                !< Pointer to vector on fine grid level
REAL(EB), POINTER, DIMENSION(:):: V1=>NULL()                !< Pointer to first vector
REAL(EB), POINTER, DIMENSION(:):: V2=>NULL()                !< Pointer to second vector

REAL(EB), POINTER, DIMENSION(:,:,:):: HP=>NULL()            !< Pointer to pressure solution
REAL(EB), POINTER, DIMENSION(:,:,:):: SIP=>NULL()           !< Pointer to structured inhomogenous Poisson MGM solution
REAL(EB), POINTER, DIMENSION(:,:,:):: UIP=>NULL()           !< Pointer to unstructured inhomogenous Poisson MGM solution
REAL(EB), POINTER, DIMENSION(:,:,:):: UHL=>NULL()           !< Pointer to unstructured homogeneous Laplace MGM solution
REAL(EB), POINTER, DIMENSION(:,:,:):: UHL2=>NULL()          !< Pointer to unstructured homogeneous Laplace MGM solution II

REAL(EB), POINTER, DIMENSION(:):: OVEL=>NULL()              !< Pointer to other velocity component
REAL(EB), POINTER, DIMENSION(:):: OUIP=>NULL()              !< Pointer to other unstructured inhomogenous Poisson solution
REAL(EB), POINTER, DIMENSION(:):: OUHL=>NULL()              !< Pointer to other unstructured homogenous Laplace solution
REAL(EB), POINTER, DIMENSION(:):: OUHL2=>NULL()             !< Pointer to other unstructured homogenous Laplace solution II

REAL(EB), POINTER, DIMENSION(:,:,:):: PRHS=>NULL()          !< Pointer to right hand side vector

REAL(EB), POINTER, DIMENSION(:,:,:):: UU=>NULL()            !< Pointer to u-velocity vector
REAL(EB), POINTER, DIMENSION(:,:,:):: VV=>NULL()            !< Pointer to v-velocity vector
REAL(EB), POINTER, DIMENSION(:,:,:):: WW=>NULL()            !< Pointer to w-velocity vector

REAL(EB), POINTER, DIMENSION(:) ::  RECV_BUFFER_REAL        !< Pointer to double precision receive vector
INTEGER,  POINTER, DIMENSION(:) ::  RECV_BUFFER_INT         !< Pointer to inter receive vector

TYPE (SCARC_CMATRIX_TYPE), POINTER:: AS=>NULL()             !< Pointer to symmetric Poisson matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: ACS=>NULL()            !< Pointer to coarse symmetric Poisson matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: AFS=>NULL()            !< Pointer to fine symmetric Poisson matrix

TYPE (SCARC_CMATRIX_TYPE), POINTER:: OAS=>NULL()            !< Pointer to neighboring symmetric Poisson matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OACS=>NULL()           !< Pointer to coarse neighboring symmetric Poisson matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER:: OAFS=>NULL()           !< Pointer to fine neighboring symmetric Poisson matrix

REAL(FB), DIMENSION(:), POINTER:: V1_FB=>NULL()             !< Pointer to first single precision vector
REAL(FB), DIMENSION(:), POINTER:: V2_FB=>NULL()             !< Pointer to second single precision vector

#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER:: MKL=>NULL()                !< Pointer to MKL type
#endif

CONTAINS

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Unset ScaRC pointers
! mainly used to test the correctness of the pointer settings in the different routines
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_NONE

M => NULL()
S => NULL()
W => NULL()
L => NULL();  LF => NULL();  LC => NULL()
G => NULL();  GF => NULL();  GC => NULL()
F => NULL();  FF => NULL();  FC => NULL()
A => NULL();  AF => NULL();  AC => NULL()
#ifdef WITH_SCARC_AMG
P => NULL();  PF => NULL();  PC => NULL()
R => NULL();  RF => NULL();  RC => NULL()
C => NULL();  CF => NULL();  CC => NULL()
Z => NULL();  ZF => NULL();  ZC => NULL()
#endif

END SUBROUTINE SCARC_POINT_TO_NONE

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified mesh
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_MESH(NM)
INTEGER, INTENT(IN):: NM

M => MESHES(NM)
S => SCARC(NM)

END SUBROUTINE SCARC_POINT_TO_MESH

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified combination of mesh and grid level
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_LEVEL(NM, NL)
INTEGER, INTENT(IN):: NM, NL

M   => MESHES(NM)
S   => SCARC(NM)
L   => S%LEVEL(NL)
MGM => L%MGM

END SUBROUTINE SCARC_POINT_TO_LEVEL

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified combination of mesh and grid level
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_MGM(NM, NL)
INTEGER, INTENT(IN):: NM, NL

M   => MESHES(NM)
S   => SCARC(NM)
L   => S%LEVEL(NL)
MGM => L%MGM

SIP   => MGM%SIP
UIP   => MGM%UIP
UHL   => MGM%UHL
UHL2  => MGM%UHL2
OUHL  => MGM%OUHL
OUHL2 => MGM%OUHL2

END SUBROUTINE SCARC_POINT_TO_MGM

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified combination of a mesh level and discretization type
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_GRID (NM, NL)
INTEGER, INTENT(IN) ::  NM, NL

CALL SCARC_POINT_TO_NONE

M => MESHES(NM)
S => SCARC(NM)
L => S%LEVEL(NL)

SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      G => L%STRUCTURED
      G%NW = L%N_WALL_CELLS_EXT                     ! TODO: set it elsewhere
   CASE (NSCARC_GRID_UNSTRUCTURED)
      G => L%UNSTRUCTURED
      G%NW = L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT
END SELECT

MGM => L%MGM
W   => G%WALL

END SUBROUTINE SCARC_POINT_TO_GRID

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified pairing of mesh levels and discretization types
! Note that pointers LF, GF, WF correspond to the fine level and LC, GC, WC to the coarsest level
! Additionally the usual L, G, W pointers are also set to the fine level for outside use
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_MULTIGRID (NM, NL1, NL2)
INTEGER, INTENT(IN) ::  NM, NL1, NL2

CALL SCARC_POINT_TO_NONE

M => MESHES(NM)
S => SCARC(NM)

L  => SCARC(NM)%LEVEL(NL1)                           ! L is set to the fines level
LF => SCARC(NM)%LEVEL(NL1)
LC => SCARC(NM)%LEVEL(NL2)

SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      G  => L%STRUCTURED                             ! G is set to the finest grid
      GF => LF%STRUCTURED
      GC => LC%STRUCTURED
      GF%NW = LF%N_WALL_CELLS_EXT                    ! TODO: set it elsewhere
      GC%NW = LC%N_WALL_CELLS_EXT
   CASE (NSCARC_GRID_UNSTRUCTURED)
      G  => L%UNSTRUCTURED
      GF => LF%UNSTRUCTURED
      GC => LC%UNSTRUCTURED
      GF%NW = LF%N_WALL_CELLS_EXT+LF%N_WALL_CELLS_INT
      GC%NW = LC%N_WALL_CELLS_EXT+LC%N_WALL_CELLS_INT
END SELECT

WC => GC%WALL
WF => GF%WALL

END SUBROUTINE SCARC_POINT_TO_MULTIGRID

! -------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified combination of a neighboring mesh level and discretization type
! Note that pointers OLF, OGF correspond to the fine level 
! Additionally the usual OL, OG pointers are also set to the fine level for outside use
! -------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
INTEGER, INTENT(IN):: NM, NOM, NL

OS  => SCARC(NM)%OSCARC(NOM)
OL  => OS%LEVEL(NL)
OLF => OS%LEVEL(NL)

SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      OG  => OL%STRUCTURED
      OGF => OLF%STRUCTURED
   CASE (NSCARC_GRID_UNSTRUCTURED)
      OG  => OL%UNSTRUCTURED
      OGF => OLF%UNSTRUCTURED
END SELECT

END SUBROUTINE SCARC_POINT_TO_OTHER_GRID

! ---------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified combination of a neighboring mesh level and a discretization type
! Note that pointers OLF, OGF correspond to the fine level and OLC, OGC to the coarsest level
! Additionally the usual OL, OG pointers are also set to the fine level for outside use
! ---------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_OTHER_MULTIGRID (NM, NOM, NL1, NL2)
INTEGER, INTENT(IN):: NM, NOM, NL1, NL2

OS  => SCARC(NM)%OSCARC(NOM)
OL  => OS%LEVEL(NL1)
OLF => OS%LEVEL(NL1)
OLC => OS%LEVEL(NL2)

SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      OG  => OL%STRUCTURED
      OGF => OLF%STRUCTURED
      OGC => OLC%STRUCTURED
   CASE (NSCARC_GRID_UNSTRUCTURED)
      OG  => OL%UNSTRUCTURED
      OGF => OLF%UNSTRUCTURED
      OGC => OLC%UNSTRUCTURED
END SELECT

END SUBROUTINE SCARC_POINT_TO_OTHER_MULTIGRID

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified matrix in compact storage technique
! --------------------------------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_CMATRIX (NTYPE)
TYPE(SCARC_CMATRIX_TYPE), POINTER:: SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN):: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_CMATRIX => G%POISSON
   CASE (NSCARC_MATRIX_LAPLACE)
      SCARC_POINT_TO_CMATRIX => G%LAPLACE
#ifdef WITH_MKL
   CASE (NSCARC_MATRIX_POISSON_SYM)
      SCARC_POINT_TO_CMATRIX => G%POISSON_SYM
   CASE (NSCARC_MATRIX_LAPLACE_SYM)
      SCARC_POINT_TO_CMATRIX => G%LAPLACE_SYM
#endif
   CASE (NSCARC_MATRIX_PROLONGATION)
      SCARC_POINT_TO_CMATRIX => G%PROLONGATION
   CASE (NSCARC_MATRIX_RESTRICTION)
      SCARC_POINT_TO_CMATRIX => G%RESTRICTION
   CASE (NSCARC_MATRIX_LOWER)
      SCARC_POINT_TO_CMATRIX => G%LOWER
   CASE (NSCARC_MATRIX_UPPER)
      SCARC_POINT_TO_CMATRIX => G%UPPER
#ifdef WITH_SCARC_AMG
   CASE (NSCARC_MATRIX_POISSON_PROL)
      SCARC_POINT_TO_CMATRIX => G%POISSON_PROL
   CASE (NSCARC_MATRIX_CONNECTION)
      SCARC_POINT_TO_CMATRIX => G%CONNECTION
   CASE (NSCARC_MATRIX_ZONES)
      SCARC_POINT_TO_CMATRIX => G%ZONES
#endif
END SELECT

END FUNCTION SCARC_POINT_TO_CMATRIX

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified matrix in bandwise storage technique
! --------------------------------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_BMATRIX (NTYPE)
TYPE(SCARC_BMATRIX_TYPE), POINTER:: SCARC_POINT_TO_BMATRIX
INTEGER, INTENT(IN):: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_BMATRIX => G%POISSONB
   CASE DEFAULT
      WRITE(*,*) 'POINT_TO_BMATRIX: Only POISSONB available in bandwise format yet!'
END SELECT

END FUNCTION SCARC_POINT_TO_BMATRIX

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified neighboring matrix in compact storage technique
! --------------------------------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_OTHER_CMATRIX (NTYPE)
TYPE(SCARC_CMATRIX_TYPE), POINTER:: SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN):: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%POISSON
#ifdef WITH_MKL
   CASE (NSCARC_MATRIX_POISSON_SYM)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%POISSON_SYM
#endif
   CASE (NSCARC_MATRIX_PROLONGATION)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%PROLONGATION
   CASE (NSCARC_MATRIX_RESTRICTION)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%RESTRICTION
#ifdef WITH_SCARC_AMG
   CASE (NSCARC_MATRIX_POISSON_PROL)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%POISSON_PROL
   CASE (NSCARC_MATRIX_CONNECTION)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%CONNECTION
   CASE (NSCARC_MATRIX_ZONES)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%ZONES
#endif
END SELECT

END FUNCTION SCARC_POINT_TO_OTHER_CMATRIX

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified neighboring matrix in bandwise storage technique
! --------------------------------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_OTHER_BMATRIX (NTYPE)
TYPE(SCARC_BMATRIX_TYPE), POINTER:: SCARC_POINT_TO_OTHER_BMATRIX
INTEGER, INTENT(IN):: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_OTHER_BMATRIX => OG%POISSONB
   CASE DEFAULT
      WRITE(*,*) 'No other bandwise matrix available yet except of POISSONB'
END SELECT

END FUNCTION SCARC_POINT_TO_OTHER_BMATRIX

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified integer receive buffer for data exchanges
! --------------------------------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_BUFFER_INT(NM, NOM, NTYPE)
INTEGER, DIMENSION(:), POINTER:: SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) ::  NM, NOM, NTYPE

SCARC_POINT_TO_BUFFER_INT => NULL()
SELECT CASE (NTYPE)
   CASE (0)
      IF (RNODE /= SNODE) THEN
         SCARC_POINT_TO_BUFFER_INT => SCARC(NM)%OSCARC(NOM)%RECV_BUFFER_INT0
      ELSE
         SCARC_POINT_TO_BUFFER_INT => SCARC(NOM)%OSCARC(NM)%SEND_BUFFER_INT0
      ENDIF
   CASE (1)
      IF (RNODE /= SNODE) THEN
         SCARC_POINT_TO_BUFFER_INT => SCARC(NM)%OSCARC(NOM)%RECV_BUFFER_INT
      ELSE
         SCARC_POINT_TO_BUFFER_INT => SCARC(NOM)%OSCARC(NM)%SEND_BUFFER_INT
      ENDIF
END SELECT

END FUNCTION SCARC_POINT_TO_BUFFER_INT

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified integer receive buffer for data exchanges
! --------------------------------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_BUFFER_REAL(NM, NOM, NTYPE)
REAL(EB), DIMENSION(:), POINTER:: SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) ::  NM, NOM, NTYPE

SCARC_POINT_TO_BUFFER_REAL => NULL()
SELECT CASE (NTYPE)
   CASE (0)
      IF (RNODE /= SNODE) THEN
         SCARC_POINT_TO_BUFFER_REAL => SCARC(NM)%OSCARC(NOM)%RECV_BUFFER_REAL0
      ELSE
         SCARC_POINT_TO_BUFFER_REAL => SCARC(NOM)%OSCARC(NM)%SEND_BUFFER_REAL0
      ENDIF
   CASE (1)
      IF (RNODE /= SNODE) THEN
         SCARC_POINT_TO_BUFFER_REAL => SCARC(NM)%OSCARC(NOM)%RECV_BUFFER_REAL
      ELSE
         SCARC_POINT_TO_BUFFER_REAL => SCARC(NOM)%OSCARC(NM)%SEND_BUFFER_REAL
      ENDIF
END SELECT

END FUNCTION SCARC_POINT_TO_BUFFER_REAL

! ----------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified vector on a given grid level
! ----------------------------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_VECTOR (NM, NL, NV)
REAL(EB), POINTER, DIMENSION(:):: SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN):: NM, NL, NV

SELECT CASE (NV)

   ! Stage one vectors (for methods on first hierarchical level)

   CASE (NSCARC_VECTOR_ONE_X)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%X
   CASE (NSCARC_VECTOR_ONE_B)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%B
   CASE (NSCARC_VECTOR_ONE_D)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%D
   CASE (NSCARC_VECTOR_ONE_R)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%R
   CASE (NSCARC_VECTOR_ONE_V)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%V
   CASE (NSCARC_VECTOR_ONE_Y)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%Y
   CASE (NSCARC_VECTOR_ONE_Z)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%Z

   ! Stage two vectors (for methods on second hierarchical level)

   CASE (NSCARC_VECTOR_TWO_X)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%X
   CASE (NSCARC_VECTOR_TWO_B)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%B
   CASE (NSCARC_VECTOR_TWO_D)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%D
   CASE (NSCARC_VECTOR_TWO_R)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%R
   CASE (NSCARC_VECTOR_TWO_V)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%V
   CASE (NSCARC_VECTOR_TWO_Y)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%Y
   CASE (NSCARC_VECTOR_TWO_Z)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%Z
END SELECT

END FUNCTION SCARC_POINT_TO_VECTOR

! ----------------------------------------------------------------------------------------------------------------------
!> \brief Point to specified vector on a given grid level (single precision version)
! ----------------------------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_VECTOR_FB(NM, NL, NV)
REAL(FB), POINTER, DIMENSION(:):: SCARC_POINT_TO_VECTOR_FB
INTEGER, INTENT(IN):: NM, NL, NV

SELECT CASE (NV)

   ! Stage one vectors (for methods on first hierarchical level)

   CASE (NSCARC_VECTOR_ONE_X)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%X_FB
   CASE (NSCARC_VECTOR_ONE_B)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%B_FB
   CASE (NSCARC_VECTOR_ONE_R)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%R_FB
   CASE (NSCARC_VECTOR_ONE_V)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%V_FB

   ! Stage two vectors (for methods on second hierarchical level)

   CASE (NSCARC_VECTOR_TWO_X)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%X_FB
   CASE (NSCARC_VECTOR_TWO_B)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%B_FB
   CASE (NSCARC_VECTOR_TWO_R)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%R_FB
   CASE (NSCARC_VECTOR_TWO_V)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%V_FB
END SELECT

END FUNCTION SCARC_POINT_TO_VECTOR_FB

! ----------------------------------------------------------------------------------------------------------------------
!> \brief Point to pressure vector in predictor or corrector
! ----------------------------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_HVECTOR(NM, NV)
REAL(EB), POINTER, DIMENSION(:,:,:):: SCARC_POINT_TO_HVECTOR
INTEGER, INTENT(IN):: NM, NV
SELECT CASE (NV)
   CASE (NSCARC_VECTOR_H)
      SCARC_POINT_TO_HVECTOR => MESHES(NM)%H
   CASE (NSCARC_VECTOR_HS)
      SCARC_POINT_TO_HVECTOR => MESHES(NM)%HS
END SELECT
END FUNCTION SCARC_POINT_TO_HVECTOR

END MODULE SCARC_POINTERS


!=======================================================================================================================
!
! MODULE SCARC_MESSAGES
!
!> \brief Organize directive-based messages services available in ScaRC/UScaRC 
!
!   This includes verbosing, debugging and dumping of different data in the course of the used methods
!
!=======================================================================================================================
MODULE SCARC_MESSAGES
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: GET_FILE_NUMBER
USE MESH_VARIABLES, ONLY: MESHES
USE SCARC_CONSTANTS
USE SCARC_TYPES
USE SCARC_VARIABLES

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief Setup debug file if requested
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MESSAGES

IF (SCARC_ERROR_FILE) HAS_CSV_DUMP = .TRUE.

! If requested, open file for CSV-information about convergence of different solvers
 
IF (HAS_CSV_DUMP) THEN
   IF (MY_RANK == 0) THEN
      WRITE (MSG%FILE_STAT, '(A,A)') TRIM(CHID),'_scarc.csv'
      MSG%LU_STAT = GET_FILE_NUMBER()
      OPEN (MSG%LU_STAT, FILE=MSG%FILE_STAT)
      WRITE(MSG%LU_STAT,*) '  #Pres,   Stack,  #ScaRC,     #CG,     #MG,   Level, #Smooth, SmoType, ', &
                           '#Coarse,     #LU,    Residual,   Cappa'
   ENDIF
ENDIF

END SUBROUTINE SCARC_SETUP_MESSAGES



END MODULE SCARC_MESSAGES


!=======================================================================================================================
!
! MODULE SCARC_TROUBLESHOOTING
!
!> \brief Organize program shutdown in case of an error based on corresponding error messages
!
!=======================================================================================================================
MODULE SCARC_TROUBLESHOOTING

USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE SCARC_CONSTANTS
USE SCARC_VARIABLES, ONLY: MSG

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief Shutdown ScaRC with error message
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ERROR(NERROR, CPARAM, NPARAM)
CHARACTER(*), INTENT(IN) :: CPARAM
INTEGER, INTENT(IN) :: NERROR, NPARAM
CHARACTER(80) :: CERROR

! Assign error message according to specified error

SELECT CASE (NERROR)
   CASE (NSCARC_ERROR_BOUNDARY_SUM)
      CERROR = 'Wrong boundary sum for IOR'
   CASE (NSCARC_ERROR_BOUNDARY_TYPE)
      CERROR = 'Wrong boundary type'
   CASE (NSCARC_ERROR_DIRECT_NOMKL)
      CERROR = 'Direct coarse grid solver is only working in combination with MKL'
   CASE (NSCARC_ERROR_EXCHANGE_RECV)
      CERROR = 'Wrong receive exchange structure'
   CASE (NSCARC_ERROR_EXCHANGE_DIAG)
      CERROR = 'Unknown type of diagonal data exchange'
   CASE (NSCARC_ERROR_EXCHANGE_SEND)
      CERROR = 'Wrong send exchange structure'
   CASE (NSCARC_ERROR_GRID_INDEX)
      CERROR = 'Wrong index for J'
   CASE (NSCARC_ERROR_GRID_NUMBER)
      CERROR = 'Number not divisable by 2'
   CASE (NSCARC_ERROR_GRID_NUMBERX)
      CERROR = 'Number of cells not divisable by 2 in x-direction, NC'
   CASE (NSCARC_ERROR_GRID_NUMBERY)
      CERROR = 'Number of cells not divisable by 2 in y-direction, NC'
   CASE (NSCARC_ERROR_GRID_NUMBERZ)
      CERROR = 'Number of cells not divisable by 2 in z-direction, NC'
   CASE (NSCARC_ERROR_GRID_RESOLUTION)
      CERROR = 'Wrong grid resolution at IOR'
   CASE (NSCARC_ERROR_NEIGHBOR_NUMBER)
      CERROR = 'More than 20 neighbors along one face not allowed'
   CASE (NSCARC_ERROR_NEIGHBOR_TYPE)
      CERROR = 'Wrong neighbor'
   CASE (NSCARC_ERROR_MATRIX_ALLOCATION)
      CERROR = 'Wrong specifier during allocation or deallocation of  matrix'
   CASE (NSCARC_ERROR_MATRIX_SUBDIAG)
      CERROR = 'Subdiagonal missing for system matrix'
   CASE (NSCARC_ERROR_MATRIX_SYMMETRY)
      CERROR = 'Matrix not symmetric for mesh'
   CASE (NSCARC_ERROR_MATRIX_SETUP)
      CERROR = 'Matrix setup failed for level type'
   CASE (NSCARC_ERROR_MATRIX_SIZE)
      CERROR = 'Matrix reducing failed because new length is too big for matrix'
   CASE (NSCARC_ERROR_MATRIX_COPY)
      CERROR = 'Matrix copy failed due to too already existing array'
   CASE (NSCARC_ERROR_MGM_PERMUTATION)
      CERROR = 'Error in MGM permutation for LU'
   CASE (NSCARC_ERROR_MGM_PARDISO)
      CERROR = 'Error in MGM - PARDISO solver for local Laplace systems not available'
   CASE (NSCARC_ERROR_MKL_CLUSTER)
      CERROR = 'MKL Library compile flag not defined, Cluster_Sparse_Solver not available'
   CASE (NSCARC_ERROR_MKL_INTERNAL)
      CERROR = 'The following MKL error was detected'
   CASE (NSCARC_ERROR_MKL_PARDISO)
      CERROR = 'MKL Library compile flag not defined, PARDISO solver not available'
   CASE (NSCARC_ERROR_MKL_STORAGE)
      CERROR = 'Wrong matrix storage scheme for MKL solvers, only COMPACT storage available'
   CASE (NSCARC_ERROR_MULTIGRID_LEVEL)
      CERROR = 'Wrong level for multigrid method'
   CASE (NSCARC_ERROR_PARSE_INPUT)
      CERROR = 'Wrong input parameter'
   CASE (NSCARC_ERROR_STENCIL)
      CERROR = 'Wrong type for matrix stencil - only constant or variable allowed'
   CASE (NSCARC_ERROR_STACK_SOLVER)
      CERROR = 'Wrong number of solvers in stack'
   CASE (NSCARC_ERROR_STACK_MESSAGE)
      CERROR = 'Too many messages in calling stack'
   CASE (NSCARC_ERROR_VECTOR_LENGTH)
      CERROR = 'Inconsistent length for vector allocation'
END SELECT

! Specify more detailed information if available

IF (CPARAM /= SCARC_NONE) THEN
   IF (MY_RANK == 0) WRITE(LU_ERR,1000)  CERROR, CPARAM, TRIM(CHID)
ELSE IF (NPARAM /= NSCARC_NONE) THEN
   IF (MY_RANK == 0) WRITE(LU_ERR,2000)  CERROR, NPARAM, TRIM(CHID)
ELSE
   IF (MY_RANK == 0) WRITE(LU_ERR,3000)  CERROR, TRIM(CHID)
ENDIF

! Also print verbose message if enabled



STOP_STATUS = SETUP_STOP
RETURN

1000 FORMAT('Error in ScaRC-solver: ', A,' : ',   A, ' (CHID: ',A,')' )
2000 FORMAT('Error in ScaRC-solver: ', A,' : ', I12, ' (CHID: ',A,')' )
3000 FORMAT('Error in ScaRC-solver: ', A, ' (CHID: ',A,')' )
END SUBROUTINE SCARC_ERROR

! --------------------------------------------------------------------------------------------------------------
!> \brief Print ScaRC warning message
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_WARNING(NWARNING, CPARAM, NPARAM)
CHARACTER(*), INTENT(IN) :: CPARAM
INTEGER, INTENT(IN) :: NWARNING, NPARAM
CHARACTER(80) :: CWARNING

! Assign warning message according to specified error

SELECT CASE (NWARNING)
   CASE (NSCARC_WARNING_NO_MKL_PRECON)
      CWARNING = 'Intel MKL library missing - only SSOR preconditioner is used'
   CASE (NSCARC_WARNING_NO_MKL_SMOOTH)
      CWARNING = 'Intel MKL library missing - only SSOR smoother is used'
   CASE (NSCARC_WARNING_NO_MKL_LU)
      CWARNING = 'Intel MKL library missing - using LU instead'
   CASE (NSCARC_WARNING_ONLY_SSOR_PRECON)
      CWARNING = 'Only SSOR preconditioner available'
   CASE (NSCARC_WARNING_ONLY_SSOR_SMOOTH)
      CWARNING = 'Only SSOR smoother available'
   CASE (NSCARC_WARNING_NO_GLOBAL_SCOPE)
      CWARNING = 'No global scope preconditioner available'
END SELECT

! Specify more detailed information if available

IF (CPARAM /= SCARC_NONE) THEN
   IF (MY_RANK == 0) WRITE(LU_ERR,1000)  CWARNING, CPARAM, TRIM(CHID)
ELSE IF (NPARAM /= NSCARC_NONE) THEN
   IF (MY_RANK == 0) WRITE(LU_ERR,2000)  CWARNING, NPARAM, TRIM(CHID)
ELSE
   IF (MY_RANK == 0) WRITE(LU_ERR,3000)  CWARNING, TRIM(CHID)
ENDIF

! Also print verbose message if enabled


RETURN

1000 FORMAT('Warning in ScaRC-solver: ', A,' : ',   A, ' (CHID: ',A,')' )
2000 FORMAT('Warning in ScaRC-solver: ', A,' : ', I12, ' (CHID: ',A,')' )
3000 FORMAT('Warning in ScaRC-solver: ', A, ' (CHID: ',A,')' )
END SUBROUTINE SCARC_WARNING

END MODULE SCARC_TROUBLESHOOTING


!=======================================================================================================================
!
! MODULE SCARC_UTILITIES
!
!> \brief Provide a set of helper routines that are needed at different points in the code.
!
!=======================================================================================================================
MODULE SCARC_UTILITIES
  
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES, ONLY: MESHES
USE SCARC_CONSTANTS
USE SCARC_TYPES
USE SCARC_VARIABLES

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS


! --------------------------------------------------------------------------------------------------------------------------
!> \brief Check if two meshes are face neighbors
! --------------------------------------------------------------------------------------------------------------------------
LOGICAL FUNCTION ARE_FACE_NEIGHBORS(NM, NOM)
USE SCARC_POINTERS, ONLY: OM
INTEGER, INTENT(IN) :: NM, NOM
ARE_FACE_NEIGHBORS = .TRUE.
OM => MESHES(NM)%OMESH(NOM)
IF (OM%NIC_R == 0 .AND. OM%NIC_S == 0) ARE_FACE_NEIGHBORS = .FALSE.
END FUNCTION ARE_FACE_NEIGHBORS


! --------------------------------------------------------------------------------------------------------------------------
!> \brief Check if two meshes are diagonal neighbors
! --------------------------------------------------------------------------------------------------------------------------
LOGICAL FUNCTION ARE_DIAG_NEIGHBORS(NM, NOM)
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: N, NOM2
ARE_DIAG_NEIGHBORS = .FALSE.
NEIGHBOR_LOOP: DO N = 1, MESHES(NM)%N_NEIGHBORING_MESHES
   NOM2 = MESHES(NM)%NEIGHBORING_MESH(N)
   IF (NOM == NOM2 .OR. .NOT. ARE_FACE_NEIGHBORS (NM, NOM)) THEN
      ARE_DIAG_NEIGHBORS = .TRUE.
      EXIT NEIGHBOR_LOOP
   ENDIF
ENDDO NEIGHBOR_LOOP
END FUNCTION ARE_DIAG_NEIGHBORS


! --------------------------------------------------------------------------------------------------------------------------
!> \brief Check if two meshes are neighbors
! --------------------------------------------------------------------------------------------------------------------------
LOGICAL FUNCTION ARE_NEIGHBORS(NM, NOM)
INTEGER, INTENT(IN) :: NM, NOM
ARE_NEIGHBORS = .FALSE.
IF (ARE_FACE_NEIGHBORS (NM, NOM) .OR. ARE_DIAG_NEIGHBORS (NM, NOM)) ARE_NEIGHBORS = .TRUE.
END FUNCTION ARE_NEIGHBORS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Assign handles to currently used grid type
!  This routine assumes, that L already points to the correct level NL of mesh NL and
!  additionally sets the requested discretization type
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SET_GRID_TYPE(NTYPE)
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE (NTYPE)
   CASE (NSCARC_GRID_STRUCTURED)
      PRES_ON_WHOLE_DOMAIN = .TRUE.
      TYPE_GRID = NSCARC_GRID_STRUCTURED
      IS_STRUCTURED   = .TRUE.
      IS_UNSTRUCTURED = .FALSE.
   CASE (NSCARC_GRID_UNSTRUCTURED)
      PRES_ON_WHOLE_DOMAIN = .FALSE.
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED
      IS_STRUCTURED   = .FALSE.
      IS_UNSTRUCTURED = .TRUE.
END SELECT

END SUBROUTINE SCARC_SET_GRID_TYPE


! ------------------------------------------------------------------------------------------------------------------
!> \brief Assign handles to currently used grid type
!  This routine assumes, that L already points to the correct level NL of mesh NL and
!  additionally sets the requested discretization type
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SET_SYSTEM_TYPE(NGRID, NMATRIX)
INTEGER, INTENT(IN) :: NGRID, NMATRIX

SELECT CASE (NMATRIX)
   CASE (NSCARC_MATRIX_POISSON)
      IS_POISSON = .TRUE.
      IS_LAPLACE = .FALSE.
   CASE (NSCARC_MATRIX_LAPLACE)
      IS_POISSON = .FALSE.
      IS_LAPLACE = .TRUE.
END SELECT

SELECT CASE (NGRID)
   CASE (NSCARC_GRID_STRUCTURED)
      PRES_ON_WHOLE_DOMAIN = .TRUE.
      TYPE_GRID = NSCARC_GRID_STRUCTURED
      IS_STRUCTURED   = .TRUE.
      IS_UNSTRUCTURED = .FALSE.
   CASE (NSCARC_GRID_UNSTRUCTURED)
      PRES_ON_WHOLE_DOMAIN = .FALSE.
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED
      IS_STRUCTURED   = .FALSE.
      IS_UNSTRUCTURED = .TRUE.
END SELECT

END SUBROUTINE SCARC_SET_SYSTEM_TYPE


! --------------------------------------------------------------------------------------------------------------
!> \brief Get full text information about the data type of the currently processed array
! --------------------------------------------------------------------------------------------------------------
CHARACTER(10) FUNCTION SET_DATA_TYPE(NDATA)
INTEGER, INTENT(IN) :: NDATA

SELECT CASE (NDATA)
   CASE (NSCARC_DATA_INTEGER)
      SET_DATA_TYPE = 'INTEGER'
   CASE (NSCARC_DATA_REAL_EB)
      SET_DATA_TYPE = 'REAL_EB'
   CASE (NSCARC_DATA_REAL_FB)
      SET_DATA_TYPE = 'REAL_FB'
   CASE (NSCARC_DATA_LOGICAL)
      SET_DATA_TYPE = 'LOGICAL'
   CASE (NSCARC_DATA_CMATRIX)
      SET_DATA_TYPE = 'CMATRIX'
   CASE (NSCARC_DATA_BMATRIX)
      SET_DATA_TYPE = 'BMATRIX'
   CASE DEFAULT
      SET_DATA_TYPE = ' '
END SELECT

END FUNCTION SET_DATA_TYPE


! --------------------------------------------------------------------------------------------------------------
!> \brief Get full text information about the dimension of the currently processed array
! --------------------------------------------------------------------------------------------------------------
CHARACTER(10) FUNCTION SET_DIMENSION(NDIM)
INTEGER, INTENT(IN) :: NDIM

SELECT CASE (NDIM)
   CASE (1)
      SET_DIMENSION = '1'
   CASE (2)
      SET_DIMENSION = '2'
   CASE (3)
      SET_DIMENSION = '3'
   CASE DEFAULT
      SET_DIMENSION = '0'
END SELECT

END FUNCTION SET_DIMENSION


! --------------------------------------------------------------------------------------------------------------
!> \brief Get full text information about the initialization type of the currently processed array
! --------------------------------------------------------------------------------------------------------------
CHARACTER(10) FUNCTION SET_INIT_TYPE(NINIT, NDATA, NSTATE)
INTEGER, INTENT(IN) :: NINIT, NDATA, NSTATE

SELECT CASE (NINIT)
   CASE (NSCARC_INIT_UNDEF)
      SET_INIT_TYPE = 'UNDEF'
   CASE (NSCARC_INIT_NONE)
      SET_INIT_TYPE = 'NONE'
   CASE (NSCARC_INIT_MINUS)
      SET_INIT_TYPE = 'MINUS'
   CASE (NSCARC_INIT_ZERO)
      SET_INIT_TYPE = 'ZERO'
   CASE (NSCARC_INIT_ONE)
      SET_INIT_TYPE = 'ONE'
   CASE (NSCARC_INIT_TRUE)
      SET_INIT_TYPE = 'TRUE'
   CASE (NSCARC_INIT_FALSE)
      SET_INIT_TYPE = 'FALSE'
   CASE (NSCARC_INIT_HUGE)
      SET_INIT_TYPE = 'HUGE'
   CASE DEFAULT
      SET_INIT_TYPE = ' '
END SELECT

IF (NDATA == NSCARC_DATA_CMATRIX .OR. NDATA == NSCARC_DATA_BMATRIX) SET_INIT_TYPE = ' '
IF (NSTATE == NSCARC_STORAGE_REMOVE) SET_INIT_TYPE = ' '

END FUNCTION SET_INIT_TYPE


! --------------------------------------------------------------------------------------------------------------
!> \brief Get type of matrix storage scheme for specified grid level
! --------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SET_MATRIX_TYPE(NL)
INTEGER, INTENT(IN) :: NL

IF (NL == NLEVEL_MAX .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
   SET_MATRIX_TYPE = NSCARC_MATRIX_COMPACT
ELSE
   SET_MATRIX_TYPE = TYPE_MATRIX
ENDIF

END FUNCTION SET_MATRIX_TYPE


! --------------------------------------------------------------------------------------------------------------
!> \brief Check if a subdiagonal entry must be computed in a specified coordinate direction
! If a structured discretization is used, then subdiagonals are built in every direction
! Else check the type of the neighboring cell in direction IOR0
! --------------------------------------------------------------------------------------------------------------
LOGICAL FUNCTION IS_VALID_DIRECTION(IX, IY, IZ, IOR0)
USE SCARC_POINTERS, ONLY: L, G
INTEGER, INTENT(IN)  :: IX, IY, IZ, IOR0
INTEGER :: IC_INDEX, IW_INDEX

IS_VALID_DIRECTION = .FALSE.
IF (TWO_D .AND. ABS(IOR0) == 2) RETURN

SELECT CASE (TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      IS_VALID_DIRECTION = .TRUE.                                           ! always build subdiagonals
      RETURN
   CASE (NSCARC_GRID_UNSTRUCTURED)
      IC_INDEX = L%CELL_INDEX_PTR(IX, IY, IZ)                               ! cell index of corresponding cell
      IW_INDEX = 0
      IF (IC_INDEX /= 0) IW_INDEX  = L%WALL_INDEX_PTR(IC_INDEX, -IOR0)      ! check its wall index

      IF (IW_INDEX == 0) THEN                                               ! if zero, build subdiagonal
         IS_VALID_DIRECTION = .TRUE.
         RETURN
      ELSE                                                                  ! if not, only build along interfaces
         IF (TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL .AND. G%WALL(IW_INDEX)%BOUNDARY_TYPE== INTERPOLATED_BOUNDARY) THEN
            IS_VALID_DIRECTION = .TRUE.
            RETURN
         ENDIF
      ENDIF
END SELECT
RETURN

END FUNCTION IS_VALID_DIRECTION


! ------------------------------------------------------------------------------------------------------------------
!> \brief Get grid permutation (MGM only)
! ------------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION GET_PERM(JC)
USE SCARC_POINTERS, ONLY : G
INTEGER, INTENT(IN) :: JC
GET_PERM = -1
IF (JC > 0 .AND. JC <= G%NC) GET_PERM = G%PERM_FW(JC)
END FUNCTION GET_PERM


! ------------------------------------------------------------------------------------------------------------------
!> \brief Filter out mean value
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE FILTER_MEANVALUE(NV, NL)
USE SCARC_POINTERS, ONLY: L, G, VC, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM, IC, I, J, K

MESH_REAL = 0.0_EB
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
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
   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
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

END SUBROUTINE FILTER_MEANVALUE


! --------------------------------------------------------------------------------------------------------------
!> \brief Restore last cell of last mesh
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE RESTORE_LAST_CELL (XX, NL)
USE SCARC_POINTERS, ONLY: S, VC, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: XX, NL

IF (UPPER_MESH_INDEX /= NMESHES .OR. TYPE_RELAX == NSCARC_RELAX_FFT) RETURN
S => SCARC(UPPER_MESH_INDEX)

VC => SCARC_POINT_TO_VECTOR (UPPER_MESH_INDEX, NL, XX)
VC(S%NC) = S%RHS_END

END SUBROUTINE RESTORE_LAST_CELL


! --------------------------------------------------------------------------------------------------------------
!> \brief Determine if cell should be considered during packing of zone numbers
! --------------------------------------------------------------------------------------------------------------
LOGICAL FUNCTION FORBIDDEN_ZONE(SEND_BUFFER_INT, IZ, ICG1, ICG2)
INTEGER, DIMENSION(:), INTENT(IN) :: SEND_BUFFER_INT
INTEGER, INTENT(IN) :: IZ, ICG1, ICG2
INTEGER :: LL, ICG

FORBIDDEN_ZONE = .FALSE.
LL = 5
DO ICG = ICG1, ICG2
   IF (SEND_BUFFER_INT(LL) == IZ) THEN
      FORBIDDEN_ZONE = .TRUE.
      RETURN
   ENDIF
   LL = LL + 4
ENDDO
END FUNCTION FORBIDDEN_ZONE


! --------------------------------------------------------------------------------------------------------------
!> \brief Control multigrid cycling (F/V/W)
! Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
! --------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION CYCLING_CONTROL(NTYPE, NL)
USE SCARC_POINTERS, ONLY: MG
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, IL, ICYCLE

SELECT CASE (NTYPE)

   ! Initialize cycle counts at beginning of multigrid method
 
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
   
   ! Reset cycle counts at beginning of each new multigrid iteration
 
   CASE (NSCARC_CYCLING_RESET)
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         DO IL = NLEVEL_MIN, NLEVEL_MAX
            MG => SCARC(NM)%LEVEL(IL)%MG
            MG%CYCLING(1)=MG%CYCLING(2)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLING_NEXT
   
   ! Determine where to proceed with cycling
 
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

CYCLING_CONTROL = ICYCLE
RETURN
END FUNCTION CYCLING_CONTROL

END MODULE SCARC_UTILITIES


!=======================================================================================================================
!
! MODULE SCARC_STORAGE
!
!> \brief Organize the allocation, deallocation and resizing of different data structures
!
!   This includes 1-, 2- or 3-dimensional vectors of different types 
!   and compactly or bandwise stored matrices
!
!=======================================================================================================================
MODULE SCARC_STORAGE

USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE SCARC_CONSTANTS
USE SCARC_TYPES
USE SCARC_VARIABLES
USE SCARC_MESSAGES
USE SCARC_TROUBLESHOOTING
USE SCARC_UTILITIES

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief Setup memory management
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STORAGE
ALLOCATE (STORAGE%ALLOCATION_LIST(NSCARC_STORAGE_MAX), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_MEMORY_MANAGMENT', 'ALLOCATION_LIST', IERROR)
END SUBROUTINE SCARC_SETUP_STORAGE

! --------------------------------------------------------------------------------------------------------------
!> \brief Update list of arrays within ScaRC memory management
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_STORAGE(NDATA, NSTATE, NDIM, NINIT, NL1, NR1, NL2, NR2, NL3, NR3, CID, CSCOPE)
USE SCARC_POINTERS, ONLY : AL
INTEGER, INTENT(IN) :: NDATA, NSTATE, NDIM, NINIT, NL1, NR1, NL2, NR2, NL3, NR3
INTEGER :: NWORK, NLEN(3), I, IP
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
CHARACTER(20) :: CTYPE, CSTATE, CINIT, CDIM

STORAGE%IP = STORAGE%IP + 1

! Extract basic name of allocated structure and name of calling routine

! Get size of requested structure

IF (NSTATE /= NSCARC_STORAGE_REMOVE) THEN

   AL => STORAGE%ALLOCATION_LIST(STORAGE%IP)

   AL%CNAME  = CID
   AL%CSCOPE = CSCOPE

   AL%LBND = NSCARC_INIT_ZERO
   AL%RBND = NSCARC_INIT_ZERO

   IF (NL1 > 0) AL%LBND(1) = NL1
   IF (NR1 > 0) AL%RBND(1) = NR1
   IF (NL2 > 0) AL%LBND(2) = NL2
   IF (NR2 > 0) AL%RBND(2) = NR2
   IF (NL3 > 0) AL%LBND(3) = NL3
   IF (NR3 > 0) AL%RBND(3) = NR3

ELSE

   DO IP = 1, STORAGE%N_ARRAYS
      AL => STORAGE%ALLOCATION_LIST(IP)
      IF (TRIM(CID) == AL%CNAME) EXIT
   ENDDO

ENDIF

NWORK = 1
DO I = 1, 3
   NLEN(I) = AL%RBND(I) - AL%LBND(I) + 1
   IF (NLEN(I) /= 0) NWORK = NWORK * NLEN(I)
ENDDO

! Get some full text information for requested structure to dump out in memory file

CTYPE = SET_DATA_TYPE(NDATA)
CDIM  = SET_DIMENSION(NDIM)
CINIT = SET_INIT_TYPE(NDATA, NINIT, NSTATE)

SELECT CASE (NSTATE)
   CASE (NSCARC_STORAGE_CREATE)
      CSTATE = 'CREATE'
      CALL SCARC_UPDATE_STORAGE_COUNTERS (NDATA, NWORK,  1)
   CASE (NSCARC_STORAGE_RESIZE)
      CSTATE = 'RESIZE'
      CALL SCARC_UPDATE_STORAGE_COUNTERS (NDATA, NWORK,  0)
   CASE (NSCARC_STORAGE_REMOVE)
      CSTATE = 'REMOVE'
      CALL SCARC_UPDATE_STORAGE_COUNTERS (NDATA, NWORK, -1)
      NWORK = -NWORK
   CASE DEFAULT
      CSTATE = ' '
END SELECT

END SUBROUTINE SCARC_UPDATE_STORAGE

! --------------------------------------------------------------------------------------------------------------
!> \brief Update memory statistics w.r.t to occupied workspace and number of allocated arrays
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_STORAGE_COUNTERS (NDATA, NWORK, NSCAL)
INTEGER, INTENT(IN) :: NDATA, NWORK, NSCAL

STORAGE%N_ARRAYS = STORAGE%N_ARRAYS + NSCAL
IF (STORAGE%N_ARRAYS > NSCARC_STORAGE_MAX) WRITE(*,*) 'ERROR in APPEND_TO_ALLOCATION_LIST: list of arrays exceeded!'

SELECT CASE (NDATA) 
   CASE (NSCARC_DATA_INTEGER)
      STORAGE%N_INT = STORAGE%N_INT + NSCAL 
      STORAGE%NWORK_INT = STORAGE%NWORK_INT + NSCAL * NWORK
   CASE (NSCARC_DATA_REAL_EB)
      STORAGE%N_INT = STORAGE%N_REAL_EB + NSCAL 
      STORAGE%NWORK_REAL_EB = STORAGE%NWORK_REAL_EB + NSCAL * NWORK
   CASE (NSCARC_DATA_REAL_FB)
      STORAGE%N_INT = STORAGE%N_REAL_FB + NSCAL 
      STORAGE%NWORK_REAL_FB = STORAGE%NWORK_REAL_FB + NSCAL * NWORK
   CASE (NSCARC_DATA_LOGICAL)
      STORAGE%N_INT = STORAGE%N_LOG + NSCAL 
      STORAGE%NWORK_LOG = STORAGE%NWORK_LOG + NSCAL * NWORK
   CASE (NSCARC_DATA_CMATRIX)
      STORAGE%N_CMATRIX = STORAGE%N_CMATRIX + NSCAL 
   CASE (NSCARC_DATA_BMATRIX)
      STORAGE%N_BMATRIX = STORAGE%N_BMATRIX + NSCAL 
END SELECT

END SUBROUTINE SCARC_UPDATE_STORAGE_COUNTERS

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize integer array of dimension 1
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT1 (WORKSPACE, NL1, NR1, NINIT, CID, CSCOPE)
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT1', CID, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
      CASE (NSCARC_INIT_HUGE)
         WORKSPACE = NSCARC_HUGE_INT
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_INTEGER, NSCARC_STORAGE_CREATE, 1, NINIT, NL1, NR1, -1, -1, -1, -1, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_INT1

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize integer array of dimension 2
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT2 (WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CID, CSCOPE)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT2', CID, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
      CASE (NSCARC_INIT_HUGE)
         WORKSPACE = NSCARC_HUGE_INT
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
   ENDIF
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_INTEGER, NSCARC_STORAGE_CREATE, 2, NINIT, NL1, NR1, NL2, NR2, -1, -1, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_INT2

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize integer array of dimension 3
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT3 (WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CID, CSCOPE)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT3', CID, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
      CASE (NSCARC_INIT_HUGE)
         WORKSPACE = NSCARC_HUGE_INT
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
      SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
   ENDIF
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_INTEGER, NSCARC_STORAGE_CREATE, 3, NINIT, NL1, NR1, NL2, NR2, NL3, NR3, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_INT3

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize Logical array of dimension 1
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_LOG1 (WORKSPACE, NL1, NR1, NINIT, CID, CSCOPE)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT1', CID, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_TRUE)
         WORKSPACE = .TRUE.
      CASE (NSCARC_INIT_FALSE)
         WORKSPACE = .FALSE.
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) THEN
      CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
   ENDIF
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_LOGICAL, NSCARC_STORAGE_CREATE, 1, NINIT, NL1, NR1, -1, -1, -1, -1, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_LOG1

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize Logical array of dimension 2
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_LOG2 (WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CID, CSCOPE)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT3', CID, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_TRUE)
         WORKSPACE = .TRUE.
      CASE (NSCARC_INIT_FALSE)
         WORKSPACE = .FALSE.
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
   ENDIF
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_LOGICAL, NSCARC_STORAGE_CREATE, 2, NINIT, NL1, NR1, NL2, NR2, -1, -1, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_LOG2

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize Logical array of dimension 3
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_LOG3 (WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CID, CSCOPE)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT3', CID, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_TRUE)
         WORKSPACE = .TRUE.
      CASE (NSCARC_INIT_FALSE)
         WORKSPACE = .FALSE.
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
      SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
   ENDIF
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_LOGICAL, NSCARC_STORAGE_CREATE, 3, NINIT, NL1, NR1, NL2, NR2, NL3, NR3, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_LOG3

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize real array of dimension 1
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL1 (WORKSPACE, NL1, NR1, NINIT, CID, CSCOPE)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN

   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL1', CID, IERROR)

   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
      CASE (NSCARC_INIT_HUGE)
         WORKSPACE = NSCARC_HUGE_REAL_EB
   END SELECT

ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_EB, NSCARC_STORAGE_CREATE, 1, NINIT, NL1, NR1, -1, -1, -1, -1, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_REAL1

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize real array of dimension 1
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL1_FB (WORKSPACE, NL1, NR1, NINIT, CID, CSCOPE)
USE PRECISION_PARAMETERS, ONLY: FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(FB), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL1_FB', CID, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_FB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_FB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_FB
      CASE (NSCARC_INIT_HUGE)
         WORKSPACE = NSCARC_HUGE_REAL_FB
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_FB, NSCARC_STORAGE_CREATE, 1, NINIT, NL1, NR1, -1, -1, -1, -1, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_REAL1_FB

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize real array of dimension 2
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL2 (WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CID, CSCOPE)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL2', CID, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
      CASE (NSCARC_INIT_HUGE)
         WORKSPACE = NSCARC_HUGE_REAL_EB
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
   ENDIF
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_EB, NSCARC_STORAGE_CREATE, 2, NINIT, NL1, NR1, NL2, NR2, -1, -1, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_REAL2

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize real array of dimension 3
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL3 (WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CID, CSCOPE)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL3', CID, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
      CASE (NSCARC_INIT_HUGE)
         WORKSPACE = NSCARC_HUGE_REAL_EB
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
      SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_ERROR(NSCARC_ERROR_VECTOR_LENGTH, CID, NSCARC_NONE)
   ENDIF
ENDIF

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_EB, NSCARC_STORAGE_CREATE, 3, NINIT, NL1, NR1, NL2, NR2, NL3, NR3, CID, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_REAL3

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate one-dimensional integer vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_INT1 (WORKSPACE, CID, CSCOPE)
INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_INTEGER, NSCARC_STORAGE_REMOVE, 1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_INT1

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate two-dimensional integer vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_INT2 (WORKSPACE, CID, CSCOPE)
INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_INTEGER, NSCARC_STORAGE_REMOVE, 2, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_INT2

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate three-dimensional integer vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_INT3 (WORKSPACE, CID, CSCOPE)
INTEGER, ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_INTEGER, NSCARC_STORAGE_REMOVE, 3, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_INT3

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate one-dimensional logical vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_LOG1 (WORKSPACE, CID, CSCOPE)
LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_LOGICAL, NSCARC_STORAGE_REMOVE, 1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_LOG1

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate two-dimensional logical vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_LOG2 (WORKSPACE, CID, CSCOPE)
LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_LOGICAL, NSCARC_STORAGE_REMOVE, 2, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_LOG2

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate three-dimensional logical vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_LOG3 (WORKSPACE, CID, CSCOPE)
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_LOGICAL, NSCARC_STORAGE_REMOVE, 3, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_LOG3

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate one-dimensional double precision vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_REAL1 (WORKSPACE, CID, CSCOPE)
REAL(EB), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_EB, NSCARC_STORAGE_REMOVE, 1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_REAL1

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate two-dimensional double precision vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_REAL2 (WORKSPACE, CID, CSCOPE)
REAL(EB), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_EB, NSCARC_STORAGE_REMOVE, 2, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_REAL2

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate three-dimensional double precision vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_REAL3 (WORKSPACE, CID, CSCOPE)
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_EB, NSCARC_STORAGE_REMOVE, 3, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_REAL3

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate one-dimensional single precision vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_REAL1_FB (WORKSPACE, CID, CSCOPE)
REAL(FB), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_EB, NSCARC_STORAGE_REMOVE, 1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_REAL1_FB

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate one-dimensional single precision vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_REAL2_FB (WORKSPACE, CID, CSCOPE)
REAL(FB), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_EB, NSCARC_STORAGE_REMOVE, 2, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_REAL2_FB

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate one-dimensional single precision vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_REAL3_FB (WORKSPACE, CID, CSCOPE)
REAL(FB), ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) :: WORKSPACE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
DEALLOCATE(WORKSPACE) 
CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_REAL_EB, NSCARC_STORAGE_REMOVE, 3, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)
END SUBROUTINE SCARC_DEALLOCATE_REAL3_FB

! --------------------------------------------------------------------------------------------------------------
!> \brief Resize one-dimensional integer vector to requested bounds
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESIZE_INT1 (WORKSPACE, NL1, NR1, CID, CSCOPE)
INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
INTEGER, ALLOCATABLE, DIMENSION(:) :: AUX
INTEGER :: NLC1, NRC1, NSC, NS

NLC1 = LBOUND(WORKSPACE, DIM=1)
NRC1 = UBOUND(WORKSPACE, DIM=1)

IF (NL1 == NLC1 .AND. NR1 == NRC1) THEN
   RETURN
ELSE

   NSC = NRC1 - NLC1 + 1
   NS  = NR1  - NL1  + 1

   ALLOCATE(AUX(1:NSC), STAT = IERROR)                        ! don't track it in memory management, only auxiliary
   AUX(1:NSC) = WORKSPACE(NLC1:NRC1)
   CALL SCARC_DEALLOCATE_INT1 (WORKSPACE, CID, CSCOPE)

   CALL SCARC_ALLOCATE_INT1 (WORKSPACE, NL1, NR1, NSCARC_INIT_NONE, CID, CSCOPE)

   IF (NS < NSC) THEN
      WORKSPACE(NL1:NL1 + NS) = AUX(1:NS)
   ELSE
      WORKSPACE(NL1:NL1 + NSC) = AUX(1:NSC)
   ENDIF
   DEALLOCATE(AUX)

ENDIF

END SUBROUTINE SCARC_RESIZE_INT1

! --------------------------------------------------------------------------------------------------------------
!> \brief Resize two-dimensional integer vector to requested bounds
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESIZE_INT2 (WORKSPACE, NL1, NR1, NL2, NR2, CID, CSCOPE)
INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: AUX
INTEGER :: NLC1, NRC1, NSC1, NS1
INTEGER :: NLC2, NRC2, NSC2, NS2

NLC1 = LBOUND(WORKSPACE, DIM=1)                ! current left  bound for dimension 1
NRC1 = UBOUND(WORKSPACE, DIM=1)                ! current right bound for dimension 1

NLC2 = LBOUND(WORKSPACE, DIM=2)                ! current left  bound for dimension 2
NRC2 = UBOUND(WORKSPACE, DIM=2)                ! current right bound for dimension 2

IF (NL1 == NLC1 .AND. NR1 == NRC1 .AND. &
    NL2 == NLC2 .AND. NR2 == NRC2) THEN
   RETURN
ELSE

   NSC1 = NRC1 - NLC1 + 1
   NSC2 = NRC2 - NLC2 + 1

   NS1 = NR1 - NL1 + 1
   NS2 = NR2 - NL2 + 1

   ALLOCATE(AUX(1: NSC1, 1: NSC2), STAT = IERROR)                     ! don't track it in memory management, only auxiliary
   AUX(1:NSC1, 1:NSC2) = WORKSPACE(NLC1:NRC1, NLC2:NRC2)
   CALL SCARC_DEALLOCATE_INT2 (WORKSPACE, CID, CSCOPE)

   CALL SCARC_ALLOCATE_INT2 (WORKSPACE, NL1, NR1, NL2, NR2, NSCARC_INIT_NONE, CID, CSCOPE)

   IF (NS1 < NSC1 .AND. NS2 < NSC2) THEN
      WORKSPACE(NL1:NL1+NS1, NL2:NL2+NS2) = AUX(1:NS1, 1:NS2)
   ELSE IF (NS1 < NSC1 .AND. NS2 > NSC2) THEN
      WORKSPACE(NL1:NL1+NS1, NL2:NL2+NSC2) = AUX(1:NS1, 1:NSC2)
   ELSE IF (NS1 > NSC1 .AND. NS2 < NSC2) THEN
      WORKSPACE(NL1:NL1+NSC1, NL2:NL2+NS2) = AUX(1:NSC1, 1:NS2)
   ELSE
      WORKSPACE(NL1:NL1+NSC1, NL2:NL2+NSC2) = AUX(1:NSC1, 1:NSC2)
   ENDIF

   DEALLOCATE(AUX)

ENDIF

END SUBROUTINE SCARC_RESIZE_INT2

! --------------------------------------------------------------------------------------------------------------
!> \brief Reduce size of integer vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_REDUCE_INT1 (WORKSPACE, NSIZE, CID, CSCOPE)
INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NSIZE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
INTEGER, ALLOCATABLE, DIMENSION(:) :: AUX

IF (NSIZE == SIZE(WORKSPACE)) THEN
   RETURN                                  ! vector already has requested size
ELSE IF (NSIZE < SIZE(WORKSPACE)) THEN

   ALLOCATE(AUX(1: NSIZE), STAT = IERROR)
   AUX(1:NSIZE) = WORKSPACE(1:NSIZE)
   CALL SCARC_DEALLOCATE_INT1 (WORKSPACE, CID, CSCOPE)

   CALL SCARC_ALLOCATE_INT1 (WORKSPACE, 1, NSIZE, NSCARC_INIT_NONE, CID, CSCOPE)
   WORKSPACE(1:NSIZE) = AUX(1:NSIZE)
   DEALLOCATE(AUX)
ENDIF

END SUBROUTINE SCARC_REDUCE_INT1

! --------------------------------------------------------------------------------------------------------------
!> \brief Reduce size of integer array with dimension 2
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_REDUCE_INT2 (WORKSPACE, NSIZE1, NSIZE2, CID, CSCOPE)
INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NSIZE1, NSIZE2
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
INTEGER :: NWORK1, NWORK2, I1, I2
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: AUX

NWORK1 = SIZE(WORKSPACE, DIM = 1)
NWORK2 = SIZE(WORKSPACE, DIM = 2)
IF (NSIZE1 == NWORK1 .AND. NSIZE2 == NWORK2 ) THEN
   RETURN                                  ! vector already has requested size
ELSE IF (NSIZE1 <= NWORK1 .AND. NSIZE2 <= NWORK2) THEN

   ALLOCATE(AUX(1: NWORK1, 1: NWORK2), STAT = IERROR) 
   DO I2 = 1, NWORK2
      DO I1 = 1, NWORK1
         AUX(I1, I2) = WORKSPACE(I1, I2)
      ENDDO
   ENDDO
   CALL SCARC_DEALLOCATE_INT2 (WORKSPACE, CID, CSCOPE)

   CALL SCARC_ALLOCATE_INT2 (WORKSPACE, 1, NSIZE1, 1, NSIZE2, NSCARC_INIT_NONE, TRIM(CID), CSCOPE)
   DO I2 = 1, NSIZE2
      DO I1 = 1, NSIZE1
         WORKSPACE(I1, I2) = AUX(I1, I2)
      ENDDO
   ENDDO
   DEALLOCATE(AUX)
ELSE
   WRITE(*,'(A,2I6,A,2I6)') 'Error in SCARC_REDUCE_INT2, orig ', NWORK1, NWORK2,': new ',NSIZE1, NSIZE2
ENDIF

END SUBROUTINE SCARC_REDUCE_INT2

! --------------------------------------------------------------------------------------------------------------
!> \brief Expand size of integer vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXPAND_INT1 (WORKSPACE, WORKSPACE_ADD, NSIZE, NSIZE_ADD, CID, CSCOPE)
INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE, WORKSPACE_ADD
INTEGER, INTENT(IN) :: NSIZE, NSIZE_ADD
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
INTEGER, ALLOCATABLE, DIMENSION(:) :: AUX

ALLOCATE(AUX(1: NSIZE + NSIZE_ADD), STAT = IERROR)
AUX(1:NSIZE) = WORKSPACE(1:NSIZE)
AUX(NSIZE+1:NSIZE+NSIZE_ADD) = WORKSPACE_ADD(NSIZE+1:NSIZE+NSIZE_ADD)
CALL SCARC_DEALLOCATE_INT1 (WORKSPACE, CID, CSCOPE)

CALL SCARC_ALLOCATE_INT1 (WORKSPACE, 1, NSIZE + NSIZE_ADD, NSCARC_INIT_NONE, TRIM(CID), CSCOPE)
WORKSPACE(1:NSIZE+NSIZE_ADD) = AUX(1:NSIZE+NSIZE_ADD)
DEALLOCATE(AUX)

END SUBROUTINE SCARC_EXPAND_INT1

! --------------------------------------------------------------------------------------------------------------
!> \brief Reduce size of integer vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_REDUCE_REAL1 (WORKSPACE, NSIZE, CID, CSCOPE)
REAL(EB), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NSIZE
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
REAL(EB), ALLOCATABLE, DIMENSION(:) :: AUX

IF (NSIZE == SIZE(WORKSPACE)) THEN
   RETURN                                  ! vector already has requested size
ELSE IF (NSIZE < SIZE(WORKSPACE)) THEN

   ALLOCATE(AUX(1: NSIZE), STAT = IERROR)
   AUX(1:NSIZE) = WORKSPACE(1:NSIZE)
   CALL SCARC_DEALLOCATE_REAL1 (WORKSPACE, CID, CSCOPE)

   CALL SCARC_ALLOCATE_REAL1 (WORKSPACE, 1, NSIZE, NSCARC_INIT_NONE, TRIM(CID), CSCOPE)
   WORKSPACE(1:NSIZE) = AUX(1:NSIZE)
   DEALLOCATE(AUX)
ENDIF

END SUBROUTINE SCARC_REDUCE_REAL1

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate matrix in compact storage format
! Allocate matrix with corresponding pointer and length structures
!    NTYPE == NSCARC_MATRIX_FULL    :  ALLOCATE VAL, COL and COLG
!    NTYPE == NSCARC_MATRIX_LIGHT   :  ALLOCATE VAL, COL 
!    NTYPE == NSCARC_MATRIX_MINIMAL :  ALLOCATE COL 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_CMATRIX (A, NL, NPREC, NTYPE, CID, CSCOPE)
TYPE (SCARC_CMATRIX_TYPE), INTENT(INOUT) :: A
INTEGER, INTENT(IN) :: NPREC, NTYPE, NL
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
INTEGER :: NDUMMY

A%CNAME = TRIM(CID)
A%NTYPE = NTYPE
A%NPREC = NPREC
NDUMMY = NL

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_CMATRIX, NSCARC_STORAGE_CREATE, -1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)

CALL SCARC_ALLOCATE_INT1 (A%ROW, 1, A%N_ROW, NSCARC_INIT_ZERO, 'A%ROW', CSCOPE)
CALL SCARC_ALLOCATE_INT1 (A%COL, 1, A%N_VAL, NSCARC_INIT_ZERO, 'A%COL', CSCOPE)

IF (TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL .AND. (NTYPE == NSCARC_MATRIX_LIGHT .OR. NTYPE == NSCARC_MATRIX_FULL)) THEN
   CALL SCARC_ALLOCATE_INT1 (A%COLG, 1, A%N_VAL, NSCARC_INIT_ZERO, 'A%COLG', CSCOPE)
ENDIF

IF (NTYPE /= NSCARC_MATRIX_MINIMAL) THEN
   IF (NPREC == NSCARC_PRECISION_SINGLE) THEN
      CALL SCARC_ALLOCATE_REAL1_FB (A%VAL_FB, 1, A%N_VAL, NSCARC_INIT_ZERO, 'A%VAL_FB', CSCOPE)
   ELSE
      CALL SCARC_ALLOCATE_REAL1 (A%VAL, 1, A%N_VAL, NSCARC_INIT_ZERO, 'A%VAL', CSCOPE)
   ENDIF
ENDIF

END SUBROUTINE SCARC_ALLOCATE_CMATRIX

! --------------------------------------------------------------------------------------------------------------
!> \brief Dellocate matrix in compact storage format
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_CMATRIX (A, CID, CSCOPE)
TYPE (SCARC_CMATRIX_TYPE), INTENT(INOUT) :: A
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

A%N_STENCIL   = 0
A%N_CONDENSED = 0
A%N_ROW       = 0
A%N_VAL       = 0
A%NTYPE       = 0
A%NPREC       = 0
A%STENCIL     = 0
A%POS         = 0

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_CMATRIX, NSCARC_STORAGE_REMOVE, -1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)

IF (ALLOCATED(A%VAL))   CALL SCARC_DEALLOCATE_REAL1 (A%VAL,   'A%VAL',   CSCOPE)
IF (ALLOCATED(A%ROW))   CALL SCARC_DEALLOCATE_INT1  (A%ROW,   'A%ROW',   CSCOPE)
IF (ALLOCATED(A%COL))   CALL SCARC_DEALLOCATE_INT1  (A%COL,   'A%COL',   CSCOPE)
IF (ALLOCATED(A%COLG))  CALL SCARC_DEALLOCATE_INT1  (A%COLG,  'A%COLG' , CSCOPE)
IF (ALLOCATED(A%RELAX)) CALL SCARC_DEALLOCATE_REAL1 (A%RELAX, 'A%RELAX', CSCOPE)

END SUBROUTINE SCARC_DEALLOCATE_CMATRIX

! --------------------------------------------------------------------------------------------------------------
!> \brief Insert value at specified position in matrix of compact storage format
! --------------------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_EVALUATE_CMATRIX(A, IC, JC) 
TYPE (SCARC_CMATRIX_TYPE), INTENT(IN) :: A
INTEGER, INTENT(IN) :: IC, JC
INTEGER :: IP

SCARC_EVALUATE_CMATRIX =  0.0_EB
DO IP = A%ROW(IC), A%ROW(IC+1)-1
   IF (A%COL(IP) == JC) THEN
      SCARC_EVALUATE_CMATRIX = A%VAL(IP)
      EXIT
   ENDIF
ENDDO

END FUNCTION SCARC_EVALUATE_CMATRIX

! --------------------------------------------------------------------------------------------------------------
!> \brief Insert value at specified position in matrix of compact storage format
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INSERT_TO_CMATRIX (A, VAL, IC, JC, NC, NP, CID) 
TYPE (SCARC_CMATRIX_TYPE), INTENT(INOUT) :: A
INTEGER, INTENT(IN) :: IC, JC, NC
INTEGER, INTENT(INOUT) :: NP
CHARACTER(*), INTENT(IN) :: CID
REAL(EB), INTENT(IN) :: VAL
INTEGER :: IP
CHARACTER(1) :: CSAVE
REAL(EB) :: TOL = 1.0E-14_EB

CSAVE = CID(1:1)              ! dummy command to prevent warning in case that DEBUG flag is not set

IF (NP == A%N_VAL) WRITE(*,*) MY_RANK+1,': SCARC_INSERT_TO_CMATRIX: Error, maximum length already reached'
IF (ABS(VAL) < TOL) RETURN

ALREADY_EXISTING_LOOP: DO IP = A%ROW(IC), A%ROW(IC+1)-1
   IF (JC /= IC .AND. JC == A%COL(IP)) THEN
      WRITE(*,*) MY_RANK+1,': SCARC_INSERT_TO_CMATRIX: Index ', JC,' already exists'
      EXIT
   ENDIF
ENDDO ALREADY_EXISTING_LOOP

IF (JC == IC) THEN
   IP = A%ROW(IC)
   A%VAL(IP) = VAL                               ! COL and ROW already correct
   A%COL(IP) = JC
ELSE
   IP = A%ROW(IC+1)
   A%VAL(IP+1:NP+1) = A%VAL(IP:NP)               ! COL and ROW must be shifted
   A%COL(IP+1:NP+1) = A%COL(IP:NP)
   A%ROW(IC+1:NC+1) = A%ROW(IC+1:NC+1) + 1
   A%VAL(IP) = VAL           
   A%COL(IP) = JC
   NP = NP + 1
ENDIF

END SUBROUTINE SCARC_INSERT_TO_CMATRIX

! --------------------------------------------------------------------------------------------------------------
!> \brief Reduce size of matrix in compact storage format
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_COPY_CMATRIX (A1, A2, CID, CSCOPE)
TYPE (SCARC_CMATRIX_TYPE), INTENT(INOUT) :: A1, A2
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

A2%CNAME       = TRIM(CID)
A2%N_STENCIL   = A1%N_STENCIL
A2%N_CONDENSED = A1%N_CONDENSED
A2%N_ROW       = A1%N_ROW
A2%N_VAL       = A1%N_VAL
A2%NTYPE       = A1%NTYPE
A2%NPREC       = A1%NPREC
A2%STENCIL     = A1%STENCIL
A2%POS         = A1%POS

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_CMATRIX, NSCARC_STORAGE_CREATE, -1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)

CALL SCARC_ALLOCATE_INT1 (A2%ROW, 1, A2%N_ROW, NSCARC_INIT_NONE, 'A2%ROW', CSCOPE)
A2%ROW = A1%ROW

CALL SCARC_ALLOCATE_INT1 (A2%COL, 1, A2%N_VAL, NSCARC_INIT_NONE, 'A2%COL', CSCOPE)
A2%COL = A1%COL

CALL SCARC_ALLOCATE_REAL1 (A2%VAL, 1, A2%N_VAL, NSCARC_INIT_ZERO, 'A2%VAL', CSCOPE)
A2%VAL = A1%VAL

END SUBROUTINE SCARC_COPY_CMATRIX

! --------------------------------------------------------------------------------------------------------------
!> \brief Reduce size of matrix in compact storage format
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_REDUCE_CMATRIX (A, CID, CSCOPE)
TYPE (SCARC_CMATRIX_TYPE), INTENT(INOUT) :: A
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
REAL(EB), ALLOCATABLE, DIMENSION(:) :: VAL
REAL(FB), ALLOCATABLE, DIMENSION(:) :: VAL_FB
INTEGER , ALLOCATABLE, DIMENSION(:) :: COL, COLG
INTEGER :: NVAL_CURRENT, NVAL_ALLOCATED

NVAL_CURRENT = A%ROW(A%N_ROW)
NVAL_ALLOCATED = SIZE(A%COL)

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_CMATRIX, NSCARC_STORAGE_RESIZE, -1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)

! If the matrix already has the desired size or specified values are to small, return or shutdown
IF (NVAL_ALLOCATED == NVAL_CURRENT) THEN
   RETURN
ELSE IF (NVAL_ALLOCATED < NVAL_CURRENT) THEN
   CALL SCARC_ERROR(NSCARC_ERROR_MATRIX_SIZE, CID, NSCARC_NONE)
ENDIF

! If the allocated size of the matrix values workspace is too large, reduce it to the real size
IF (NVAL_CURRENT < SIZE(A%COL)) THEN

   ALLOCATE(COL(1: NVAL_CURRENT), STAT = IERROR)
   COL(1:NVAL_CURRENT) = A%COL(1:NVAL_CURRENT)
   CALL SCARC_DEALLOCATE_INT1 (A%COL, 'A%COL', CSCOPE)
   CALL SCARC_ALLOCATE_INT1 (A%COL, 1, NVAL_CURRENT, NSCARC_INIT_NONE, 'A%COL', CSCOPE)
   A%COL(1:NVAL_CURRENT) = COL(1:NVAL_CURRENT)
   DEALLOCATE(COL)

   IF (ALLOCATED(A%COLG)) THEN
      ALLOCATE(COLG(1: NVAL_CURRENT), STAT = IERROR)
      COLG(1:NVAL_CURRENT) = A%COLG(1:NVAL_CURRENT)
      CALL SCARC_DEALLOCATE_INT1 (A%COLG, 'A%COLG',CSCOPE)
      CALL SCARC_ALLOCATE_INT1 (A%COLG, 1, NVAL_CURRENT, NSCARC_INIT_NONE, 'A%COLG', CSCOPE)
      A%COLG(1:NVAL_CURRENT) = COLG(1:NVAL_CURRENT)
      DEALLOCATE(COLG)
   ENDIF

   IF (A%NTYPE /= NSCARC_MATRIX_MINIMAL) THEN
      SELECT CASE (A%NPREC)
         CASE (NSCARC_PRECISION_SINGLE)
            ALLOCATE(VAL_FB(1: NVAL_CURRENT), STAT = IERROR)
            VAL_FB(1:NVAL_CURRENT) = A%VAL_FB(1:NVAL_CURRENT)
            CALL SCARC_DEALLOCATE_REAL1_FB (A%VAL_FB, 'A%VAL_FB', CSCOPE)
            CALL SCARC_ALLOCATE_REAL1_FB (A%VAL_FB, 1, NVAL_CURRENT, NSCARC_INIT_NONE, 'A%VAL_FB', CSCOPE)
            A%VAL_FB(1:NVAL_CURRENT) = VAL_FB(1:NVAL_CURRENT)
            DEALLOCATE(VAL_FB)
         CASE (NSCARC_PRECISION_DOUBLE)
            ALLOCATE(VAL(1: NVAL_CURRENT), STAT = IERROR)
            VAL(1:NVAL_CURRENT) = A%VAL(1:NVAL_CURRENT)
            CALL SCARC_DEALLOCATE_REAL1 (A%VAL, 'A%VAL', CSCOPE)
            CALL SCARC_ALLOCATE_REAL1 (A%VAL, 1, NVAL_CURRENT, NSCARC_INIT_NONE, 'A%VAL', CSCOPE)
            A%VAL(1:NVAL_CURRENT) = VAL(1:NVAL_CURRENT)
            DEALLOCATE(VAL)
         END SELECT
   ENDIF
   A%N_VAL = NVAL_CURRENT

ENDIF

END SUBROUTINE SCARC_REDUCE_CMATRIX

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate matrix in bandwise storage format
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_BMATRIX (A, NL, CID, CSCOPE)
TYPE (SCARC_BMATRIX_TYPE), INTENT(INOUT) :: A
CHARACTER(*), INTENT(IN) :: CID, CSCOPE
INTEGER, INTENT(IN) :: NL
CHARACTER(40) :: CINFO

A%CNAME = CID

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_BMATRIX, NSCARC_STORAGE_CREATE, -1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)

WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CID),'_LEV',NL,'.AUX'
CALL SCARC_ALLOCATE_REAL1 (A%AUX, 1, A%N_DIAG, NSCARC_INIT_ZERO, CINFO, CSCOPE)
WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CID),'_LEV',NL,'.VAL'
CALL SCARC_ALLOCATE_REAL2 (A%VAL, 1, A%N_DIAG, 1, A%N_STENCIL, NSCARC_INIT_ZERO, CINFO, CSCOPE)

END SUBROUTINE SCARC_ALLOCATE_BMATRIX

! --------------------------------------------------------------------------------------------------------------
!> \brief Deallocate matrix in bandwise storage format
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_BMATRIX (A, CID, CSCOPE)
TYPE (SCARC_BMATRIX_TYPE), INTENT(INOUT) :: A
CHARACTER(*), INTENT(IN) :: CID, CSCOPE

A%N_STENCIL   = 0
A%N_CONDENSED = 0
A%N_VAL       = 0
A%N_DIAG      = 0

CALL SCARC_UPDATE_STORAGE(NSCARC_DATA_BMATRIX, NSCARC_STORAGE_REMOVE, -1, -1, -1, -1, -1, -1, -1, -1, CID, CSCOPE)

IF (ALLOCATED(A%AUX))    CALL SCARC_DEALLOCATE_REAL1 (A%AUX, 'A%AUX', CSCOPE)
IF (ALLOCATED(A%VAL))    CALL SCARC_DEALLOCATE_REAL2 (A%VAL, 'A%VAL', CSCOPE)
IF (ALLOCATED(A%RELAX))  CALL SCARC_DEALLOCATE_REAL2 (A%RELAX, 'A%RELAX', CSCOPE)
IF (ALLOCATED(A%RELAXD)) CALL SCARC_DEALLOCATE_REAL1 (A%RELAXD, 'A%RELAXD', CSCOPE)
A%STENCIL = 0

A%OFFSET = 0
A%LENGTH = 0
A%SOURCE = 0
A%TARGET = 0

END SUBROUTINE SCARC_DEALLOCATE_BMATRIX

END MODULE SCARC_STORAGE


!=======================================================================================================================
!
! MODULE SCARC_CONVERGENCE
!
!> \brief Manage iteration parameters of the currently used ScaRC/UscaRC solver
!
!=======================================================================================================================
MODULE SCARC_CONVERGENCE

USE PRECISION_PARAMETERS, ONLY: EB
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_MESSAGES

IMPLICIT NONE (TYPE,EXTERNAL)
  
REAL(EB) :: DT                                  !< TS width 
REAL(EB) :: DTI                                 !< Inverse of TS width 
REAL(EB) :: OMEGA                               !< Relaxation parameter for current solver
REAL(EB) :: EPS                                 !< Requested accuracy for current solver
REAL(EB) :: RES                                 !< Current residual of current solver
REAL(EB) :: RESIN = -1.0_EB                     !< Initial residual of current solver
REAL(EB) :: CAPPA = -1.0_EB                     !< Convergence rate of current solver

REAL(EB) :: VELOCITY_ERROR_MGM                  !< Requested threshold for velocity error in MGM method
REAL(EB) :: VELOCITY_ERROR_GLOBAL               !< Global velocity error achieved in MGM pass

INTEGER :: NIT        = 0                       !< Maximum number of iterations in current solver
INTEGER :: NIT_MGM    = 0                       !< Maximum number of iterations in MGM solver
INTEGER :: ITE        = 0                       !< Current number of iterations in current solver
INTEGER :: ITE_CG     = 0                       !< Statistical information about number of Krylov iterations
INTEGER :: ITE_MG     = 0                       !< Statistical information about number of multigrid iterations
INTEGER :: ITE_LU     = 0                       !< Statistical information about number of LU iterations
INTEGER :: ITE_PRES   = 0                       !< Statistical information about number of pressure iterations
INTEGER :: ITE_TOTAL  = 0                       !< Statistical information about number of total iterations
INTEGER :: ITE_SMOOTH = 0                       !< Statistical information about number of smoothing iterations
INTEGER :: ITE_COARSE = 0                       !< Statistical information about number of coarse grid iterations
INTEGER :: ITE_GLOBAL = 0                       !< Statistical information about number of global iterations

INTEGER  :: X                                   !< Handle for solution 1D-vector 
INTEGER  :: B                                   !< Handle for right hand side one dimensional-vector 
INTEGER  :: D                                   !< Handle for auxiliary one-dimensional vector
INTEGER  :: R                                   !< Handle for auxiliary one-dimensional vector
INTEGER  :: V                                   !< Handle for auxiliary one-dimensional vector
INTEGER  :: Y                                   !< Handle for auxiliary one-dimensional vector
INTEGER  :: Z                                   !< Handle for auxiliary one-dimensional vector


CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief Set current iteration state
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SET_ITERATION_STATE (DT_CURRENT)
REAL(EB), INTENT(IN) :: DT_CURRENT

DT  = DT_CURRENT
DTI = 1.0_EB/DT_CURRENT

ITE_PRES = ITE_PRES + 1
ITE_GLOBAL = ICYC

END SUBROUTINE SCARC_SET_ITERATION_STATE

! --------------------------------------------------------------------------------------------------------------
!> \brief Check if solver converges or diverges and print out residual information
! --------------------------------------------------------------------------------------------------------------
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
            NSTATE = NSCARC_STATE_CONV_INITIAL
         ELSE
            NSTATE = NSCARC_STATE_CONV
         ENDIF
         NIT = 0
      ENDIF
END SELECT
IF (RES > NSCARC_THRESHOLD_DIVGERGENCE) NSTATE = NSCARC_STATE_DIVG

SCARC_CONVERGENCE_STATE = NSTATE

IF (HAS_CSV_DUMP) CALL SCARC_DUMP_CSV(ISM, NS, NL)

END FUNCTION SCARC_CONVERGENCE_STATE

! --------------------------------------------------------------------------------------------------------------
!> \brief Compute convergence rate and print out residual information for final loop
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)
INTEGER, INTENT(IN) :: NSTATE, NS, NL

IF (NSTATE == NSCARC_STATE_DIVG) THEN
   ITE   = - 1
   CAPPA = 1.0_EB
ELSE
   IF (NSTATE == NSCARC_STATE_CONV_INITIAL) THEN
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
         IF (NSTATE == NSCARC_STATE_CONV_INITIAL) THEN
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

END SUBROUTINE SCARC_CONVERGENCE_RATE

! --------------------------------------------------------------------------------------------------------------
!> \brief Increase corresponding iteration count (just for visualization of convergence behavior)
! --------------------------------------------------------------------------------------------------------------
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

! ------------------------------------------------------------------------------------------------------------------
!> \brief Dump residual information
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_CSV(ISM, NS, NL)
INTEGER, INTENT(IN) :: ISM, NS, NL

IF (.NOT.HAS_CSV_DUMP .OR. MY_RANK /= 0) RETURN
IF (ITE_TOTAL == 0 .AND. TYPE_SOLVER /= NSCARC_SOLVER_MAIN) RETURN
IF (TYPE_SOLVER == NSCARC_SOLVER_COARSE) RETURN
WRITE(MSG%LU_STAT,1000) ITE_PRES, NS, ITE_TOTAL, ITE_CG, ITE_MG, NL, ITE_SMOOTH, ISM, ITE_COARSE, ITE_LU, RES, CAPPA

1000 FORMAT(10(I8,','), E14.6,',',E14.6)
END SUBROUTINE SCARC_DUMP_CSV

END MODULE SCARC_CONVERGENCE


!=======================================================================================================================
!
! MODULE SCARC_CPU
!
!> \brief Measure and dump CPU timings for different parts of the ScaRC/UScaRC solvers
!
!=======================================================================================================================
MODULE SCARC_CPU
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE SCARC_CONSTANTS
USE SCARC_MESSAGES

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief Setup time measurements
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CPU
ALLOCATE (CPU(0:N_MPI_PROCESSES-1), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_CPU', 'CPU', IERROR)
END SUBROUTINE SCARC_SETUP_CPU

! ------------------------------------------------------------------------------------------------------------------
!> \brief Dump CPU times of several routines
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_TIMERS
INTEGER, PARAMETER :: LINE_LENGTH = 5 + 12*11
TYPE (MPI_STATUS) :: STATUS
INTEGER :: N
CHARACTER(LEN=LINE_LENGTH) :: LINE
CHARACTER(LEN=LINE_LENGTH), DIMENSION(0:N_MPI_PROCESSES-1) :: LINE_ARRAY

! All MPI processes except root send their timings to the root process. The root process then writes them out to a file.
WRITE(LINE,'(I5,12(",",ES10.3))') MY_RANK,                       &
                                  CPU(MY_RANK)%OVERALL,          &
                                  CPU(MY_RANK)%SETUP,            &
                                  CPU(MY_RANK)%SOLVER,           &
                                  CPU(MY_RANK)%ITERATION,        &
                                  CPU(MY_RANK)%MATVEC_PRODUCT,   &
                                  CPU(MY_RANK)%SCALAR_PRODUCT,   &
                                  CPU(MY_RANK)%RELAXATION,       &
                                  CPU(MY_RANK)%SMOOTHER,         &
                                  CPU(MY_RANK)%COARSE,           &
                                  CPU(MY_RANK)%EXCHANGE,         &
                                  CPU(MY_RANK)%BUFFER_PACKING,   &
                                  CPU(MY_RANK)%BUFFER_UNPACKING

IF (MY_RANK>0) THEN
   CALL MPI_SEND(LINE,LINE_LENGTH,MPI_CHARACTER,0,MY_RANK,MPI_COMM_WORLD,IERROR)
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

END MODULE SCARC_CPU


!=======================================================================================================================
!
! MODULE SCARC_STACK
!
!> \brief Introduce stack hierarchy for the different consecutive solution methods and
!   organize their alternate calls
!
!=======================================================================================================================
MODULE SCARC_STACK
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_MESSAGES
USE SCARC_STORAGE
USE SCARC_CONVERGENCE

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup environent on specified stack level
! ------------------------------------------------------------------------------------------------------------------
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

! ------------------------------------------------------------------------------------------------------------------
!> \brief Restore environent on specified stack level
! ------------------------------------------------------------------------------------------------------------------
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

! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for Krylov method
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_VECTORS()
USE SCARC_POINTERS, ONLY: G, SV, ST, SCARC_POINT_TO_GRID
INTEGER :: NM, NSTACK, NL

CROUTINE = 'SCARC_SETUP_VECTORS'

! If multiple grid types are used (currently only in MGM method) take care that the vectors are
! allocated in the longest necessary length which corresponds to the structured discretization.
! The related workspaces are also used for possible shorter instances in the unstructured discretization
IF (HAS_MULTIPLE_GRIDS) CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_STRUCTURED)   

DO NSTACK = 1, N_STACK_TOTAL

   SV  => STACK(NSTACK)%SOLVER

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      DO NL = SV%TYPE_LEVEL(1), SV%TYPE_LEVEL(2)

         CALL SCARC_POINT_TO_GRID (NM, NL)                                    
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)

         IF (SV%X /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1 (ST%X, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%X', CROUTINE)
         IF (SV%B /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1 (ST%B, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%B', CROUTINE)
         IF (SV%D /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1 (ST%D, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%D', CROUTINE)
         IF (SV%R /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1 (ST%R, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%R', CROUTINE)
         IF (SV%V /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1 (ST%V, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%V', CROUTINE)
         IF (SV%Y /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1 (ST%Y, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Y', CROUTINE)
         IF (SV%Z /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1 (ST%Z, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Z', CROUTINE)


         IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN 
            IF (SV%X_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB (ST%X_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%X_FB', CROUTINE)
            IF (SV%B_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB (ST%B_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%B_FB', CROUTINE)
            IF (SV%R_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB (ST%R_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%R_FB', CROUTINE)
            IF (SV%V_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB (ST%V_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%V_FB', CROUTINE)
         ENDIF

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_VECTORS

! --------------------------------------------------------------------------------------------------------------
!> \brief  Setup scope for current solver 
! i.e. set pointers to used vectors related to current position in stack
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SCOPE(NS, NP)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NS, NP                          ! references to current stack and parent
 
! If not first solver in stack, store last iteration parameters of parent solver NP
 
IF (NP > 0) THEN
   SVP => STACK(NP)%SOLVER
   SVP%ITE   = ITE
   SVP%RES   = RES
   SVP%RESIN = RESIN
   SVP%CAPPA = CAPPA
ENDIF
 
! Set new environment for solver on stack position NS
 
SV => STACK(NS)%SOLVER

CNAME = SV%CNAME
ITE   = 0
NIT   = SV%NIT
EPS   = SV%EPS
OMEGA = SV%OMEGA
RESIN = SV%RESIN
CAPPA = -1.0

TYPE_PARENT = NP

TYPE_ACCURACY       = SV%TYPE_ACCURACY
TYPE_COARSE         = SV%TYPE_COARSE
TYPE_COARSENING     = SV%TYPE_COARSENING
TYPE_CYCLING        = SV%TYPE_CYCLING
TYPE_GRID           = SV%TYPE_GRID
TYPE_EXCHANGE       = SV%TYPE_EXCHANGE
TYPE_INTERPOL       = SV%TYPE_INTERPOL
TYPE_LEVEL          = SV%TYPE_LEVEL
TYPE_MATRIX         = SV%TYPE_MATRIX
TYPE_METHOD         = SV%TYPE_METHOD
TYPE_METHOD         = SV%TYPE_METHOD
TYPE_MKL            = SV%TYPE_MKL
TYPE_MKL_PRECISION  = SV%TYPE_MKL_PRECISION
TYPE_MULTIGRID      = SV%TYPE_MULTIGRID
TYPE_PARENT         = SV%TYPE_PARENT
TYPE_PRECON         = SV%TYPE_PRECON
TYPE_RELAX          = SV%TYPE_RELAX
TYPE_SCOPE          = SV%TYPE_SCOPE
TYPE_SMOOTH         = SV%TYPE_SMOOTH
TYPE_SOLVER         = SV%TYPE_SOLVER
TYPE_STAGE          = SV%TYPE_STAGE
TYPE_STENCIL        = SV%TYPE_STENCIL
TYPE_TWOLEVEL       = SV%TYPE_TWOLEVEL
TYPE_VECTOR         = SV%TYPE_VECTOR

X = SV%X
B = SV%B
D = SV%D
R = SV%R
V = SV%V
Y = SV%Y
Z = SV%Z

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) ITE_TOTAL = 0

END SUBROUTINE SCARC_SETUP_SCOPE

! --------------------------------------------------------------------------------------------------------------
!> \brief Reset environment of calling routine when leaving solver
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELEASE_SCOPE(NS, NP)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN)  :: NS, NP                            ! references to current stack and parent

SV  => STACK(NS)%SOLVER

! Store convergence information of preceding solver for FDS dump routine
 
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   SCARC_CAPPA      = CAPPA
   SCARC_RESIDUAL   = RES
   SCARC_ITERATIONS = ITE
   !CALL SCARC_CLOSE_CSV_FILE()
ENDIF

SV%RESIN = RESIN
SV%RES   = RES
SV%ITE   = ITE
SV%EPS   = EPS
SV%CAPPA = CAPPA

! If not first solver in stack, reset environment of parent (calling) routine
 
IF (NP > 0) THEN

   SVP => STACK(NP)%SOLVER

   ITE   = SVP%ITE
   NIT   = SVP%NIT
   EPS   = SVP%EPS
   RESIN = SVP%RESIN
   RES   = SVP%RES
   OMEGA = SVP%OMEGA
   CAPPA = SVP%CAPPA

   TYPE_ACCURACY       = SVP%TYPE_ACCURACY
   TYPE_COARSE         = SVP%TYPE_COARSE
   TYPE_COARSENING     = SVP%TYPE_COARSENING
   TYPE_CYCLING        = SVP%TYPE_CYCLING
   TYPE_GRID           = SVP%TYPE_GRID
   TYPE_EXCHANGE       = SVP%TYPE_EXCHANGE
   TYPE_INTERPOL       = SVP%TYPE_INTERPOL
   TYPE_LEVEL          = SVP%TYPE_LEVEL
   TYPE_MATRIX         = SVP%TYPE_MATRIX
   TYPE_METHOD         = SVP%TYPE_METHOD
   TYPE_MKL            = SVP%TYPE_MKL
   TYPE_MKL_PRECISION  = SVP%TYPE_MKL_PRECISION
   TYPE_MULTIGRID      = SVP%TYPE_MULTIGRID
   TYPE_PARENT         = SVP%TYPE_PARENT
   TYPE_PRECON         = SVP%TYPE_PRECON
   TYPE_RELAX          = SVP%TYPE_RELAX
   TYPE_SCOPE          = SVP%TYPE_SCOPE
   TYPE_SMOOTH         = SVP%TYPE_SMOOTH
   TYPE_SOLVER         = SVP%TYPE_SOLVER
   TYPE_STAGE          = SVP%TYPE_STAGE
   TYPE_STENCIL        = SVP%TYPE_STENCIL
   TYPE_TWOLEVEL       = SVP%TYPE_TWOLEVEL
   TYPE_VECTOR         = SVP%TYPE_VECTOR

   X = SVP%X
   B = SVP%B
   D = SVP%D
   R = SVP%R
   V = SVP%V
   Y = SVP%Y
   Z = SVP%Z

ENDIF

END SUBROUTINE SCARC_RELEASE_SCOPE

! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup references to solution vectors related to used stack position (main solver or preconditioner)
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STACK_VECTORS(BX, BB, BD, BR, BV, BY, BZ, NSTACK)
USE SCARC_POINTERS, ONLY: SV
LOGICAL, INTENT(IN) :: BX, BB, BD, BR, BV, BY, BZ
INTEGER, INTENT(IN) :: NSTACK

SV  => STACK(NSTACK)%SOLVER

SELECT CASE (SV%TYPE_STAGE)

   ! Solver from working stage ONE, e.g. main Krylov or MG solver
 
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
 
   ! Solver from working stage TWO, e.g. MG solver as preconditioner
 
   CASE (NSCARC_STAGE_TWO)
      IF (BX) SV%X = NSCARC_VECTOR_TWO_X
      IF (BB) SV%B = NSCARC_VECTOR_TWO_B
      IF (BD) SV%D = NSCARC_VECTOR_TWO_D
      IF (BR) SV%R = NSCARC_VECTOR_TWO_R
      IF (BV) SV%V = NSCARC_VECTOR_TWO_V
      IF (BY) SV%Y = NSCARC_VECTOR_TWO_Y
      IF (BZ) SV%Z = NSCARC_VECTOR_TWO_Z
#ifdef WITH_MKL
     IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
        IF (BX) SV%X_FB = NSCARC_VECTOR_TWO_X
        IF (BB) SV%B_FB = NSCARC_VECTOR_TWO_B
        IF (BR) SV%R_FB = NSCARC_VECTOR_TWO_R
        IF (BV) SV%V_FB = NSCARC_VECTOR_TWO_V
     ENDIF
#endif

END SELECT

END SUBROUTINE SCARC_SETUP_STACK_VECTORS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup convergence environment for stack position of requested preconditioner
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STACK_PRECON(NSTACK, NSCOPE)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NSTACK, NSCOPE
 
! Point to stack entry of called preconditioner and types for preconditioner
 
SV => STACK(NSTACK)%SOLVER

SV%TYPE_SCOPE(0) = NSCOPE

SV%EPS   =  SCARC_PRECON_ACCURACY
SV%NIT   =  SCARC_PRECON_ITERATIONS
SV%OMEGA =  SCARC_PRECON_OMEGA

! Preset name of preconditioning method
 
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
   CASE (NSCARC_RELAX_LU)
      SV%CNAME = 'SCARC_PRECON_LU'
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
   CASE (NSCARC_RELAX_OPTIMIZED)
      SV%CNAME = 'SCARC_PRECON_OPTIMIZED'
      SV%NIT   = 1
      SV%OMEGA = 1.0_EB
#endif
   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_PRECON)
END SELECT

! Preset types for preconditioner (use same as for calling solver)
 
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

! Preset pointers for preconditioner (use same as for alling solver)
 
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

END SUBROUTINE SCARC_SETUP_STACK_PRECON

! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup convergence environment for stack position of requested smoother
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCOPE)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NSTACK, NSCOPE
 
! Basic setup of stack information and types/names for smoother
 
CALL SCARC_SETUP_STACK(NSTACK)

SV => STACK(NSTACK)%SOLVER
SV%TYPE_SCOPE(0) = NSCOPE
SV%OMEGA =  SCARC_SMOOTH_OMEGA

SELECT CASE(TYPE_SMOOTH)
   CASE (NSCARC_RELAX_JAC)
      SV%CNAME = 'SCARC_SMOOTHER_JAC'
   CASE (NSCARC_RELAX_SSOR)
      SV%CNAME = 'SCARC_SMOOTHER_SSOR'
   CASE (NSCARC_RELAX_MJAC)
      SV%CNAME = 'SCARC_SMOOTHER_MJAC'
   CASE (NSCARC_RELAX_MGS)
      SV%CNAME = 'SCARC_SMOOTHER_MGS'
   CASE (NSCARC_RELAX_MSGS)
      SV%CNAME = 'SCARC_SMOOTHER_MSGS'
   CASE (NSCARC_RELAX_MSOR)
      SV%CNAME = 'SCARC_SMOOTHER_MSOR'
   CASE (NSCARC_RELAX_MSSOR)
      SV%CNAME = 'SCARC_SMOOTHER_MSSOR'
   CASE (NSCARC_RELAX_ILU)
      SV%CNAME = 'SCARC_SMOOTHER_ILU'
   CASE (NSCARC_RELAX_FFT)
      SV%CNAME = 'SCARC_SMOOTHER_FFT'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_FFTO)
      SV%CNAME = 'SCARC_SMOOTHER_FFTO'
      SV%OMEGA = 1.0_EB
#ifdef WITH_MKL
   CASE (NSCARC_RELAX_MKL)
      SV%OMEGA = 1.0_EB
      IF (NSCOPE == NSCARC_SCOPE_GLOBAL) THEN
         SV%CNAME = 'SCARC_SMOOTHER_CLUSTER'
      ELSE
         SV%CNAME = 'SCARC_SMOOTHER_PARDISO'
      ENDIF
   CASE (NSCARC_RELAX_OPTIMIZED)
      SV%CNAME = 'SCARC_SMOOTHER_OPTIMIZED'
      SV%OMEGA = 1.0_EB
#endif
   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_SMOOTH)
END SELECT

! Preset iteration parameters for Multigrid method
 
SV%EPS = SCARC_SMOOTH_ACCURACY
SV%NIT = SCARC_SMOOTH_ITERATIONS

! Preset types for preconditioner (use same descriptors as calling solver)
 
SVP => STACK(NSTACK-1)%SOLVER

SV%TYPE_SOLVER        = NSCARC_SOLVER_SMOOTH
SV%TYPE_STAGE         = SVP%TYPE_STAGE
SV%TYPE_LEVEL(0:2)    = SVP%TYPE_LEVEL(0:2)
SV%TYPE_RELAX         = SVP%TYPE_RELAX
SV%TYPE_GRID          = SVP%TYPE_GRID
SV%TYPE_INTERPOL      = SVP%TYPE_INTERPOL
SV%TYPE_ACCURACY      = SVP%TYPE_ACCURACY
SV%TYPE_MKL_PRECISION = SVP%TYPE_MKL_PRECISION
SV%TYPE_MATRIX        = SVP%TYPE_MATRIX

! Preset references for preconditioner (use same pointers as calling solver)

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

END SUBROUTINE SCARC_SETUP_STACK_SMOOTH

! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup convergence environment for stack position of MGM method 
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STACK_MGM(NSTACK, NSCOPE)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NSTACK, NSCOPE
 
! Point to stack entry of called preconditioner and types for preconditioner
 
SV => STACK(NSTACK)%SOLVER

SV%TYPE_SCOPE(0) = NSCOPE
SV%TYPE_SOLVER   = NSCARC_SOLVER_MGM
SV%TYPE_RELAX    = TYPE_PRECON
SV%TYPE_GRID     = TYPE_GRID

SV%EPS   =  SCARC_PRECON_ACCURACY
SV%NIT   =  SCARC_PRECON_ITERATIONS
SV%OMEGA =  SCARC_PRECON_OMEGA

! Preset name of preconditioning method
 
SELECT CASE(TYPE_MGM_LAPLACE)
   CASE (NSCARC_MGM_LAPLACE_CG)
      SV%CNAME = 'SCARC_MGM_LAPLACE_CG'
   CASE (NSCARC_MGM_LAPLACE_LU)
      SV%CNAME = 'SCARC_MGM_LAPLACE_LU'
   CASE (NSCARC_MGM_LAPLACE_LUPERM)
      SV%CNAME = 'SCARC_MGM_LAPLACE_LUPERM'
#ifdef WITH_MKL
   CASE (NSCARC_MGM_LAPLACE_PARDISO)
      SV%CNAME = 'SCARC_MGM_LAPLACE_PARDISO'
   CASE (NSCARC_MGM_LAPLACE_OPTIMIZED)
      SV%CNAME = 'SCARC_MGM_LAPLACE_OPTIMIZED'
#endif
   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_MGM_LAPLACE)
END SELECT

! Preset types for preconditioner (use same as for calling solver)
 
SVP => STACK(NSTACK-1)%SOLVER

SV%TYPE_ACCURACY      = SVP%TYPE_ACCURACY
SV%TYPE_LEVEL(0:2)    = SVP%TYPE_LEVEL(0:2)
SV%TYPE_MATRIX        = SVP%TYPE_MATRIX
SV%TYPE_MKL_PRECISION = SVP%TYPE_MKL_PRECISION
SV%TYPE_STAGE         = SVP%TYPE_STAGE

! Preset pointers for preconditioner (use same as for alling solver)
 
SV%X = SVP%X
SV%B = SVP%B

IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
   SV%X_FB = SVP%X_FB
   SV%B_FB = SVP%B_FB
ENDIF

END SUBROUTINE SCARC_SETUP_STACK_MGM

! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup convergence environment for stack position of Krylov method
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STACK_KRYLOV(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX

! Basic setup of stack information and types for Krylov method
 
CALL SCARC_SETUP_STACK(NSTACK)

SV  => STACK(NSTACK)%SOLVER
SV%TYPE_METHOD   = NSCARC_METHOD_KRYLOV
SV%TYPE_SOLVER   = NSOLVER
SV%TYPE_SCOPE(0) = NSCOPE
SV%TYPE_STAGE    = NSTAGE
SV%TYPE_LEVEL(1) = NLMIN
SV%TYPE_LEVEL(2) = NLMAX
SV%TYPE_MATRIX   = TYPE_MATRIX

! Preset iteration parameters for Krylov method

SELECT CASE(NSOLVER)

   ! Krylov method is used as main solver

   CASE (NSCARC_SOLVER_MAIN)
   
      SV%CNAME = 'SCARC_MAIN_KRYLOV'
   
      SV%EPS = SCARC_KRYLOV_ACCURACY
      SV%NIT = SCARC_KRYLOV_ITERATIONS
   
      SV%TYPE_RELAX    = TYPE_PRECON
      SV%TYPE_TWOLEVEL = TYPE_TWOLEVEL
   
   ! Krylov method is used as local MGM Laplace solver solver

   CASE (NSCARC_SOLVER_MGM)
   
      SV%CNAME = 'SCARC_MGM_KRYLOV'
   
      SV%EPS = SCARC_KRYLOV_ACCURACY
      SV%NIT = SCARC_KRYLOV_ITERATIONS
   
      SV%TYPE_RELAX    = TYPE_PRECON
      SV%TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE
   
   ! Krylov method is used as coarse grid solver solver

   CASE (NSCARC_SOLVER_COARSE)
   
      SV%CNAME = 'SCARC_COARSE_KRYLOV'
   
      SV%EPS = SCARC_COARSE_ACCURACY
      SV%NIT = SCARC_COARSE_ITERATIONS
   
      SV%TYPE_RELAX    = NSCARC_RELAX_SSOR             ! only use SSOR-preconditioning for coarse solver
      SV%TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE          ! only use one level for coarse solver
   
   ! Otherwise: print error message

   CASE DEFAULT
   
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)
   
END SELECT

 
! Point to solution vectors (in corresponding scope)
 
CALL SCARC_SETUP_STACK_VECTORS(.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_STACK_KRYLOV

! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup convergence environment for stack position of Multigrid method
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STACK_MULTIGRID(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX
 
! Basic setup of stack information and types for multigrid method
 
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

   ! Multigrid method is used as main solver

   CASE (NSCARC_SOLVER_MAIN)
      SV%CNAME = 'SCARC_MAIN_MG'

   ! Multigrid method is only used as preconditioner for global Krylov method

   CASE (NSCARC_SOLVER_PRECON)
      SV%CNAME = 'SCARC_PRECONDITIONER'

   ! Otherwise: print error message

   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)

END SELECT

! Preset iteration parameters for Multigrid method
 
SV%EPS = SCARC_MULTIGRID_ACCURACY
SV%NIT = SCARC_MULTIGRID_ITERATIONS
 
! Point to solution vectors (in corresponding scope)
 
CALL SCARC_SETUP_STACK_VECTORS(.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_STACK_MULTIGRID

END MODULE SCARC_STACK


!=======================================================================================================================
! 
! MODULE SCARC_PARSER
! 
!> \brief Parse ScaRC related input parameters and initialize basic ScaRC structures
!
!=======================================================================================================================
MODULE SCARC_PARSER
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_MESSAGES, ONLY: MSG
USE SCARC_TROUBLESHOOTING
USE SCARC_STACK, ONLY: STACK

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! ------------------------------------------------------------------------------------------------------------------
!> \brief Determine input parameters needed for setting up the discretization structure
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PARSE_INPUT

ITERATE_PRESSURE = .TRUE.  ! Although there is no need to do pressure iterations to drive down 
                           ! velocity error leave it .TRUE. to write out velocity error diagnostics
 
! ------------- Set type of discretization
 
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
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_GRID, NSCARC_NONE)
END SELECT
 
! ------------ Set type of matrix storage (COMPACT/BANDWISE)
 
SELECT CASE (TRIM(SCARC_MATRIX))
   CASE ('COMPACT')
      TYPE_MATRIX = NSCARC_MATRIX_COMPACT
   CASE ('BANDWISE')
      TYPE_MATRIX = NSCARC_MATRIX_BANDWISE
   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MATRIX, NSCARC_NONE)
END SELECT

! ------------ Set type of matrix stencil (CONSTANT/VARIABLE)
 
SELECT CASE (TRIM(SCARC_STENCIL))
   CASE ('CONSTANT')
      TYPE_STENCIL = NSCARC_STENCIL_CONSTANT
   CASE ('VARIABLE')
      TYPE_STENCIL = NSCARC_STENCIL_VARIABLE
   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_STENCIL, NSCARC_NONE)
END SELECT
 
! ------------ Set type of global solver
 
SELECT CASE (TRIM(SCARC_METHOD))

   ! ------------------------- Global Krylov solver ----------------------------------
   CASE ('KRYLOV')

      TYPE_METHOD = NSCARC_METHOD_KRYLOV

      ! Set type of two-level method
      SELECT CASE (TRIM(SCARC_TWOLEVEL))
         CASE ('NONE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE
         CASE ('ADDITIVE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_COARSE_ADD
         CASE ('MULTIPLICATIVE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_COARSE_MUL1
         CASE ('MULTIPLICATIVE2')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_COARSE_MUL2
         CASE ('COARSE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_COARSE_ONLY
         CASE ('MACRO')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_MACRO
         CASE ('XMEAN_ADD')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_XMEAN_ADD
         CASE ('XMEAN_MUL1')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_XMEAN_MUL1
         CASE ('XMEAN_MUL2')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_XMEAN_MUL2
         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_TWOLEVEL, NSCARC_NONE)
      END SELECT

      ! Set type of interpolation for two-level Krylov method
      SELECT CASE (TRIM(SCARC_KRYLOV_INTERPOL))
         CASE ('NONE')
            TYPE_INTERPOL = NSCARC_UNDEF_INT
         CASE ('CONSTANT')
            TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
         CASE ('BILINEAR')
            TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
         CASE ('BILINEAR2')
            TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR2
         CASE ('AMG')
            TYPE_INTERPOL = NSCARC_INTERPOL_AMG
         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_KRYLOV_INTERPOL, NSCARC_NONE)
      END SELECT

      ! Set type of preconditioner (JACOBI/SSOR/MGS/MSGS/MSOR/MSSOR/ILU/LU/FFT/GMG/PARDISO/CLUSTER/OPTIMIZED)
      SELECT CASE (TRIM(SCARC_PRECON))
         CASE ('JACOBI')                                    ! Jacobi preconditioner
            TYPE_PRECON = NSCARC_RELAX_JAC
         CASE ('SSOR')                                      ! Symmetric SOR preconditioner
            TYPE_PRECON = NSCARC_RELAX_SSOR
         CASE ('MJAC')                                      ! Jacobi preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MJAC
         CASE ('MGS')                                       ! Gauss-Seidel preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MGS
         CASE ('MSGS')                                      ! Symmetric Gauss-Seidel preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MSGS
         CASE ('MSOR')                                      ! SOR preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MSOR
         CASE ('MSSOR')                                     ! Symmetric SOR preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MSSOR
         CASE ('ILU')                                       ! ILU preconditioner
            TYPE_PRECON = NSCARC_RELAX_ILU
         CASE ('LU')                                        ! LU preconditioner
            TYPE_PRECON = NSCARC_RELAX_LU
         CASE ('MULTIGRID')                                 ! Multigrid preconditioner
            TYPE_PRECON = NSCARC_RELAX_GMG
            SELECT CASE (TRIM(SCARC_SMOOTH))
               CASE ('JACOBI')
                  TYPE_SMOOTH = NSCARC_RELAX_JAC
               CASE ('SSOR')
                  TYPE_SMOOTH = NSCARC_RELAX_SSOR
               CASE ('ILU')
                  TYPE_SMOOTH = NSCARC_RELAX_ILU
               CASE ('FFT')
                  IF (IS_UNSTRUCTURED) CALL SCARC_ERROR(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
                  TYPE_SMOOTH = NSCARC_RELAX_FFT
               CASE ('PARDISO')
#ifdef WITH_MKL
                  TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
                  TYPE_SMOOTH = NSCARC_RELAX_SSOR
                  CALL SCARC_WARNING(NSCARC_WARNING_NO_MKL_SMOOTH, SCARC_NONE, NSCARC_NONE)
#endif

               CASE ('CLUSTER')
#ifdef WITH_MKL
                  TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
                  TYPE_SMOOTH = NSCARC_RELAX_SSOR
                  CALL SCARC_ERROR(NSCARC_WARNING_NO_MKL_SMOOTH, SCARC_NONE, NSCARC_NONE)
#endif
               CASE DEFAULT
                  CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_SMOOTH, NSCARC_NONE)
            END SELECT
         CASE ('FFT')                                                ! FFT preconditioner
            IF (IS_UNSTRUCTURED) CALL SCARC_ERROR(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            TYPE_PRECON = NSCARC_RELAX_FFT
         CASE ('FFTO')                                               ! FFT with overlap preconditioner
            IF (IS_UNSTRUCTURED) CALL SCARC_ERROR(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            IF (NMESHES == 1) THEN
               TYPE_PRECON = NSCARC_RELAX_FFT
            ELSE
               TYPE_PRECON = NSCARC_RELAX_FFTO
            ENDIF
         CASE ('OPTIMIZED')                                          ! LU preconditioner based on either FFT or PARDISO
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_OPTIMIZED
            TYPE_MKL(0)   = NSCARC_MKL_LOCAL
            TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
#else
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_PRECON, NSCARC_NONE)
#endif
         CASE ('PARDISO')                                            ! LU preconditioner based on MKL-PARDISO
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_MKL
            TYPE_MKL(0)   = NSCARC_MKL_LOCAL
            TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
#else
            TYPE_PRECON   = NSCARC_RELAX_SSOR
            CALL SCARC_WARNING(NSCARC_WARNING_NO_MKL_PRECON, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('CLUSTER')                            !  LU-preconditioner based on MKL Cluster_Sparse_Solver
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_MKL
            TYPE_MKL(0)   = NSCARC_MKL_GLOBAL
            TYPE_SCOPE(1) = NSCARC_SCOPE_GLOBAL
#else
            TYPE_PRECON   = NSCARC_RELAX_SSOR
            CALL SCARC_WARNING(NSCARC_WARNING_NO_MKL_PRECON, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            TYPE_PRECON   = NSCARC_RELAX_SSOR
            TYPE_MKL(0)   = NSCARC_MKL_NONE
            TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
            CALL SCARC_WARNING(NSCARC_WARNING_ONLY_SSOR_PRECON, SCARC_NONE, NSCARC_NONE)
      END SELECT

      ! set type scope for preconditioner (GLOBAL/LOCAL)
      SELECT CASE (TRIM(SCARC_PRECON_SCOPE))
         CASE ('GLOBAL')
            TYPE_SCOPE(1) = NSCARC_SCOPE_GLOBAL
         CASE ('LOCAL')
            TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
         CASE DEFAULT
            TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
            CALL SCARC_WARNING(NSCARC_WARNING_NO_GLOBAL_SCOPE, SCARC_PRECON_SCOPE, NSCARC_NONE)
      END SELECT

   ! ------------------------- Global geometric multigrid solver -------------------------------
   CASE ('MULTIGRID')

      TYPE_METHOD = NSCARC_METHOD_MULTIGRID

      ! Set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_MULTIGRID))
         CASE ('GEOMETRIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
            TYPE_COARSENING = NSCARC_COARSENING_GMG         ! GMG-default, may be overwritten by SCARC_COARSENING 
         CASE ('ALGEBRAIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
            TYPE_COARSENING = NSCARC_COARSENING_CUBIC       ! AMG-default, may be overwritten by SCARC_COARSENING 
         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID, NSCARC_NONE)
      END SELECT

      ! Set type of smoother (JACOBI/SGS/SSOR/MSSOR/ILU/PARDISO/CLUSTER)
      SELECT CASE (TRIM(SCARC_SMOOTH))                      ! use same parameters as for preconditioner
         CASE ('JACOBI')                                    ! Jacobi preconditioner
            TYPE_SMOOTH = NSCARC_RELAX_JAC
         CASE ('SSOR')                                      ! SSOR preconditioner
            TYPE_SMOOTH = NSCARC_RELAX_SSOR
         CASE ('MJAC')                                      ! Jacobi preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MJAC
         CASE ('MGS')                                       ! Gauss-Seidel preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MGS
         CASE ('MSGS')                                      ! Symmetric Gauss-Seidel preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MSGS
         CASE ('MSOR')                                      ! SOR preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MSOR
         CASE ('MSSOR')                                     ! Symmetric SOR preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MSSOR
         CASE ('ILU')                                       ! ILU preconditioner
            TYPE_SMOOTH = NSCARC_RELAX_ILU
         CASE ('FFT')
            IF (IS_UNSTRUCTURED) CALL SCARC_ERROR(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            TYPE_SMOOTH = NSCARC_RELAX_FFT
         CASE ('FFTO')
            IF (IS_UNSTRUCTURED) CALL SCARC_ERROR(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            IF (NMESHES == 1) THEN
               TYPE_SMOOTH = NSCARC_RELAX_FFT
            ELSE
               TYPE_SMOOTH = NSCARC_RELAX_FFTO
            ENDIF
         CASE ('OPTIMIZED')
#ifdef WITH_MKL
            TYPE_SMOOTH = NSCARC_RELAX_OPTIMIZED
#else
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_SMOOTH, NSCARC_NONE)
#endif
         CASE ('PARDISO')
#ifdef WITH_MKL
            TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
            TYPE_SMOOTH = NSCARC_RELAX_SSOR
            CALL SCARC_WARNING(NSCARC_WARNING_NO_MKL_SMOOTH, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('CLUSTER')
#ifdef WITH_MKL
            TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
            TYPE_SMOOTH = NSCARC_RELAX_SSOR
            CALL SCARC_WARNING(NSCARC_WARNING_NO_MKL_SMOOTH, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            TYPE_SMOOTH = NSCARC_RELAX_SSOR
            CALL SCARC_WARNING(NSCARC_WARNING_ONLY_SSOR_SMOOTH, SCARC_NONE, NSCARC_NONE)
      END SELECT

      ! set type scope for smoother (GLOBAL/LOCAL)
      SELECT CASE (TRIM(SCARC_SMOOTH_SCOPE))
         CASE ('GLOBAL')
            TYPE_SCOPE(2) = NSCARC_SCOPE_GLOBAL
         CASE ('LOCAL')
            TYPE_SCOPE(2) = NSCARC_SCOPE_LOCAL
         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_SMOOTH_SCOPE, NSCARC_NONE)
      END SELECT

   ! ------------------------- Global LU-decomposition solver -------------------------------
   CASE ('MKL')

      TYPE_METHOD  = NSCARC_METHOD_LU

      ! Set type of MKL method (global/local)
      SELECT CASE (TRIM(SCARC_MKL_SCOPE))                  
         CASE ('GLOBAL')
#ifdef WITH_MKL
            TYPE_MKL(0)   = NSCARC_MKL_GLOBAL
#else
            CALL SCARC_ERROR(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('LOCAL')
#ifdef WITH_MKL
            TYPE_MKL(0)   = NSCARC_MKL_LOCAL
#else
            CALL SCARC_ERROR(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MKL_SCOPE, NSCARC_NONE)
      END SELECT

   ! ------------------------- McKenny-Greengard-Mayo solver -------------------------
   CASE ('MGM')

      !  Both structured and unstructured discretization are required
      !  The second pass of this method is purely locally based (the Laplace solutions) 

      HAS_MULTIPLE_GRIDS = .TRUE.

      TYPE_METHOD = NSCARC_METHOD_MGM

      ! set type of MGM interface BCs of Laplace problems
      SELECT CASE (TRIM(SCARC_MGM_BC))
         CASE ('MEAN')
            TYPE_MGM_BC = NSCARC_MGM_BC_MEAN
         CASE ('TRUE')
            TYPE_MGM_BC = NSCARC_MGM_BC_TRUE
         CASE ('EXTRAPOLATION')
            TYPE_MGM_BC = NSCARC_MGM_BC_EXPOL
         CASE ('TAYLOR')
            TYPE_MGM_BC = NSCARC_MGM_BC_TAYLOR
         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MGM_BC, NSCARC_NONE)
      END SELECT

      ! set type of MGM solver for local Laplace problems 
      SELECT CASE (TRIM(SCARC_MGM_LAPLACE_SOLVER))
         CASE ('AMG')
            TYPE_MGM_LAPLACE = NSCARC_MGM_LAPLACE_AMG
         CASE ('CG','KRYLOV')
            TYPE_MGM_LAPLACE = NSCARC_MGM_LAPLACE_CG
         CASE ('FFT')
            TYPE_MGM_LAPLACE = NSCARC_MGM_LAPLACE_FFT
         CASE ('LU')
            TYPE_MGM_LAPLACE = NSCARC_MGM_LAPLACE_LU
         CASE ('LUPERM')
            TYPE_MGM_LAPLACE = NSCARC_MGM_LAPLACE_LUPERM
         CASE ('OPTIMIZED')
            TYPE_MGM_LAPLACE = NSCARC_MGM_LAPLACE_OPTIMIZED
         CASE ('PARDISO')
            TYPE_MGM_LAPLACE = NSCARC_MGM_LAPLACE_PARDISO
         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MGM_LAPLACE_SOLVER, NSCARC_NONE)
      END SELECT

      ! set type of MGM interpolation for interface BCs of Laplace problems
      SELECT CASE (TRIM(SCARC_MGM_INTERPOLATION))
         CASE ('LINEAR')
            TYPE_MGM_INTERPOL = NSCARC_MGM_INTERPOL_LINEAR
         CASE ('SQUARE')
            TYPE_MGM_INTERPOL = NSCARC_MGM_INTERPOL_SQUARE
         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MGM_INTERPOLATION, NSCARC_NONE)
      END SELECT

   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_METHOD, NSCARC_NONE)

END SELECT

! If a multigrid solver is used (either as main solver or as preconditioner)
! set types for multigrid, coarse grid solver and cycling pattern
 
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_RELAX_GMG) THEN

   ! Set type of multigrid (GEOMETRIC/ALGEBRAIC with corresponding coarsening strategy)
   SELECT CASE (TRIM(SCARC_MULTIGRID))
      CASE ('GEOMETRIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
      CASE ('ALGEBRAIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
      CASE DEFAULT
         CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID, NSCARC_NONE)
   END SELECT

   ! Set type of cycling pattern (F/V/W)
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
         CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_CYCLE, NSCARC_NONE)
   END SELECT

   ! Set type of interpolation 
   SELECT CASE (TRIM(SCARC_MULTIGRID_INTERPOL))
      CASE ('STANDARD')
         TYPE_INTERPOL = NSCARC_INTERPOL_STANDARD
      CASE ('CONSTANT')
         TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
      CASE ('BILINEAR')
         TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
      CASE ('BILINEAR2')
         TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR2
      CASE ('CLASSICAL')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL
      CASE ('DIRECT')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT
      CASE ('AMG')
         TYPE_INTERPOL = NSCARC_INTERPOL_AMG
      CASE DEFAULT
         CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_INTERPOL, NSCARC_NONE)
   END SELECT

ENDIF

! ------------ Set type of coarsening strategy in case of multi-level methods

SELECT CASE (TRIM(SCARC_COARSENING))
   CASE ('AGGREGATED')
      TYPE_COARSENING = NSCARC_COARSENING_AGGREGATED
   CASE ('AGGREGATEDS')
      TYPE_COARSENING = NSCARC_COARSENING_AGGREGATED_S
   CASE ('CUBIC')
      TYPE_COARSENING = NSCARC_COARSENING_CUBIC
   CASE ('GMG')
      TYPE_COARSENING = NSCARC_COARSENING_GMG
   CASE DEFAULT
      !CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_COARSENING, NSCARC_NONE)
      TYPE_COARSENING = NSCARC_COARSENING_AGGREGATED
END SELECT

! Set type of accuracy (ABSOLUTE/RELATIVE)
 
SELECT CASE (TRIM(SCARC_ACCURACY))
   CASE ('ABSOLUTE')
      TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
   CASE ('RELATIVE')
      TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_ACCURACY, NSCARC_NONE)
END SELECT

! Set type of precision for MKL solver (SINGLE/DOUBLE)
 
SELECT CASE (TRIM(SCARC_MKL_PRECISION))
   CASE ('SINGLE')
      TYPE_MKL_PRECISION = NSCARC_PRECISION_SINGLE
   CASE ('DOUBLE')
      TYPE_MKL_PRECISION = NSCARC_PRECISION_DOUBLE
   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_MKL_PRECISION, NSCARC_NONE)
END SELECT

! -------- Define some logical variables - just for notational convenience
 
IS_STRUCTURED   = (TYPE_GRID == NSCARC_GRID_STRUCTURED)
IS_UNSTRUCTURED = (TYPE_GRID == NSCARC_GRID_UNSTRUCTURED)

IS_CG     = (TYPE_METHOD == NSCARC_METHOD_KRYLOV)
IS_CG_MG  = IS_CG .AND. (TYPE_PRECON == NSCARC_RELAX_GMG) 
IS_CG_GMG = IS_CG .AND. (TYPE_PRECON == NSCARC_RELAX_GMG) .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC)
IS_CG_AMG = IS_CG .AND. (TYPE_PRECON == NSCARC_RELAX_GMG) .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC)

IS_MG  = (TYPE_METHOD == NSCARC_METHOD_MULTIGRID)
IS_GMG = IS_MG .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC)
IS_AMG = IS_MG .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC)

IS_FFT = (TYPE_PRECON == NSCARC_RELAX_FFT)  .OR. (TYPE_SMOOTH == NSCARC_RELAX_FFT)
IS_FFTO= (TYPE_PRECON == NSCARC_RELAX_FFTO) .OR. (TYPE_SMOOTH == NSCARC_RELAX_FFTO)
IS_MKL = (TYPE_PRECON >= NSCARC_RELAX_MKL)  .OR. (TYPE_SMOOTH >= NSCARC_RELAX_MKL) 

IF (IS_CG .AND. (TYPE_PRECON /= NSCARC_RELAX_GMG)) THEN
   IF (TYPE_TWOLEVEL >= NSCARC_TWOLEVEL_XMEAN_ADD) THEN 
      HAS_XMEAN_LEVELS = .TRUE.
   ELSE IF (TYPE_TWOLEVEL > NSCARC_TWOLEVEL_NONE) THEN
      HAS_TWO_LEVELS = .TRUE.
   ENDIF
ENDIF
HAS_MULTIPLE_LEVELS = IS_MG .OR. IS_CG_MG .OR. HAS_TWO_LEVELS 

IS_CG_ADD    = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_COARSE_ADD)
IS_CG_MUL1   = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_COARSE_MUL1)
IS_CG_MUL2   = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_COARSE_MUL1)
IS_CG_COARSE = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_COARSE_ONLY)
IS_CG_MACRO  = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_MACRO)

HAS_GMG_LEVELS = IS_GMG .OR. IS_CG_GMG .OR. HAS_TWO_LEVELS
HAS_AMG_LEVELS = IS_AMG .OR. IS_CG_AMG 

IS_MGM = TYPE_METHOD == NSCARC_METHOD_MGM

! If two or more grid levels are used, also set type of coarse grid solver

IF (HAS_TWO_LEVELS .OR. HAS_MULTIPLE_LEVELS) THEN
   SELECT CASE (TRIM(SCARC_COARSE))
      CASE ('ITERATIVE')
         TYPE_COARSE = NSCARC_COARSE_ITERATIVE
      CASE ('DIRECT')
#ifdef WITH_MKL
         TYPE_COARSE   = NSCARC_COARSE_DIRECT
         TYPE_MKL(0)   = NSCARC_MKL_COARSE
#else
         CALL SCARC_ERROR(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
      CASE DEFAULT
         CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_COARSE, NSCARC_NONE)
   END SELECT
ENDIF

END SUBROUTINE SCARC_PARSE_INPUT

END MODULE SCARC_PARSER


!=======================================================================================================================
!
! MODULE SCARC_MPI
!
!> \brief Organize MPI communication structures needed for the different data exchanges in ScaRC/UScaRC
!
!=======================================================================================================================
MODULE SCARC_MPI
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MESH_VARIABLES, ONLY: MESHES
USE SCARC_CONSTANTS
USE SCARC_TYPES, ONLY: SCARC_GRID_TYPE
USE SCARC_VARIABLES
USE SCARC_UTILITIES, ONLY: ARE_FACE_NEIGHBORS, ARE_FACE_NEIGHBORS, SET_MATRIX_TYPE, SCARC_SET_GRID_TYPE
USE SCARC_TROUBLESHOOTING, ONLY: SCARC_ERROR
USE SCARC_STORAGE, ONLY: SCARC_ALLOCATE_INT1, SCARC_ALLOCATE_REAL1
USE SCARC_MESSAGES, ONLY: MSG

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS


! ---------------------------------------------------------------------------------------------------------------
!> \brief Setup dimensions for data exchanges
! ---------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS(G, OG, NOM, IREFINE)
USE SCARC_POINTERS, ONLY: GWC
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G, OG
INTEGER, INTENT(IN) :: NOM, IREFINE
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IW
LOGICAL :: FOUND

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

   GWC => G%WALL(IW)
   IF (GWC%NOM/=NOM) CYCLE SEARCH_LOOP
   FOUND = .TRUE.

   SELECT CASE (GWC%IOR)
      CASE ( 1) 
         IMIN=MAX(IMIN,GWC%IXN(1)-1)
      CASE (-1) 
         IMAX=MIN(IMAX,GWC%IXN(2)+1)
      CASE ( 2) 
         JMIN=MAX(JMIN,GWC%IYN(1)-1)
      CASE (-2) 
         JMAX=MIN(JMAX,GWC%IYN(2)+1)
      CASE ( 3) 
         KMIN=MAX(KMIN,GWC%IZN(1)-1)
      CASE (-3) 
         KMAX=MIN(KMAX,GWC%IZN(2)+1)
   END SELECT
ENDDO SEARCH_LOOP

N_EXCHANGES = N_EXCHANGES+1

END SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate several global structures for data exchange
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NP

CROUTINE = 'SCARC_SETUP_GLOBALS'

IF (N_MPI_PROCESSES > 1) THEN

   ! Allocate and preset counter and displacement vector for global data exchanges
   CALL SCARC_ALLOCATE_INT1 (COUNTS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'COUNTS', CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (DISPLS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'DISPLS', CROUTINE)

   ! Get number of data to send per process
   DO NP = 0, N_MPI_PROCESSES-1
      DO NM = 1, NMESHES
         IF (PROCESS(NM)==NP) COUNTS(NP) = COUNTS(NP) + 1
      ENDDO
   ENDDO

   ! Get displacements on communication vector for all meshes
   DO NP = 1, N_MPI_PROCESSES-1
      DISPLS(NP) = COUNTS(NP-1) + DISPLS(NP-1)
   ENDDO

ENDIF

CALL SCARC_ALLOCATE_INT1 (MESH_INT , 1, NMESHES, NSCARC_INIT_ZERO, 'MESH_INT', CROUTINE)
CALL SCARC_ALLOCATE_REAL1 (MESH_REAL, 1, NMESHES, NSCARC_INIT_ZERO, 'MESH_REAL', CROUTINE)

END SUBROUTINE SCARC_SETUP_GLOBALS


! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate workspace for data exchanges of different data types and sizes and perform basic exchanges
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGES
USE SCARC_POINTERS, ONLY: S, OS, OL, OG, OGF, SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID
INTEGER :: NL, NM, NOM, NLEN
INTEGER :: INBR

CROUTINE = 'SCARC_SETUP_EXCHANGES'
 
! Allocate request array for data exchanges
! Exchange basic information about wall sizes (needed for the dimensioning of the exchange buffers)
 
IF (N_MPI_PROCESSES>1) THEN
   ALLOCATE (REQ(N_EXCHANGES*40), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_EXCHANGES', 'REQ', IERROR)
   REQ = MPI_REQUEST_NULL
ENDIF
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_BASIC_SIZES, NSCARC_NONE, NLEVEL_MIN)

! Allocate send and receive buffers (real and integer) in correct lengths
! These are allocated with sizes according to the requirements of the finest grid level
! In case of a multi-level method, they are also used for the coarser levels (with shorter exchange sizes)

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)                          

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NLEVEL_MIN)

      ! Allocate buffers in maximum length for finest grid level (same buffers are used in shorter length on coarse levels, too)

      IF (.NOT.IS_MGM) THEN
   
         OS%NLEN_MAX_BUFFER_LAYER1  = MAX(OS%NLEN_MAX_BUFFER_LAYER1, (OG%NCG+1) * 1)
         OS%NLEN_MAX_BUFFER_LAYER2  = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OG%NCG+1) * 2)
         OS%NLEN_MAX_BUFFER_LAYER4  = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OG%NCG+1) * 4)
         OS%NLEN_MAX_BUFFER_STENCIL = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OG%NCG+1) * NSCARC_MAX_STENCIL)
         OS%NLEN_MAX_BUFFER_FULL    = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OG%NCG+1) * NSCARC_MAX_STENCIL * 4)

         OG%NLEN_BUFFER_LAYER1  = OS%NLEN_MAX_BUFFER_LAYER1
         OG%NLEN_BUFFER_LAYER2  = OS%NLEN_MAX_BUFFER_LAYER2
         OG%NLEN_BUFFER_LAYER4  = OS%NLEN_MAX_BUFFER_LAYER4
         OG%NLEN_BUFFER_STENCIL = OS%NLEN_MAX_BUFFER_STENCIL
         OG%NLEN_BUFFER_FULL    = OS%NLEN_MAX_BUFFER_FULL

         NLEN = 2*OS%NLEN_MAX_BUFFER_FULL 
   
         CALL SCARC_ALLOCATE_INT1 (OS%SEND_BUFFER_INT , 1, NLEN, NSCARC_INIT_HUGE, 'OS%SEND_BUFFER_INT', CROUTINE)
         CALL SCARC_ALLOCATE_INT1 (OS%RECV_BUFFER_INT , 1, NLEN, NSCARC_INIT_HUGE, 'OS%RECV_BUFFER_INT', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1 (OS%SEND_BUFFER_REAL, 1, NLEN, NSCARC_INIT_HUGE, 'OS%SEND_BUFFER_REAL', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1 (OS%RECV_BUFFER_REAL, 1, NLEN, NSCARC_INIT_HUGE, 'OS%RECV_BUFFER_REAL', CROUTINE)

      ELSE

         OGF => OL%UNSTRUCTURED

         OS%NLEN_MAX_BUFFER_LAYER1  = MAX(OS%NLEN_MAX_BUFFER_LAYER1, (OGF%NCG+1) * 1)
         OS%NLEN_MAX_BUFFER_LAYER2  = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OGF%NCG+1) * 2)
         OS%NLEN_MAX_BUFFER_LAYER4  = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OGF%NCG+1) * 4)
         OS%NLEN_MAX_BUFFER_STENCIL = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OGF%NCG+1) * NSCARC_MAX_STENCIL)
         OS%NLEN_MAX_BUFFER_FULL    = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OGF%NCG+1) * NSCARC_MAX_STENCIL * 4)

         OGF%NLEN_BUFFER_LAYER1  = OS%NLEN_MAX_BUFFER_LAYER1
         OGF%NLEN_BUFFER_LAYER2  = OS%NLEN_MAX_BUFFER_LAYER2
         OGF%NLEN_BUFFER_LAYER4  = OS%NLEN_MAX_BUFFER_LAYER4
         OGF%NLEN_BUFFER_STENCIL = OS%NLEN_MAX_BUFFER_STENCIL
         OGF%NLEN_BUFFER_FULL    = OS%NLEN_MAX_BUFFER_FULL

         OGF => OL%STRUCTURED

         OS%NLEN_MAX_BUFFER_LAYER1  = MAX(OS%NLEN_MAX_BUFFER_LAYER1, (OGF%NCG+1) * 1)
         OS%NLEN_MAX_BUFFER_LAYER2  = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OGF%NCG+1) * 2)
         OS%NLEN_MAX_BUFFER_LAYER4  = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OGF%NCG+1) * 4)
         OS%NLEN_MAX_BUFFER_STENCIL = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OGF%NCG+1) * NSCARC_MAX_STENCIL)
         OS%NLEN_MAX_BUFFER_FULL    = MAX(OS%NLEN_MAX_BUFFER_LAYER2, (OGF%NCG+1) * NSCARC_MAX_STENCIL * 4)

         OGF%NLEN_BUFFER_LAYER1  = OS%NLEN_MAX_BUFFER_LAYER1
         OGF%NLEN_BUFFER_LAYER2  = OS%NLEN_MAX_BUFFER_LAYER2
         OGF%NLEN_BUFFER_LAYER4  = OS%NLEN_MAX_BUFFER_LAYER4
         OGF%NLEN_BUFFER_STENCIL = OS%NLEN_MAX_BUFFER_STENCIL
         OGF%NLEN_BUFFER_FULL    = OS%NLEN_MAX_BUFFER_FULL

         NLEN = 2*OS%NLEN_MAX_BUFFER_FULL 
   
         CALL SCARC_ALLOCATE_INT1 (OS%SEND_BUFFER_INT , 1, NLEN, NSCARC_INIT_HUGE, 'OS%SEND_BUFFER_INT', CROUTINE)
         CALL SCARC_ALLOCATE_INT1 (OS%RECV_BUFFER_INT , 1, NLEN, NSCARC_INIT_HUGE, 'OS%RECV_BUFFER_INT', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1 (OS%SEND_BUFFER_REAL, 1, NLEN, NSCARC_INIT_HUGE, 'OS%SEND_BUFFER_REAL', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1 (OS%RECV_BUFFER_REAL, 1, NLEN, NSCARC_INIT_HUGE, 'OS%RECV_BUFFER_REAL', CROUTINE)
      ENDIF

      ! Neighboring wall structures for common wall cells
      ALLOCATE (OG%WALL(OG%NCG), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OG%WALL', IERROR)

      ! In case of GMG with a predefined grid hierarchy allocate corresponding level-structures

      IF (NLEVEL_MAX > NLEVEL_MIN) THEN
         DO NL=NLEVEL_MIN+1,NLEVEL_MAX
            CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
            ALLOCATE (OG%WALL(OG%NCG), STAT=IERROR)
            CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OG%WALL', IERROR)
         ENDDO
      ENDIF

   ENDDO
ENDDO

! If there is more than 1 mesh, initialize communication structures on finest level 
! and setup mapping from local to global numbering
 
IF (NMESHES > 1) THEN

   IF (IS_MGM) CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_STRUCTURED)

   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_NUMBERS, NSCARC_NONE, NLEVEL_MIN)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_SIZES,   NSCARC_NONE, NLEVEL_MIN)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_SOLIDS,       NSCARC_NONE, NLEVEL_MIN)

   IF (HAS_MULTIPLE_LEVELS .AND. .NOT.HAS_AMG_LEVELS) THEN
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX

         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            DO INBR = 1, SCARC(NM)%N_NEIGHBORS
      
               CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)                          
               NOM = S%NEIGHBORS(INBR)
               CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
   
               OG%NLEN_BUFFER_LAYER1  = MAX(OG%NLEN_BUFFER_LAYER1,  (OG%NCG+1) * 1)
               OG%NLEN_BUFFER_LAYER2  = MAX(OG%NLEN_BUFFER_LAYER2,  (OG%NCG+1) * 2)
               OG%NLEN_BUFFER_LAYER4  = MAX(OG%NLEN_BUFFER_LAYER4,  (OG%NCG+1) * 4)
               OG%NLEN_BUFFER_STENCIL = MAX(OG%NLEN_BUFFER_STENCIL, (OG%NCG+1) * NSCARC_MAX_STENCIL)
               OG%NLEN_BUFFER_FULL    = MAX(OG%NLEN_BUFFER_FULL,    (OG%NCG+1) * NSCARC_MAX_STENCIL * 4)
   
            ENDDO
          ENDDO

         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_NUMBERS, NSCARC_NONE, NL)
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_SIZES,   NSCARC_NONE, NL)
      ENDDO
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_EXCHANGES


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform data exchange corresponding to requested exchange type 
! 
! NSCARC_EXCHANGE_AUXILIARY       :  exchange various auxiliary information (AMG only)
! NSCARC_EXCHANGE_BASIC_SIZES     :  exchange initial information about interface sizes
! NSCARC_EXCHANGE_CELL_NEIGHBORS  :  exchange neighboring cell numbers 
! NSCARC_EXCHANGE_CELL_NUMBERS    :  exchange neighboring grid widths 
! NSCARC_EXCHANGE_CELL_SIZES      :  exchange neighboring grid widths 
! NSCARC_EXCHANGE_LAYER2_NUMS     :  exchange number of second layer cells (AMG only)
! NSCARC_EXCHANGE_LAYER2_VALS     :  exchange values of second layer cells (AMG only)
! NSCARC_EXCHANGE_MATRIX_COLS     :  exchange columns of neighboring matrices 
! NSCARC_EXCHANGE_MATRIX_COLSG    :  exchange global columns of neighboring matrices 
! NSCARC_EXCHANGE_MATRIX_DIAGS    :  exchange diagonal entries of neighboring matrices
! NSCARC_EXCHANGE_MATRIX_SIZES    :  exchange sizes of neighboring matrices 
! NSCARC_EXCHANGE_MATRIX_VALS     :  exchange values of neighboring matrices 
! NSCARC_EXCHANGE_MGM_DOUBLE      :  exchange double neighbouring layers for true BC setting of MGM method
! NSCARC_EXCHANGE_MGM_SINGLE      :  exchange single neighbouring layer for all other BC settings of MGM method
! NSCARC_EXCHANGE_MGM_VELO        :  exchange neighboring velocity entries in MGM method
! NSCARC_EXCHANGE_NULLSPACE       :  exchange sum of nullspace entries (AMG only)
! NSCARC_EXCHANGE_PRESSURE        :  exchange pressure values (3d-vector) along internal boundaries 
! NSCARC_EXCHANGE_SOLIDS          :  exchange information about adjacent solids
! NSCARC_EXCHANGE_VECTOR_MEAN     :  exchange vector and build mean of neighboring and own data
! NSCARC_EXCHANGE_VECTOR_PLAIN    :  exchange plain vector and only use neighboring data
! NSCARC_EXCHANGE_ZONE_NEIGHBORS  :  exchange number of aggregation zones (AMG only)
! NSCARC_EXCHANGE_ZONE_TYPES      :  exchange aggregation zone types (AMG only)
! 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE (NTYPE, NPARAM, NL)
USE SCARC_POINTERS, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID
INTEGER, INTENT(IN) :: NTYPE, NPARAM, NL
REAL(EB) :: TNOW
INTEGER :: NM, NOM

N_REQ = 0
TYPE_EXCHANGE = NTYPE

! ---------- Receive data from neighbors 
 
RECEIVE_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)             

   RECEIVE_OMESHES_LOOP: DO NOM = 1, NMESHES

      RNODE = PROCESS(NM)
      SNODE = PROCESS(NOM)

      IF (RNODE==SNODE .OR. .NOT.ARE_FACE_NEIGHBORS(NM, NOM)) CYCLE RECEIVE_OMESHES_LOOP
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

      SELECT_EXCHANGE_TYPE: SELECT CASE (NTYPE)

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_AUXILIARY)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'AUXILIARY')

#endif
         CASE (NSCARC_EXCHANGE_BASIC_SIZES)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'BASIC SIZES')

         CASE (NSCARC_EXCHANGE_CELL_NEIGHBORS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER4, 'CELL NEIGHBORS')

         CASE (NSCARC_EXCHANGE_CELL_NUMBERS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'CELL NUMBERS')

         CASE (NSCARC_EXCHANGE_CELL_SIZES)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'CELL SIZES')

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_LAYER2_NUMS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'LAYER2_NUMS')

         CASE (NSCARC_EXCHANGE_LAYER2_VALS)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'LAYER2_VALS')
#endif

         CASE (NSCARC_EXCHANGE_MATRIX_COLS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_FULL, 'POISSON COLS')

         CASE (NSCARC_EXCHANGE_MATRIX_COLSG)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_FULL, 'POISSON COLSG')

         CASE (NSCARC_EXCHANGE_MATRIX_DIAGS)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'POISSON DIAGS')

         CASE (NSCARC_EXCHANGE_MATRIX_SIZES)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'POISSON SIZES')

         CASE (NSCARC_EXCHANGE_MATRIX_VALS)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_FULL, 'POISSON VALS')

         CASE (NSCARC_EXCHANGE_MGM_DOUBLE)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'MGM_TRUE')

         CASE (NSCARC_EXCHANGE_MGM_SINGLE)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'MGM_MEAN')

         CASE (NSCARC_EXCHANGE_MGM_VELO)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'MGM_VELO')

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_NULLSPACE)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'NULLSPACE')
#endif

         CASE (NSCARC_EXCHANGE_PRESSURE)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'PRESSURE')

         CASE (NSCARC_EXCHANGE_SOLIDS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'SOLIDS')

         CASE (NSCARC_EXCHANGE_VECTOR_MEAN)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'VECTOR MEAN')

         CASE (NSCARC_EXCHANGE_VECTOR_PLAIN)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'VECTOR PLAIN')

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_ZONE_NEIGHBORS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER4, 'ZONE NEIGHBORS')

         CASE (NSCARC_EXCHANGE_ZONE_TYPES)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'ZONE TYPES')
#endif

         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_EXCHANGE_RECV, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SELECT_EXCHANGE_TYPE

   ENDDO RECEIVE_OMESHES_LOOP
ENDDO RECEIVE_MESHES_LOOP

! ---------- Pack data for requested exchange type in corresponding SEND-buffer
  
TNOW = CURRENT_TIME()

SEND_PACK_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)      

   SEND_PACK_OMESHES_LOOP: DO NOM = 1, NMESHES

      IF (.NOT. ARE_FACE_NEIGHBORS(NM, NOM)) CYCLE SEND_PACK_OMESHES_LOOP
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      SEND_PACK_OMESHES_SELECT: SELECT CASE (NTYPE)

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_AUXILIARY)
            CALL SCARC_PACK_AUXILIARY(NL)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'AUXILIARY')
#endif

         CASE (NSCARC_EXCHANGE_BASIC_SIZES)
            CALL SCARC_PACK_BASIC_SIZES
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'BASIC SIZES')

         CASE (NSCARC_EXCHANGE_CELL_NEIGHBORS)
            CALL SCARC_PACK_CELL_NEIGHBORS (NL)
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER4, 'CELL NEIGHBORS')

         CASE (NSCARC_EXCHANGE_CELL_NUMBERS)
            CALL SCARC_PACK_CELL_NUMBERS
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'CELL NUMBERS')

         CASE (NSCARC_EXCHANGE_CELL_SIZES)
            CALL SCARC_PACK_CELL_SIZES
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'CELL SIZES')

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_LAYER2_NUMS)
            CALL SCARC_PACK_LAYER2_NUMS
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'LAYER2_NUMS')

         CASE (NSCARC_EXCHANGE_LAYER2_VALS)
            CALL SCARC_PACK_LAYER2_VALS
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'LAYER2_VALS')
#endif

         CASE (NSCARC_EXCHANGE_MATRIX_SIZES)
            CALL SCARC_PACK_MATRIX_SIZES(NL)
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'POISSON SIZES')

         CASE (NSCARC_EXCHANGE_MATRIX_DIAGS)
            CALL SCARC_PACK_MATRIX_DIAGS(NSCARC_MATRIX_POISSON)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'POISSON DIAGS')

         CASE (NSCARC_EXCHANGE_MATRIX_COLS)
            CALL SCARC_PACK_MATRIX_COLS(NPARAM)
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_FULL, 'POISSON COLS')

         CASE (NSCARC_EXCHANGE_MATRIX_COLSG)
            CALL SCARC_PACK_MATRIX_COLSG(NPARAM)
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_FULL, 'POISSON COLSG')

         CASE (NSCARC_EXCHANGE_MATRIX_VALS)
            CALL SCARC_PACK_MATRIX_VALS(NPARAM, NL)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_FULL, 'POISSON VALS')

         CASE (NSCARC_EXCHANGE_MGM_DOUBLE)
            CALL SCARC_PACK_MGM_DOUBLE(NM)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'MGM_TRUE')

         CASE (NSCARC_EXCHANGE_MGM_SINGLE)
            CALL SCARC_PACK_MGM_SINGLE(NM)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'MGM_MEAN')

         CASE (NSCARC_EXCHANGE_MGM_VELO)
            CALL SCARC_PACK_MGM_VELO(NM)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'MGM_VELO')

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_NULLSPACE)
            CALL SCARC_PACK_NULLSPACE(NL)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'NULLSPACE')
#endif

         CASE (NSCARC_EXCHANGE_PRESSURE)
            CALL SCARC_PACK_PRESSURE(NM)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'PRESSURE')

         CASE (NSCARC_EXCHANGE_SOLIDS)
            CALL SCARC_PACK_SOLIDS
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'SOLIDS')

         CASE (NSCARC_EXCHANGE_VECTOR_MEAN)
            CALL SCARC_PACK_VECTOR_MEAN(NM, NL, NPARAM)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'VECTOR MEAN')

         CASE (NSCARC_EXCHANGE_VECTOR_PLAIN)
            CALL SCARC_PACK_VECTOR_PLAIN(NM, NL, NPARAM)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'VECTOR PLAIN')

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_ZONE_NEIGHBORS)
            CALL SCARC_PACK_ZONE_NEIGHBORS(NL)
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER4, 'ZONE NEIGHBORS')

         CASE (NSCARC_EXCHANGE_ZONE_TYPES)
            CALL SCARC_PACK_ZONE_TYPES
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'ZONE TYPES')
#endif

         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_EXCHANGE_SEND, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SEND_PACK_OMESHES_SELECT
   ENDDO SEND_PACK_OMESHES_LOOP
ENDDO SEND_PACK_MESHES_LOOP


! ---------- Wait for all meshes to have sent and received their data
 
IF (N_MPI_PROCESSES > 1 .AND. N_REQ /= 0) CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),MPI_STATUSES_IGNORE,IERROR)

! ---------- Unpack received data from corresponding RECEIVE-buffers
 
TNOW = CURRENT_TIME()
SEND_UNPACK_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   

   SEND_UNPACK_OMESHES_LOOP: DO NOM = 1, NMESHES

      SNODE  = PROCESS(NM)
      RNODE  = PROCESS(NOM)

      IF (.NOT. ARE_FACE_NEIGHBORS(NM, NOM)) CYCLE SEND_UNPACK_OMESHES_LOOP
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

      SEND_UNPACK_OMESHES_SELECT: SELECT CASE (NTYPE)

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_AUXILIARY)
            CALL SCARC_UNPACK_AUXILIARY (NM, NOM, NL)
#endif

         CASE (NSCARC_EXCHANGE_BASIC_SIZES)
            CALL SCARC_UNPACK_BASIC_SIZES (NM, NOM)

         CASE (NSCARC_EXCHANGE_CELL_NEIGHBORS)
            CALL SCARC_UNPACK_CELL_NEIGHBORS (NM, NOM, NL)

         CASE (NSCARC_EXCHANGE_CELL_NUMBERS)
            CALL SCARC_UNPACK_CELL_NUMBERS (NM, NOM)

         CASE (NSCARC_EXCHANGE_CELL_SIZES)
            CALL SCARC_UNPACK_CELL_SIZES(NM, NOM)

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_LAYER2_NUMS)
            CALL SCARC_UNPACK_LAYER2_NUMS (NM, NOM)

         CASE (NSCARC_EXCHANGE_LAYER2_VALS)
            CALL SCARC_UNPACK_LAYER2_VALS (NM, NOM)
#endif

         CASE (NSCARC_EXCHANGE_MATRIX_COLS)
            CALL SCARC_UNPACK_MATRIX_COLS (NM, NOM, NPARAM)

         CASE (NSCARC_EXCHANGE_MATRIX_COLSG)
            CALL SCARC_UNPACK_MATRIX_COLSG (NM, NOM, NPARAM)

         CASE (NSCARC_EXCHANGE_MATRIX_SIZES)
            CALL SCARC_UNPACK_MATRIX_SIZES (NM, NOM, NL)

         CASE (NSCARC_EXCHANGE_MATRIX_DIAGS)
            CALL SCARC_UNPACK_MATRIX_DIAGS (NM, NOM)

         CASE (NSCARC_EXCHANGE_MATRIX_VALS)
            CALL SCARC_UNPACK_MATRIX_VALS (NM, NOM, NL, NPARAM)

         CASE (NSCARC_EXCHANGE_MGM_SINGLE)
            CALL SCARC_UNPACK_MGM_SINGLE (NM, NOM)

         CASE (NSCARC_EXCHANGE_MGM_DOUBLE)
            CALL SCARC_UNPACK_MGM_DOUBLE (NM, NOM)

         CASE (NSCARC_EXCHANGE_MGM_VELO)
            CALL SCARC_UNPACK_MGM_VELO (NM, NOM)

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_NULLSPACE)
            CALL SCARC_UNPACK_NULLSPACE (NM, NOM, NL)
#endif
         CASE (NSCARC_EXCHANGE_PRESSURE)
            CALL SCARC_UNPACK_PRESSURE (NM, NOM)

         CASE (NSCARC_EXCHANGE_SOLIDS)
            CALL SCARC_UNPACK_SOLIDS (NM, NOM)

         CASE (NSCARC_EXCHANGE_VECTOR_MEAN)
            CALL SCARC_UNPACK_VECTOR_MEAN (NM, NOM, NL, NPARAM)

         CASE (NSCARC_EXCHANGE_VECTOR_PLAIN)
            CALL SCARC_UNPACK_VECTOR_PLAIN (NM, NOM, NL, NPARAM)

#ifdef WITH_SCARC_AMG
         CASE (NSCARC_EXCHANGE_ZONE_NEIGHBORS)
            CALL SCARC_UNPACK_ZONE_NEIGHBORS (NM, NOM, NL)

         CASE (NSCARC_EXCHANGE_ZONE_TYPES)
            CALL SCARC_UNPACK_ZONE_TYPES (NM, NOM)
#endif

         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_EXCHANGE_SEND, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SEND_UNPACK_OMESHES_SELECT
   ENDDO SEND_UNPACK_OMESHES_LOOP
ENDDO SEND_UNPACK_MESHES_LOOP

END SUBROUTINE SCARC_EXCHANGE


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform data exchange including diagonal neighbors
! 
! NSCARC_EXCHANGE_MGM_TRUE   :  exchange information for MGM true BC settings 
! 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_WITH_DIAG (NTYPE, NL)
USE SCARC_POINTERS, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID
INTEGER, INTENT(IN) :: NTYPE, NL
REAL(EB) :: TNOW
INTEGER :: NM, NOM

N_REQ = 0
TYPE_EXCHANGE = NTYPE

! ---------- Receive data from neighbors 
 
RECEIVE_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)             

   RECEIVE_OMESHES_LOOP: DO NOM = 1, NMESHES

      RNODE = PROCESS(NM)
      SNODE = PROCESS(NOM)

      IF (RNODE==SNODE .OR. .NOT. ARE_FACE_NEIGHBORS(NM, NOM)) CYCLE RECEIVE_OMESHES_LOOP
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

      SELECT_EXCHANGE_TYPE: SELECT CASE (NTYPE)

         CASE (NSCARC_EXCHANGE_MGM_TRUE)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'BASIC SIZES')

         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_EXCHANGE_DIAG, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SELECT_EXCHANGE_TYPE

   ENDDO RECEIVE_OMESHES_LOOP
ENDDO RECEIVE_MESHES_LOOP

! ---------- Pack data for requested exchange type in corresponding SEND-buffer
  
TNOW = CURRENT_TIME()

SEND_PACK_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)      

   SEND_PACK_OMESHES_LOOP: DO NOM = 1, NMESHES

      IF (.NOT. ARE_FACE_NEIGHBORS(NM, NOM)) CYCLE SEND_PACK_OMESHES_LOOP
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      SEND_PACK_OMESHES_SELECT: SELECT CASE (NTYPE)

         CASE (NSCARC_EXCHANGE_MGM_TRUE)
            CALL SCARC_PACK_MGM_TRUE (NM)
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'BASIC SIZES')

         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_EXCHANGE_DIAG, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SEND_PACK_OMESHES_SELECT
   ENDDO SEND_PACK_OMESHES_LOOP
ENDDO SEND_PACK_MESHES_LOOP


! ---------- Wait for all meshes to have sent and received their data
 
IF (N_MPI_PROCESSES > 1 .AND. N_REQ /= 0) CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),MPI_STATUSES_IGNORE,IERROR)

! ---------- Unpack received data from corresponding RECEIVE-buffers
 
TNOW = CURRENT_TIME()
SEND_UNPACK_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   

   SEND_UNPACK_OMESHES_LOOP: DO NOM = 1, NMESHES

      SNODE  = PROCESS(NM)
      RNODE  = PROCESS(NOM)

      IF (.NOT. ARE_FACE_NEIGHBORS(NM, NOM)) CYCLE SEND_UNPACK_OMESHES_LOOP
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

      SEND_UNPACK_OMESHES_SELECT: SELECT CASE (NTYPE)

         CASE (NSCARC_EXCHANGE_MGM_TRUE)
            CALL SCARC_UNPACK_MGM_TRUE (NM, NOM)

         CASE DEFAULT
            CALL SCARC_ERROR(NSCARC_ERROR_EXCHANGE_DIAG, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SEND_UNPACK_OMESHES_SELECT
   ENDDO SEND_UNPACK_OMESHES_LOOP
ENDDO SEND_UNPACK_MESHES_LOOP

END SUBROUTINE SCARC_EXCHANGE_WITH_DIAG


! --------------------------------------------------------------------------------------------------------------
!> \brief Receive data of type integer
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RECV_MESSAGE_INT(NM, NOM, NL, NTYPE, CTEXT)
USE SCARC_POINTERS, ONLY: OS, OG
INTEGER, INTENT(IN) :: NM, NOM, NL, NTYPE
CHARACTER(*), INTENT(IN) :: CTEXT
!INTEGER,  DIMENSION(:), POINTER :: RECV_BUFFER_INT
INTEGER,  POINTER :: RECV_BUFFER_INT
INTEGER :: NLEN, NDUMMY
CHARACTER(40) :: CDUMMY

NDUMMY = NM; NDUMMY = NOM; NDUMMY = NL; CDUMMY=CTEXT     ! prevent compilation warning if VERBOSE flag is not set (temporarily)

IF (RNODE == SNODE) RETURN
N_REQ = N_REQ+1

SELECT CASE(NTYPE)
   CASE (NSCARC_BUFFER_BASIC)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT0(1);  NLEN = NSCARC_MAX_BUFFER0
   CASE (NSCARC_BUFFER_FULL)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_FULL
   CASE (NSCARC_BUFFER_LAYER1)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER1
   CASE (NSCARC_BUFFER_LAYER2)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER2
   CASE (NSCARC_BUFFER_LAYER4)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER4
   CASE (NSCARC_BUFFER_STENCIL)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_STENCIL
END SELECT


RECV_BUFFER_INT = NSCARC_HUGE_INT
CALL MPI_IRECV(RECV_BUFFER_INT, NLEN, MPI_INTEGER, SNODE, TAG, MPI_COMM_WORLD, REQ(N_REQ),IERROR)

END SUBROUTINE SCARC_RECV_MESSAGE_INT


! --------------------------------------------------------------------------------------------------------------
!> \brief Receive data of type real
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RECV_MESSAGE_REAL(NM, NOM, NL, NTYPE, CTEXT)
USE SCARC_POINTERS, ONLY: OS, OG
INTEGER, INTENT(IN) :: NM, NOM, NL, NTYPE
CHARACTER(*), INTENT(IN) :: CTEXT
!REAL(EB), DIMENSION(:), POINTER :: RECV_BUFFER_REAL
REAL(EB), POINTER :: RECV_BUFFER_REAL
INTEGER :: NLEN, NDUMMY
CHARACTER(40) :: CDUMMY

NDUMMY = NM; NDUMMY = NOM; NDUMMY = NL; CDUMMY=CTEXT     ! prevent compilation warning if VERBOSE flag is not set (temporarily)

IF (RNODE == SNODE) RETURN
N_REQ = N_REQ+1

SELECT CASE(NTYPE)
   CASE (NSCARC_BUFFER_BASIC)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL0(1);  NLEN = NSCARC_MAX_BUFFER0
   CASE (NSCARC_BUFFER_FULL)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_FULL
   CASE (NSCARC_BUFFER_LAYER1)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER1
   CASE (NSCARC_BUFFER_LAYER2)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER2
   CASE (NSCARC_BUFFER_LAYER4)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER4
   CASE (NSCARC_BUFFER_STENCIL)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_STENCIL
END SELECT


RECV_BUFFER_REAL = NSCARC_INIT_UNDEF
CALL MPI_IRECV(RECV_BUFFER_REAL, NLEN, MPI_DOUBLE_PRECISION, SNODE, TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)

END SUBROUTINE SCARC_RECV_MESSAGE_REAL


! --------------------------------------------------------------------------------------------------------------
!> \brief Send data of integer type
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SEND_MESSAGE_INT(NM, NOM, NL, NTYPE, CTEXT)
USE SCARC_POINTERS, ONLY: OS, OG
CHARACTER(*), INTENT(IN) :: CTEXT
INTEGER, INTENT(IN) :: NM, NOM, NL, NTYPE
!INTEGER,  DIMENSION(:), POINTER :: SEND_BUFFER_INT
INTEGER,  POINTER :: SEND_BUFFER_INT
INTEGER :: NLEN, NDUMMY
CHARACTER(40) :: CDUMMY

NDUMMY = NM; NDUMMY = NOM; NDUMMY = NL; CDUMMY = CTEXT    ! prevent compilation warning if VERBOSE flag is not set (temporarily)

IF (RNODE == SNODE) RETURN
N_REQ = N_REQ+1

SELECT CASE(NTYPE)
   CASE (NSCARC_BUFFER_BASIC)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT0(1);  NLEN = NSCARC_MAX_BUFFER0
   CASE (NSCARC_BUFFER_FULL)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_FULL
   CASE (NSCARC_BUFFER_LAYER1)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER1
   CASE (NSCARC_BUFFER_LAYER2)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER2
   CASE (NSCARC_BUFFER_LAYER4)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER4
   CASE (NSCARC_BUFFER_STENCIL)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_STENCIL
END SELECT

CALL MPI_ISEND(SEND_BUFFER_INT, NLEN, MPI_INTEGER, SNODE, TAG, MPI_COMM_WORLD, REQ(N_REQ),IERROR)

END SUBROUTINE SCARC_SEND_MESSAGE_INT


! --------------------------------------------------------------------------------------------------------------
!> \brief Send data of real type
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SEND_MESSAGE_REAL(NM, NOM, NL, NTYPE, CTEXT)
USE SCARC_POINTERS, ONLY: OS, OG
CHARACTER(*), INTENT(IN) :: CTEXT
INTEGER, INTENT(IN) :: NM, NOM, NL, NTYPE
!REAL(EB), DIMENSION(:), POINTER :: SEND_BUFFER_REAL
REAL(EB), POINTER :: SEND_BUFFER_REAL
INTEGER :: NLEN, NDUMMY
CHARACTER(40) :: CDUMMY

NDUMMY = NM; NDUMMY = NOM; NDUMMY = NL; CDUMMY = CTEXT    ! prevent compilation warning if VERBOSE flag is not set (temporarily)

IF (RNODE == SNODE) RETURN
N_REQ = N_REQ+1

SELECT CASE(NTYPE)
   CASE (NSCARC_BUFFER_BASIC)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL0(1);  NLEN = NSCARC_MAX_BUFFER0
   CASE (NSCARC_BUFFER_FULL)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_FULL
   CASE (NSCARC_BUFFER_LAYER1)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER1
   CASE (NSCARC_BUFFER_LAYER2)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER2
   CASE (NSCARC_BUFFER_LAYER4)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER4
   CASE (NSCARC_BUFFER_STENCIL)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_STENCIL
END SELECT

CALL MPI_ISEND(SEND_BUFFER_REAL, NLEN, MPI_DOUBLE_PRECISION, SNODE, TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)

END SUBROUTINE SCARC_SEND_MESSAGE_REAL


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack initial exchange sizes along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_BASIC_SIZES
USE SCARC_POINTERS, ONLY: OS, OG
OS%SEND_BUFFER_INT0(1)=OG%NCG
OS%SEND_BUFFER_INT0(2)=OG%NZG
END SUBROUTINE SCARC_PACK_BASIC_SIZES


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack initial exchange sizes along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_BASIC_SIZES (NM, NOM)
USE SCARC_POINTERS, ONLY: OG, RECV_BUFFER_INT, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM
RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 0)
OG%NCG = RECV_BUFFER_INT(1)
OG%NZG = RECV_BUFFER_INT(2)
END SUBROUTINE SCARC_UNPACK_BASIC_SIZES


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack numbers of cells which are overlapped by neighbor
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_CELL_NUMBERS
USE SCARC_POINTERS, ONLY: G, OS, OL, OG
INTEGER :: IOR0, ICG, IWG, IXW, IYW, IZW

OS%SEND_BUFFER_INT = NSCARC_HUGE_INT
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      IXW = G%WALL(IWG)%IXW
      IYW = G%WALL(IWG)%IYW
      IZW = G%WALL(IWG)%IZW
      OS%SEND_BUFFER_INT(ICG) = G%CELL_NUMBER(IXW, IYW, IZW)
   ENDDO
ENDDO
END SUBROUTINE SCARC_PACK_CELL_NUMBERS


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack numbers of cells which are overlapped by neighbor
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_CELL_NUMBERS (NM, NOM)
USE SCARC_POINTERS, ONLY: L, G, OL, OG, RECV_BUFFER_INT, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: LL, IOR0, ICG, ICE, IWG, IXG, IYG, IZG

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      G%ICE_TO_ICN(ICE) = RECV_BUFFER_INT(LL)
      IWG = OG%ICG_TO_IWG(ICG)
      IXG = G%WALL(IWG)%IXG
      IYG = G%WALL(IWG)%IYG
      IZG = G%WALL(IWG)%IZG
      IF (RECV_BUFFER_INT(ICG) < 0) THEN
         L%IS_SOLID(IXG, IYG, IZG) = .TRUE.
         G%CELL_NUMBER(IXG, IYG, IZG) = -G%CELL_NUMBER(IXG, IYG, IZG)     ! mark solid cell with negative sign
      ELSE
         L%IS_SOLID(IXG, IYG, IZG) = .FALSE.
      ENDIF
      LL = LL + 1
   ENDDO
ENDDO
END SUBROUTINE SCARC_UNPACK_CELL_NUMBERS


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack cell width information 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_CELL_SIZES
USE SCARC_POINTERS, ONLY: L, OS

OS%SEND_BUFFER_REAL0(1) = L%DXL(0)
OS%SEND_BUFFER_REAL0(2) = L%DXL(L%NX)
OS%SEND_BUFFER_REAL0(3) = L%DYL(0)
OS%SEND_BUFFER_REAL0(4) = L%DYL(L%NY)
OS%SEND_BUFFER_REAL0(5) = L%DZL(0)
OS%SEND_BUFFER_REAL0(6) = L%DZL(L%NZ)

END SUBROUTINE SCARC_PACK_CELL_SIZES


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack cell width information 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_CELL_SIZES (NM, NOM)
USE SCARC_POINTERS, ONLY: L, OL, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 0)

IF (OL%GHOST_LASTW(-1) /= 0) L%DXL(0)    = 0.5_EB*(RECV_BUFFER_REAL(1) + L%DXL(0))
IF (OL%GHOST_LASTW( 1) /= 0) L%DXL(L%NX) = 0.5_EB*(RECV_BUFFER_REAL(2) + L%DXL(L%NX))
IF (OL%GHOST_LASTW(-2) /= 0) L%DYL(0)    = 0.5_EB*(RECV_BUFFER_REAL(3) + L%DYL(0))
IF (OL%GHOST_LASTW( 2) /= 0) L%DYL(L%NY) = 0.5_EB*(RECV_BUFFER_REAL(4) + L%DYL(L%NY))
IF (OL%GHOST_LASTW(-3) /= 0) L%DZL(0)    = 0.5_EB*(RECV_BUFFER_REAL(5) + L%DZL(0))
IF (OL%GHOST_LASTW( 3) /= 0) L%DZL(L%NZ) = 0.5_EB*(RECV_BUFFER_REAL(6) + L%DZL(L%NZ))

END SUBROUTINE SCARC_UNPACK_CELL_SIZES


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack zones numbers along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_CELL_NEIGHBORS(NL)
USE SCARC_POINTERS, ONLY: L, G, OS, OL, OG, A, F, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: IOR0, ICG, ICW, ICWG, LL, IWG, IXW, IYW, IZW

A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
OS%SEND_BUFFER_INT = NSCARC_HUGE_INT
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW  = OG%ICG_TO_ICW(ICG, 1)                                 ! local cell number adjacent to interface
      ICWG = G%LOCAL_TO_GLOBAL(ICW)                                ! global cell number adjacent to interface
      OS%SEND_BUFFER_INT(LL)   = ICW
      OS%SEND_BUFFER_INT(LL+1) = ICWG
      LL = LL + 2
   ENDDO
   IF (NL /= NLEVEL_MIN) CYCLE
   F => L%FACE(IOR0)
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      IXW = G%WALL(IWG)%IXW + F%INCRX
      IYW = G%WALL(IWG)%IYW + F%INCRY
      IZW = G%WALL(IWG)%IZW + F%INCRZ
      ICW = G%CELL_NUMBER(IXW, IYW, IZW)
      OG%ICG_TO_ICW(ICG, 2) = ICW                              
      ICWG = G%LOCAL_TO_GLOBAL(ICW)
      OS%SEND_BUFFER_INT(LL)   = ICW
      OS%SEND_BUFFER_INT(LL+1) = ICWG
      LL = LL + 2
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_CELL_NEIGHBORS


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack zones numbers along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_CELL_NEIGHBORS(NM, NOM, NL)
USE SCARC_POINTERS, ONLY: G, OL, OG, RECV_BUFFER_INT, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM, NL
INTEGER :: ICG, IOR0, LL, IOFF

CROUTINE = 'SCARC_UNPACK_CELL_NEIGHBORS'

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)

CALL SCARC_ALLOCATE_INT1 (OG%ICG_TO_OCELL, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'G%ICG_TO_OCELL', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OG%ICG_TO_GCELL, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'G%ICG_TO_GCELL', CROUTINE)

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      OG%ICG_TO_OCELL(ICG) = RECV_BUFFER_INT(LL)
      OG%ICG_TO_GCELL(ICG) = RECV_BUFFER_INT(LL+1)
      LL = LL + 2
   ENDDO
   IF (NL /= NLEVEL_MIN) CYCLE
   IOFF = OL%GHOST_LASTW(IOR0)
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      G%ICE2 = G%ICE2 + 1
      OG%ICG_TO_OCELL(ICG + IOFF) = RECV_BUFFER_INT(LL)
      OG%ICG_TO_GCELL(ICG + IOFF) = RECV_BUFFER_INT(LL+1)
      OG%ICG_TO_ICE(ICG, 2) = G%ICE2
      LL = LL + 2
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_CELL_NEIGHBORS


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping information about matrix columns (compact storage technique only)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MATRIX_COLS(NM, NOM, NMATRIX)
USE SCARC_POINTERS, ONLY: OL, OG, OAC, RECV_BUFFER_INT, SCARC_POINT_TO_OTHER_CMATRIX, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM, NMATRIX
INTEGER :: IOR0, ICG, LL, ICP

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)
OAC => SCARC_POINT_TO_OTHER_CMATRIX (NMATRIX)

LL = 1                                 
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTE(IOR0) == 0) CYCLE
   ICP = OL%GHOST_FIRSTE(IOR0)
   OAC%ROW(ICP) = LL
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IF (OG%ICG_TO_ICE(ICG,1) < 0) CYCLE                          ! skip solid cells
      OAC%COL(LL) = ABS(RECV_BUFFER_INT(LL))
      DO WHILE (RECV_BUFFER_INT(LL+1) >= 0)
         LL = LL + 1
         OAC%COL(LL) = ABS(RECV_BUFFER_INT(LL))
      ENDDO
      LL = LL + 1
      ICP = ICP + 1
      OAC%ROW(ICP) = LL
   ENDDO
   OAC%N_ROW = ICP  
   OAC%N_VAL = LL - 1
ENDDO

END SUBROUTINE SCARC_UNPACK_MATRIX_COLS


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping information about matrix columns (compact storage technique only)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MATRIX_COLSG(NMATRIX)                
USE SCARC_POINTERS, ONLY: OS, OL, OG, AC, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NMATRIX
INTEGER :: IOR0, ICG, ICW, LL, ICOL
INTEGER, POINTER, DIMENSION(:) :: COLG

AC => SCARC_POINT_TO_CMATRIX (NMATRIX)
IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) THEN
   COLG => AC%COL
ELSE
   COLG => AC%COLG
ENDIF

LL = 1
OS%SEND_BUFFER_INT = NSCARC_HUGE_INT

DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE

   DO ICG= OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW < 0) CYCLE                               ! skip solid cells
      ICOL = AC%ROW(ICW)
      OS%SEND_BUFFER_INT(LL) = -COLG(ICOL)          ! send first element with negative sign (thus, mark beginning)
      LL = LL + 1                              
      DO ICOL = AC%ROW(ICW)+1, AC%ROW(ICW+1)-1
         OS%SEND_BUFFER_INT(LL) = COLG(ICOL)   
         LL = LL + 1
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_MATRIX_COLSG


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping information about matrix columns (compact storage technique only)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MATRIX_COLSG(NM, NOM, NMATRIX)
USE SCARC_POINTERS, ONLY: OL, OG, OAC, RECV_BUFFER_INT, SCARC_POINT_TO_OTHER_CMATRIX, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM, NMATRIX
INTEGER :: IOR0, ICG, LL, ICP
INTEGER, POINTER, DIMENSION(:) :: COLG

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)
OAC => SCARC_POINT_TO_OTHER_CMATRIX (NMATRIX)
IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) THEN
   COLG => OAC%COL
ELSE
   COLG => OAC%COLG
ENDIF

LL = 1                                 
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTE(IOR0) == 0) CYCLE
   ICP = OL%GHOST_FIRSTE(IOR0)
   OAC%ROW(ICP) = LL
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IF (OG%ICG_TO_ICE(ICG,1) < 0) CYCLE                     ! skip solid cells
      COLG(LL) = ABS(RECV_BUFFER_INT(LL))
      DO WHILE (RECV_BUFFER_INT(LL+1) > 0)
         LL = LL + 1
         COLG(LL) = ABS(RECV_BUFFER_INT(LL))
      ENDDO
      LL = LL + 1
      ICP = ICP + 1
      OAC%ROW(ICP) = LL
   ENDDO
   OAC%N_ROW = ICP  
   OAC%N_VAL = LL - 1
ENDDO

END SUBROUTINE SCARC_UNPACK_MATRIX_COLSG


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping information about matrix values (both storage techniques)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MATRIX_VALS(NMATRIX, NL)
USE SCARC_POINTERS, ONLY: AB, AC, OS, OL, OG, SCARC_POINT_TO_BMATRIX, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NMATRIX, NL
INTEGER :: IOR0, ICG, ICW, LL, ID, ICOL

LL = 1
OS%SEND_BUFFER_INT = NSCARC_ZERO_REAL_EB
SELECT CASE (SET_MATRIX_TYPE(NL))

   CASE (NSCARC_MATRIX_BANDWISE)                

      AB => SCARC_POINT_TO_BMATRIX (NMATRIX)
      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            DO ID = 1, AB%N_STENCIL
               OS%SEND_BUFFER_REAL(LL) = AB%VAL(ICW, ID)
               LL = LL + 1
            ENDDO
         ENDDO
      ENDDO

   CASE (NSCARC_MATRIX_COMPACT)

      AC => SCARC_POINT_TO_CMATRIX (NMATRIX)
      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         DO ICG= OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            IF (ICW < 0) CYCLE                                ! skip solid cells
            DO ICOL = AC%ROW(ICW), AC%ROW(ICW+1)-1
               OS%SEND_BUFFER_REAL(LL) = AC%VAL(ICOL)
               LL = LL + 1
            ENDDO
         ENDDO
      ENDDO

END SELECT

END SUBROUTINE SCARC_PACK_MATRIX_VALS


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping information about matrix values (both storage techniques)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MATRIX_VALS(NM, NOM, NL, NMATRIX)
USE SCARC_POINTERS, ONLY: OL, OG, OAB, OAC, RECV_BUFFER_REAL, &
                          SCARC_POINT_TO_OTHER_BMATRIX, SCARC_POINT_TO_OTHER_CMATRIX, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM, NL, NMATRIX
INTEGER :: IOR0, ICG, ICOL, ID, LL

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)

SELECT CASE (SET_MATRIX_TYPE(NL))

   ! Bandwise matrix on level NL
   CASE (NSCARC_MATRIX_BANDWISE)              

      OAB => SCARC_POINT_TO_OTHER_BMATRIX (NMATRIX)
      LL = 1
      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
            DO ID = 1, OAB%N_STENCIL
               OAB%VAL(ICG, ID) = RECV_BUFFER_REAL(LL)
               LL = LL + 1
            ENDDO
         ENDDO
      ENDDO

   ! Compact matrix on level NL
   CASE (NSCARC_MATRIX_COMPACT)

      OAC => SCARC_POINT_TO_OTHER_CMATRIX (NMATRIX)
      LL = 1
      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTE(IOR0) == 0) CYCLE
         DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
            IF (OG%ICG_TO_ICE(ICG,1) < 0) CYCLE                     ! skip solid cells
            DO ICOL = OAC%ROW(ICG), OAC%ROW(ICG+1)-1
               OAC%VAL(ICOL) = RECV_BUFFER_REAL(LL)
               LL = LL + 1
            ENDDO
         ENDDO
      ENDDO

END SELECT

END SUBROUTINE SCARC_UNPACK_MATRIX_VALS


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack information about matrix sizes into send vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MATRIX_SIZES(NL)
USE SCARC_POINTERS, ONLY: OS, OG
INTEGER, INTENT(IN) :: NL

SELECT CASE (SET_MATRIX_TYPE(NL))
   CASE (NSCARC_MATRIX_BANDWISE)
      OS%SEND_BUFFER_INT0(1) = OG%POISSONB%N_VAL
      OS%SEND_BUFFER_INT0(2) = OG%POISSONB%N_DIAG
      OS%SEND_BUFFER_INT0(3) = OG%POISSONB%N_STENCIL
   CASE (NSCARC_MATRIX_COMPACT)
      OS%SEND_BUFFER_INT0(1) = OG%POISSON%N_VAL
      OS%SEND_BUFFER_INT0(2) = OG%POISSON%N_ROW
      OS%SEND_BUFFER_INT0(3) = OG%POISSON%N_STENCIL
END SELECT

END SUBROUTINE SCARC_PACK_MATRIX_SIZES
   

! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack information about matrix sizes into send vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MATRIX_SIZES(NM, NOM, NL)
USE SCARC_POINTERS, ONLY: OG, RECV_BUFFER_INT, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM, NL

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 0)
SELECT CASE (SET_MATRIX_TYPE(NL))
   CASE (NSCARC_MATRIX_BANDWISE)
      OG%POISSONB%N_VAL     = RECV_BUFFER_INT(1)
      OG%POISSONB%N_DIAG    = RECV_BUFFER_INT(2)
      OG%POISSONB%N_STENCIL = RECV_BUFFER_INT(3)
   CASE (NSCARC_MATRIX_COMPACT)
      OG%POISSON%N_VAL     = RECV_BUFFER_INT(1)
      OG%POISSON%N_ROW     = RECV_BUFFER_INT(2)
      OG%POISSON%N_STENCIL = RECV_BUFFER_INT(3)
END SELECT

END SUBROUTINE SCARC_UNPACK_MATRIX_SIZES


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping information about matrix diagonals (compact storage technique only)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MATRIX_DIAGS(NTYPE)
USE SCARC_POINTERS, ONLY: AC, OS, OL, OG, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NTYPE
INTEGER :: IOR0, ICG, ICW, ICOL

AC => SCARC_POINT_TO_CMATRIX (NTYPE)

OS%SEND_BUFFER_REAL = NSCARC_ZERO_REAL_EB
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW < 0) CYCLE                                   ! skip solid cells
      ICOL = AC%ROW(ICW)
      OS%SEND_BUFFER_REAL(ICG) = AC%VAL(ICOL)
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_MATRIX_DIAGS
      

! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping information about matrix diagonals (compact storage technique only)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MATRIX_DIAGS(NM, NOM)
USE SCARC_POINTERS, ONLY: G, OL, OG, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: IOR0, ICG, ICE, LL

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      IF (ICE < 0) CYCLE                                   ! skip solid cells
      G%DIAG(ICE) = RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_MATRIX_DIAGS


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping parts of specified pressure vector (predictor/corrector)
! Note: Vector VC is numbered via I, J, K values
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_PRESSURE(NM)
USE SCARC_POINTERS, ONLY: G, OS, OL, OG, OS, SCARC_POINT_TO_HVECTOR
INTEGER, INTENT(IN) :: NM
REAL(EB), DIMENSION(:,:,:), POINTER :: VC
INTEGER :: IOR0, ICG, IWG

IF (PREDICTOR) THEN
   VC => SCARC_POINT_TO_HVECTOR (NM, NSCARC_VECTOR_H)
ELSE
   VC => SCARC_POINT_TO_HVECTOR (NM, NSCARC_VECTOR_HS)
ENDIF
OS%SEND_BUFFER_REAL = NSCARC_HUGE_REAL_EB
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      OS%SEND_BUFFER_REAL(ICG) = VC(G%WALL(IWG)%IXW, G%WALL(IWG)%IYW, G%WALL(IWG)%IZW)
   ENDDO
ENDDO
END SUBROUTINE SCARC_PACK_PRESSURE


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping parts of specified pressure vector (predictor/corrector)
! Note: Vector VC is numbered via I, J, K values
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_PRESSURE(NM, NOM)
USE SCARC_POINTERS, ONLY: G, OL, OG, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL, SCARC_POINT_TO_HVECTOR
INTEGER, INTENT(IN) :: NM, NOM
REAL(EB), DIMENSION(:,:,:), POINTER :: VC
INTEGER :: LL, IOR0, IWG, IXG, IYG, IZG, ICG

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
IF (PREDICTOR) THEN
   VC => SCARC_POINT_TO_HVECTOR (NM, NSCARC_VECTOR_H)
ELSE
   VC => SCARC_POINT_TO_HVECTOR (NM, NSCARC_VECTOR_HS)
ENDIF
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   UNPACK_PRESSURE: DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      IXG=G%WALL(IWG)%IXG
      IYG=G%WALL(IWG)%IYG
      IZG=G%WALL(IWG)%IZG
      VC(IXG, IYG, IZG) = RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO UNPACK_PRESSURE
ENDDO
END SUBROUTINE SCARC_UNPACK_PRESSURE


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack single layer of interface boundary cells for MGM mean and extrapolation boundary settings
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MGM_SINGLE(NM)
USE SCARC_POINTERS, ONLY: G, OL, OG, OS, UHL
INTEGER, INTENT(IN) :: NM
INTEGER :: IOR0, ICG, ICW, IWG, IXW, IYW, IZW

UHL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%UHL
OS%SEND_BUFFER_REAL = NSCARC_ZERO_REAL_EB

DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   PACK_MGM_SINGLE: DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW < 0) CYCLE                                  ! skip solid cells
      IWG = OG%ICG_TO_IWG(ICG)
      IXW=G%WALL(IWG)%IXW
      IYW=G%WALL(IWG)%IYW
      IZW=G%WALL(IWG)%IZW
      OS%SEND_BUFFER_REAL(ICG) = UHL(IXW, IYW, IZW)
   ENDDO PACK_MGM_SINGLE
ENDDO

END SUBROUTINE SCARC_PACK_MGM_SINGLE


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack single layer of interface boundary cells for MGM mean and extrapolation boundary settings
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MGM_SINGLE(NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, OUHL, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: IOR0, ICG, IWG, LL

OUHL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OUHL
RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   UNPACK_MGM_SINGLE: DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      OUHL(IWG) = RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO UNPACK_MGM_SINGLE
ENDDO

END SUBROUTINE SCARC_UNPACK_MGM_SINGLE

! --------------------------------------------------------------------------------------------------------------
!> \brief Pack single layer of interface boundary cells for MGM mean and extrapolation boundary settings
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MGM_TRUE (NM)
USE SCARC_POINTERS, ONLY: G, OL, OG, OS, UHL
INTEGER, INTENT(IN) :: NM
INTEGER :: IOR0, ICG, ICW, IWG, IXW, IYW, IZW

UHL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%UHL
OS%SEND_BUFFER_REAL = NSCARC_ZERO_REAL_EB

DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   PACK_MGM_TRUE: DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW < 0) CYCLE                                  ! skip solid cells
      IWG = OG%ICG_TO_IWG(ICG)
      IXW=G%WALL(IWG)%IXW
      IYW=G%WALL(IWG)%IYW
      IZW=G%WALL(IWG)%IZW
      OS%SEND_BUFFER_REAL(ICG) = UHL(IXW, IYW, IZW)
   ENDDO PACK_MGM_TRUE
ENDDO

END SUBROUTINE SCARC_PACK_MGM_TRUE


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack single layer of interface boundary cells for MGM mean and extrapolation boundary settings
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MGM_TRUE(NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, OUHL, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: IOR0, ICG, IWG, LL

OUHL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OUHL
RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   UNPACK_MGM_TRUE: DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      OUHL(IWG) = RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO UNPACK_MGM_TRUE
ENDDO

END SUBROUTINE SCARC_UNPACK_MGM_TRUE


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack double layer of interface boundary cells for MGM true approximate boundary settings
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MGM_DOUBLE(NM)
USE SCARC_POINTERS, ONLY: G, OL, OG, OS
INTEGER, INTENT(IN) :: NM
REAL(EB), DIMENSION(:,:,:), POINTER :: UL
INTEGER :: IOR0, ICG, ICW, IWG, IXW, IYW, IZW, LL

UL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%UHL
OS%SEND_BUFFER_REAL = NSCARC_ZERO_REAL_EB

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   PACK_MGM_DOUBLE: DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW < 0) CYCLE                                  ! skip solid cells
      IWG = OG%ICG_TO_IWG(ICG)
      IXW=G%WALL(IWG)%IXW
      IYW=G%WALL(IWG)%IYW
      IZW=G%WALL(IWG)%IZW
      OS%SEND_BUFFER_REAL(LL) = UL(IXW, IYW, IZW)
      SELECT CASE (IOR0)
         CASE (1)
            OS%SEND_BUFFER_REAL(LL+1) = UL(IXW+1, IYW, IZW)
         CASE (-1)
            OS%SEND_BUFFER_REAL(LL+1) = UL(IXW-1, IYW, IZW)
         CASE ( 2)
            OS%SEND_BUFFER_REAL(LL+1) = UL(IXW, IYW+1, IZW)
         CASE (-2)
            OS%SEND_BUFFER_REAL(LL+1) = UL(IXW, IYW-1, IZW)
         CASE ( 3)
            OS%SEND_BUFFER_REAL(LL+1) = UL(IXW, IYW, IZW+1)
         CASE (-3)
            OS%SEND_BUFFER_REAL(LL+1) = UL(IXW, IYW, IZW-1)
      END SELECT
      LL = LL + 2
   ENDDO PACK_MGM_DOUBLE
ENDDO

END SUBROUTINE SCARC_PACK_MGM_DOUBLE


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack double layer of interface boundary cells for MGM true approximate boundary setting
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MGM_DOUBLE(NM, NOM)
USE SCARC_POINTERS, ONLY: 
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_REAL, OUHL, OUHL2, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: LL, IOR0, ICG, IWG

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
OUHL  => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OUHL
OUHL2 => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OUHL2
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   UNPACK_MGM_DOUBLE: DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      OUHL(IWG)  = RECV_BUFFER_REAL(LL)
      OUHL2(IWG) = RECV_BUFFER_REAL(LL+1)
      LL = LL + 2
   ENDDO UNPACK_MGM_DOUBLE
ENDDO
END SUBROUTINE SCARC_UNPACK_MGM_DOUBLE


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping parts of specified pressure vector (predictor/corrector)
! Note: Vector VC is numbered via I, J, K values
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MGM_VELO(NM)
USE SCARC_POINTERS, ONLY: G, OL, OG, OS, UU, VV, WW
INTEGER, INTENT(IN) :: NM
INTEGER :: IOR0, ICG, IWG, IXW, IYW, IZW, LL

OS%SEND_BUFFER_REAL = NSCARC_HUGE_REAL_EB

UU => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%UVEL
VV => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%VVEL
WW => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%WVEL

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      IXW=G%WALL(IWG)%IXW
      IYW=G%WALL(IWG)%IYW
      IZW=G%WALL(IWG)%IZW
      SELECT CASE (IOR0)
         CASE ( 1)
            OS%SEND_BUFFER_REAL(LL) = UU(IXW-1, IYW, IZW) 
         CASE (-1)
            OS%SEND_BUFFER_REAL(LL) = UU(IXW, IYW, IZW) 
         CASE ( 2)
            OS%SEND_BUFFER_REAL(LL) = VV(IXW, IYW-1, IZW) 
         CASE (-2)
            OS%SEND_BUFFER_REAL(LL) = VV(IXW, IYW, IZW) 
         CASE ( 3)
            OS%SEND_BUFFER_REAL(LL) = WW(IXW, IYW, IZW-1) 
         CASE (-3)
            OS%SEND_BUFFER_REAL(LL) = WW(IXW, IYW, IZW) 
      END SELECT
      LL = LL + 1
   ENDDO
ENDDO
END SUBROUTINE SCARC_PACK_MGM_VELO


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping parts of specified pressure vector (predictor/corrector)
! Note: Vector VC is numbered via I, J, K values
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MGM_VELO2(NM)
USE SCARC_POINTERS, ONLY: G, OL, OG, OS, UU, VV, WW
INTEGER, INTENT(IN) :: NM
REAL(EB), DIMENSION(:,:,:), POINTER :: H3
INTEGER :: IOR0, ICG, IWG, IXW, IYW, IZW, LL

OS%SEND_BUFFER_REAL = NSCARC_HUGE_REAL_EB

H3 => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%UIP
UU => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%UVEL
VV => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%VVEL
WW => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%WVEL

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      IXW=G%WALL(IWG)%IXW
      IYW=G%WALL(IWG)%IYW
      IZW=G%WALL(IWG)%IZW
      SELECT CASE (IOR0)
         CASE ( 1)
            OS%SEND_BUFFER_REAL(LL)   = H3(IXW, IYW, IZW) - H3(IXW-1, IYW, IZW)
            OS%SEND_BUFFER_REAL(LL+1) = UU(IXW-1, IYW, IZW) 
         CASE (-1)
            OS%SEND_BUFFER_REAL(LL)   = H3(IXW+1, IYW, IZW) - H3(IXW, IYW, IZW)
            OS%SEND_BUFFER_REAL(LL+1) = UU(IXW, IYW, IZW) 
         CASE ( 2)
            OS%SEND_BUFFER_REAL(LL)   = H3(IXW, IYW, IZW) - H3(IXW, IYW-1, IZW)
            OS%SEND_BUFFER_REAL(LL+1) = VV(IXW, IYW-1, IZW) 
         CASE (-2)
            OS%SEND_BUFFER_REAL(LL)   = H3(IXW, IYW+1, IZW) - H3(IXW, IYW, IZW)
            OS%SEND_BUFFER_REAL(LL+1) = VV(IXW, IYW, IZW) 
         CASE ( 3)
            OS%SEND_BUFFER_REAL(LL)   = H3(IXW, IYW, IZW) - H3(IXW, IYW, IZW-1)
            OS%SEND_BUFFER_REAL(LL+1) = WW(IXW, IYW, IZW-1) 
         CASE (-3)
            OS%SEND_BUFFER_REAL(LL)   = H3(IXW, IYW, IZW+1) - H3(IXW, IYW, IZW)
            OS%SEND_BUFFER_REAL(LL+1) = WW(IXW, IYW, IZW) 
      END SELECT
      LL = LL + 2
   ENDDO
ENDDO
END SUBROUTINE SCARC_PACK_MGM_VELO2


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping parts of specified pressure vector (predictor/corrector)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MGM_VELO(NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, OVEL, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: LL, IOR0, ICG, IWG

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   UNPACK_MGM_VELO: DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      SELECT CASE (ABS(IOR0))
         CASE (1)
            OVEL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OUVEL
         CASE (2)
            OVEL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OVVEL
         CASE (3)
            OVEL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OWVEL
      END SELECT
      OVEL(IWG) = RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO UNPACK_MGM_VELO
ENDDO
END SUBROUTINE SCARC_UNPACK_MGM_VELO


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping parts of specified pressure vector (predictor/corrector)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_MGM_VELO2(NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, OUIP, OVEL, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: LL, IOR0, ICG, IWG

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
OUIP => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OUIP

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   UNPACK_MGM_VELO: DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      OUIP(IWG) = RECV_BUFFER_REAL(LL)
      SELECT CASE (ABS(IOR0))
         CASE (1)
            OVEL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OUVEL
         CASE (2)
            OVEL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OVVEL
         CASE (3)
            OVEL => SCARC(NM)%LEVEL(NLEVEL_MIN)%MGM%OWVEL
      END SELECT
      OVEL(IWG) = RECV_BUFFER_REAL(LL+1)
      LL = LL + 2
   ENDDO UNPACK_MGM_VELO
ENDDO
END SUBROUTINE SCARC_UNPACK_MGM_VELO2


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping parts of specified vector VC (numbered via IC values)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_VECTOR_PLAIN(NM, NL, NV)
USE SCARC_POINTERS, ONLY: OL, OG, OS, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NM, NL, NV
REAL(EB), DIMENSION(:), POINTER :: VC
INTEGER :: IOR0, ICG, ICW

VC => SCARC_POINT_TO_VECTOR(NM, NL, NV)
OS%SEND_BUFFER_REAL = NSCARC_ZERO_REAL_EB

DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW < 0) CYCLE                                  ! skip solid cells
      OS%SEND_BUFFER_REAL(ICG) = VC(ICW)
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_VECTOR_PLAIN


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping parts of specified vector VC (numbered via IC values)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_VECTOR_PLAIN(NM, NOM, NL, NVECTOR)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NM, NOM, NL, NVECTOR
REAL(EB), DIMENSION(:), POINTER :: VC
INTEGER :: IOR0, LL, ICG, ICE

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
VC => SCARC_POINT_TO_VECTOR(NM, NL, NVECTOR)

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG,1)
      IF (ICE < 0) CYCLE                            ! skip solid cells
      VC(ICE) = RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_VECTOR_PLAIN


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping parts of solid vector IS_SOLID
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_SOLIDS
USE SCARC_POINTERS, ONLY: L, G, OL, OG, OS
INTEGER :: IOR0, ICG, IWG, IXW, IYW, IZW

OS%SEND_BUFFER_INT = NSCARC_ZERO_INT        

DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      IXW = G%WALL(IWG)%IXW 
      IYW = G%WALL(IWG)%IYW 
      IZW = G%WALL(IWG)%IZW 
      IF (L%IS_SOLID(IXW, IYW, IZW)) OS%SEND_BUFFER_INT(ICG) = 1
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_SOLIDS


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping parts of specified vector VC (numbered via IC values)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_SOLIDS(NM, NOM)
USE SCARC_POINTERS, ONLY: L, G, OL, OG, RECV_BUFFER_INT, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: IOR0, LL, ICG, IWG, IXG, IYG, IZG

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      IXG = G%WALL(IWG)%IXG
      IYG = G%WALL(IWG)%IYG
      IZG = G%WALL(IWG)%IZG
      IF (RECV_BUFFER_INT(LL) == 1) THEN
         L%IS_SOLID(IXG, IYG, IZG) = .TRUE.
      ELSE
         L%IS_SOLID(IXG, IYG, IZG) = .FALSE.
      ENDIF
      LL = LL + 1
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_SOLIDS


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping and internal parts of specified vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_VECTOR_MEAN(NM, NL, NVECTOR)
USE SCARC_POINTERS, ONLY: OL, OG, OS, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NM, NL, NVECTOR
REAL(EB), DIMENSION(:), POINTER :: VC
INTEGER :: IOR0, ICG, ICW, ICE, LL

VC => SCARC_POINT_TO_VECTOR(NM, NL, NVECTOR)
LL = 1
OS%SEND_BUFFER_REAL = NSCARC_HUGE_REAL_EB
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW < 0) CYCLE                                          ! skip solid internal cells
      OS%SEND_BUFFER_REAL(LL) = VC(ICW)
      LL = LL + 1
   ENDDO
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      IF (ICE < 0) CYCLE                                          ! skip solid external cells
      OS%SEND_BUFFER_REAL(LL) = VC(ICE)
      LL = LL + 1
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_VECTOR_MEAN


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping and internal parts of specified vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_VECTOR_MEAN(NM, NOM, NL, NVECTOR)
USE SCARC_POINTERS, ONLY: G, OG, OL, RECV_BUFFER_REAL, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM, NL, NVECTOR
REAL(EB), DIMENSION(:), POINTER :: VC
INTEGER :: IOR0, LL, ICG, ICW, ICE

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
VC => SCARC_POINT_TO_VECTOR(NM, NL, NVECTOR)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICW = G%ICG_TO_ICW(ICG,1)
      IF (ICW < 0) CYCLE                                                       ! skip internal solid cells
      VC(ICW) = 2.0_EB/3.0_EB * VC(ICW) + 1.0_EB/3.0_EB * RECV_BUFFER_REAL(LL)
      !VC(ICW) = 0.5_EB * VC(ICW) + 0.5_EB * RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG,1)
      IF (ICE < 0) CYCLE                                                       ! skip external solid cells
      VC(ICE) = 1.0_EB/3.0_EB * VC(ICE) + 2.0_EB/3.0_EB * RECV_BUFFER_REAL(LL)
      !VC(ICE) = 0.5_EB * VC(ICE) + 0.5_EB * RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_VECTOR_MEAN


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping information about matrix columns (compact storage technique only)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_MATRIX_COLS(NMATRIX)                
USE SCARC_POINTERS, ONLY: OS, OL, OG, AC, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NMATRIX
INTEGER :: IOR0, LL, ICOL, ICG, ICW

AC => SCARC_POINT_TO_CMATRIX (NMATRIX)

LL = 1
OS%SEND_BUFFER_INT = NSCARC_HUGE_INT

DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE

   DO ICG= OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW < 0) CYCLE                               ! skip solid cells
      ICOL = AC%ROW(ICW)
      OS%SEND_BUFFER_INT(LL) = -AC%COL(ICOL)           ! send first element with negative sign (thus, mark beginning)
      LL = LL + 1                            
      DO ICOL = AC%ROW(ICW)+1, AC%ROW(ICW+1)-1
         OS%SEND_BUFFER_INT(LL) = AC%COL(ICOL)
         LL = LL + 1
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_MATRIX_COLS


#ifdef WITH_SCARC_AMG
! --------------------------------------------------------------------------------------------------------------
!> \brief Pack overlapping auxiliary vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_AUXILIARY(NL)
USE SCARC_POINTERS, ONLY: L, G, OL, OG, OS, F
INTEGER, INTENT(IN) :: NL
INTEGER :: IOR0, ICG, ICW1, ICW2, IWG, IXW, IYW, IZW, LL

LL = 1
OS%SEND_BUFFER_REAL = NSCARC_HUGE_REAL_EB
DO IOR0 = -3, 3

   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   F => L%FACE(IOR0)
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

      ICW1 = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW1 > 0) OS%SEND_BUFFER_REAL(LL) = G%AUX1(ICW1)
      LL = LL + 1

      IF (NL /= NLEVEL_MIN) CYCLE

      IWG  = OG%ICG_TO_IWG(ICG)
      IXW  = G%WALL(IWG)%IXW + F%INCRX
      IYW  = G%WALL(IWG)%IYW + F%INCRY
      IZW  = G%WALL(IWG)%IZW + F%INCRZ

      ICW2 = G%CELL_NUMBER(IXW, IYW, IZW)
      IF (ICW2 > 0) OS%SEND_BUFFER_REAL(LL+1) = G%AUX1(ICW2)
      LL = LL + 1

   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_AUXILIARY


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping auxiliary vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_AUXILIARY (NM, NOM, NL)
USE SCARC_POINTERS, ONLY: G, OL, OG, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM, NL
INTEGER :: IOR0, ICG, ICE1, ICE2, LL

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)

      ICE1 = OG%ICG_TO_ICE(ICG,1)
      IF (ICE1 > 0) G%AUX1(ICE1) = RECV_BUFFER_REAL(LL)
      LL = LL + 1

      IF (NL /= NLEVEL_MIN) CYCLE

      ICE2 = OG%ICG_TO_ICE(ICG,2)
      IF (ICE2 > 0) G%AUX1(ICE2) = RECV_BUFFER_REAL(LL+1)
      LL = LL + 1

   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_AUXILIARY


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack and unpack overlapping nullspace vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_NULLSPACE(NL)
USE SCARC_POINTERS, ONLY: L, G, F, OL, OG, OS
INTEGER, INTENT(IN) :: NL
INTEGER :: IOR0, ICG, ICW1, ICW2, IWG, IXW, IYW, IZW, LL

LL = 1
OS%SEND_BUFFER_REAL = NSCARC_HUGE_REAL_EB
DO IOR0 = -3, 3

   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   F => L%FACE(IOR0)
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

      ICW1 = OG%ICG_TO_ICW(ICG, 1)
      IF (ICW1 > 0) OS%SEND_BUFFER_REAL(LL)   = G%NULLSPACE(ICW1)
      LL = LL + 1

      IF (NL /= NLEVEL_MIN) CYCLE

      IWG  = OG%ICG_TO_IWG(ICG)
      IXW  = G%WALL(IWG)%IXW + F%INCRX
      IYW  = G%WALL(IWG)%IYW + F%INCRY
      IZW  = G%WALL(IWG)%IZW + F%INCRZ

      ICW2 = G%CELL_NUMBER(IXW, IYW, IZW)
      IF (ICW2 > 0) OS%SEND_BUFFER_REAL(LL+1) = G%NULLSPACE(ICW2)
      LL = LL + 1

   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_NULLSPACE


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack overlapping nullspace vector 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_NULLSPACE (NM, NOM, NL)
USE SCARC_POINTERS, ONLY: G, OL, OG, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM, NL
INTEGER :: IOR0, ICG, ICE1, ICE2, LL

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)

      ICE1 = OG%ICG_TO_ICE(ICG,1)
      IF (ICE1 > 0) G%NULLSPACE(ICE1) = RECV_BUFFER_REAL(LL)
      LL = LL + 1

      IF (NL /= NLEVEL_MIN) CYCLE

      ICE2 = OG%ICG_TO_ICE(ICG,2)
      IF (ICE2 > 0) G%NULLSPACE(ICE2) = RECV_BUFFER_REAL(LL+1)
      LL = LL + 1

   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_NULLSPACE

! --------------------------------------------------------------------------------------------------------------
!> \brief Pack zones numbers along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_ZONE_NEIGHBORS(NL)
USE SCARC_POINTERS, ONLY: L, G, OS, OL, OG, A, F, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: IOR0, ICG, ICW, LL, IWG, IXW, IYW, IZW, IZWG

A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
OS%SEND_BUFFER_INT = NSCARC_HUGE_INT
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW  = OG%ICG_TO_ICW(ICG, 1)
      IZW  = G%ZONES_LOCAL(ICW)                                ! local zone number adjacent to interface
      IZWG = G%ZONES_GLOBAL(ICW)                               ! global zone number adjacent to interface
      OS%SEND_BUFFER_INT(LL  ) = IZW
      OS%SEND_BUFFER_INT(LL+1) = IZWG
      LL = LL + 2
   ENDDO
   IF (NL /= NLEVEL_MIN) CYCLE
   F => L%FACE(IOR0)
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG  = OG%ICG_TO_IWG(ICG)
      IXW  = G%WALL(IWG)%IXW + F%INCRX
      IYW  = G%WALL(IWG)%IYW + F%INCRY
      IZW  = G%WALL(IWG)%IZW + F%INCRZ
      ICW  = G%CELL_NUMBER(IXW, IYW, IZW)
      IZW  = G%ZONES_LOCAL(ICW)                       
      IZWG = G%ZONES_GLOBAL(ICW)                      
      OS%SEND_BUFFER_INT(LL  ) = IZW
      OS%SEND_BUFFER_INT(LL+1) = IZWG
      LL = LL + 2
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_ZONE_NEIGHBORS


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack zones numbers along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_ZONE_NEIGHBORS(NM, NOM, NL)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_INT, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM, NL
INTEGER :: ICG, IOR0, LL, KK

CROUTINE = 'SCARC_UNPACK_ZONE_NEIGHBORS'
RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)

CALL SCARC_ALLOCATE_INT1 (OG%ICG_TO_OZONE, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'G%ICG_TO_OZONE', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OG%ICG_TO_GZONE, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'G%ICG_TO_GZONE', CROUTINE)

LL = 1
KK = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      OG%ICG_TO_OZONE(KK) = RECV_BUFFER_INT(LL)
      OG%ICG_TO_GZONE(KK) = RECV_BUFFER_INT(LL+1)
      KK = KK + 1
      LL = LL + 2
   ENDDO
   IF (NL /= NLEVEL_MIN) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      OG%ICG_TO_OZONE(KK) = RECV_BUFFER_INT(LL)
      OG%ICG_TO_GZONE(KK) = RECV_BUFFER_INT(LL+1)
      KK = KK + 1
      LL = LL + 2
   ENDDO
ENDDO
END SUBROUTINE SCARC_UNPACK_ZONE_NEIGHBORS


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack zones numbers along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_LAYER2_NUMS
USE SCARC_POINTERS, ONLY: OS, OL, OG, L, G
INTEGER :: IOR0, ICG, INUM

CALL SCARC_ALLOCATE_INT1 (G%ELAYER2_NUMS, 1, G%NCE2 - G%NC + 1, NSCARC_INIT_ZERO, 'G%ELAYER2_NUMS', CROUTINE)

OS%SEND_BUFFER_INT = NSCARC_HUGE_INT
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      INUM = OG%ICG_TO_ELAYER2(ICG)
      OS%SEND_BUFFER_INT(ICG) = INUM
      IF (INUM /= 0) THEN 
         L%N_LAYER2_TOTAL = L%N_LAYER2_TOTAL + 1
         G%ELAYER2_NUMS(L%N_LAYER2_TOTAL) = OG%ICG_TO_ELAYER2(ICG)
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_LAYER2_NUMS


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack zones numbers along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_LAYER2_NUMS(NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_INT, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: ICG, IOR0

CROUTINE = 'SCARC_UNPACK_LAYER2_NUMS'
RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)

CALL SCARC_ALLOCATE_INT1 (OG%ICG_TO_ILAYER2, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'OG%ICG_TO_ILAYER2', CROUTINE)

DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      OG%ICG_TO_ILAYER2(ICG) = RECV_BUFFER_INT(ICG)
   ENDDO
ENDDO
END SUBROUTINE SCARC_UNPACK_LAYER2_NUMS


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack zones numbers along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_LAYER2_VALS
USE SCARC_POINTERS, ONLY: G, OS, OL, OG
INTEGER :: IOR0, ICG, IC

OS%SEND_BUFFER_REAL = NSCARC_HUGE_INT
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IC = FINDLOC (G%LOCAL_TO_GLOBAL(1:G%NC), VALUE = OG%ICG_TO_ILAYER2(ICG), DIM = 1)
      IF (IC /= 0) OS%SEND_BUFFER_REAL(ICG) = G%NULLSPACE(IC)
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_LAYER2_VALS


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack zones numbers along interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_LAYER2_VALS(NM, NOM)
USE SCARC_POINTERS, ONLY: L, OL, G, OG, RECV_BUFFER_REAL, SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: ICG, IOR0

CROUTINE = 'SCARC_UNPACK_LAYER2_VALS'
RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)

CALL SCARC_ALLOCATE_REAL1 (G%ELAYER2_VALS, 1, L%N_LAYER2_TOTAL, NSCARC_INIT_ZERO, 'G%ELAYER2_VALS', CROUTINE)

DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IF (OG%ICG_TO_ELAYER2(ICG) /= 0) THEN
         L%L2PTR = L%L2PTR + 1
         G%ELAYER2_VALS(L%L2PTR) = RECV_BUFFER_REAL(ICG)
      ENDIF
   ENDDO
ENDDO
END SUBROUTINE SCARC_UNPACK_LAYER2_VALS


! --------------------------------------------------------------------------------------------------------------
!> \brief Pack zones information into send vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_ZONE_TYPES
USE SCARC_POINTERS, ONLY: G, OS, OL, OG
INTEGER :: IOR0, ICG, ICW, ICE, LL

LL = 1
OS%SEND_BUFFER_INT = NSCARC_HUGE_INT
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE

   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (G%ZONES_LOCAL(ICW) /= 0) THEN
         OS%SEND_BUFFER_INT(LL)   = G%ZONES_LOCAL(ICW)
      ELSE
         OS%SEND_BUFFER_INT(LL)   = 0
      ENDIF
      LL = LL + 1
   ENDDO
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      IF (G%ZONES_LOCAL(ICW) /= 0) THEN
         OS%SEND_BUFFER_INT(LL) = G%ZONES_LOCAL(ICE)
      ELSE
         OS%SEND_BUFFER_INT(LL) = 0
      ENDIF
      LL = LL + 1
   ENDDO
ENDDO
END SUBROUTINE SCARC_PACK_ZONE_TYPES


! --------------------------------------------------------------------------------------------------------------
!> \brief Unpack zones information into send vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_ZONE_TYPES(NM, NOM)
USE SCARC_POINTERS, ONLY: G, OL, OG, RECV_BUFFER_INT, SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: IOR0, LL, ICG, ICW, ICE

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      IF (G%ZONES_LOCAL(ICE) == 0) G%ZONES_LOCAL(ICE) = RECV_BUFFER_INT(LL)
      LL = LL + 1
   ENDDO
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IF (G%ZONES_LOCAL(ICW) == 0) G%ZONES_LOCAL(ICW) = RECV_BUFFER_INT(LL)
      LL = LL + 1
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_ZONE_TYPES
#endif

END MODULE SCARC_MPI


#ifdef WITH_MKL
!=======================================================================================================================
!
! Module SCARC_MKL
!
!> \brief Setup environment necessary for the call of the IntelMKL local and global LU-solvers
!   PARDISO and CLUSTER_SPARSE_SOLVER
!
!=======================================================================================================================
MODULE SCARC_MKL
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MKL_PARDISO
USE MKL_CLUSTER_SPARSE_SOLVER
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_TROUBLESHOOTING, ONLY: SCARC_ERROR
USE SCARC_STACK, ONLY: STACK, SCARC_SETUP_STACK, SCARC_SETUP_STACK_VECTORS
USE SCARC_STORAGE, ONLY: SCARC_ALLOCATE_INT1, SCARC_ALLOCATE_INT2, SCARC_ALLOCATE_INT3

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for LU-solvers (based on MKL)
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MKL(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX
 
! Basic setup of stack information and types for MKL
 
CALL SCARC_SETUP_STACK(NSTACK)

SV  => STACK(NSTACK)%SOLVER
SELECT CASE (NSOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      SV%CNAME = 'SCARC_MAIN_SOLVER'
   CASE (NSCARC_SOLVER_COARSE)
      SV%CNAME = 'SCARC_COARSE_LU'
   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)
END SELECT

! Reset types for LU-decomposition method
 
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

! Point to solution vectors (in corresponding scope)
 
CALL SCARC_SETUP_STACK_VECTORS(.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_MKL


! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for LU-solvers (based on MKL)
!  In the multi-mesh case use CLUSTER_SPARSE_SOLVER, else PARDISO solver (only on finest grid level)
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MKL_ENVIRONMENT

STACK(NSCARC_STACK_ROOT)%SOLVER => POISSON_SOLVER
CALL SCARC_SETUP_MKL(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSCARC_STACK_ROOT, NLEVEL_MIN, NLEVEL_MIN)

IF (NMESHES > 1) THEN 
   CALL SCARC_SETUP_CLUSTER (NLEVEL_MIN, NLEVEL_MAX)
ELSE 
   CALL SCARC_SETUP_PARDISO (NLEVEL_MIN, NLEVEL_MAX)
ENDIF

END SUBROUTINE SCARC_SETUP_MKL_ENVIRONMENT


! --------------------------------------------------------------------------------------------------------------
!> \brief Initialize CLUSTER_SPARSE_SOLVER from MKL-library
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CLUSTER(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO NL = NLMIN, NLMAX
      CALL SCARC_SETUP_CLUSTER_MESH(NM, NL)
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_CLUSTER


! --------------------------------------------------------------------------------------------------------------
!> \brief Initialize PARDISO from MKL-library
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PARDISO(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO NL = NLMIN, NLMAX
      CALL SCARC_SETUP_PARDISO_MESH(NM, NL)
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_PARDISO


! --------------------------------------------------------------------------------------------------------------
!> \brief Initialize CLUSTER_SPARSE_SOLVER from MKL-library on local mesh
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CLUSTER_MESH(NM, NL)
USE SCARC_POINTERS, ONLY: L, G, AS, MKL, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I 
REAL (EB) :: TNOW
REAL (EB) :: DUMMY(1)=0.0_EB
REAL (FB) :: DUMMY_FB(1)=0.0_FB

TNOW = CURRENT_TIME()
CROUTINE = 'SCARC_SETUP_CLUSTER'

CALL SCARC_POINT_TO_GRID (NM, NL)                                    

MKL => L%MKL
AS  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_SYM)

! Allocate workspace for parameters and pointers needed in MKL-routine
 
CALL SCARC_ALLOCATE_INT1 (MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL%IPARM', CROUTINE)

IF (.NOT.ALLOCATED(MKL%CT)) THEN
   ALLOCATE(MKL%CT(64), STAT=IERROR)
   CALL CHKMEMERR ('SCARC', 'CT', IERROR)
   DO I=1,64
      MKL%CT(I)%DUMMY = 0
   ENDDO
ENDIF

! Define corresponding parameters
! Note: IPARM-vectory is allocate from 1:64, not from 0:63
 
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
!MKL%IPARM(41) = G%NC_OFFSET(NM)+G%NC_LOCAL(NM)          ! last global cell number for mesh NM
 
! First perform only reordering and symbolic factorization
! Then perform only factorization
 
IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

   MKL%IPARM(28) = 1         ! single precision
   MKL%PHASE = 11
   CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MPI_COMM_WORLD, MKL%ERROR)
   MKL%PHASE = 22
   CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MPI_COMM_WORLD, MKL%ERROR)

   IF (MKL%ERROR /= 0) THEN
      WRITE(*,*) 'ERROR in MKL SETUP, MKL%ERROR=', MKL%ERROR
      CALL MPI_FINALIZE(IERROR)
      STOP
   ENDIF

ELSE

   MKL%IPARM(28) = 0         ! double precision
   MKL%PHASE = 11
   CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
   MKL%PHASE = 22
   CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
   IF (MKL%ERROR /= 0) THEN
      WRITE(*,*) 'ERROR in MKL SETUP, MKL%ERROR=', MKL%ERROR
      CALL MPI_FINALIZE(IERROR)
      STOP
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_CLUSTER_MESH


! --------------------------------------------------------------------------------------------------------------
!> \brief Initialize PARDISO solver from MKL-library on local mesh
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PARDISO_MESH(NM, NL)
USE SCARC_POINTERS, ONLY: L, G, AS, MKL, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, IDUMMY(1)=0
REAL (EB) :: TNOW
REAL (EB) :: DUMMY(1)=0.0_EB
REAL (FB) :: DUMMY_FB(1)=0.0_FB

TNOW = CURRENT_TIME()
CROUTINE = 'SCARC_SETUP_PARDISO'

CALL SCARC_POINT_TO_GRID (NM, NL)                                    

MKL => L%MKL
AS  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_SYM)

! Allocate workspace for parameters and pointers needed in MKL-routine
 
CALL SCARC_ALLOCATE_INT1 (MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL%IPARM', CROUTINE)

IF (.NOT.ALLOCATED(MKL%PT)) THEN
   ALLOCATE(MKL%PT(64), STAT=IERROR)
   CALL CHKMEMERR ('SCARC', 'PT', IERROR)
   DO I=1,64
      MKL%PT(I)%DUMMY = 0
   ENDDO
ENDIF

! Define corresponding parameters
! Note: IPARM-vectory is allocate from 1:64, not from 0:63
 
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

! First perform only reordering and symbolic factorization
! Then perform only factorization

IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
   MKL%IPARM(28) = 1         ! single precision
   MKL%PHASE = 11
   CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                  AS%VAL_FB, AS%ROW, AS%COL, IDUMMY, &
                  MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)
   MKL%PHASE = 22
   CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                  AS%VAL_FB, AS%ROW, AS%COL, IDUMMY, &
                  MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)
ELSE
   MKL%IPARM(28) = 0         ! double precision
   MKL%PHASE = 11
   CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                  AS%VAL, AS%ROW, AS%COL, IDUMMY, &
                  MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
   MKL%PHASE = 22
   CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                  AS%VAL, AS%ROW, AS%COL, IDUMMY, &
                  MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
ENDIF

END SUBROUTINE SCARC_SETUP_PARDISO_MESH


! --------------------------------------------------------------------------------------------------------------
!> \brief Initialize PARDISO solver from MKL-library
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_PARDISO(NM, NL)
USE SCARC_POINTERS, ONLY: L, G, AS, MKL, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, IDUMMY(1)=0
REAL (EB) :: TNOW
REAL (EB) :: DUMMY(1)=0.0_EB
REAL (FB) :: DUMMY_FB(1)=0.0_FB

TNOW = CURRENT_TIME()
CROUTINE = 'SCARC_SETUP_PARDISO'

CALL SCARC_POINT_TO_GRID (NM, NL)                                    

MKL => L%MKL
AS  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE_SYM)

! Allocate workspace for parameters and pointers needed in MKL-routine
 
CALL SCARC_ALLOCATE_INT1 (MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL%IPARM', CROUTINE)

IF (.NOT.ALLOCATED(MKL%PT)) THEN
   ALLOCATE(MKL%PT(64), STAT=IERROR)
   CALL CHKMEMERR ('SCARC', 'PT', IERROR)
   DO I=1,64
      MKL%PT(I)%DUMMY = 0
   ENDDO
ENDIF

! Define corresponding parameters
! Note: IPARM-vectory is allocate from 1:64, not from 0:63
 
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

! First perform only reordering and symbolic factorization
! Then perform only factorization

IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
   MKL%IPARM(28) = 1         ! single precision
   MKL%PHASE = 11
   CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                  AS%VAL_FB, AS%ROW, AS%COL, IDUMMY, &
                  MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)
   MKL%PHASE = 22
   CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                  AS%VAL_FB, AS%ROW, AS%COL, IDUMMY, &
                  MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)
ELSE
   MKL%IPARM(28) = 0         ! double precision
   MKL%PHASE = 11
   CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                  AS%VAL, AS%ROW, AS%COL, IDUMMY, &
                  MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
   MKL%PHASE = 22
   CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                  AS%VAL, AS%ROW, AS%COL, IDUMMY, &
                  MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
ENDIF
END SUBROUTINE SCARC_SETUP_MGM_PARDISO

END MODULE SCARC_MKL

#endif
!=======================================================================================================================
!
! MODULE SCARC_VECTORS
!
!> \brief Define a set of Linear Algebra operations based on vectors and matrices
!
!=======================================================================================================================
MODULE SCARC_VECTORS
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE MESH_VARIABLES, ONLY: MESHES
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_UTILITIES, ONLY: SET_MATRIX_TYPE
USE SCARC_CPU, ONLY: SCARC_SETUP_CPU
USE SCARC_MPI, ONLY: SCARC_SETUP_GLOBALS, SCARC_SETUP_EXCHANGES, SCARC_EXCHANGE

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief Vector multiplied with a constant scalar is added to another vector 
!     DY(I) = DA * DX(I) + DY(I) 
! --------------------------------------------------------------------------------------------------------------
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

! --------------------------------------------------------------------------------------------------------------
!> \brief Vector multiplied with a constant scalar is added to vector multiplied with another scalar
!     DY(I) = DA1 * DX(I) + DA2 * DY(I) 
! --------------------------------------------------------------------------------------------------------------
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

! --------------------------------------------------------------------------------------------------------------
!> \brief Vector multiplied with variable scalars (componentwise) is added to another vector 
!     DY(I) = DA(I)*DX(I) + DY(I)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DAXPY_VARIABLE(N, DA, DX, DY)
REAL(EB), INTENT(IN), DIMENSION(:):: DA, DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) = DY(I) + DA(I) * DX(I)
ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE SCARC_DAXPY_VARIABLE

! --------------------------------------------------------------------------------------------------------------
!> \brief Vector is multiplied with a constant scalar 
!     DY(I) = DA(I)*DX(I) 
! --------------------------------------------------------------------------------------------------------------
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

! --------------------------------------------------------------------------------------------------------------
!> \brief Vector is multiplied with variable scalars (componentwise)
!     DY(I) = DA(I)*DX(I) 
! --------------------------------------------------------------------------------------------------------------
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

! --------------------------------------------------------------------------------------------------------------
!> \brief Compute global scalar-product including global data exchange
! --------------------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_SCALAR_PRODUCT(NV1, NV2, NL)
USE SCARC_POINTERS, ONLY: G, V1, V2, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
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

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

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

! Compute global scalar product as sum of local scalar products
 
IF (N_MPI_PROCESSES>1) & 
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_REAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERROR)

SCARC_SCALAR_PRODUCT = RANK_REAL

CPU(MY_RANK)%SCALAR_PRODUCT = CPU(MY_RANK)%SCALAR_PRODUCT + CURRENT_TIME()-TNOW
END FUNCTION SCARC_SCALAR_PRODUCT

! --------------------------------------------------------------------------------------------------------------
!> \brief Compute global L2-norm including global data exchange
! --------------------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_L2NORM(NV1, NL)
INTEGER, INTENT(IN) :: NV1, NL
REAL(EB) :: TNOW
TNOW = CURRENT_TIME()

GLOBAL_REAL = SCARC_SCALAR_PRODUCT(NV1, NV1, NL)
GLOBAL_REAL = SQRT (GLOBAL_REAL)

SCARC_L2NORM = GLOBAL_REAL

CPU(MY_RANK)%L2NORM =CPU(MY_RANK)%L2NORM+CURRENT_TIME()-TNOW
END FUNCTION SCARC_L2NORM

! --------------------------------------------------------------------------------------------------------------
!> \brief Compute linear combination of two vectors for bandwise storage technique
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_SUM(NV1, NV2, SCAL1, SCAL2, NL)
USE SCARC_POINTERS, ONLY: G, V1, V2, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB), INTENT(IN) :: SCAL1, SCAL2
INTEGER :: NM
#ifdef WITH_MKL
EXTERNAL :: DAXPBY
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

   V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
   V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

#ifdef WITH_MKL
   CALL DAXPBY(G%NCE, SCAL1, V1, 1, SCAL2, V2, 1)
#else
   CALL SCARC_DAXPY_CONSTANT_DOUBLE(G%NCE, SCAL1, V1, SCAL2, V2)
#endif

ENDDO

END SUBROUTINE SCARC_VECTOR_SUM

! --------------------------------------------------------------------------------------------------------------
!> \brief Define vector2 to be a scaled copy of vector 1
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_COPY(NV1, NV2, SCAL1, NL)
USE SCARC_POINTERS, ONLY: G, V1, V2, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB), INTENT(IN) :: SCAL1
INTEGER :: NM
#ifdef WITH_MKL
EXTERNAL :: DCOPY, DSCAL
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

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

! --------------------------------------------------------------------------------------------------------------
!> \brief Clear vector
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_CLEAR(NV, NL)
USE SCARC_POINTERS, ONLY: VC, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   VC => SCARC_POINT_TO_VECTOR(NM, NL, NV)
   VC =  0.0_EB
ENDDO

END SUBROUTINE SCARC_VECTOR_CLEAR

! --------------------------------------------------------------------------------------------------------------
!> \brief Preset vector with specified value
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_RANDOM_INIT (NV, NL)
USE SCARC_POINTERS, ONLY: L, G, VC, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: IC, NM, I, J, K
REAL (EB) :: VAL

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

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

! --------------------------------------------------------------------------------------------------------------
!> \brief Preset vector with specified value
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_INIT (NV, VAL, NL)
USE SCARC_POINTERS, ONLY: L, G, VC, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV, NL
REAL (EB), INTENT(IN) :: VAL
INTEGER :: IC, NM, I, J, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

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

! --------------------------------------------------------------------------------------------------------------
!> \brief Set exact solution according to specified function
! --------------------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION EXACT(X,Z)
REAL (EB), INTENT(IN) :: X, Z
!EXACT = (X**2 - X**4) * (Z**4 - Z**2)                              
!EXACT = (X**2 - 1) * (Z**2 - 1)                                    
!EXACT =  625.0_EB/16.0_EB * X * (0.8_EB - X) * Z * (0.8_EB - Z)    
EXACT = - X * (0.8_EB - X) * Z * (0.8_EB - Z)        ! FUNCTION 3
END FUNCTION EXACT

! --------------------------------------------------------------------------------------------------------------
!> \brief Set right hand side according to specified function
! --------------------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION RHS(X,Z)
REAL (EB), INTENT(IN) :: X, Z
!RHS = 2.0_EB*((1.0_EB - 6.0_EB*X**2)*Z**2*(1.0_EB-Z**2)+(1.0_EB-6.0_EB*Z**2)*X**2*(1.0_EB-X**2))
!RHS = -X**2 - Z**2 +2
!RHS = 625.0_EB/8.0_EB * (X * (0.8_EB - X) + Z * (0.8_EB - Z))
RHS = 2.0_EB * (X * (0.8_EB - X) + Z * (0.8_EB - Z))
END FUNCTION RHS

! --------------------------------------------------------------------------------------------------------------
!> \brief Preset right hand side in such a way that exact solution is known
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_EXACT (NE, NL)
USE SCARC_POINTERS, ONLY: M, L, G, VC, XMID, ZMID, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
USE SCARC_CONVERGENCE, ONLY: ITE_TOTAL
INTEGER, INTENT(IN) :: NE, NL
INTEGER :: IC, NM, I, K

IF (ITE_TOTAL == 0) WRITE(*,*) 'TODO: PRESET_EXACT is active !!!'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                
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
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_EXACT

! --------------------------------------------------------------------------------------------------------------
!> \brief Preset vector with specific values
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_VECTOR (NV, NL)
USE SCARC_POINTERS, ONLY: M, L, G, VC, XMID, ZMID, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: IC, NM, I, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   VC = 0.0_EB
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
         IF (XMID(I) < 0.1_EB) THEN
            VC(IC) = 1000.0_EB
         ELSE
            VC(IC) = 0.0_EB
         ENDIF
 ! WRITE(MSG%LU_DEBUG,*) 'IC=',IC,':X=',XMID,':Z=',ZMID,': RHS=',VC(IC)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_VECTOR

! --------------------------------------------------------------------------------------------------------------
!> \brief Preset right hand side in such a way that exact solution is known
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_RHS (NV, NL)
USE SCARC_POINTERS, ONLY: M, L, G, VC, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
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

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I,1,K)) CYCLE
         IC = G%CELL_NUMBER(I,1,K)
         X  = M%XC(I)
         Z  = M%ZC(K)
         !WRITE(MSG%LU_DEBUG,'(A,i3,a,e10.2,a,e10.2,a,E14.6)') 'IC=',IC,':X=',X,':Z=',Z,': RHS=',VC(IC)
         VC(IC) = RHS(X,Z)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_RHS

! --------------------------------------------------------------------------------------------------------------
!> \brief Compute global matrix-vector product A*x = y on grid level NL
! where NV1 is a reference to X and NV2 is a reference to Y
! including data exchange along internal boundaries
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATVEC_PRODUCT(NV1, NV2, NL)
USE SCARC_POINTERS, ONLY: L, OL, G, F, OG, GWC, A, AB, V1, V2, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                          SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_BMATRIX
INTEGER, INTENT(IN) :: NV1, NV2, NL           
REAL(EB) :: TNOW
INTEGER :: NM, NOM, IC, JC, IOR0, ICOL, INBR, ICE, ICW, ICG
INTEGER :: I, J, K, IW, IS=0, IT=0, IL=0, INUM1, INUM2
REAL(EB) :: TMP, VSAVE
#ifdef WITH_MKL
EXTERNAL :: DAXPBY, DAXPY
#endif

TNOW = CURRENT_TIME()

! If this call is related to a globally acting solver, exchange internal boundary values of
! vector1 such that the ghost values contain the corresponding overlapped values of adjacent neighbor
IF (TYPE_MATVEC == NSCARC_MATVEC_GLOBAL) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NV1, NL)

! Perform global matrix-vector product:
! Note: - matrix already contains subdiagonal values from neighbor along internal boundaries
!       - if vector1 contains neighboring values, then correct values of global matvec are achieved
 
SELECT_MATRIX_TYPE: SELECT CASE (SET_MATRIX_TYPE(NL))

   ! COMPACT storage technique
 
   CASE (NSCARC_MATRIX_COMPACT)
   

      MESHES_COMPACT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_GRID (NM, NL)                                    
         IF (IS_LAPLACE) THEN
            A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE)
         ELSE
            A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
         ENDIF
         
         V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
         V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)
   
         IF (NL == NLEVEL_MIN) THEN
            DO IC = 1, G%NC
               ICOL = A%ROW(IC)                                                ! diagonal entry
               JC   = A%COL(ICOL)
               V2(IC) = A%VAL(ICOL)* V1(JC)
               DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1                            ! subdiagonal entries
                  JC = A%COL(ICOL)
                  V2(IC) =  V2(IC) + A%VAL(ICOL) * V1(JC)
               ENDDO
            ENDDO
         ELSE
            DO IC = 1, G%NC
               ICOL = A%ROW(IC)                                                ! diagonal entry
               JC   = A%COL(ICOL)
               V2(IC) = A%VAL(ICOL)* V1(JC)
               DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1                            ! subdiagonal entries
                  JC = A%COL(ICOL)
                  IF (JC == 0) CYCLE
                  V2(IC) =  V2(IC) + A%VAL(ICOL) * V1(JC)
               ENDDO
            ENDDO
         ENDIF
   
      ENDDO MESHES_COMPACT_LOOP
   
   ! Bandwise storage technique:
   ! matrix diagonals are supposed to be constant
   ! matrix-vector multiplication is based on daxpy-routines using the constant matrix stencil
   ! the 'wrong' entries due to boundary conditions and zero entries in subdiagonals are explicitly corrected 
 
   CASE (NSCARC_MATRIX_BANDWISE)
   
      SELECT_STENCIL_TYPE: SELECT CASE (TYPE_STENCIL)

         ! Variable entries with own implementation of daxpyv 
         ! matrix-vector multiplication is based on variable matrix stencil
 
         CASE (NSCARC_STENCIL_VARIABLE)
         
            MESHES_BANDWISE_VARIABLE_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID (NM, NL)                                    
               AB => SCARC_POINT_TO_BMATRIX (NSCARC_MATRIX_POISSON)
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               ! point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               ! point to Y-vector
               V2 = 0.0_EB
         
               !!$OMP PARALLEL default(none) PRIVATE(IOR0, IL, IS, IT) SHARED(AB, V1, V2) 
               DO IOR0 = 3, -3, -1

                  IF (AB%POS(IOR0) == 0) CYCLE

                  IL = AB%LENGTH(IOR0) 
                  IS = AB%SOURCE(IOR0)
                  IT = AB%TARGET(IOR0)

                  CALL SCARC_DAXPY_VARIABLE(IL, AB%VAL(IT:IT+IL-1,AB%POS(IOR0)), V1(IS:IS+IL-1), V2(IT:IT+IL-1))

               ENDDO
               !!$OMP END PARALLEL

               IF (NM==NMESHES) VSAVE = V2(G%NC)                        ! save value in case of pure Neumann bdry

               DO IOR0 = 3, -3, -1
                  F => L%FACE(IOR0)
                  DO INBR = 1, F%N_NEIGHBORS
                     NOM = F%NEIGHBORS(INBR)
                     CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
                     DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
                        ICW = OG%ICG_TO_ICW(ICG, 1)
                        ICE = OG%ICG_TO_ICE(ICG, 1)
                        V2(ICW) = V2(ICW) + F%INCR_FACE * V1(ICE)
                     ENDDO
                  ENDDO
               ENDDO

               IF (IS_PURE_NEUMANN.AND.NM==NMESHES) V2(G%NC) = VSAVE   ! restore value in last cell

            ENDDO MESHES_BANDWISE_VARIABLE_LOOP

 
         ! Storage of constant matrix entries - with corrections at subdiagonals and diagonal (BC's)
 
         CASE (NSCARC_STENCIL_CONSTANT)

            MESHES_BANDWISE_CONSTANT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID (NM, NL)                                    
               AB => G%POISSONB
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               ! point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               ! point to Y-vector
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
               WALL_CELLS_BANDWISE_LOOP: DO IW = 1, G%NW 

                  IOR0 = G%WALL(IW)%IOR
                  IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE              ! cycle boundaries in y-direction for 2D-cases

                  GWC => G%WALL(IW)
                  F  => L%FACE(IOR0)
            
                  I = GWC%IXW
                  J = GWC%IYW
                  K = GWC%IZW
            
                  IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE 
            
                  IC  = G%CELL_NUMBER(I, J, K) 
            
                  ! SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE

                  IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN
                     TMP = V2(IC)
                     SELECT CASE (GWC%BTYPE)
                        CASE (DIRICHLET)
                           V2(IC) = V2(IC) - F%INCR_BOUNDARY * V1(IC)
                        CASE (NEUMANN)
                           V2(IC) = V2(IC) + F%INCR_BOUNDARY * V1(IC)
                     END SELECT
                  ENDIF 
            
               ENDDO WALL_CELLS_BANDWISE_LOOP
         
            ENDDO MESHES_BANDWISE_CONSTANT_LOOP

      END SELECT SELECT_STENCIL_TYPE
END SELECT SELECT_MATRIX_TYPE

CPU(MY_RANK)%MATVEC_PRODUCT =CPU(MY_RANK)%MATVEC_PRODUCT+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_MATVEC_PRODUCT

END MODULE SCARC_VECTORS


!=======================================================================================================================
!
! MODULE SCARC_GRIDS
!
!> \brief Setup all structures related to the different grid types (structured/unstructured) and
!   the different grid resolution levels
!
!   This includes information w.r.t the mesh faces and the wall cells along external, interface and
!   internal boundaries
!
!=======================================================================================================================
MODULE SCARC_GRIDS
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE MESH_VARIABLES, ONLY: MESHES, MESH_TYPE
USE SCARC_CONSTANTS
USE SCARC_TYPES, ONLY: SCARC_LEVEL_TYPE, SCARC_GRID_TYPE
USE SCARC_VARIABLES
USE SCARC_MESSAGES, ONLY: MSG
USE SCARC_STORAGE, ONLY: SCARC_ALLOCATE_INT1, SCARC_DEALLOCATE_INT1, &
                         SCARC_ALLOCATE_INT2, SCARC_ALLOCATE_INT3, &
                         SCARC_ALLOCATE_REAL1, SCARC_ALLOCATE_LOG3
USE SCARC_UTILITIES, ONLY: ARE_FACE_NEIGHBORS, SCARC_SET_GRID_TYPE
USE SCARC_CPU, ONLY: CPU
USE SCARC_TROUBLESHOOTING, ONLY: SCARC_ERROR
USE SCARC_MPI, ONLY: SCARC_SETUP_EXCHANGE_DIMENSIONS

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate basic ScaRC-structures for all needed levels
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BASICS
USE SCARC_POINTERS, ONLY: S, SCARC_POINT_TO_MESH
INTEGER :: NM

ALLOCATE (SCARC(NMESHES), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERROR)

ALLOCATE (STACK(NSCARC_STACK_MAX), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'STACK', IERROR)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MESH(NM)

   ! Needed information about other neighboring meshes
   ALLOCATE (S%OSCARC(NMESHES), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_BASICS', 'OSCARC', IERROR)

   ! Information for single grid levels
   ALLOCATE (S%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_BASICS', 'LEVEL', IERROR)

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_BASICS

! --------------------------------------------------------------------------------------------------------------
!> \brief Determine number of grid levels 
! NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LEVELS
#ifdef WITH_MKL
INTEGER :: NL
#endif

SELECT_SCARC_METHOD: SELECT CASE (TYPE_METHOD)
 
   ! Global data-parallel Krylov method 
 
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

#ifdef WITH_MKL
         ! Preconditioning by defect correction based on LU-decomposition
         ! If two-level method, also use coarse grid level, otherwise only use single (finest) grid level
         ! Either using a global CLUSTER_SPARSE_SOLVER or local PARDISO solvers from MKL
 
         CASE (NSCARC_RELAX_MKL)

            IF (HAS_TWO_LEVELS) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF

            IF (TYPE_SCOPE(1) == NSCARC_SCOPE_GLOBAL) THEN
               TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
            ELSE IF (TYPE_SCOPE(1) == NSCARC_SCOPE_LOCAL) THEN
               TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
            ENDIF

         CASE (NSCARC_RELAX_OPTIMIZED)

            IF (HAS_TWO_LEVELS) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF
            TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
#endif
         ! Preconditioning by defect correction based on geometric multigrid method,
         ! use specified hierarchy of grid levels
 
         CASE (NSCARC_RELAX_GMG)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
            IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif

         ! Preconditioning by defect correction based on local basic iterations (JACOBI/SSOR),
         ! if two-level method, also use coarse grid, otherwise only use single (finest) grid level
 
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

   ! Global data-parallel Multigrid method 
 
   CASE (NSCARC_METHOD_MULTIGRID)

         ! If not specified by user, determine number of possible grid levels

         CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)

#ifdef WITH_MKL

         ! In case of smoothing by different MKL solvers, mark levels for the use of MKL,
         ! either by locally acting PARDISO solvers or Globally acting CLUSTER_SPARSE_SOLVER

         IF (TYPE_SMOOTH == NSCARC_RELAX_MKL) THEN

            IF (TYPE_SCOPE(2) == NSCARC_SCOPE_LOCAL) THEN
               TYPE_MKL(NLEVEL_MIN:NLEVEL_MAX-1) = NSCARC_MKL_LOCAL
            ELSE IF (TYPE_SCOPE(2) == NSCARC_SCOPE_GLOBAL) THEN
               TYPE_MKL(NLEVEL_MIN:NLEVEL_MAX-1) = NSCARC_MKL_GLOBAL
            ENDIF

         ENDIF

         IF (TYPE_MKL(0) == NSCARC_MKL_COARSE) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif
 
   ! Global LU-decomposition 
 
   CASE (NSCARC_METHOD_LU)

      ! Only use single (finest) grid level and mark this level for the use of MKL methods

      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
      TYPE_MKL(NLEVEL_MIN) = TYPE_MKL(0)
 
   ! Global McKenney-Greengard-Mayo method - only finest level 
 
   CASE (NSCARC_METHOD_MGM)

      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)

END SELECT SELECT_SCARC_METHOD

#ifdef WITH_MKL
! Define MKL related logical short names based on number of levels
 
DO NL = NLEVEL_MIN, NLEVEL_MAX
   IS_MKL_LEVEL(NL) = (TYPE_MKL(0)  == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
                      (TYPE_MKL(0)  == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
                      (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL)
ENDDO
#endif
END SUBROUTINE SCARC_SETUP_LEVELS

! --------------------------------------------------------------------------------------------------------------
!> \brief Setup single level in case of default Krylov method
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS(NTYPE)
USE SCARC_POINTERS, ONLY: M
INTEGER, INTENT(IN) :: NTYPE
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM, NLEVEL

SELECT_LEVEL_TYPE: SELECT CASE (NTYPE)

   ! only use finest grid level
   CASE(NSCARC_LEVEL_SINGLE)
   
      NLEVEL     = 1
      NLEVEL_MIN = 1
      NLEVEL_MAX = 1
   
   ! determine maximum number of possible levels based on number of grid cells (based on doubling)

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
      IF (IS_MG .OR. IS_CG_MG) THEN
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
   
END SELECT SELECT_LEVEL_TYPE

END SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS

! --------------------------------------------------------------------------------------------------------------
!> \brief Determine maximum number of possible levels 
! In case of GMG- or 2-Level-method, NC must be divisable by 2 at least one time
! --------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_GET_MAX_LEVEL(NC, IOR0)
INTEGER, INTENT(IN) :: NC, IOR0
INTEGER :: NC0, NL

! Print error message if not divisable by 2

IF (IS_GMG .AND. SCARC_MULTIGRID_LEVEL > 1 .AND. MOD(NC,2)/=0) THEN
   SELECT CASE (ABS(IOR0))
      CASE (1)
         CALL SCARC_ERROR(NSCARC_ERROR_GRID_NUMBERX, SCARC_NONE, NC)
      CASE (2)
         CALL SCARC_ERROR(NSCARC_ERROR_GRID_NUMBERY, SCARC_NONE, NC)
      CASE (3)
         CALL SCARC_ERROR(NSCARC_ERROR_GRID_NUMBERZ, SCARC_NONE, NC)
   END SELECT
ENDIF

! Divide by 2 as often as possible or until user defined max-level is reached

IF (SCARC_MULTIGRID_LEVEL > 1) THEN
   NC0=NC
   DO NL=1,NSCARC_LEVEL_MAX
      NC0=NC0/2
      IF (MOD(NC0,2)/=0) EXIT                ! if no longer divisable by two, leave loop ...
      IF (NL==SCARC_MULTIGRID_LEVEL) EXIT    ! if max possible number of levels reached, leave loop ...
      IF (NC0==1) EXIT                       ! if corresponding power of two has been found, leave loop ...
   ENDDO
   SCARC_GET_MAX_LEVEL=NL
ELSE
   SCARC_GET_MAX_LEVEL=NLEVEL_MIN
ENDIF

RETURN
END FUNCTION SCARC_GET_MAX_LEVEL

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Setup discretization information
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRIDS
USE SCARC_POINTERS, ONLY: M, S, L, G, XCOR, YCOR, ZCOR, XMID, YMID, ZMID, &
                          SCARC_POINT_TO_MESH, SCARC_POINT_TO_GRID
INTEGER :: NL, NM, NC, IX, IY, IZ, IO
INTEGER :: IBAR, JBAR, KBAR

CROUTINE = 'SCARC_SETUP_GRIDS'

! On all grid levels:
! Specify general mesh related geometry information

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MESH(NM)

   ! Store bounds of mesh in SCARC-structure

   S%XS = M%XS;  S%XF = M%XF
   S%YS = M%YS;  S%YF = M%YF
   S%ZS = M%ZS;  S%ZF = M%ZF

   S%IBAR = M%IBAR;  S%JBAR = M%JBAR;  S%KBAR = M%KBAR
   IBAR   = M%IBAR;  JBAR   = M%JBAR;  KBAR   = M%KBAR

   LEVEL_LOOP1: DO NL = NLEVEL_MIN, NLEVEL_MAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    

      L%NX = IBAR;  L%NY = JBAR;  L%NZ = KBAR

      L%N_CELLS = L%NX*L%NY*L%NZ

      L%N_WALL_CELLS_EXT = M%N_EXTERNAL_WALL_CELLS
      L%N_WALL_CELLS_INT = M%N_INTERNAL_WALL_CELLS
      L%N_WALL_CELLS     = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT

      ! Get coordination information

      L%DX = (S%XF-S%XS)/REAL(L%NX,EB) ; L%DXH = 0.5_EB*L%DX ; L%DX2 = L%DX
      L%DY = (S%YF-S%YS)/REAL(L%NY,EB) ; L%DYH = 0.5_EB*L%DY ; L%DY2 = L%DY
      L%DZ = (S%ZF-S%ZS)/REAL(L%NZ,EB) ; L%DZH = 0.5_EB*L%DZ ; L%DZ2 = L%DZ

      L%DXI  = 1.0_EB/L%DX;  L%DYI  = 1.0_EB/L%DY;  L%DZI =  1.0_EB/L%DZ
      L%DXI2 = L%DXI**2;     L%DYI2 = L%DYI**2;     L%DZI2 = L%DZI**2

      ! Needed in case of GMG with multiple grid levels

      IBAR=IBAR/2
      IF (.NOT.TWO_D) JBAR=JBAR/2
      KBAR=KBAR/2

      IF (NL == NLEVEL_MIN) THEN

         ! On finest level store information about obstructions

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

         ! Point to already existing arrays from main FDS program

         XCOR => M%X ;  YCOR => M%Y ;  ZCOR => M%Z
         XMID => M%XC;  YMID => M%YC;  ZMID => M%ZC

      ELSE

         ! Allocate and compute coordinate information for coarser levels

         CALL SCARC_ALLOCATE_REAL1 (L%XCOR, 0, L%NX, NSCARC_INIT_NONE, 'L%XCOR', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1 (L%YCOR, 0, L%NY, NSCARC_INIT_NONE, 'L%YCOR', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1 (L%ZCOR, 0, L%NZ, NSCARC_INIT_NONE, 'L%ZCOR', CROUTINE)

         DO IX = 0, L%NX
            L%XCOR(IX) = S%XS + IX*L%DX
         ENDDO
         DO IY = 0, L%NY
            L%YCOR(IY) = S%YS + IY*L%DY
         ENDDO
         DO IZ = 0, L%NZ
            L%ZCOR(IZ) = S%ZS + IZ*L%DZ
         ENDDO

         ! Allocate and compute midpoint information for coarser levels

         CALL SCARC_ALLOCATE_REAL1 (L%XMID, 0, L%NX+1, NSCARC_INIT_NONE, 'L%XMID', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1 (L%YMID, 0, L%NY+1, NSCARC_INIT_NONE, 'L%YMID', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1 (L%ZMID, 0, L%NZ+1, NSCARC_INIT_NONE, 'L%ZMID', CROUTINE)

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

         XCOR => L%XCOR;  YCOR => L%YCOR;  ZCOR => L%ZCOR
         XMID => L%XMID;  YMID => L%YMID;  ZMID => L%ZMID

      ENDIF

      ! Allocate vectors for step sizes in different directions

      CALL SCARC_ALLOCATE_REAL1 (L%DXL, 0, L%NX, NSCARC_INIT_ZERO, 'L%DXL', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (L%DYL, 0, L%NY, NSCARC_INIT_ZERO, 'L%DYL', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (L%DZL, 0, L%NZ, NSCARC_INIT_ZERO, 'L%DZL', CROUTINE)

      ! Set step sizes between cell midpoints, use interior step sizes for ghost cells as initial values
      ! correct sizes for ghost cells are exchanged later

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

! On finest grid level:
! Allocate several arrays for the administration of discretization related data
 
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)               ! sets pointers M, L and G

   ! Set pointers to already existing cell and wall index arrays from main program (on finest level)
   CALL SCARC_SETUP_CELL_INDEX(L, M, NLEVEL_MIN)
   CALL SCARC_SETUP_WALL_INDEX(L, G, M, NLEVEL_MIN)

   ! Allocate and initialize IS_SOLID array which indicates the state of a cell (gasphase/solid)

   CALL SCARC_ALLOCATE_LOG3 (L%IS_SOLID, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_TRUE, 'L%IS_SOLID', CROUTINE)
   L%IS_SOLID (1:L%NX, 1:L%NY, 1:L%NZ) = .FALSE.

   ! Identify and mark solid obstruction cells in IS_SOLID-part of the discretization

   DO IZ = 1, L%NZ
      DO IY = 1, L%NY
         DO IX = 1, L%NX
            IF (M%SOLID(M%CELL_INDEX(IX, IY, IZ))) L%IS_SOLID(IX, IY, IZ) = .TRUE.
         ENDDO
      ENDDO
   ENDDO
 
   ! If both discretization types (structured/unstructured) must be administrated (MGM method only):
   ! allocate all arrays which are related to a specific discretization type
 
   IF (HAS_MULTIPLE_GRIDS) THEN
 
      ! First process structured discretization
 
      CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_STRUCTURED)
      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

      ! Number of local cells per mesh

      CALL SCARC_ALLOCATE_INT1 (G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_LOCAL', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_OFFSET', CROUTINE)

      ! Allocate wall information array

      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      ! Allocate and preset cell numbers array

      CALL SCARC_ALLOCATE_INT3 (G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'CELL_NUMBER', CROUTINE)

      ! Allocate index array which specifies I, J, K components for all degrees of freedom

      NC = L%NX * L%NY * L%NZ
      CALL SCARC_ALLOCATE_INT1 (G%ICX, 1, NC, NSCARC_INIT_UNDEF, 'G%ICX', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%ICY, 1, NC, NSCARC_INIT_UNDEF, 'G%ICY', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%ICZ, 1, NC, NSCARC_INIT_UNDEF, 'G%ICZ', CROUTINE)

      ! Define local cell numbers for Poisson equation

      DO IZ=1,L%NZ
         DO IY=1,L%NY
            DO IX=1,L%NX
               G%NC_LOCAL(NM) = G%NC_LOCAL(NM) + 1
               G%CELL_NUMBER(IX,IY,IZ) = G%NC_LOCAL(NM)
               G%ICX(G%NC_LOCAL(NM)) = IX
               G%ICY(G%NC_LOCAL(NM)) = IY
               G%ICZ(G%NC_LOCAL(NM)) = IZ
            ENDDO
         ENDDO
      ENDDO
      G%NC   = G%NC_LOCAL(NM)
      G%NCE  = G%NC_LOCAL(NM)
      G%NCE2 = G%NC_LOCAL(NM)
 
      ! Then process unstructured discretization
 
      CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)
      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

      ! Also allocate and preset cell numbers and state arrays for unstructured discretization

      CALL SCARC_ALLOCATE_INT3 (G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'G%CELL_NUMBER', CROUTINE)

      ! Allocate index array which specifies I, J, K components for all degrees of freedom

      NC = L%NX * L%NY * L%NZ
      CALL SCARC_ALLOCATE_INT1 (G%ICX, 1, NC, NSCARC_INIT_UNDEF, 'G%ICX', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%ICY, 1, NC, NSCARC_INIT_UNDEF, 'G%ICY', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%ICZ, 1, NC, NSCARC_INIT_UNDEF, 'G%ICZ', CROUTINE)

      ! Number of local cells per mesh

      CALL SCARC_ALLOCATE_INT1 (G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_LOCAL', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_OFFSET', CROUTINE)

      ! Allocate wall information array

      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      ! Define local cell numbers for Poisson equation

      DO IZ=1,L%NZ
         DO IY=1,L%NY
            DO IX=1,L%NX
               IF (.NOT.L%IS_SOLID(IX,IY,IZ)) THEN
                  G%NC_LOCAL(NM) = G%NC_LOCAL(NM) + 1
                  G%CELL_NUMBER(IX,IY,IZ) = G%NC_LOCAL(NM)
                  G%ICX(G%NC_LOCAL(NM)) = IX
                  G%ICY(G%NC_LOCAL(NM)) = IY
                  G%ICZ(G%NC_LOCAL(NM)) = IZ
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      G%NC   = G%NC_LOCAL(NM)
      G%NCE  = G%NC_LOCAL(NM)
      G%NCE2 = G%NC_LOCAL(NM)
      
      ! Check whether the current meshes contains obstructions or not
      ! This information can be used to replace local preconditioning by the faster FFT for structured meshes
      IF (G%NC < L%NX*L%NY*L%NZ) L%HAS_OBSTRUCTIONS = .TRUE.

   ! If only one specified type of discretization must be admistrated:
   ! allocate and preset cell numbers and state arrays for requested type of discretization
 
   ELSE

      ! Only process specified type of discretization

      CALL SCARC_SET_GRID_TYPE(TYPE_GRID)
      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

      ! Also allocate and preset cell numbers and state arrays for unstructured discretization

      CALL SCARC_ALLOCATE_INT3 (G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'CELL_NUMBER', CROUTINE)

      ! Allocate index array which specifies I, J, K components for all degrees of freedom

      NC = L%NX * L%NY * L%NZ
      CALL SCARC_ALLOCATE_INT1 (G%ICX, 1, NC, NSCARC_INIT_UNDEF, 'G%ICX', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%ICY, 1, NC, NSCARC_INIT_UNDEF, 'G%ICY', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%ICZ, 1, NC, NSCARC_INIT_UNDEF, 'G%ICZ', CROUTINE)

      ! Number of local cells per mesh

      CALL SCARC_ALLOCATE_INT1 (G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_LOCAL', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_OFFSET', CROUTINE)

      ! Allocate wall information array

      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      ! Define local cell numbers for Poisson equation

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
      G%NC   = G%NC_LOCAL(NM)
      G%NCE  = G%NC_LOCAL(NM)
      G%NCE2 = G%NC_LOCAL(NM)

      ! Check whether the current meshes contains obstructions or not
      ! This information can be used to replace local preconditioning by the faster FFT for structured meshes
      IF (TYPE_GRID == NSCARC_GRID_UNSTRUCTURED .AND. G%NC < L%NX*L%NY*L%NZ) L%HAS_OBSTRUCTIONS = .TRUE.

   ENDIF

ENDDO MESHES_LOOP2

END SUBROUTINE SCARC_SETUP_GRIDS

! --------------------------------------------------------------------------------------------------------------------------
!> \brief Setup discretization information on coarser levels
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRID_LEVEL(NL)
USE SCARC_POINTERS, ONLY: LF, LC, GC, SCARC_POINT_TO_MULTIGRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NC, IXF, IYF, IZF, IX, IY, IZ, NSTEP

CROUTINE = 'SCARC_SETUP_GRID_LEVEL'

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NLEVEL_MIN, NL)

   CALL SCARC_ALLOCATE_LOG3 (LC%IS_SOLID, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_FALSE, 'LC%IS_SOLID', CROUTINE)
   LC%IS_SOLID (1:LC%NX, 1:LC%NY, 1:LC%NZ)  = .FALSE.

   NSTEP = 2**(NL - NLEVEL_MIN)

   SELECT CASE(TYPE_GRID)

 
      ! Get cell numberings for coarser grid in case of structured discretization
 
      CASE (NSCARC_GRID_STRUCTURED)

         CALL SCARC_ALLOCATE_INT1 (GC%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_LOCAL', CROUTINE)
         CALL SCARC_ALLOCATE_INT1 (GC%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_OFFSET', CROUTINE)

         CALL SCARC_ALLOCATE_INT3 (GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, &
                                  NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER', CROUTINE)

         NC = LC%NX * LC%NY * LC%NZ
         CALL SCARC_ALLOCATE_INT1 (GC%ICX , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICX', CROUTINE)
         CALL SCARC_ALLOCATE_INT1 (GC%ICY , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICY', CROUTINE)
         CALL SCARC_ALLOCATE_INT1 (GC%ICZ , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICZ', CROUTINE)

         GC%CELL_NUMBER = NSCARC_UNDEF_INT

         DO IZ = 1, LC%NZ
            IZF = (IZ-1)*NSTEP + 1
            DO IY = 1, LC%NY
               IYF = (IY-1)*NSTEP + 1
               DO IX = 1, LC%NX

                  IXF = (IX-1)*NSTEP + 1
                  LC%IS_SOLID(IX,IY,IZ) = LF%IS_SOLID(IXF, IYF, IZF)

                  GC%NC_LOCAL(NM) = GC%NC_LOCAL(NM) + 1
                  GC%CELL_NUMBER(IX,IY,IZ) = GC%NC_LOCAL(NM)

                  GC%ICX(GC%NC_LOCAL(NM)) = IX
                  GC%ICY(GC%NC_LOCAL(NM)) = IY
                  GC%ICZ(GC%NC_LOCAL(NM)) = IZ

               ENDDO
            ENDDO
         ENDDO

         GC%NC = GC%NC_LOCAL(NM)

 
      ! Get cell numberings for coarser grid in case of unstructured discretization
 
      CASE (NSCARC_GRID_UNSTRUCTURED)

         CALL SCARC_ALLOCATE_INT1 (GC%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_LOCAL', CROUTINE)
         CALL SCARC_ALLOCATE_INT1 (GC%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_OFFSET', CROUTINE)

         CALL SCARC_ALLOCATE_INT3 (GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, &
                                  NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER', CROUTINE)

         NC = LC%NX * LC%NY * LC%NZ
         CALL SCARC_ALLOCATE_INT1 (GC%ICX , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICX', CROUTINE)
         CALL SCARC_ALLOCATE_INT1 (GC%ICY , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICY', CROUTINE)
         CALL SCARC_ALLOCATE_INT1 (GC%ICZ , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICZ', CROUTINE)

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


! ------------------------------------------------------------------------------------------------------------------------
!> \brief Get information about global numbers of unknowns for unstructured discretization
! ------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DIMENSIONS(NL)
USE SCARC_POINTERS, ONLY: G, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2

! Preset communication array MESH_INT with local numbers of cells for all meshes depending on type of discretization
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   MESH_INT(NM) = G%NC_LOCAL(NM)

!   DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!      MESH_INT(NM2) = G%NC_LOCAL(NM2)
!   ENDDO

ENDDO MESHES_LOOP1


! Broadcast number of local mesh cells on level NL to all and build global sum
IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,MESH_INT,COUNTS,DISPLS, MPI_INTEGER,MPI_COMM_WORLD,IERROR)
NC_GLOBAL(NL) = SUM(MESH_INT(1:NMESHES))

! Store information on local and global cells numbers on data structure of corresponding discretization type
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

   G%NC_LOCAL(1:NMESHES) = MESH_INT(1:NMESHES)
   G%NC_GLOBAL = SUM(MESH_INT(1:NMESHES))

   ! compute offset between local grid numberings
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

! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup structures related to mesh faces on finest grid level
!   - get dimensions for each of the 6 faces of a mesh
!   - get grid width vector along face
!   - get information for adjacent neighbors
!   - allocate pointer arrays for data exchanges with neighbors
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACES
USE SCARC_POINTERS, ONLY: M, S, L, LC, F, OL, SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID
INTEGER :: NL, NM, NOM
INTEGER :: IFACE, IOR0, JOR0, INBR, IWG, ICW
LOGICAL :: IS_KNOWN(-3:3)
INTEGER :: FACE_NEIGHBORS(-3:3, NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: MESH_NEIGHBORS(6*NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: N_FACE_NEIGHBORS(-3:3)
INTEGER :: N_MESH_NEIGHBORS

CROUTINE = 'SCARC_SETUP_FACES'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)             ! consider only finest grid level

   ! Allocate FACE arrays on different grid levels

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

   CALL SCARC_SETUP_FACE_BASICS(L)
 
   ! Store first wall cell number for each face
 
   ICW = 1
   FACE_ORDER_LOOP: DO IFACE = 1, 6
      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)
      F%NCW0 = ICW
      ICW = ICW + F%NCW
   ENDDO FACE_ORDER_LOOP

   ! Loop over external wall cells:
   ! store basic data and determine number of adajacent neighbors to each face
 
   EXTERNAL_WALL_CELLS_LOOP: DO IWG = 1, L%N_WALL_CELLS_EXT

      NOM  = M%EXTERNAL_WALL(IWG)%NOM
      IOR0 = M%WALL(IWG)%ONE_D%IOR

      IF (NOM /= 0) THEN
         IS_KNOWN = .FALSE.
         DO JOR0 = -3, 3                                              ! neighbor already known?
            IF (JOR0 == 0) CYCLE
            DO INBR = 1, N_FACE_NEIGHBORS(JOR0)
               IF (FACE_NEIGHBORS(JOR0, INBR) == NOM) THEN
                  IS_KNOWN(JOR0) = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDDO
         IF (.NOT.IS_KNOWN(IOR0)) THEN
            N_FACE_NEIGHBORS(IOR0) = N_FACE_NEIGHBORS(IOR0) + 1       ! increase neighbor counter for face
            FACE_NEIGHBORS(IOR0, N_FACE_NEIGHBORS(IOR0)) = NOM        ! store number of neighbor for face
         ENDIF
         IF (.NOT.ANY(IS_KNOWN)) THEN
            N_MESH_NEIGHBORS = N_MESH_NEIGHBORS + 1                   ! increase neighbor counter for mesh
            MESH_NEIGHBORS(N_FACE_NEIGHBORS(IOR0)) = NOM              ! store number of neighbor for mesh
         ENDIF
      ENDIF

   ENDDO EXTERNAL_WALL_CELLS_LOOP

 
   ! Allocate array which stores numbers of all neighboring meshes
 
   IF (N_MESH_NEIGHBORS /= 0) &
      CALL SCARC_ALLOCATE_INT1 (S%NEIGHBORS, 1, N_MESH_NEIGHBORS, NSCARC_INIT_UNDEF, 'S%NEIGHBORS', CROUTINE)
   S%N_NEIGHBORS = N_MESH_NEIGHBORS

   NEIGHBORS_OF_FACE_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE NEIGHBORS_OF_FACE_LOOP

      ! If there are neighbors at face IOR0 store information about them

      F => L%FACE(IOR0)
      IF (N_FACE_NEIGHBORS(IOR0) /= 0) THEN

         ! Allocate array for storing the numbers of the single neighbors

         F%N_NEIGHBORS = N_FACE_NEIGHBORS(IOR0)
         CALL SCARC_ALLOCATE_INT1 (F%NEIGHBORS, 1, N_FACE_NEIGHBORS(IOR0), NSCARC_INIT_NONE, 'F%NEIGHBORS', CROUTINE)

         ! Store every neighbor and allocate corresponding administration arrays on finest level

         DO INBR = 1, N_FACE_NEIGHBORS(IOR0)

            NOM = FACE_NEIGHBORS(IOR0, INBR)
            F%NEIGHBORS(INBR) = NOM                          ! store NOM as a neighbor of that face and if
            CALL SCARC_STORE_NEIGHBOR(NM, NOM)               ! not already done also as mesh neighbor itself

            CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NLEVEL_MIN)
            IF (.NOT.ALLOCATED(OL%FACE)) THEN
               ALLOCATE(OL%FACE(-3:3), STAT=IERROR)
               CALL ChkMemErr('SCARC_SETUP_FACES','OL%FACE',IERROR)
            ENDIF

         ENDDO

      ENDIF
   ENDDO NEIGHBORS_OF_FACE_LOOP

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FACES


! --------------------------------------------------------------------------------------------------------------
!> \brief Setup subdivision information 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIVISION
USE SCARC_POINTERS, ONLY: SUB
INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFFER_INT
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS_NBR   
INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLS_NBR    
INTEGER :: N, NM, INBR, IP, MAX_NBR

CROUTINE = 'SCARC_SETUP_SUBDIVISION'

! Determine number of neighbors for each mesh and make them available in a global array
SUB => SUBDIVISION

CALL SCARC_ALLOCATE_INT1 (SUB%N_NEIGHBORS, 1, NMESHES, NSCARC_INIT_ZERO, 'SUB%N_NEIGHBORS', CROUTINE)
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   SUB%N_NEIGHBORS(NM) = SCARC(NM)%N_NEIGHBORS
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,SUB%N_NEIGHBORS,COUNTS,DISPLS, MPI_INTEGER,MPI_COMM_WORLD,IERROR)
SUB%N_NEIGHBORS_TOTAL = SUM(SUB%N_NEIGHBORS)

CALL SCARC_ALLOCATE_INT1 (COUNTS_NBR, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'COUNTS_NBR', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (DISPLS_NBR, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'DISPLS_NBR', CROUTINE)

DO N = 0, N_MPI_PROCESSES - 1
   DO NM = 1, NMESHES
      IF (PROCESS(NM) == N) COUNTS_NBR(N) = COUNTS_NBR(N) + SUB%N_NEIGHBORS(NM)
   ENDDO
ENDDO
DO N = 1, N_MPI_PROCESSES -1
   DISPLS_NBR(N) = COUNTS_NBR(N-1) + DISPLS_NBR(N-1)
ENDDO

CALL SCARC_ALLOCATE_INT1 (BUFFER_INT, 1, SUB%N_NEIGHBORS_TOTAL, NSCARC_INIT_ZERO, 'BUFFER_INT', CROUTINE)
IP = 1
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO INBR = 1, SCARC(NM)%N_NEIGHBORS
      BUFFER_INT(DISPLS_NBR(PROCESS(NM)) + IP) = SCARC(NM)%NEIGHBORS(INBR)
      IP = IP + 1
   ENDDO
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,BUFFER_INT,COUNTS_NBR,DISPLS_NBR, MPI_INTEGER,MPI_COMM_WORLD,IERROR)
MAX_NBR = MAXVAL(SUB%N_NEIGHBORS)

CALL SCARC_ALLOCATE_INT2 (SUB%NEIGHBORS, 1, MAX_NBR, 1, NMESHES,  NSCARC_INIT_ZERO, 'SUB%NEIGHBORS', CROUTINE)

DO NM = 1, NMESHES
   DO INBR = 1, SUB%N_NEIGHBORS(NM)
      SUB%NEIGHBORS(INBR, NM) = BUFFER_INT(DISPLS_NBR(PROCESS(NM)) + INBR)
   ENDDO
ENDDO

CALL SCARC_DEALLOCATE_INT1 (COUNTS_NBR, 'COUNTS_NBR', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (DISPLS_NBR, 'DISPLS_NBR', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (BUFFER_INT, 'BUFFER_INT', CROUTINE)

END SUBROUTINE SCARC_SETUP_SUBDIVISION


! --------------------------------------------------------------------------------------------------------------
!> \brief Setup neighborship structure for data exchanges along mesh interfaces
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_NEIGHBORS
USE SCARC_POINTERS, ONLY: OS, OLF, OLC
INTEGER :: NM, NOM, NL

!> Setup information about global numbers of unknowns 
CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

! Initialize level structures on neighboring meshes
LEVEL_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LEVEL_NEIGHBOR_LOOP: DO NOM = 1, NMESHES



      ! On finest level point to exchange structures from surrounding FDS 
 
      IF (.NOT. ARE_FACE_NEIGHBORS(NM, NOM)) CYCLE LEVEL_NEIGHBOR_LOOP

      N_EXCHANGES = N_EXCHANGES+1                                         ! count number of exchanges

      OS => SCARC(NM)%OSCARC(NOM)
      ALLOCATE (OS%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)             ! allocate neighboring structures
      CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OS%LEVEL', IERROR)

      OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)                      ! point to neighbor on finest grid level

      OLF%NX = MESHES(NOM)%IBAR                                           ! number of cells in x-direction on neighbor
      OLF%NY = MESHES(NOM)%JBAR                                           ! number of cells in y-direction on neighbor
      OLF%NZ = MESHES(NOM)%KBAR                                           ! number of cells in z-direction on neighbor

      OLF%N_WALL_CELLS_EXT = MESHES(NOM)%N_EXTERNAL_WALL_CELLS            ! number of external wall cells on neighbor
      OLF%N_WALL_CELLS_INT = MESHES(NOM)%N_INTERNAL_WALL_CELLS            ! number of external wall cells on neighbor
      OLF%N_WALL_CELLS     = OLF%N_WALL_CELLS_EXT + OLF%N_WALL_CELLS_INT  ! number of walls cell on neighbor

      OLF%N_CELLS = OLF%NX*OLF%NY*OLF%NZ                                  ! number of cells on neighbor (structured)

 
      ! In case of GMG with a predefined grid hierarchy define corresponding level-structures
 
      IF (NLEVEL_MAX > NLEVEL_MIN) THEN     

         DO NL=NLEVEL_MIN+1,NLEVEL_MAX

            OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)                        ! OLF points to finer, OLC to coarser level

            OLC%NX = OLF%NX/2                                             ! use double grid width
            IF (TWO_D) THEN 
               OLC%NY = 1
            ELSE 
               OLC%NY = OLF%NY/2
            ENDIF
            OLC%NZ = OLF%NZ/2

            OLC%N_CELLS          = OLC%NX * OLC%NY * OLC%NZ               ! set new number of cells
            OLC%N_WALL_CELLS     = OLC%N_WALL_CELLS_EXT                   ! set new number of wall cells
            OLC%N_WALL_CELLS_EXT = 2 * (OLC%NX*OLC%NZ + OLC%NX*OLC%NY + OLC%NY*OLC%NZ)    ! TODO: CHECK!

         ENDDO
      ENDIF

   ENDDO LEVEL_NEIGHBOR_LOOP
ENDDO LEVEL_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_NEIGHBORS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Determine basic data for single faces (orientation, dimensions, numbers)
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE_BASICS(L)
USE SCARC_POINTERS, ONLY: F
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
INTEGER:: IOR0

FACES_OF_MESH_LOOP: DO IOR0 = -3, 3

   IF (IOR0 == 0) CYCLE FACES_OF_MESH_LOOP

   F => L%FACE(IOR0)
   
   SELECT CASE (ABS(IOR0))

      ! ---------- Faces in x-direction
 
      CASE (1)

         F%NOP =  L%NX                           ! number of cells between opposite mesh faces

         F%NX  =  1                              ! number of cells in x-direction
         F%NY  =  L%NY                           ! number of cells in y-direction
         F%NZ  =  L%NZ                           ! number of cells in z-direction

         F%NCW =  L%NY*L%NZ                      ! number of wall cells at that face
         F%DH  => L%DXL                          ! step size vector between opposite mesh faces

         F%INCR_BOUNDARY  = L%DXI2               ! contribution due to boundary condition 
         F%SCAL_DIRICHLET = -2.0_EB * L%DXI2

         F%INCRY =  0
         F%INCRZ =  0
         IF (IOR0 > 0) THEN
            F%SCAL_NEUMANN = L%DXI
            F%INCR_FACE    = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
            F%INCR_INSIDE  = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
            F%INCRX =  1                           ! offset to next internal cell in that direction
         ELSE
            F%SCAL_NEUMANN = -L%DXI
            F%INCR_FACE    = 2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
            F%INCR_INSIDE  = 2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
            F%INCRX = -1
         ENDIF     
         IF (TWO_D) THEN
            F%INCR_STENCIL = (/ F%NY, 0, 0, 0, 0,  0, -F%NY /)
         ELSE
            F%INCR_STENCIL = (/ F%NY, 1, 0, 0, 0, -1, -F%NY /)
         ENDIF
 
      ! ---------- Faces in y-direction
 
      CASE (2)

         F%NOP =  L%NY                   ! dito

         F%NX  =  L%NX
         F%NY  =  1
         F%NZ  =  L%NZ

         F%NCW =  L%NX*L%NZ
         F%DH  => L%DYL

         F%INCR_BOUNDARY  = L%DYI2
         F%SCAL_DIRICHLET = -2.0_EB * L%DYI2

         F%INCRX =  0                           ! offset to next internal cell in that direction
         F%INCRY =  0
         F%INCRZ =  0
         IF (IOR0>0) THEN
            F%SCAL_NEUMANN = L%DYI
            IF (.NOT.TWO_D) THEN
               F%INCR_FACE   = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
               F%INCR_INSIDE = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
               F%INCRY = 1
            ENDIF
         ELSE
            F%SCAL_NEUMANN = -L%DYI
            IF (.NOT.TWO_D) THEN
               F%INCR_FACE   =  2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
               F%INCR_INSIDE =  2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
               F%INCRY = -1
            ENDIF
         ENDIF
         IF (TWO_D) THEN
            F%INCR_STENCIL = (/ 0   , 0, 0, 0,  0, 0, 0     /)             ! special case, not used
         ELSE
            F%INCR_STENCIL = (/ F%NX, 0, 1, 0, -1, 0, -F%NX /)
         ENDIF

      ! ---------- Faces in z-direction
 
      CASE (3)

         F%NOP =  L%NZ                   ! dito

         F%NX  =  L%NX
         F%NY  =  L%NY
         F%NZ  =  1

         F%NCW =  L%NX*L%NY
         F%DH  => L%DZL

         F%NX  = L%NX
         F%INCR_BOUNDARY  = L%DZI2
         F%SCAL_DIRICHLET = -2.0_EB * L%DZI2

         F%INCRX =  0
         F%INCRY =  0
         IF (IOR0>0) THEN
            F%SCAL_NEUMANN = L%DZI
            F%INCR_FACE    = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
            F%INCR_INSIDE  = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
            F%INCRZ =  1
         ELSE
            F%SCAL_NEUMANN = -L%DZI
            F%INCR_FACE    =  2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
            F%INCR_INSIDE  =  2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
            F%INCRZ = -1
         ENDIF
         IF (TWO_D) THEN
            F%INCR_STENCIL = (/ 0, 0   , 1, 0, -1,     0, 0 /)
         ELSE
            F%INCR_STENCIL = (/ 0, F%NX, 1, 0, -1, -F%NX, 0 /)
         ENDIF

   END SELECT

ENDDO FACES_OF_MESH_LOOP

END SUBROUTINE SCARC_SETUP_FACE_BASICS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup wall related structures and boundary conditions
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLS(NGRID_TYPE)
USE SCARC_POINTERS, ONLY: M, L, LF, LC, FF, FC, OL, OLF, OLC, G, GC, GF, OG, OGC, OGF, GWC, MWC, EWC, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_MULTIGRID, &  
                          SCARC_POINT_TO_OTHER_MULTIGRID
INTEGER, INTENT(IN) :: NGRID_TYPE
INTEGER :: NL, NM, NOM
INTEGER :: IREFINE, IFACE, IOR0, JOR0, INBR, IWG, IWC, ICW, IW
LOGICAL :: IS_KNOWN(-3:3), IS_DIRIC, IS_OPEN

CROUTINE = 'SCARC_SETUP_WALLS'
 
! -------- Setup pointers for chosen grid type (structured/unstructured)

CALL SCARC_SET_GRID_TYPE (NGRID_TYPE)

! -------- Get dimensionings for wall cells
 
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

 
   ! First loop over external wall cells:
   ! Determine number of adajacent neighbors to each face with corresponding number of IW's
   ! Store neighbors, orientation and number of couplings for a single wall cell
 
   EXTERNAL_WALL_CELLS_LOOP1: DO IWG = 1, L%N_WALL_CELLS_EXT

      MWC => M%WALL(IWG)
      EWC => M%EXTERNAL_WALL(IWG)

      NOM  =  EWC%NOM
      IOR0 =  MWC%ONE_D%IOR

      GWC => G%WALL(IWG)
      GWC%NOM  = NOM                                    ! store number of neighbor in wall cell
      GWC%IOR  = IOR0                                   ! store orientation of that cell

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

         G%NCE  = G%NCE  + 1                                                ! increase number of extended grid cells
         IF (HAS_AMG_LEVELS) G%NCE2 = G%NCE2 + 2                            ! increase number of extended grid cells type2
         CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NLEVEL_MIN)

         IF (ANY(IS_KNOWN)) OG%NCG = OG%NCG + 1                             ! increase counter for local ghost cells
         IF (OL%GHOST_FIRSTW(IOR0) == 0) OL%GHOST_FIRSTW(IOR0) = OG%NCG     ! save first ghost cell for -IOR0
         IF (OL%GHOST_FIRSTE(IOR0) == 0) OL%GHOST_FIRSTE(IOR0) = OG%NCG     ! save first extended cell for -IOR0
         OL%GHOST_LASTW(IOR0) = OG%NCG                                     
         OL%GHOST_LASTE(IOR0) = OG%NCG                                     

      ENDIF

   ENDDO EXTERNAL_WALL_CELLS_LOOP1
   IF (HAS_AMG_LEVELS) G%ICE2 = G%NCE                                       ! initialize counter for second layer ghost cells

 
   ! Then process internal wall cells
 
   INTERNAL_WALL_CELLS_LOOP1: DO IWG = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT

      MWC => M%WALL(IWG)
      GWC => G%WALL(IWG)

      GWC%IOR  = MWC%ONE_D%IOR
      GWC%NOM  = 0

      GWC%BTYPE = NEUMANN
      GWC%BOUNDARY_TYPE = M%WALL(IWG)%BOUNDARY_TYPE

      GWC%IXG =  MWC%ONE_D%II                        ! ghost cell indices
      GWC%IYG =  MWC%ONE_D%JJ
      GWC%IZG =  MWC%ONE_D%KK

      GWC%IXW =  MWC%ONE_D%IIG                       ! (internal) wall cell indices
      GWC%IYW =  MWC%ONE_D%JJG
      GWC%IZW =  MWC%ONE_D%KKG

   ENDDO INTERNAL_WALL_CELLS_LOOP1

 
   ! Allocate corresponding pointer arrays for data exchanges with neighbors
 
   IF (G%NCE > G%NC) THEN
      G%N_FINE = G%NCE
      CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_IWG, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'G%ICE_TO_IWG', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_ICN, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'G%ICE_TO_IWG', CROUTINE)
   ENDIF

   FACE_NEIGHBORS_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE FACE_NEIGHBORS_LOOP
      DO INBR = 1, L%FACE(IOR0)%N_NEIGHBORS

         NOM = L%FACE(IOR0)%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NLEVEL_MIN)

         CALL SCARC_ALLOCATE_INT1 (OG%ICG_TO_IWG, 1, OG%NCG, NSCARC_INIT_ZERO, 'OG%ICG_TO_IWG', CROUTINE)

         IF (HAS_AMG_LEVELS) THEN
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICW, 1, OG%NCG, 1, 2, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICW', CROUTINE)
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICE, 1, OG%NCG, 1, 2, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICE', CROUTINE)
         ELSE
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICW, 1, OG%NCG, 1, 1, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICW', CROUTINE)
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICE, 1, OG%NCG, 1, 1, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICE', CROUTINE)
         ENDIF

      ENDDO

   ENDDO FACE_NEIGHBORS_LOOP

 
   ! Second loop over external wall cells:
   ! Store detailed coordinate and cell data and get type of boundary condition
 
   G%ICE = G%NC
   WALL_CELLS_LOOP2: DO IWG = 1, L%N_WALL_CELLS_EXT

      IOR0 = G%WALL(IWG)%IOR
      NOM  = G%WALL(IWG)%NOM

      MWC => M%WALL(IWG)
      EWC => M%EXTERNAL_WALL(IWG)

 
      ! Preset ScaRC's boundary type indicator BTYPE
      ! INTERNAL  : the global Poisson problem is solved, so no BC's along mesh interfaces are needed
      ! DIRICHLET : - in the structured case face-wise BC-settings are used ccording to original FFT-solver
      !               (this also allows to use FFT as local preconditioner)
      !             - in the unstructured case Dirichlet BC's are only used for open boundary cells
      ! NEUMANN   : is used for the rest
 
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

      GWC%IXG = MWC%ONE_D%II                                 ! ghost cell indices
      GWC%IYG = MWC%ONE_D%JJ
      GWC%IZG = MWC%ONE_D%KK

      GWC%IXW = MWC%ONE_D%IIG                                ! (internal) wall cell indices
      GWC%IYW = MWC%ONE_D%JJG
      GWC%IZW = MWC%ONE_D%KKG

      ! If there exists a neighbor for that wall cell, setup corresponding neighborship information
      IF (NOM /= 0) CALL SCARC_SETUP_WALL_NEIGHBOR(G, OG, &
                                                   EWC%IIO_MIN, EWC%IIO_MAX, &
                                                   EWC%JJO_MIN, EWC%JJO_MAX, &
                                                   EWC%KKO_MIN, EWC%KKO_MAX, &
                                                   IWG, NM, NOM, NLEVEL_MIN)

   ENDDO WALL_CELLS_LOOP2

ENDDO MESHES_LOOP1

 
! Set dimensions on finest level for requested type(s) of discretization
! and mapping from local to global cell numbering
 
CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

 
! -------- For multi-level variants get discretization information and dimensions on coarser levels
 
DO NL = NLEVEL_MIN+1, NLEVEL_MAX
   CALL SCARC_SETUP_GRID_LEVEL(NL)
   CALL SCARC_SETUP_DIMENSIONS(NL)
ENDDO

 
! -------- Check whether there are no Dirichlet BC's available - TODO: Check !!!
 
MESH_INT = 0                            
RANK_INT = 0

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)
   MESH_INT(NM) = G%N_DIRIC   
   RANK_INT = RANK_INT + MESH_INT(NM)
ENDDO

IF (N_MPI_PROCESSES>1) &
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_INT, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERROR)
N_DIRIC_GLOBAL(NLEVEL_MIN) = RANK_INT

IS_PURE_NEUMANN = N_DIRIC_GLOBAL(NLEVEL_MIN) == 0 .AND. &
                  (TYPE_PRECON /= NSCARC_RELAX_FFT .OR. TYPE_PRECON /= NSCARC_RELAX_FFTO)


 
! -------- Only for multi-level variants 
! (twolevel-CG or GMG method as main solver or preconditioner):
! Determine WALL, FACE and OSCARC types for coarser levels
 
MULTI_LEVEL_IF: IF (HAS_GMG_LEVELS) THEN

   MESHES_LOOP3: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      IREFINE=1
      MULTI_LEVELS_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX

         CALL SCARC_POINT_TO_MULTIGRID (NM, NL-1, NL)

         IREFINE=IREFINE*2
         IF (IS_GMG) CALL SCARC_CHECK_DIVISIBILITY(GF%NCE-GF%NC, 'GF%NCE')

         ! Initialize counts for overlapping and wall cells
 
         GC%NCE = GC%NC + (GF%NCE-GF%NC)/2
         GC%ICE = GC%NC

         LC%N_WALL_CELLS_EXT = SCARC_COUNT_EXTERNAL_WALL_CELLS(LF, LC, GF)
         LC%N_WALL_CELLS_INT = SCARC_COUNT_INTERNAL_WALL_CELLS(LF, LC, GC)

         LC%N_WALL_CELLS = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT

         SELECT CASE(TYPE_GRID)
            CASE (NSCARC_GRID_STRUCTURED)
               GC%NW = LC%N_WALL_CELLS_EXT
            CASE (NSCARC_GRID_UNSTRUCTURED)
               GC%NW = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT
         END SELECT

         ALLOCATE(GC%WALL(LC%N_WALL_CELLS), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_NEIGHBORS','WALL',IERROR)

         ! First allocate administrative mapping arrays for own mesh if there is an overlap
 
         IF (GC%NCE > GC%NC) THEN
            CALL SCARC_ALLOCATE_INT1 (GC%ICE_TO_IWG , GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'GC%ICE_TO_IWG', CROUTINE)
            CALL SCARC_ALLOCATE_INT1 (GC%ICE_TO_ICN , GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'GC%ICE_TO_ICN', CROUTINE)
         ENDIF

         ! Setup basic face information for coarser mesh
 
         CALL SCARC_SETUP_FACE_BASICS(LC)

         IWC = 1
         IWG = 1
         ICW = 1
         FACES_LOOP: DO IFACE = 1, 6

            IOR0 = FACE_ORIENTATION(IFACE)

            FF => LF%FACE(IOR0)
            FC => LC%FACE(IOR0)

            ! initialize FACE type for coarser mesh
            FC%NCW0 = ICW
            FC%N_NEIGHBORS = FF%N_NEIGHBORS

            IF (FC%N_NEIGHBORS /= 0) &
               CALL SCARC_ALLOCATE_INT1 (FC%NEIGHBORS, 1, FC%N_NEIGHBORS, NSCARC_INIT_NONE, 'FC%FACE_NEIGHBORS', CROUTINE)
            DO INBR= 1, FC%N_NEIGHBORS
               FC%NEIGHBORS(INBR) = FF%NEIGHBORS(INBR)
            ENDDO
            
            FC%NCW = FC%NX * FC%NY * FC%NZ                                ! get number of wall cells for that face
            ICW = ICW + FC%NCW                                            ! increase global wall cell counter

            ! Get related data and pointer structures for every mesh neighbor

            IF (LF%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
               DO INBR = 1, LF%FACE(IOR0)%N_NEIGHBORS

                  NOM = LF%FACE(IOR0)%NEIGHBORS(INBR)
                  CALL SCARC_POINT_TO_OTHER_MULTIGRID (NM, NOM, NL-1, NL)

                  IF (IS_GMG) THEN
                     CALL SCARC_CHECK_DIVISIBILITY(OLF%N_WALL_CELLS_LOCAL, 'OLF%N_WALL_CELLS_LOCAL')
                     CALL SCARC_CHECK_DIVISIBILITY(OGF%NCG, 'OGF%NCG')
                  ENDIF

                  IF (.NOT.TWO_D) THEN
                     OLC%N_WALL_CELLS_LOCAL = OLF%N_WALL_CELLS_LOCAL/4
                     OGC%NCG = OGF%NCG/4
                  ELSE
                     OLC%N_WALL_CELLS_LOCAL = OLF%N_WALL_CELLS_LOCAL/2
                     OGC%NCG = OGF%NCG/2
                  ENDIF

                  CALL SCARC_SETUP_EXCHANGE_DIMENSIONS(GC, OGC, NOM, IREFINE)

                  CALL SCARC_ALLOCATE_INT1 (OGC%ICG_TO_IWG, 1, OGC%NCG, NSCARC_INIT_ZERO, 'OGC%ICG_TO_IWG', CROUTINE)
                  CALL SCARC_ALLOCATE_INT2 (OGC%ICG_TO_ICW, 1, OGC%NCG, 1, 1, NSCARC_INIT_ZERO, 'OGC%ICG_TO_ICW', CROUTINE)
                  CALL SCARC_ALLOCATE_INT2 (OGC%ICG_TO_ICE, 1, OGC%NCG, 1, 1, NSCARC_INIT_ZERO, 'OGC%ICG_TO_ICE', CROUTINE)

               ENDDO
            ENDIF

            ! Setup complete wall information on coarser mesh level

            CALL SCARC_SETUP_WALL_LEVEL(LF, LC, GF, GC, IOR0, IWC, IREFINE, NM, NL)

         ENDDO FACES_LOOP

         CALL SCARC_SETUP_CELL_INDEX (LC, M, NL)
         CALL SCARC_SETUP_WALL_COORDS(LC, GC)
         CALL SCARC_SETUP_WALL_INDEX (LC, GC, M, NL)

         ! Setup order in which ghost cells are processed during data exchanges

         WALLCELLS_LOOP: DO IW = 1, LC%N_WALL_CELLS
         
           NOM = GC%WALL(IW)%NOM
           IF (NOM /= 0) THEN
              IOR0 = GC%WALL(IW)%IOR
              CALL SCARC_POINT_TO_OTHER_MULTIGRID (NM, NOM, NL-1, NL)
              OGC%ICG2 = OGC%ICG2 + 1
              IF (OLC%GHOST_FIRSTW(IOR0) == 0) OLC%GHOST_FIRSTW(IOR0) = OGC%ICG2
              IF (OLC%GHOST_FIRSTE(IOR0) == 0) OLC%GHOST_FIRSTE(IOR0) = OGC%ICG2
              OLC%GHOST_LASTW(IOR0) = OGC%ICG2 
              OLC%GHOST_LASTE(IOR0) = OGC%ICG2 
           ENDIF
         
         ENDDO WALLCELLS_LOOP

      ENDDO MULTI_LEVELS_LOOP
   ENDDO MESHES_LOOP3
ENDIF MULTI_LEVEL_IF


! Correct boundary types for cells adjacent to obstructions on ghost cells

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)                  ! sets level and grid pointers L and G
   IF (IS_UNSTRUCTURED) THEN
      CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(L, G)
      IF (.NOT.HAS_AMG_LEVELS) THEN
         DO NL = NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(L, G)
         ENDDO
      ENDIF
   ENDIF
ENDDO

! Debug FACE, WALL and DISCRET structures - only if directive SCARC_DEBUG is set


END SUBROUTINE SCARC_SETUP_WALLS


! -------------------------------------------------------------------------------------------------------------
!> \brief Store all neighbors of a mesh
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_STORE_NEIGHBOR(NM, NOM)
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: INBR
DO INBR = 1, SCARC(NM)%N_NEIGHBORS
   IF (SCARC(NM)%NEIGHBORS(INBR) == NSCARC_UNDEF_INT) EXIT      ! not found, to be stored
   IF (SCARC(NM)%NEIGHBORS(INBR) == NOM) RETURN                 ! nothing to do, already stored
ENDDO
SCARC(NM)%NEIGHBORS(INBR) = NOM
RETURN
END SUBROUTINE SCARC_STORE_NEIGHBOR


! -------------------------------------------------------------------------------------------------------------
!> \brief Setup cells indexing array on coarser grid levels in case of MG method
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CELL_INDEX(L, M, NL)
USE SCARC_POINTERS, ONLY: OB
TYPE (MESH_TYPE), POINTER, INTENT(IN) :: M
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
INTEGER, INTENT(IN) :: NL
INTEGER :: I, J, K, NOBST

CROUTINE = 'SCARC_SETUP_CELL_INDEX'

! If finest level, the corresponding CELL_INDEX array is already available by surrounding routines
! on coarser levels, it must still be computed

IF (NL == NLEVEL_MIN) THEN

   L%CELL_INDEX_PTR => M%CELL_INDEX

ELSE

   CALL SCARC_ALLOCATE_INT3 (L%CELL_INDEX, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'L%CELL_INDEX', CROUTINE)
   L%CELL_INDEX_PTR => L%CELL_INDEX
   L%N_CELL_INDEX = 0

   ! Preset it for all grid cells
 
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

 
   ! Consider cells in obstructions
 
   DO NOBST=1,L%N_OBST
      OB => L%OBST(NOBST)
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


! -------------------------------------------------------------------------------------------------------------
!> \brief Setup wall cells indexing array on coarser grid levels
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_INDEX(L, G, M, NL)
TYPE (MESH_TYPE), POINTER, INTENT(IN) :: M
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: NL
INTEGER :: I, J, K, ICG, IW, IOR0

CROUTINE = 'SCARC_SETUP_WALL_INDEX'

! if on finest level, the array WALL_INDEX is already available by surrounding routines
IF (NL == NLEVEL_MIN) THEN

   L%WALL_INDEX_PTR => M%WALL_INDEX

   ! if on coarser levels, it must still be computed
ELSE

   CALL SCARC_ALLOCATE_INT2 (L%WALL_INDEX, 1, L%N_CELL_INDEX, -3, 3, NSCARC_INIT_ZERO, 'L%WALL_INDEX', CROUTINE)
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


! ---------------------------------------------------------------------------------------------------------------
!> \brief Setup all necessary information for a wall cell with neighbor in case of MG method
! Number of obstructions on coarse level is the same as on fine level
! TODO: Only works for special cases which run for GMG, must still be extended!!
! ---------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_COORDS(L, G)
USE SCARC_POINTERS, ONLY: OB, GWC
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN)  :: G
INTEGER :: IC, IO, IWC
INTEGER :: I, J, K

IWC = L%N_WALL_CELLS_EXT + 1
DO IO = 1, L%N_OBST

   OB => L%OBST(IO)

   ! Analyze IOR = 1

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

   ! Analyze IOR = -1
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

   ! Analyze IOR = 2

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

   ! Analyze IOR = -2

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

   ! Analyze IOR = 3

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

   ! Analyze IOR = -3

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


! ---------------------------------------------------------------------------------------------------------------
!> \brief Correct boundary type array related to internal obstructions on ghost cells
! ---------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS(L, G) 
USE SCARC_POINTERS, ONLY: GWC
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: G
INTEGER :: IWG
INTEGER :: IX, IY, IZ, IOR0, BTYPE0

DO IWG = 1, L%N_WALL_CELLS_EXT

   GWC => G%WALL(IWG)
   IF (GWC%NOM == 0) CYCLE                    ! TODO: equal or not equal ??

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


! ---------------------------------------------------------------------------------------------------------------
!> \brief Setup all necessary information for a wall cell with neighbor
! ---------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS(LF, LC, GF)
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: LF, LC
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: GF
INTEGER :: IXC, IYC, IZC
INTEGER :: IXF, IYF, IZF
INTEGER :: IWC, ICF(4)=0, IWF(4)=0, IOR0

ICF = 0
IWC = 0
IWF = 0

IF (TWO_D) THEN
   IYC = 1
   IYF = 1

   ! IOR = 1

   IOR0 = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      ICF(1) = LF%CELL_INDEX_PTR(1  , IYF  , IZF  )
      ICF(2) = LF%CELL_INDEX_PTR(1  , IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 2)) IWC = IWC + 1
   ENDDO

   ! IOR = -1

   IOR0 = -1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      ICF(1) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF  )
      ICF(2) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 2)) IWC = IWC + 1
   ENDDO

   ! IOR = 2

   IOR0 = 2
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF  , IYF, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1, IYF, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = -2

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
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = 3

   IOR0 = 3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF  , 1)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF  , 1)
      IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 2)) IWC = IWC + 1
   ENDDO

   ! IOR = -3

   IOR0 = -3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF  , LF%NZ)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF  , LF%NZ)
      IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 2)) IWC = IWC + 1
   ENDDO

ELSE

   ! IOR = 1

   IOR0 = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IYC = 1, LC%NY
         IYF = 2*IYC - 1
         ICF(1) = LF%CELL_INDEX_PTR(1  , IYF  , IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(1  , IYF+1, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(1  , IYF  , IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(1  , IYF+1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = -1

   IOR0 = -1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IYC = 1, LC%NY
         IYF = 2*IYC - 1
         ICF(1) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(LF%NX, IYF+1, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(LF%NX, IYF+1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = 2

   IOR0 = 2
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF  , 1, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1, 1, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF  , 1, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1, 1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = -2

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
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = 3

   IOR0 = 3
   DO IYC = 1, LC%NY
      IYF = 2*IYC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF  , 1)
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF  , 1)
         ICF(3) = LF%CELL_INDEX_PTR(IXF    , IYF+1, 1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1  , IYF+1, 1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = -3

   IOR0 = -3
   DO IYC = 1, LC%NY
      IYF = 2*IYC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF  , LF%NZ)
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF  , LF%NZ)
         ICF(3) = LF%CELL_INDEX_PTR(IXF  , IYF+1, LF%NZ)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1, IYF+1, LF%NZ)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

ENDIF

SCARC_COUNT_EXTERNAL_WALL_CELLS = IWC
END FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS


! ---------------------------------------------------------------------------------------------------------------
!> \brief Setup all necessary information for a wall cell with neighbor
! ---------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS(LF, LC, GC)
USE SCARC_POINTERS, ONLY: OBF, OBC
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: LF, LC
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: GC
INTEGER :: IWC, NW_INT
INTEGER :: IC, IO
INTEGER :: I, J, K

LC%N_OBST = LF%N_OBST                   ! Number of obstructions is the same on all levels

ALLOCATE(LC%OBST(LC%N_OBST), STAT=IERROR)
CALL ChkMemErr('SCARC_COUNT_INTERNAL_WALL_CELLS','OBST',IERROR)

NW_INT = 0
IWC = LC%N_WALL_CELLS_EXT + 1

DO IO = 1, LF%N_OBST

   OBF => LF%OBST(IO)
   OBC => LC%OBST(IO)

   OBC%I1 = (OBF%I1+1)/2
   OBC%I2 =  OBF%I2/2

   IF (TWO_D) THEN
      OBC%J1 = 0
      OBC%J2 = 1
   ELSE
      OBC%J1 = (OBF%J1+1)/2
      OBC%J2 =  OBF%J2/2
   ENDIF

   OBC%K1 = (OBF%K1+1)/2
   OBC%K2 =  OBF%K2/2

   ! Analyze IOR = 1

   I = OBC%I1
   DO K = OBC%K1+1, OBC%K2
      DO J = OBC%J1+1, OBC%J2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -1

   I = OBC%I2
   DO K = OBC%K1+1, OBC%K2
      DO J = OBC%J1+1, OBC%J2
         IC = GC%CELL_NUMBER(I+1, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = 2

   J = OBC%J1
   DO K = OBC%K1+1, OBC%K2
      DO I = OBC%I1+1, OBC%I2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -2

   J = OBC%J2
   DO K = OBC%K1+1, OBC%K2
      DO I = OBC%I1+1, OBC%I2
         IC = GC%CELL_NUMBER(I, J+1, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = 3

   K = OBC%K1
   DO J = OBC%J1+1, OBC%J2
      DO I = OBC%I1+1, OBC%I2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -3

   K = OBC%K2
   DO J = OBC%J1+1, OBC%J2
      DO I = OBC%I1+1, OBC%I2
         IC = GC%CELL_NUMBER(I, J, K+1)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

ENDDO

SCARC_COUNT_INTERNAL_WALL_CELLS = NW_INT
END FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS


! ---------------------------------------------------------------------------------------------------------------
!> \brief Count external wall cells on specified face if mesh
! ---------------------------------------------------------------------------------------------------------------
LOGICAL FUNCTION IS_EXTERNAL_WALLCELL(L, G, IOR0, ICF, NCNT)
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: IOR0, NCNT
INTEGER, DIMENSION(:), INTENT(IN) :: ICF
INTEGER :: I, IWF_LAST, IWF(4)=0
REAL(EB) :: BSUM

IS_EXTERNAL_WALLCELL = .FALSE.

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
   CALL SCARC_ERROR(NSCARC_ERROR_BOUNDARY_SUM, SCARC_NONE, IOR0)
ENDIF

END FUNCTION IS_EXTERNAL_WALLCELL


! ---------------------------------------------------------------------------------------------------------------
!> \brief Setup all necessary information for a wall cell with neighbor
! ---------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_NEIGHBOR(G, OG, NX1, NX2, NY1, NY2, NZ1, NZ2, IWG, NM, NOM, NL)
USE SCARC_POINTERS, ONLY: GWC, SCARC_POINT_TO_OTHER_GRID
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: G, OG
INTEGER, INTENT(IN) :: NX1, NX2, NY1, NY2, NZ1, NZ2
INTEGER, INTENT(IN) :: IWG, NM, NOM, NL
INTEGER :: NOMX, NOMY, NOMZ
INTEGER :: ICG, ICE, IX, IY, IZ, JL, IXW, IYW, IZW, IXG, IYG, IZG

CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

ICE  = G%ICE
ICG  = OG%ICG

! set neighboring coordinates
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

! store information about overlapped cells and set mapping arrays
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ICG  = ICG  + 1
         ICE  = ICE  + 1

         GWC%ICE = ICE                                         ! number of extended grid cell
         GWC%ICG = ICG                                         ! number of ghost grid cell

         G%ICE_TO_IWG(ICE) = IWG                               ! map extended cell to global wall cell
         
         IXG = G%WALL(IWG)%IXG
         IYG = G%WALL(IWG)%IYG
         IZG = G%WALL(IWG)%IZG

         IXW = G%WALL(IWG)%IXW
         IYW = G%WALL(IWG)%IYW
         IZW = G%WALL(IWG)%IZW

         G%CELL_NUMBER(IXG, IYG, IZG) = ICE

         OG%ICG_TO_IWG(ICG)    = IWG                              ! map ghost cell to global wall cell
         OG%ICG_TO_ICW(ICG, 1) = G%CELL_NUMBER(IXW, IYW, IZW)     ! get cell number of adjacent internal cell

      ENDDO
   ENDDO
ENDDO

G%ICE  = ICE                                                   ! store extended cell counter
OG%ICG = ICG                                                   ! store ghost cell counter
OG%ICG_TO_IWG(ICG) = IWG                                       ! map local wall cell to global wall cell
OG%ICG_TO_ICE(ICG, 1) = ICE                                    ! map local wall cell to global wall cell

END SUBROUTINE SCARC_SETUP_WALL_NEIGHBOR


! ------------------------------------------------------------------------------------------------------------------
!> \brief Check divisibility by 2 of a given number of elements (in one grid direction)
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_DIVISIBILITY(NN, CDIR)
INTEGER, INTENT(IN) :: NN
CHARACTER(*) , INTENT(IN) :: CDIR
IF (MOD(NN,2) /= 0) CALL SCARC_ERROR(NSCARC_ERROR_GRID_NUMBER, CDIR, NSCARC_NONE)
END SUBROUTINE SCARC_CHECK_DIVISIBILITY


! ------------------------------------------------------------------------------------------------------------------
!> \brief Set wall cell information on coarse level
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_LEVEL(LF, LC, GF, GC, IOR0, IWC, IREFINE, NM, NL)
USE SCARC_POINTERS, ONLY: FF, FC, WF, WC, OGC, SCARC_POINT_TO_OTHER_MULTIGRID
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: LF, LC
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: GF, GC
INTEGER, INTENT(INOUT) :: IWC
INTEGER, INTENT(IN) :: IOR0, IREFINE, NM, NL
INTEGER :: IWF(4) , IBCF(4), NOMF(4)
INTEGER :: IX,  IY,  IZ, I
INTEGER :: NX1, NY1, NZ1
INTEGER :: NX2, NY2, NZ2
INTEGER :: IX1, IY1, IZ1
INTEGER :: IX2, IY2, IZ2
INTEGER :: IDIFF, JDIFF, KDIFF

WC => GC%WALL
WF => GF%WALL

! set coordinate dimensions for correspoding face

SELECT CASE (ABS(IOR0))
   CASE (1)
      IF (IOR0 > 0) THEN                                           ! set dimensions for wall cell counting
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

 
! Loop over all wall cells of face IOR0
 
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ! Set orientation of neiboring face, indices of ghost and adjacent cell for coarse IW

         WC(IWC)%IOR = IOR0

         FF => LF%FACE(IOR0)
         FC => LC%FACE(IOR0)

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

         ! ------------------------------------------------------------
         !  2D-version
         ! ------------------------------------------------------------
         IF (TWO_D) THEN

            ! determine fine IW's, which must be merged to one coarse IW

            SELECT CASE (ABS(IOR0))
               CASE ( 1)
                  IWF(1) = FF%NCW0 + 2*(IZ-1)
               CASE ( 2)
                  IWF(1) = FF%NCW0 + 2*(IZ-1)*LF%NX + 2*(IX - 1)
               CASE ( 3)
                  IWF(1) = FF%NCW0 + 2*(IX-1)
            END SELECT
            IWF(2) = IWF(1)+1

            ! set fine cell neighbors (they must be the same for all fine IW's)

            NOMF(1) = WF(IWF(1))%NOM
            NOMF(2) = WF(IWF(2))%NOM
            IF (NOMF(1) /= NOMF(2)) CALL SCARC_ERROR(NSCARC_ERROR_NEIGHBOR_TYPE, SCARC_NONE, NOMF(1))

            WC(IWC)%NOM = NOMF(1)

            ! set corresponding pressure_bc_index on coarser level

            IBCF(1) = WF(IWF(1))%BTYPE
            IBCF(2) = WF(IWF(2))%BTYPE
            IF (IBCF(1) == INTERNAL .OR. IBCF(2) == INTERNAL) THEN
               WC(IWC)%BTYPE = INTERNAL
            ELSE IF (IBCF(1) == DIRICHLET .OR. IBCF(2) == DIRICHLET) THEN
               WC(IWC)%BTYPE = DIRICHLET
            ELSE
               WC(IWC)%BTYPE = NEUMANN
            ENDIF

            ! set corresponding pressure_bc_index on coarser level

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
               CALL SCARC_ERROR(NSCARC_ERROR_BOUNDARY_TYPE, SCARC_NONE, IBCF(1))
            ENDIF

            ! in case of an internal boundary set neighboring wall cells

            IF (NOMF(1) > 0) THEN

               CALL SCARC_POINT_TO_OTHER_MULTIGRID (NM, NOMF(1), NL-1, NL)

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
                        CALL SCARC_ERROR(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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
                        CALL SCARC_ERROR(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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

                
               ! Allocate and specify ICN and ICE arrays for OC
                
               CALL SCARC_SETUP_WALL_NEIGHBOR(GC, OGC, IX1, IX2, 1, 1, IZ1, IZ2, IWC, NM, NOMF(1), NL)

            ENDIF

         ! ------------------------------------------------------------
         ! 3D-version
         ! ------------------------------------------------------------
         ELSE

            ! determine fine IW's, which must be merged to one coarse IW

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

            ! set fine cell neighbors (they must be the same for all fine IW's)

            DO I=1,4
               NOMF(I) = WF(IWF(I))%NOM
            ENDDO

            IF (NOMF(1) /= NOMF(2) .OR. NOMF(1) /= NOMF(3) .OR. NOMF(1) /= NOMF(4)) &
               CALL SCARC_ERROR(NSCARC_ERROR_NEIGHBOR_TYPE, SCARC_NONE, IOR0)
            WC(IWC)%NOM = NOMF(1)

            ! set corresponding pressure_bc_index on coarser level

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

            ! set corresponding pressure_bc_index on coarser level

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
               CALL SCARC_ERROR(NSCARC_ERROR_BOUNDARY_TYPE, SCARC_NONE, -999)
            ENDIF

            ! in case of an internal boundary set WALL(10:15,IWC)

            IF (NOMF(1) > 0) THEN

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
                        CALL SCARC_ERROR(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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
                        CALL SCARC_ERROR(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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
                        CALL SCARC_ERROR(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
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

               CALL SCARC_SETUP_WALL_NEIGHBOR(GC, OGC, IX1, IX2, IY1, IY2, IZ1, IZ2, IWC, NM, NOMF(1), NL)

            ENDIF
         ENDIF
         IWC = IWC + 1
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_WALL_LEVEL

! -------------------------------------------------------------------------------------------------------------
!> \brief Setup mapping from local to global cell numbering
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBAL_CELL_MAPPING(NL)
USE SCARC_POINTERS, ONLY: G, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, IC, IW, ICE, ICN

!IF (NMESHES == 1) RETURN                                     ! TODO: omit LOCAL_TO_GLOBAL vector for single-mesh cases
CROUTINE = 'SCARC_SETUP_GLOBAL_CELL_MAPPING'

!IF (IS_MGM .AND. IS_UNSTRUCTURED) RETURN                     ! TODO: omit LOCAL_TO_GLOBAL vector for local pass

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)      

   IF (HAS_AMG_LEVELS) THEN
      CALL SCARC_ALLOCATE_INT1 (G%LOCAL_TO_GLOBAL, 1, G%NCE2, NSCARC_INIT_ZERO, 'G%LOCAL_TO_GLOBAL', CROUTINE)  
   ELSE
      CALL SCARC_ALLOCATE_INT1 (G%LOCAL_TO_GLOBAL, 1, G%NCE , NSCARC_INIT_ZERO, 'G%LOCAL_TO_GLOBAL', CROUTINE) 
   ENDIF

   DO IC = 1, G%NC
      G%LOCAL_TO_GLOBAL(IC) = IC + G%NC_OFFSET(NM)
   ENDDO

   DO IW = 1, G%NW
      NOM = G%WALL(IW)%NOM
      IF (NOM == 0) CYCLE
      ICE = G%WALL(IW)%ICE
      ICN = G%ICE_TO_ICN(ICE)
      G%LOCAL_TO_GLOBAL(ICE) = ICN + G%NC_OFFSET(NOM)
   ENDDO

ENDDO

END SUBROUTINE SCARC_SETUP_GLOBAL_CELL_MAPPING

END MODULE SCARC_GRIDS


!=======================================================================================================================
!
! MODULE SCARC_MATRICES
!
!> \brief Setup and organize the matrix types needed for the different ScaRC/UscaRC solvers
!
!   This inlcudes local/global Poisson and Laplace matrices, their boundary conditions and 
!   a corresponding condensing in the purely Neumann case
!
!=======================================================================================================================
MODULE SCARC_MATRICES
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE SCARC_CONSTANTS
USE SCARC_UTILITIES
USE SCARC_STORAGE
USE SCARC_MPI, ONLY: SCARC_EXCHANGE
USE SCARC_GRIDS, ONLY: SCARC_SETUP_GLOBAL_CELL_MAPPING

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup system of equations (Poisson matrix + BC's) for different variants of ScaRC
! Define matrix stencils and initialize matrices and boundary conditions on all needed levels
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEMS
USE SCARC_POINTERS, ONLY: SCARC_POINT_TO_GRID
#ifdef WITH_MKL
USE SCARC_POINTERS, ONLY: L
INTEGER :: TYPE_MKL_SAVE(0:1), TYPE_SCOPE_SAVE(0:1)
#endif
INTEGER :: NM, NL
  
CROUTINE = 'SCARC_SETUP_SYSTEMS'

! ------ Setup sizes for system matrices
  
SELECT_SCARC_METHOD_SIZES: SELECT CASE (TYPE_METHOD)

   ! -------- Global Krylov method

   CASE (NSCARC_METHOD_KRYLOV)
   
      CALL SCARC_SET_GRID_TYPE (TYPE_GRID)                      ! process specified discretization type
      CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)               ! setup sizes on finest level
   
      IF (HAS_TWO_LEVELS .AND. .NOT.HAS_AMG_LEVELS) &
         CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MAX)            ! twolevel-precon: also setup size for coarse level
   
      IF (IS_CG_GMG) THEN                                                   
         DO NL=NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_SETUP_POISSON_SIZES (NL)                 ! GMG-precon: also setup size for all other levels
         ENDDO
      ENDIF
   
   ! -------- Global Multigrid method

   CASE (NSCARC_METHOD_MULTIGRID)
   
      CALL SCARC_SET_GRID_TYPE (TYPE_GRID)                      ! process specified discretization type
      SELECT CASE (TYPE_MULTIGRID)
         CASE (NSCARC_MULTIGRID_GEOMETRIC)                                   
            DO NL=NLEVEL_MIN, NLEVEL_MAX
               CALL SCARC_SETUP_POISSON_SIZES (NL)              ! GMG: setup size for all levels
            ENDDO
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)                                   
            CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)         ! AMG: setup sizes only on finest level
      END SELECT
   
   ! -------- Global MGM method - currently just proof of concept

   CASE (NSCARC_METHOD_MGM)
   
      CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_STRUCTURED)         ! First process structured discretization
      CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)        
   
      CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)       ! Then process unstructured discretization
      IF (SCARC_MGM_CHECK_LAPLACE .OR. SCARC_MGM_EXACT_INITIAL) &
         CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)            ! ... for global Poisson matrix (only if requested)
      CALL SCARC_SETUP_LOCAL_LAPLACE_SIZES (NLEVEL_MIN)         ! ... for local Laplace matrices
   
END SELECT SELECT_SCARC_METHOD_SIZES

! ------ Assemble system matrices on requested grid levels and set boundary conditions
  
MESHES_POISSON_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   SELECT_SCARC_METHOD: SELECT CASE (TYPE_METHOD)

      ! ---------- Krylov method (CG) as main solver, different preconditioners possible

      CASE (NSCARC_METHOD_KRYLOV)

         ! For all different possible Krylov variants, first setup Poisson matrix on finest level including BC's 

         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         ! Depending on the requested preconditioner, also assemble the Poisson matrix with BC's on specific coarser levels

         SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

            ! In case of multigrid as preconditioner:
            ! only build higher level structures in case of geometric multigrid (algebraic variant is done elsewhere)

            CASE (NSCARC_RELAX_GMG)

               IF (IS_CG_GMG) THEN
                  DO NL = NLEVEL_MIN+1, NLEVEL_MAX
                     CALL SCARC_SETUP_POISSON (NM, NL)
                     CALL SCARC_SETUP_BOUNDARY(NM, NL)
                  ENDDO
               ENDIF

#ifdef WITH_MKL
            ! In case of LU-decomposition as preconditioner
            ! locally acting: PARDISO from MKL as preconditioners on fine level with possible coarse grid correction

            CASE (NSCARC_RELAX_MKL)

               IF (TYPE_SCOPE(1) == NSCARC_SCOPE_LOCAL .AND. HAS_TWO_LEVELS .AND. .NOT.HAS_AMG_LEVELS) THEN
                  CALL SCARC_SETUP_POISSON (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
               ENDIF
#endif

            ! in case of default preconditioners (JACOBI/SSOR/FFT/...):
            ! if there is an additional coarse grid correction which is NOT AMG-based, 
            ! then also assemble matrix on coarse grid level

            CASE DEFAULT
   
               IF (HAS_TWO_LEVELS .AND. .NOT.HAS_AMG_LEVELS) THEN
                  CALL SCARC_SETUP_POISSON (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
               ENDIF

         END SELECT SELECT_KRYLOV_PRECON

      ! ---------- Multigrid as main solver

      CASE (NSCARC_METHOD_MULTIGRID)

         ! For all different possible multigrid-variants, first setup Poisson matrix on finest level including BC's 

         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         ! On case of a  geometric multigrid, assemble standard n-point-matrix hierarchy on all coarser levels, too
         ! Note: in case of an algebraic multigrid, this will be done in a separate routine later

         IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
            DO NL = NLEVEL_MIN + 1, NLEVEL_MAX
               CALL SCARC_SETUP_POISSON (NM, NL)
               CALL SCARC_SETUP_BOUNDARY(NM, NL)
            ENDDO
         ENDIF


      ! ---------- McKenny-Greengard-Mayo method:
      ! Solving for the structured and unstructured Poisson matrix
      ! Assemble both, the structured and unstructured Poisson matrix
      ! temporarily they will be stored separately in matrices AC and ACU due to the different
      ! settings along internal boundary cells,
      ! in the medium term, a toggle mechanism will be implemented which only switches the corresponding
      ! entries while keeping the entries which are the same for both discretization types

      CASE (NSCARC_METHOD_MGM)
   
         ! First assemble structured matrix with inhomogeneous boundary conditions

         TYPE_SCOPE(0) = NSCARC_SCOPE_GLOBAL
         CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_STRUCTURED)
         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         ! Then assemble unstructured matrix with homogeneous Dirichlet boundary conditions along
         ! external boundaries and special MGM BC-settings along mesh interfaces

         CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
         IF (SCARC_MGM_CHECK_LAPLACE .OR. SCARC_MGM_EXACT_INITIAL) THEN
            CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
            CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
         ENDIF

         TYPE_SCOPE(0) = NSCARC_SCOPE_LOCAL
         CALL SCARC_SETUP_LAPLACE (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY_WITH_INTERFACES(NM, NLEVEL_MIN) 
         CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_STRUCTURED)

   END SELECT SELECT_SCARC_METHOD

ENDDO MESHES_POISSON_LOOP

! Setup mappings for the global numbering of vectors and the Poisson matrix (compact storage technique only)
 
IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
   IF (IS_MGM) THEN

      TYPE_SCOPE = NSCARC_SCOPE_GLOBAL
      CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_STRUCTURED)
      CALL SCARC_SETUP_GLOBAL_CELL_MAPPING(NLEVEL_MIN)
      CALL SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NLEVEL_MIN)

      IF (SCARC_MGM_CHECK_LAPLACE .OR. SCARC_MGM_EXACT_INITIAL) THEN
         CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)
         CALL SCARC_SETUP_GLOBAL_CELL_MAPPING(NLEVEL_MIN)
         CALL SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NLEVEL_MIN)
      ENDIF

   ELSE

      CALL SCARC_SETUP_GLOBAL_CELL_MAPPING(NLEVEL_MIN)
      CALL SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NLEVEL_MIN)

   ENDIF
ENDIF
 
! If there is more than one mesh, exchange matrix values in overlapping parts
! This must be done for all multilevel methods at least at the finest grid level
! Furthermore also at all higher levels except for the AMG method,
! in this case it will be done later in routine SETUP_ALGEBRAIC_MULTIGRID

IF (SET_MATRIX_TYPE(NLEVEL_MIN) == NSCARC_MATRIX_COMPACT) CALL SCARC_SETUP_GLOBAL_POISSON_OVERLAPS(NLEVEL_MIN)

MULTI_LEVEL_IF: IF (HAS_MULTIPLE_LEVELS .AND. .NOT.HAS_AMG_LEVELS) THEN
   DO NL = NLEVEL_MIN+1, NLEVEL_MAX
      IF (SET_MATRIX_TYPE(NL) /= NSCARC_MATRIX_COMPACT) CYCLE
      CALL SCARC_SETUP_GLOBAL_CELL_MAPPING(NL)
      CALL SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NL)
      CALL SCARC_SETUP_GLOBAL_POISSON_OVERLAPS(NL)
   ENDDO 
ENDIF MULTI_LEVEL_IF

! ------ If MKL-solver is used on specific levels, then setup symmetric Poisson matrix there
#ifdef WITH_MKL
IF (.NOT. IS_MGM) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)
      IF (TYPE_PRECON == NSCARC_RELAX_OPTIMIZED .OR. TYPE_SMOOTH == NSCARC_RELAX_OPTIMIZED) THEN
         IF (.NOT.L%HAS_OBSTRUCTIONS) TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_NONE
      ENDIF
      IF (TYPE_MKL(NLEVEL_MIN) /= NSCARC_MKL_NONE) THEN
         CALL SCARC_SETUP_MATRIX_MKL(NSCARC_MATRIX_POISSON, NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY_MKL(NSCARC_MATRIX_POISSON, NM, NLEVEL_MIN)
      ENDIF

      IF (HAS_GMG_LEVELS) THEN
         DO NL = NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_POINT_TO_GRID (NM, NL)
            IF (TYPE_PRECON == NSCARC_RELAX_OPTIMIZED .OR. TYPE_SMOOTH == NSCARC_RELAX_OPTIMIZED) THEN
               IF (.NOT.L%HAS_OBSTRUCTIONS) TYPE_MKL(NL) = NSCARC_MKL_NONE
            ENDIF
            IF (TYPE_MKL(NL) /= NSCARC_MKL_NONE)  CALL SCARC_SETUP_MATRIX_MKL  (NSCARC_MATRIX_POISSON, NM, NL)
            IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL) CALL SCARC_SETUP_BOUNDARY_MKL(NSCARC_MATRIX_POISSON, NM, NL)
         ENDDO
      ENDIF

   ENDDO 

ELSE 

   CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)

   TYPE_SCOPE_SAVE(0:1) = TYPE_SCOPE(0:1)
   TYPE_MKL_SAVE(0:1) = TYPE_MKL(0:1)
   IF (SCARC_MGM_CHECK_LAPLACE .OR. SCARC_MGM_EXACT_INITIAL) THEN
      TYPE_SCOPE(0) = NSCARC_SCOPE_GLOBAL
      IF (TRIM(SCARC_PRECON) == 'CLUSTER') THEN
         TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
      ELSE
         TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
      ENDIF
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_MATRIX_MKL(NSCARC_MATRIX_POISSON, NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY_MKL(NSCARC_MATRIX_POISSON, NM, NLEVEL_MIN)
      ENDDO 
   ENDIF

   TYPE_SCOPE(0:1) = NSCARC_SCOPE_LOCAL
   TYPE_MKL(0:1) = NSCARC_MKL_LOCAL
   SELECT CASE (TYPE_MGM_LAPLACE)
      CASE (NSCARC_MGM_LAPLACE_CG, NSCARC_MGM_LAPLACE_PARDISO) 
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)
            CALL SCARC_SETUP_MATRIX_MKL (NSCARC_MATRIX_LAPLACE, NM, NLEVEL_MIN)
         ENDDO 
      CASE (NSCARC_MGM_LAPLACE_OPTIMIZED) 
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)
            IF (L%STRUCTURED%NC > L%UNSTRUCTURED%NC) CALL SCARC_SETUP_MATRIX_MKL (NSCARC_MATRIX_LAPLACE, NM, NLEVEL_MIN)
         ENDDO 
   END SELECT
   
   TYPE_SCOPE(0:1) = TYPE_SCOPE_SAVE(0:1)
   TYPE_MKL(0:1) = TYPE_MKL_SAVE(0:1)
   CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_STRUCTURED)

ENDIF
#endif

END SUBROUTINE SCARC_SETUP_SYSTEMS


! --------------------------------------------------------------------------------------------------------------
!> \brief Define sizes for system matrix A (including extended regions related to overlaps)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSON_SIZES(NL)
USE SCARC_POINTERS, ONLY: S, L, G, OG, A, OA, AB, OAB, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX, &
                          SCARC_POINT_TO_BMATRIX, SCARC_POINT_TO_OTHER_BMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, INBR

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   
   SELECT_MATRIX_TYPE: SELECT CASE (SET_MATRIX_TYPE(NL))
 
      ! -------- Matrix in compact storage technique
 
      CASE (NSCARC_MATRIX_COMPACT)

         A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)

         ! Assign IOR settings to corresponding positions in stencil

         IF (TWO_D) THEN
            A%N_STENCIL = 5
            A%POS(-3:3) = (/1,0,2,3,4,0,5/)     
         ELSE
            A%N_STENCIL = 7
            A%POS(-3:3) = (/1,2,3,4,5,6,7/)
         ENDIF

         A%N_VAL  = G%NCE * A%N_STENCIL
         A%N_ROW  = G%NCE + 1

         ! Allocate matrices on overlapping parts for later data exchanges with neighbors

         DO INBR = 1, SCARC(NM)%N_NEIGHBORS
            NOM = S%NEIGHBORS(INBR)
            CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
            OA => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_POISSON)
            OA%N_STENCIL = A%N_STENCIL
            OA%N_VAL = 4 * OG%NCG * A%N_STENCIL            ! TODO: CHECK LENGTH
            OA%N_ROW = OG%NCG + 1
         ENDDO
 
      ! -------- Matrix in bandwise storage technique
 
      CASE (NSCARC_MATRIX_BANDWISE)

         AB => SCARC_POINT_TO_BMATRIX (NSCARC_MATRIX_POISSON)

         IF (TWO_D) THEN
   
            AB%N_STENCIL   = 5                      ! 5-point Laplacian
            AB%POS(-3:3)   = (/5,0,4,3,2,0,1/)      ! assignment of IOR settings to columns in matrix array
   
            AB%OFFSET( 3)  = -L%NX                  ! lower z
            AB%OFFSET( 1)  = -1                     ! lower x
            AB%OFFSET( 0)  =  0                     ! diag
            AB%OFFSET(-1)  =  1                     ! upper x
            AB%OFFSET(-3)  =  L%NX                  ! upper z
   
            AB%SOURCE( 3)   =  1                    ! lower z
            AB%SOURCE( 1)   =  1                    ! lower x
            AB%SOURCE( 0)   =  1                    ! diag
            AB%SOURCE(-1)   =  2                    ! upper x
            AB%SOURCE(-3)   =  L%NX+1               ! upper z
   
            AB%TARGET( 3)   =  L%NX+1               ! lower z
            AB%TARGET( 1)   =  2                    ! lower x
            AB%TARGET( 0)   =  1                    ! diag
            AB%TARGET(-1)   =  1                    ! upper x
            AB%TARGET(-3)   =  1                    ! upper z
   
            AB%LENGTH( 3)  =  G%NC - L%NX           ! lower z
            AB%LENGTH( 1)  =  G%NC - 1              ! lower x
            AB%LENGTH( 0)  =  G%NC                  ! diag
            AB%LENGTH(-1)  =  G%NC - 1              ! upper x
            AB%LENGTH(-3)  =  G%NC - L%NX           ! upper z
   
         ELSE
   
            AB%N_STENCIL   = 7                      ! 7-point Laplacian
            AB%POS(-3:3)   = (/7,6,5,4,3,2,1/)      ! assignment of IOR settings to columns in matrix array
   
            AB%OFFSET( 3)  = -L%NX*L%NY             ! lower z
            AB%OFFSET( 2)  = -L%NX                  ! lower y
            AB%OFFSET( 1)  = -1                     ! lower x
            AB%OFFSET( 0)  =  0                     ! diag
            AB%OFFSET(-1)  =  1                     ! upper x
            AB%OFFSET(-2)  =  L%NX                  ! upper y
            AB%OFFSET(-3)  =  L%NX*L%NY             ! upper z
   
            AB%SOURCE( 3)  =  1                     ! lower z
            AB%SOURCE( 2)  =  1                     ! lower y
            AB%SOURCE( 1)  =  1                     ! lower x
            AB%SOURCE( 0)  =  1                     ! diag
            AB%SOURCE(-1)  =  2                     ! upper x
            AB%SOURCE(-2)  =  L%NX+1                ! upper y
            AB%SOURCE(-3)  =  L%NX*L%NY+1           ! upper z
   
            AB%TARGET( 3)  =  L%NX*L%NY+1           ! lower z
            AB%TARGET( 2)  =  L%NX+1                ! lower y
            AB%TARGET( 1)  =  2                     ! lower x
            AB%TARGET( 0)  =  1                     ! diag
            AB%TARGET(-1)  =  1                     ! upper x
            AB%TARGET(-2)  =  1                     ! upper y
            AB%TARGET(-3)  =  1                     ! upper z
   
            AB%LENGTH( 3)  =  G%NC - L%NX*L%NY      ! lower z
            AB%LENGTH( 2)  =  G%NC - L%NX           ! lower y
            AB%LENGTH( 1)  =  G%NC - 1              ! lower x
            AB%LENGTH( 0)  =  G%NC                  ! diag
            AB%LENGTH(-1)  =  G%NC - 1              ! upper x
            AB%LENGTH(-2)  =  G%NC - L%NX           ! upper y
            AB%LENGTH(-3)  =  G%NC - L%NX*L%NY      ! upper z
   
         ENDIF

         AB%N_VAL  = G%NC * AB%N_STENCIL
         AB%N_DIAG = G%NC

         ! Determine sizes of overlapping parts for later communication with corresponding neighbors

         DO INBR = 1, SCARC(NM)%N_NEIGHBORS
            NOM = S%NEIGHBORS(INBR)
            CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
            OAB => SCARC_POINT_TO_OTHER_BMATRIX (NSCARC_MATRIX_POISSON)
            OAB%N_STENCIL = AB%N_STENCIL
            OAB%N_VAL     = OG%NCG * AB%N_STENCIL
            OAB%N_DIAG    = OG%NCG 
         ENDDO

   END SELECT SELECT_MATRIX_TYPE
   
ENDDO MESHES_LOOP
   
! -------- Exchange matrix sizes in case of a multi-mesh geometry
 
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SIZES, NSCARC_MATRIX_POISSON, NL)

END SUBROUTINE SCARC_SETUP_POISSON_SIZES


! --------------------------------------------------------------------------------------------------------------
!> \brief Define sizes for local unstructured Laplace matrices
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LOCAL_LAPLACE_SIZES(NL)
USE SCARC_POINTERS, ONLY: G, A, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
  
   CALL SCARC_POINT_TO_GRID (NM, NL)
   A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE)

   IF (TWO_D) THEN
      A%N_STENCIL = 5
      A%POS(-3:3) = (/1,0,2,3,4,0,5/)
   ELSE
      A%N_STENCIL = 7
      A%POS(-3:3) = (/1,2,3,4,5,6,7/)
   ENDIF

   A%N_VAL  = G%NC * A%N_STENCIL
   A%N_ROW  = G%NC + 1

ENDDO MESHES_LOOP
  
END SUBROUTINE SCARC_SETUP_LOCAL_LAPLACE_SIZES


! -------------------------------------------------------------------------------------------------------------
!> \brief Get global numberings for compact column vector of Poisson matrix 
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NL)
USE SCARC_POINTERS, ONLY: G, A, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC, ICOL, JC

IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      CALL SCARC_POINT_TO_GRID (NM, NL)                                    
      A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
      A%COLG = A%COL
   ENDDO
ELSE
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      CALL SCARC_POINT_TO_GRID (NM, NL)                                    
      A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
      DO IC = 1, G%NC
         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
            JC = A%COL(ICOL)
            A%COLG(ICOL) = G%LOCAL_TO_GLOBAL(JC)
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_GLOBAL_POISSON_COLUMNS


! -------------------------------------------------------------------------------------------------------------
!> \brief Make Poisson matrix global by exchanging adjacent overlaps
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBAL_POISSON_OVERLAPS(NL)
USE SCARC_POINTERS, ONLY: S, A, OA, SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, INBR, NOM

IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) RETURN

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_POISSON, NL)
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLSG, NSCARC_MATRIX_POISSON, NL)
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_POISSON, NL)

CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_POISSON, 1, NL)

MESHES_FINE_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)    
   A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
   CALL SCARC_REDUCE_CMATRIX (A, 'G%POISSON', CROUTINE)

   OMESHES_FINE_LOOP: DO INBR = 1, S%N_NEIGHBORS
      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
      OA => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_POISSON)
      CALL SCARC_REDUCE_CMATRIX (OA, 'OG%POISSON', CROUTINE)
   ENDDO OMESHES_FINE_LOOP

ENDDO MESHES_FINE_LOOP
    
END SUBROUTINE SCARC_SETUP_GLOBAL_POISSON_OVERLAPS


! --------------------------------------------------------------------------------------------------------------
!> \brief Check if specified cell is within a given mesh
! --------------------------------------------------------------------------------------------------------------
LOGICAL FUNCTION SCARC_CELL_WITHIN_MESH(G, NM, IC)
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: NM, IC
INTEGER :: IC_START, IC_STOP

SCARC_CELL_WITHIN_MESH = .FALSE.
IC_START = G%NC_OFFSET(NM) + 1
IF (NM < NMESHES) THEN
   IC_STOP  = G%NC_OFFSET(NM+1)
ELSE
   IC_STOP  = G%NC_GLOBAL
ENDIF
IF (IC_START <=  IC .AND. IC <= IC_STOP) SCARC_CELL_WITHIN_MESH = .TRUE.
RETURN

END FUNCTION SCARC_CELL_WITHIN_MESH


! --------------------------------------------------------------------------------------------------------------
!> \brief Allocate Poisson matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
! Compact storage technique (POISSON)
!    Compression technique to store sparse matrices, non-zero entries are stored
!    in a 1D-vector B(.), row after row,
!    Each row starts with its diagonal entry followed by the other non-zero entries
!    In order to identify each element, pointer arrays ROW and COL are needed,
!    ROW points to the several diagonal entries in vector B(.),
!    COL points to the columns which non-zero entries in the matrix stencil
! Bandwise storage technique (POISSONB)
!    explanation to come ...
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSON (NM, NL)
USE SCARC_POINTERS, ONLY: S, L, G, A, AB, OA, OAB, &
                          SCARC_POINT_TO_GRID,    SCARC_POINT_TO_OTHER_GRID, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX, &
                          SCARC_POINT_TO_BMATRIX, SCARC_POINT_TO_OTHER_BMATRIX
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IX, IY, IZ, IC, IP, INBR, NOM

CROUTINE = 'SCARC_SETUP_POISSON'
 
! Compute single matrix entries and corresponding row and column pointers
! Along internal boundaries use placeholders for the neighboring matrix entries
! which will be communicated in a following step
 
SELECT_STORAGE_TYPE: SELECT CASE (SET_MATRIX_TYPE(NL))

 
   ! ---------- COMPACT Storage technique
 
   CASE (NSCARC_MATRIX_COMPACT)
   
      ! Allocate main matrix on non-overlapping part of mesh

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    
      A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
      CALL SCARC_ALLOCATE_CMATRIX (A, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'G%POISSON', CROUTINE)

      ! For every neighbor allocate small matrix on overlapping part of mesh

      DO INBR = 1, SCARC(NM)%N_NEIGHBORS
         NOM = S%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
         OA => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_POISSON)
         CALL SCARC_ALLOCATE_CMATRIX (OA, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OG%POISSON', CROUTINE)
      ENDDO

      IP = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
               IC = G%CELL_NUMBER(IX, IY, IZ)

               ! Main diagonal 

               CALL SCARC_SETUP_MAINDIAG (IC, IX, IY, IZ, IP)
   
               ! Lower subdiagonals

               IF (IS_VALID_DIRECTION(IX, IY, IZ,  3)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ-1, IP,  3)
               IF (IS_VALID_DIRECTION(IX, IY, IZ,  2)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY-1, IZ  , IP,  2)
               IF (IS_VALID_DIRECTION(IX, IY, IZ,  1)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX-1, IY  , IZ  , IP,  1)
   
               ! Upper subdiagonals

               IF (IS_VALID_DIRECTION(IX, IY, IZ, -1)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX+1, IY  , IZ  , IP, -1)
               IF (IS_VALID_DIRECTION(IX, IY, IZ, -2)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY+1, IZ  , IP, -2)
               IF (IS_VALID_DIRECTION(IX, IY, IZ, -3)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ+1, IP, -3)
   
            ENDDO
         ENDDO
      ENDDO
   
      A%ROW(G%NC+1) = IP
      A%N_VAL = IP
   
      CALL SCARC_GET_MATRIX_STENCIL_MAX(A, G%NC)

 
   ! ---------- bandwise storage technique
 
   CASE (NSCARC_MATRIX_BANDWISE)
   
      ! Allocate main matrix on non-overlapping part of mesh

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    
      AB => SCARC_POINT_TO_BMATRIX (NSCARC_MATRIX_POISSON)
      CALL SCARC_ALLOCATE_BMATRIX (AB, NL, 'G%POISSONB', CROUTINE)
   
      ! For every neighbor allocate little matrix on overlapping part of mesh

      DO INBR = 1, SCARC(NM)%N_NEIGHBORS
         NOM = S%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
         OAB => SCARC_POINT_TO_BMATRIX (NSCARC_MATRIX_POISSON)
         CALL SCARC_ALLOCATE_BMATRIX (OAB, NL, 'OG%POISSONB', CROUTINE)
      ENDDO
   
      IP  = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
               IC = G%CELL_NUMBER(IX, IY, IZ)
   
               ! Lower subdiagonals

               IF (IS_VALID_DIRECTION(IX, IY, IZ,  3)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX  , IY  , IZ-1,  3)
               IF (IS_VALID_DIRECTION(IX, IY, IZ,  2)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX  , IY-1, IZ  ,  2)
               IF (IS_VALID_DIRECTION(IX, IY, IZ,  1)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX-1, IY  , IZ  ,  1)
   
               ! Main diagonal

               CALL SCARC_SETUP_MAINDIAGB (IC, IX, IY, IZ)

               ! Upper subdiagonals

               IF (IS_VALID_DIRECTION(IX, IY, IZ, -1)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX+1, IY  , IZ  , -1)
               IF (IS_VALID_DIRECTION(IX, IY, IZ, -2)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX  , IY+1, IZ  , -2)
               IF (IS_VALID_DIRECTION(IX, IY, IZ, -3)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX  , IY  , IZ+1, -3)
   
            ENDDO
         ENDDO
      ENDDO
   
END SELECT SELECT_STORAGE_TYPE

END SUBROUTINE SCARC_SETUP_POISSON


! --------------------------------------------------------------------------------------------------------------
!> \brief Assemble local unstructured Laplace matrices
! The grid numbering is permuted in such a way that all the nonzero entries of the RHS 
! are located of the end of the corresponding vector
! this concerns the entries along internal obstructions and in case of a multi-mesh computation
! also the entries along the internal interfaces
! All other entries of the RHS are zero for the local Laplace problems, such that the
! forward substitution process Ly=b only must start from the nonzero entries on
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LAPLACE (NM, NL)
USE SCARC_POINTERS, ONLY: L, G, A, GWC, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IX, IY, IZ, IC, JC, KC, IOR0, IW, IP, KKC(-3:3), JJC(-3:3)
INTEGER :: TYPE_SCOPE_SAVE

CROUTINE = 'SCARC_SETUP_LAPLACE'
TYPE_SCOPE_SAVE = TYPE_SCOPE(0)
TYPE_SCOPE(0) = NSCARC_SCOPE_LOCAL
 

! Point to unstructured grid

CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)
CALL SCARC_POINT_TO_GRID (NM, NL)              

! Allocate permutation vectors 

CALL SCARC_ALLOCATE_INT1 (G%PERM_FW , 1, G%NC, NSCARC_INIT_ZERO, 'G%PERM_FW', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (G%PERM_BW , 1, G%NC, NSCARC_INIT_ZERO, 'G%PERM_BW', CROUTINE)
   
! Obstruction cells are numbered last such that they appear at the end of a vector

IF (TYPE_MGM_LAPLACE == NSCARC_MGM_LAPLACE_LUPERM) THEN

   JC = G%NC
   DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
      
      GWC => G%WALL(IW)
      
      IX = GWC%IXW ;  IY = GWC%IYW ;  IZ = GWC%IZW
      
      IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
      
      IOR0 = GWC%IOR
      IC   = G%CELL_NUMBER(IX, IY, IZ)

      G%PERM_FW(IC) = JC
      G%PERM_BW(JC) = IC
      JC = JC - 1

   ENDDO

   ! Interface cells are numbered second last

   DO IW = 1, L%N_WALL_CELLS_EXT
      
      GWC => G%WALL(IW)
      IF (GWC%BTYPE /= INTERNAL) CYCLE
      
      IX = GWC%IXW ;  IY = GWC%IYW ;  IZ = GWC%IZW
      
      IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
      
      IOR0 = GWC%IOR
      IC   = G%CELL_NUMBER(IX,IY,IZ)

      IF (G%PERM_FW(IC) == 0) THEN
         G%PERM_FW(IC) = JC
         G%PERM_BW(JC) = IC
         JC = JC - 1
      ENDIF

   ENDDO

   ! The rest is used from beginning to first interface cell

   KC = 1
   DO IC = 1, G%NC
      IF (G%PERM_FW(IC) /= 0) CYCLE
      G%PERM_BW(KC) = IC
      G%PERM_FW(IC) = KC
      KC = KC + 1
   ENDDO
   IF (KC /= JC + 1) THEN
      WRITE(*,*) 'KC =', KC,': JC =', JC
      CALL SCARC_ERROR(NSCARC_ERROR_MGM_PERMUTATION, SCARC_NONE, NSCARC_NONE)
   ENDIF

   G%NONZERO = KC

ELSE

   DO IC = 1, G%NC
      G%PERM_BW(IC) = IC
      G%PERM_FW(IC) = IC
   ENDDO

   G%NONZERO = 1

ENDIF

! Allocate Laplace matrix on non-overlapping part of mesh

A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE)
CALL SCARC_ALLOCATE_CMATRIX (A, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'G%LAPLACE', CROUTINE)

! Assemble Laplace matrix with grid permutation based on MGM-method 

IP = 1
IF (TYPE_MGM_LAPLACE == NSCARC_MGM_LAPLACE_LUPERM) THEN

   DO IC = 1, G%NC

      JJC = -1
      KKC = -1

      JJC(0) = G%PERM_BW(IC);  KKC(0) = G%PERM_FW(JJC(0))
      
      IX = G%ICX(JJC(0)); IY = G%ICY(JJC(0)); IZ = G%ICZ(JJC(0))

      JJC(-3) = G%CELL_NUMBER(IX  , IY, IZ+1)     ; KKC(-3) = GET_PERM(JJC(-3))  
      JJC(-1) = G%CELL_NUMBER(IX+1, IY, IZ  )     ; KKC(-1) = GET_PERM(JJC(-1))   
      JJC( 1) = G%CELL_NUMBER(IX-1, IY, IZ  )     ; KKC( 1) = GET_PERM(JJC( 1))    
      JJC( 3) = G%CELL_NUMBER(IX  , IY, IZ-1)     ; KKC( 3) = GET_PERM(JJC( 3))     
      IF (.NOT.TWO_D) THEN
        JJC(-2) = G%CELL_NUMBER(IX, IY+1, IZ)     ; KKC(-2) = GET_PERM(JJC(-2))     
        JJC( 2) = G%CELL_NUMBER(IX, IY-1, IZ)     ; KKC( 2) = GET_PERM(JJC( 2))     
      ENDIF

      ! Main diagonal 
      CALL SCARC_SETUP_MAINDIAG (IC, IX, IY, IZ, IP)
         
      ! Lower subdiagonals

      IF (IS_VALID_DIRECTION(IX, IY, IZ,  3)) CALL SCARC_SETUP_SUBDIAG_PERM(IX, IY, IZ, IX  , IY  , IZ-1, KKC( 3), IP,  3)
      IF (IS_VALID_DIRECTION(IX, IY, IZ,  2)) CALL SCARC_SETUP_SUBDIAG_PERM(IX, IY, IZ, IX  , IY-1, IZ  , KKC( 2), IP,  2)
      IF (IS_VALID_DIRECTION(IX, IY, IZ,  1)) CALL SCARC_SETUP_SUBDIAG_PERM(IX, IY, IZ, IX-1, IY  , IZ  , KKC( 1), IP,  1)

      ! Upper subdiagonals

      IF (IS_VALID_DIRECTION(IX, IY, IZ, -1)) CALL SCARC_SETUP_SUBDIAG_PERM(IX, IY, IZ, IX+1, IY  , IZ  , KKC(-1), IP, -1)
      IF (IS_VALID_DIRECTION(IX, IY, IZ, -2)) CALL SCARC_SETUP_SUBDIAG_PERM(IX, IY, IZ, IX  , IY+1, IZ  , KKC(-2), IP, -2)
      IF (IS_VALID_DIRECTION(IX, IY, IZ, -3)) CALL SCARC_SETUP_SUBDIAG_PERM(IX, IY, IZ, IX  , IY  , IZ+1, KKC(-3), IP, -3)

   ENDDO
      
! Assemble Laplace matrix without grid permutation 

ELSE

   DO IZ = 1, L%NZ
      DO IY = 1, L%NY
         DO IX = 1, L%NX

            IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
            IC = G%CELL_NUMBER(IX, IY, IZ)

            ! Main diagonal 

            CALL SCARC_SETUP_MAINDIAG (IC, IX, IY, IZ, IP)

            ! Lower subdiagonals

            IF (IS_VALID_DIRECTION(IX, IY, IZ,  3)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ-1, IP,  3)
            IF (IS_VALID_DIRECTION(IX, IY, IZ,  2)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY-1, IZ  , IP,  2)
            IF (IS_VALID_DIRECTION(IX, IY, IZ,  1)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX-1, IY  , IZ  , IP,  1)

            ! Upper subdiagonals

            IF (IS_VALID_DIRECTION(IX, IY, IZ, -1)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX+1, IY  , IZ  , IP, -1)
            IF (IS_VALID_DIRECTION(IX, IY, IZ, -2)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY+1, IZ  , IP, -2)
            IF (IS_VALID_DIRECTION(IX, IY, IZ, -3)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ+1, IP, -3)

         ENDDO
      ENDDO
   ENDDO
      
ENDIF

A%ROW(G%NC+1) = IP
A%N_VAL = IP
      
CALL SCARC_GET_MATRIX_STENCIL_MAX(A, G%NC)

TYPE_SCOPE(0) = TYPE_SCOPE_SAVE

 
END SUBROUTINE SCARC_SETUP_LAPLACE


! --------------------------------------------------------------------------------------------------------------
!> \brief Set main diagonal entry for Poisson matrix in compact storage technique
! These values correspond to the full matrix of the global problem
! In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
! If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MAINDIAG (IC, IX, IY, IZ, IP)
USE SCARC_POINTERS, ONLY: L, A
INTEGER, INTENT(IN) :: IC, IX, IY, IZ
INTEGER, INTENT(INOUT) :: IP

A%VAL(IP) = - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D) A%VAL(IP) = A%VAL(IP) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
A%VAL(IP) = A%VAL(IP) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

A%ROW(IC) = IP
A%COL(IP) = IC

A%STENCIL(0) = A%VAL(IP)

IP = IP + 1
END SUBROUTINE SCARC_SETUP_MAINDIAG


! --------------------------------------------------------------------------------------------------------------
!> \brief Determine if cell has a neighbor and, if yes, return corresponding wall cell index
! --------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)
USE SCARC_POINTERS, ONLY: L, G, F, GWC
INTEGER, INTENT(IN) :: IC, IOR0
INTEGER :: IXW, IYW, IZW
INTEGER :: IXG, IYG, IZG
INTEGER :: IW

SCARC_ASSIGN_SUBDIAG_TYPE = -1

F => L%FACE(IOR0)
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


! --------------------------------------------------------------------------------------------------------------
!> \brief Set subdigonal entries for Poisson matrix in compact storage technique on specified face
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIAG (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IP, IOR0)
USE SCARC_POINTERS, ONLY: L, F, G, A
INTEGER, INTENT(IN) :: IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0
INTEGER, INTENT(INOUT) :: IP
INTEGER :: IW
LOGICAL :: IS_INTERNAL_CELL

! Decide wheter cell is interior or exterior cell

F => L%FACE(IOR0)
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

! If IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
IF (IS_INTERNAL_CELL) THEN

   IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX2, IY2, IZ2)) THEN
      A%VAL(IP) = A%VAL(IP) + F%INCR_INSIDE
      A%COL(IP) = G%CELL_NUMBER(IX2, IY2, IZ2)
      A%STENCIL(-IOR0) = A%VAL(IP)
      IP = IP + 1
   ELSE
      CALL SCARC_ERROR(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

! If IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell

ELSE IF (TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL .AND. L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN

   IW = SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)           ! get IW of a possibly suitable neighbor at face IOR0
   IF (IW > 0) then                                    ! if available, build corresponding subdiagonal entry
      A%VAL(IP) = A%VAL(IP) + F%INCR_FACE
      A%COL(IP) = G%WALL(IW)%ICE                       ! store its extended number in matrix column pointers
      A%STENCIL(-IOR0) = A%VAL(IP)
      IP = IP + 1
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_SUBDIAG


! --------------------------------------------------------------------------------------------------------------
!> \brief Set subdigonal entries for permuted matrix in compact storage technique on specified face
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIAG_PERM (IX1, IY1, IZ1, IX2, IY2, IZ2, ICOL, IP, IOR0)
USE SCARC_POINTERS, ONLY: L, F, A
INTEGER, INTENT(IN) :: IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0, ICOL
INTEGER, INTENT(INOUT) :: IP
LOGICAL :: IS_INTERNAL_CELL

!A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)

! Decide wheter cell is interior or exterior cell

F => L%FACE(IOR0)
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

! If IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
IF (IS_INTERNAL_CELL .AND. .NOT.L%IS_SOLID(IX2, IY2, IZ2)) THEN
   A%VAL(IP) = A%VAL(IP) + F%INCR_INSIDE
   A%COL(IP) = ICOL
   A%STENCIL(-IOR0) = A%VAL(IP)
   IP = IP + 1
ENDIF

END SUBROUTINE SCARC_SETUP_SUBDIAG_PERM


! --------------------------------------------------------------------------------------------------------------
!> \brief Set boundary conditions including the interfaces between the meshes
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BOUNDARY_WITH_INTERFACES (NM, NL)
USE SCARC_POINTERS, ONLY: L, G, F, GWC, A, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP
INTEGER :: ICXM, ICXP, ICYM, ICYP, ICZM, ICZP

CALL SCARC_POINT_TO_GRID (NM, NL)       

A => G%LAPLACE

SELECT CASE (TYPE_MGM_BC)

   ! --------------------------------------------------------------------------
   CASE (NSCARC_MGM_BC_TAYLOR)

      !DO IW = MGM%NW1, MGM%NW2
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

         ICXM  = G%CELL_NUMBER(I-1, J, K)
         ICXP  = G%CELL_NUMBER(I+1, J, K)
         IF (.NOT. TWO_D) THEN
            ICYM  = G%CELL_NUMBER(I, J-1, K)
            ICYP  = G%CELL_NUMBER(I, J+1, K)
         ENDIF
         ICZM  = G%CELL_NUMBER(I, J, K-1)
         ICZP  = G%CELL_NUMBER(I, J, K+1)

         GWC%ICW = IC

         IP = A%ROW(IC)
         A%VAL(IP) = 0.0_EB
         SELECT CASE(ABS(IOR0))
            CASE (1)
               A%VAL(IP) = A%VAL(IP) - L%DXI2 - 5.0_EB/2.0_EB*L%DZI2
               DO IP = A%ROW(IC)+1, A%ROW(IC+1)-1
                  IF (ICXM <= G%NC .AND. A%COL(IP) == ICXM) THEN
                     A%VAL(IP) = A%VAL(IP) + L%DXI2 
                  ENDIF
                  IF (ICXP <= G%NC .AND. A%COL(IP) == ICXP) THEN
                     A%VAL(IP) = A%VAL(IP) + L%DXI2 
                  ENDIF
                  IF (ICZM <= G%NC .AND. A%COL(IP) == ICZM) THEN
                     A%VAL(IP) = A%VAL(IP) + 5.0_EB/4.0_EB*L%DZI2
                  ENDIF
                  IF (ICZP <= G%NC .AND. A%COL(IP) == ICZP) THEN
                        A%VAL(IP) = A%VAL(IP) + 5.0_EB/4.0_EB*L%DZI2
                  ENDIF
               ENDDO
            CASE (2)
               IF (.NOT. TWO_D) THEN
                  WRITE(*,*) 'TAYLOR-3D: Not yet finished!'
               ENDIF
            CASE (3)
               A%VAL(IP) = A%VAL(IP) - L%DZI2 - 5.0_EB/2.0_EB*L%DXI2
               DO IP = A%ROW(IC)+1, A%ROW(IC+1)-1
                  IF (ICZM <= G%NC .AND. A%COL(IP) == ICZM) THEN
                     A%VAL(IP) = A%VAL(IP) + L%DZI2 
                  ENDIF
                  IF (ICZP <= G%NC .AND. A%COL(IP) == ICZP) THEN
                     A%VAL(IP) = A%VAL(IP) + L%DZI2 
                  ENDIF
                  IF (ICXM <= G%NC .AND. A%COL(IP) == ICXM) THEN
                     A%VAL(IP) = A%VAL(IP) + 5.0_EB/4.0_EB*L%DXI2
                  ENDIF
                  IF (ICXP <= G%NC .AND. A%COL(IP) == ICXP) THEN
                     A%VAL(IP) = A%VAL(IP) + 5.0_EB/4.0_EB*L%DXI2
                  ENDIF
               ENDDO
         END SELECT

      ENDDO 

   ! --------------------------------------------------------------------------
   CASE DEFAULT

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
         IF (TYPE_MGM_LAPLACE == NSCARC_MGM_LAPLACE_LUPERM) THEN
            IC  = G%PERM_FW(G%CELL_NUMBER(I, J, K))
         ELSE
            IC  = G%CELL_NUMBER(I, J, K)
         ENDIF
         !GWC%ICW = IC

         IP = A%ROW(IC)
         SELECT CASE (GWC%BTYPE)
            CASE (DIRICHLET, INTERNAL)
               A%VAL(IP) = A%VAL(IP) - F%INCR_BOUNDARY
            CASE (NEUMANN)
               A%VAL(IP) = A%VAL(IP) + F%INCR_BOUNDARY
         END SELECT

      ENDDO 

END SELECT 


END SUBROUTINE SCARC_SETUP_BOUNDARY_WITH_INTERFACES


! --------------------------------------------------------------------------------------------------------------
!> \brief Set main diagonal entry for Poisson matrix in bandwise storage technique
! These values correspond to the full matrix of the global problem
! In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
! If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MAINDIAGB (IC, IX, IY, IZ)
USE SCARC_POINTERS, ONLY: L, G, AB
INTEGER, INTENT(IN)  :: IC, IX, IY, IZ
INTEGER :: ID

AB => G%POISSONB
ID = AB%POS(0)               ! get column vector corresponding to matrix diagonal

AB%VAL(IC, ID) = AB%VAL(IC, ID) - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D)  AB%VAL(IC, ID) = AB%VAL(IC, ID) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
AB%VAL(IC, ID) = AB%VAL(IC, ID) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

AB%STENCIL(0) = AB%VAL(IC, ID)

END SUBROUTINE SCARC_SETUP_MAINDIAGB


! --------------------------------------------------------------------------------------------------------------
!> \brief Set subdigonal entries for Poisson matrix in bandwise storage technique on specified face
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIAGB (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0)
USE SCARC_POINTERS, ONLY: L, F, G, AB
INTEGER, INTENT(IN) :: IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0
INTEGER :: IW, ID
LOGICAL  :: IS_INTERNAL_CELL

F => L%FACE(IOR0)

! Decide wheter cell is interior or exterior cell
 
AB => G%POISSONB
ID = AB%POS(IOR0)                                
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

! If IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
 
IF (IS_INTERNAL_CELL) THEN

   IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX2, IY2, IZ2)) THEN
      AB%VAL(IC, ID)   = AB%VAL(IC, ID) + F%INCR_INSIDE
      AB%STENCIL(IOR0) = AB%VAL(IC, ID)
   ELSE
      CALL SCARC_ERROR(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

! If IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
 
!ELSE IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
ELSE IF (L%FACE(IOR0)%N_NEIGHBORS == 123456) THEN       ! CAUTION: TO FIX AGAIN, ONLY FOR TESTING, IMPOSSIBLE CONDITION

   IW = SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)            ! get IW of a possibly suitable neighbor at face IOR0
   IF (IW /= 0) THEN
      AB%VAL(IC, ID)   = AB%VAL(IC, ID) + F%INCR_FACE
      AB%STENCIL(IOR0) = AB%VAL(IC, ID)
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_SUBDIAGB


! --------------------------------------------------------------------------------------------------------------
!> \brief Get maximum stencil size in specified matrix 
! This is known to be 7 for the 3D-Poisson matrix on finest level
! In algebraic multigrid-method this size results only in the course and can be much larger
! (required for dimensioning the coarse-level matrices)
! If NTYPE == 0, only internal matrix part is considered, if NTYPE == 1, also the overlap
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_MATRIX_STENCIL_MAX (A, NLEN)
TYPE (SCARC_CMATRIX_TYPE), INTENT(INOUT) :: A
INTEGER, INTENT(IN) :: NLEN
INTEGER :: IC

A%N_STENCIL_MAX = 0
DO IC = 1, NLEN
   A%N_STENCIL_MAX = MAX(A%N_STENCIL_MAX, A%ROW(IC+1)-A%ROW(IC)+1)
ENDDO

END SUBROUTINE SCARC_GET_MATRIX_STENCIL_MAX


#ifdef WITH_MKL
! --------------------------------------------------------------------------------------------------------------
!> \brief Setup symmetric version of Poisson matrix for MKL solver in double precision
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MKL (NMATRIX, NM, NL)
USE SCARC_POINTERS, ONLY: G, A, AS, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NMATRIX, NM, NL
INTEGER :: IC, JC, JC0, ICS, JCS, JCG
INTEGER :: ICOL, JCOL, IAS
INTEGER :: ISYM, JSYM, NSYM
REAL(EB) :: VAL = 0.0_EB, VALS = 0.0_EB, DIFF
LOGICAL  :: BSYM, BCHECK_SYMMETRY = .FALSE.
INTEGER, DIMENSION(:), ALLOCATABLE :: ICOL_AUX, IC_AUX
INTEGER, POINTER, DIMENSION(:) :: ACOLG, ASCOLG

CROUTINE = 'SCARC_SETUP_MATRIX_MKL'

CALL SCARC_POINT_TO_GRID (NM, NL)                                    

SELECT CASE (NMATRIX)
   CASE (NSCARC_MATRIX_POISSON)
      A  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
      AS => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_SYM)
   CASE (NSCARC_MATRIX_LAPLACE)
      A  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE)
      AS => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE_SYM)
END SELECT

IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) THEN
   ACOLG  => A%COL
ELSE
   ACOLG  => A%COLG
ENDIF
  
! ---------- Store only symmetric parts of matrix (diagonal and upper part)
  
IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN

   IF (BCHECK_SYMMETRY) THEN
      ! First check whether symmetry of system matrix is guaranteed
      DO IC = 1, G%NC
         COLUMN_LOOP: DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
            ICS = ACOLG(ICOL)
            VAL = A%VAL(ICOL)
            IF (ICS > IC .AND. ICS <= G%NC) THEN
               BSYM = .FALSE.
               DO JCOL = A%ROW(ICS)+1, A%ROW(ICS+1)-1
                  JCS = ACOLG(JCOL)
                  IF (JCS == IC) THEN
                     VALS = A%VAL(JCOL)
                     DIFF = ABS(VAL-VALS)
                     IF (ABS(VAL - VALS) < 1E-6) THEN
                        BSYM=.TRUE.
                        CYCLE COLUMN_LOOP
                     ENDIF
                  ENDIF
               ENDDO
               IF (.NOT.BSYM) CALL SCARC_ERROR(NSCARC_ERROR_MATRIX_SYMMETRY, SCARC_NONE, NM)
            ENDIF
         ENDDO COLUMN_LOOP
      ENDDO
   ENDIF

 
   ! Compute number of entries in symmetric matrix
 
   AS%N_VAL = 0
   DO IC = 1, G%NC
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL) THEN
            JC = ACOLG(ICOL)
            IF (JC >= IC .AND. JC <= G%NC) AS%N_VAL = AS%N_VAL+1
         ELSE IF (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) THEN
            IF (NL == NLEVEL_MIN) THEN
               JCG = G%LOCAL_TO_GLOBAL(ACOLG(ICOL))
            ELSE
               JCG = ACOLG(ICOL)
            ENDIF
            IF (JCG >= IC + G%NC_OFFSET(NM)) AS%N_VAL = AS%N_VAL+1
         ELSE
            CALL SCARC_ERROR(NSCARC_ERROR_MATRIX_SETUP, SCARC_NONE, TYPE_MKL(NL))
         ENDIF
      ENDDO
   ENDDO

ELSE
   AS%N_VAL = A%N_VAL
ENDIF

! Allocate storage for symmetric matrix and its column and row pointers
  
CALL SCARC_GET_MATRIX_STENCIL_MAX(A, G%NC)
AS%N_ROW = G%NC + 1
AS%N_VAL = A%N_STENCIL_MAX * G%NC
CALL SCARC_ALLOCATE_CMATRIX (AS, NL, TYPE_MKL_PRECISION, NSCARC_MATRIX_FULL, 'G%AS', CROUTINE)

IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) THEN
   ASCOLG  => AS%COL
ELSE
   ASCOLG  => AS%COLG
ENDIF

! If global MKL method is used, also allocate auxiliary space for computation of global numbering

IF (IS_MKL_LEVEL(NL)) THEN
   CALL SCARC_ALLOCATE_INT1 (ICOL_AUX, 1, A%N_STENCIL_MAX, NSCARC_HUGE_INT, 'ICOL_AUX', CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (IC_AUX  , 1, A%N_STENCIL_MAX, NSCARC_HUGE_INT, 'IC_AUX', CROUTINE)
ENDIF
  
! Subtract symmetric matrix part from usual system matrix
  
IAS = 1
DO IC = 1, AS%N_ROW - 1
   AS%ROW(IC) = IAS

   TYPE_MKL_SELECT: SELECT CASE (TYPE_MKL(NL)) 

      ! Blockwise use of local MKL solvers - no global numbering required

      CASE(NSCARC_MKL_LOCAL) 

         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
            JC = A%COL(ICOL)
            IF (JC >= IC .AND. JC <= G%NC) THEN
               AS%COL(IAS) = A%COL(ICOL)
               ASCOLG(IAS) = A%COL(ICOL)
               SELECT CASE (TYPE_MKL_PRECISION)
                  CASE (NSCARC_PRECISION_DOUBLE)
                     AS%VAL(IAS) = A%VAL(ICOL)
                  CASE (NSCARC_PRECISION_SINGLE)
                     AS%VAL_FB(IAS) = REAL(A%VAL(ICOL),FB)
                  END SELECT
               IAS = IAS + 1
            ENDIF
         ENDDO
         AS%ROW(IC+1) = IAS

      ! Global use of MKL solver - get global numbering of matrix elements

      CASE(NSCARC_MKL_GLOBAL) 

         ! Store indices of all diagonal and upper-diagonal entries

         ICOL_AUX = 0
         IC_AUX   = NSCARC_HUGE_INT
         ISYM = 1
         JC0 = ACOLG(A%ROW(IC))
         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
            JC = ACOLG(ICOL)
            IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
               IF (JC >= JC0) THEN
                  ICOL_AUX(ISYM) = ICOL
                  IC_AUX(ISYM) = JC
                  ISYM  = ISYM  + 1
               ENDIF
            ELSE
               ICOL_AUX(ISYM) = ICOL
               IC_AUX(ISYM) = JC
               ISYM  = ISYM  + 1
            ENDIF
         ENDDO
         AS%ROW(IC+1) = IAS

         NSYM = ISYM - 1
         JSYM = 1

         ! Sort them in increasing order (for the use of Cluster_Sparse_Solver and PARDISO functionality)

         SORT_LOOP: DO WHILE (JSYM <= NSYM)
            DO ISYM = 1, NSYM
               JC = IC_AUX(ISYM)
               IF (JC == NSCARC_HUGE_INT) CYCLE
               IF (JC <= MINVAL(ABS(IC_AUX(1:NSYM)))) THEN
                  ICOL = ICOL_AUX(ISYM)
                  SELECT CASE (TYPE_MKL_PRECISION)
                     CASE (NSCARC_PRECISION_DOUBLE)
                        AS%VAL(IAS) = A%VAL(ICOL)
                     CASE (NSCARC_PRECISION_SINGLE)
                        AS%VAL_FB(IAS) = REAL(A%VAL(ICOL), FB)
                  END SELECT
                  AS%COL(IAS) = ACOLG(ICOL)
                  IC_AUX(ISYM) = NSCARC_HUGE_INT            ! mark entry as already used
                  IAS  = IAS  + 1
               ENDIF
            ENDDO
            JSYM = JSYM + 1
         ENDDO SORT_LOOP

   END SELECT TYPE_MKL_SELECT
ENDDO

AS%ROW(AS%N_ROW) = IAS

IF (IS_MKL_LEVEL(NL)) THEN
   CALL SCARC_DEALLOCATE_INT1 (ICOL_AUX, 'COL_AUX', CROUTINE)
   CALL SCARC_DEALLOCATE_INT1 (IC_AUX,  'IC_AUX', CROUTINE)
ENDIF

CALL SCARC_REDUCE_CMATRIX (AS, 'AS', CROUTINE)

END SUBROUTINE SCARC_SETUP_MATRIX_MKL
#endif


! --------------------------------------------------------------------------------------------------------------
!> \brief Insert correct boundary conditions into system matrix

! If A is a pure Neumann matrix, get neighboring cell indices of communicated stencil legs for 
! condensed system, also save values and column indices of last matrix row of last mesh

! Set correct boundary conditions for system matrix
! Take care of whether the structured or unstructured discretization is used

! If there are no Dirichlet BC's transform sytem into condensed one by replacing the
! matrix entries in last column and row by the stored ones (zeros and one at diaonal position)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
USE SCARC_POINTERS, ONLY: L, G, F, GWC, A, AB, ACO, ABCO, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP, ICO, ICOL

CALL SCARC_POINT_TO_GRID (NM, NL)                                    

SELECT CASE (SET_MATRIX_TYPE(NL))

   ! ---------- Matrix in compact storage technique
 
   CASE (NSCARC_MATRIX_COMPACT)

      A => G%POISSON

      ! Setup condensing if there are no Dirichlet BC's 

      IF (IS_PURE_NEUMANN) CALL SCARC_SETUP_CMATRIX_CONDENSED(NM)

      ! Set correct boundary conditions 

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

         ! SPD-matrix with mixture of Dirichlet and Neumann BC's according to BTYPE

         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN

            IP = A%ROW(IC)
            SELECT CASE (GWC%BTYPE)
               CASE (DIRICHLET)
                  A%VAL(IP) = A%VAL(IP) - F%INCR_BOUNDARY
               CASE (NEUMANN)
                  A%VAL(IP) = A%VAL(IP) + F%INCR_BOUNDARY
            END SELECT

         ! Purely Neumann matrix

         ELSE IF (GWC%BTYPE == NEUMANN) THEN
            IP = A%ROW(IC)
            A%VAL(IP) = A%VAL(IP) + F%INCR_BOUNDARY
         ENDIF

      ENDDO 

      ! Transform into condensed system, if there are no Dirichlet BC's 

      IF (IS_PURE_NEUMANN) THEN
         DO ICO = 1, A%N_CONDENSED
            ACO => A%CONDENSED(ICO)
            DO ICOL = 1, ACO%N_COL
               IP = ACO%PTR(ICOL)
               A%VAL(IP) = ACO%VAL2(ICOL)
            ENDDO
         ENDDO
      ENDIF 

 
   ! ---------- Matrix in bandwise storage technique
 
   CASE (NSCARC_MATRIX_BANDWISE)

      ! Preset matrix switch if no Dirichlet BC's available

      AB => G%POISSONB
      IF (IS_PURE_NEUMANN) CALL SCARC_SETUP_BMATRIX_CONDENSED(NM)

      ! Set right boundary conditions 

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

         ! SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE

         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN

            SELECT CASE (GWC%BTYPE)
               CASE (DIRICHLET)
                  AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) - F%INCR_BOUNDARY
               CASE (NEUMANN)
                  AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) + F%INCR_BOUNDARY
            END SELECT

         ! Purely Neumann matrix

         ELSE
            IF (GWC%BTYPE == NEUMANN) AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) + F%INCR_BOUNDARY
         ENDIF

      ENDDO 
   
      ! Transform into condensed system, if there are no Dirichlet BC's 

      IF (IS_PURE_NEUMANN) THEN
         DO ICO = 1, AB%N_CONDENSED
            ABCO => AB%CONDENSED(ICO)
            IF (ICO == 1) THEN
               AB%VAL(ABCO%ICO, 1:AB%N_STENCIL) = ABCO%VAL2(1:AB%N_STENCIL)
            ELSE
               IP = AB%POS(ABCO%IOR0)
               AB%VAL(ABCO%ICO, IP) = ABCO%VAL2(IP)
            ENDIF
         ENDDO
      ENDIF 
 
END SELECT 

END SUBROUTINE SCARC_SETUP_BOUNDARY


#ifdef WITH_MKL
! --------------------------------------------------------------------------------------------------------------
!> \brief Insert internal Dirichlet boundary conditions to local MKL preconditioning matrices
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BOUNDARY_MKL (NMATRIX, NM, NL)
USE SCARC_POINTERS, ONLY: L, G, F, GWC, AS, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NMATRIX, NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP

CALL SCARC_POINT_TO_GRID (NM, NL)                                    

IF (NMATRIX == NSCARC_MATRIX_POISSON) THEN
   AS => G%POISSON_SYM
ELSE
   AS => G%LAPLACE_SYM
ENDIF

SELECT CASE (TYPE_MKL_PRECISION) 

   CASE (NSCARC_PRECISION_DOUBLE)

      DO IW = 1, G%NW
      
         GWC => G%WALL(IW)
         IOR0 = GWC%IOR
         IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE       
         IF (GWC%BTYPE /= INTERNAL) CYCLE
      
         F  => L%FACE(IOR0)
      
         I = GWC%IXW
         J = GWC%IYW
         K = GWC%IZW
      
         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
      
         NOM = GWC%NOM
         IC  = G%CELL_NUMBER(I, J, K)
      
         IP = AS%ROW(IC)
         AS%VAL(IP) = AS%VAL(IP) - F%INCR_BOUNDARY
      
      ENDDO 
      
   CASE (NSCARC_PRECISION_SINGLE)

      DO IW = 1, G%NW
      
         GWC => G%WALL(IW)
         IOR0 = GWC%IOR
         IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE       
         IF (GWC%BTYPE /= INTERNAL) CYCLE
      
         F  => L%FACE(IOR0)
      
         I = GWC%IXW
         J = GWC%IYW
         K = GWC%IZW
      
         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
      
         NOM = GWC%NOM
         IC  = G%CELL_NUMBER(I, J, K)
      
         IP = AS%ROW(IC)
         AS%VAL_FB(IP) = AS%VAL_FB(IP) - REAL(F%INCR_BOUNDARY, FB)
      
      ENDDO 
END SELECT

 
END SUBROUTINE SCARC_SETUP_BOUNDARY_MKL
#endif


! --------------------------------------------------------------------------------------------------------------
!> \brief Setup condensed system for compact matrix storage technique
! Define switch entries for toggle between original and condensed values
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CMATRIX_CONDENSED (NM)
USE SCARC_POINTERS, ONLY: L, G, A, ACO, GWC
INTEGER, INTENT(IN) :: NM
INTEGER :: ICO = 0, NC, NOM, IP, IC, JC, ICE, ICN, ICOL, IOR0, IW, I, J, K

A => G%POISSON
LAST_CELL_IN_LAST_MESH_IF: IF (NM == NMESHES) THEN

   NC = G%NC_LOCAL(NMESHES)
   IP = A%ROW(NC)

   ! Store column indices and values of diagonal and all off-diagonal entries in last row
   ! index '1' corresponds to main diagonal entry

   ICO = ICO + 1
   ACO => A%CONDENSED(ICO)

   ICOL = 1
   ACO%PTR(ICOL)  = IP
   ACO%COL(ICOL)  = A%COL(IP)
   ACO%VAL1(ICOL) = A%VAL(IP)
   ACO%VAL2(ICOL) = 1.0_EB

   DO IP = A%ROW(NC)+1, A%ROW(NC+1)-1
      ICOL = ICOL + 1
      ACO%PTR(ICOL)  = IP
      ACO%COL(ICOL)  = A%COL(IP)
      ACO%VAL1(ICOL) = A%VAL(IP)
      ACO%VAL2(ICOL) = 0.0_EB
   ENDDO
   ACO%N_COL = ICOL                                ! number of stored columns

   ! Within last mesh: check which other cells have a connection to the last cell;
   ! in each corresponding matrix row store the column index and value of just that matrix entry
   ! for each direction only one value has to be stored
 
   JC = NC - 1
   DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
      IF (A%COL(IP) == NC) THEN
         ICO = ICO + 1
         ACO => A%CONDENSED(ICO)
         ACO%PTR(1)  = IP
         ACO%COL(1)  = JC
         ACO%VAL1(1) = A%VAL(IP)                     ! store original value of system matrix
         ACO%VAL2(1) = 0.0_EB                        ! store new value of condensed system matrix
         ACO%N_COL   = 1
         EXIT
      ENDIF
   ENDDO

   JC = NC - L%NX
   DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
      IF (A%COL(IP) == NC) THEN
         ICO = ICO + 1
         ACO => A%CONDENSED(ICO)
         ACO%PTR(1)  = IP
         ACO%COL(1)  = JC
         ACO%VAL1(1) = A%VAL(IP)                     ! store original value of system matrix
         ACO%VAL2(1) = 0.0_EB                        ! store new value of condensed system matrix
         ACO%N_COL   = 1
         EXIT
      ENDIF
   ENDDO

   IF (.NOT.TWO_D) THEN
      JC = NC - L%NX * L%NY
      DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
         IF (A%COL(IP) == NC) THEN
            ICO = ICO + 1
            ACO => A%CONDENSED(ICO)
            ACO%PTR(1)  = IP
            ACO%COL(1)  = JC
            ACO%VAL1(1) = A%VAL(IP)                  ! store original value of system matrix
            ACO%VAL2(1) = 0.0_EB                     ! store new value of condensed system matrix
            ACO%N_COL   = 1
            EXIT
         ENDIF
      ENDDO
   ENDIF

ENDIF LAST_CELL_IN_LAST_MESH_IF

! Cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled
! with the last cell of last mesh;
! this can be a cell on the opposite side of the own mesh or on a different mesh
! if such a cell exists, store corresponding matrix entry
 
DO IW = 1, G%NW

   GWC => G%WALL(IW)

   IOR0 = GWC%IOR
   IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE

   I = GWC%IXW
   J = GWC%IYW
   K = GWC%IZW

   IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

   NOM = GWC%NOM
   IC  = G%CELL_NUMBER(I, J, K)
   GWC%ICW = IC

   IF (NOM == NMESHES) THEN

      ICE = GWC%ICE                               ! adjacent ghost cell number
      ICN = G%ICE_TO_ICN(ICE)                     ! get column index of neighboring offdiagonal matrix entry
      IF (ICN /= SCARC(NMESHES)%NC) CYCLE         ! if no relation to last cell in last mesh, cycle

      DO IP = A%ROW(IC)+1, A%ROW(IC+1)-1
         IF (A%COL(IP) == ICE) THEN
            ICO = ICO + 1
            ACO => A%CONDENSED(ICO)
            ACO%PTR(1)  = IP
            ACO%COL(1)  = ICN
            ACO%VAL1(1) = A%VAL(IP)
            ACO%VAL2(1) = 0.0_EB
            ACO%N_COL   = 1
            EXIT
         ENDIF
      ENDDO

   ENDIF 
ENDDO 

A%N_CONDENSED = ICO

END SUBROUTINE SCARC_SETUP_CMATRIX_CONDENSED


! --------------------------------------------------------------------------------------------------------------
!> \brief Setup condensed system for bandwise matrix storage technique
! Define switch entries for toggle between original and condensed values
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BMATRIX_CONDENSED (NM)
USE SCARC_POINTERS, ONLY: L, G, AB, ABCO, GWC
INTEGER, INTENT(IN) :: NM
INTEGER :: ICO = 0, NC, NOM, IOR0, IC, JC, ICE, ICN, IW, I, J, K

AB => G%POISSONB
LAST_CELL_IN_LAST_MESH_BANDWISE_IF: IF (NM == NMESHES) THEN

   NC = G%NC_LOCAL(NMESHES)

   ! Store column indices and values of diagonal and all off-diagonal entries in last row
   ! index '1' corresponds to main diagonal entry

   ICO = ICO + 1
   ABCO => AB%CONDENSED(ICO)

   ABCO%IOR0 = 0
   ABCO%ICO  = NC
   ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(NC, 1:AB%N_STENCIL)
   ABCO%VAL2(1:AB%N_STENCIL) = 0.0_EB
   ABCO%VAL2(AB%POS(0)) = 1.0_EB

   ! Within last mesh: check which other cells have a connection to the last cell;
   ! in each corresponding matrix row store the column index and value of just that matrix entry
   ! for each direction only one value has to be stored
 
   JC = NC - 1
   DO IOR0 = -3, 3
      IF (JC + AB%OFFSET(IOR0) == NC) THEN
         ICO = ICO + 1
         ABCO => AB%CONDENSED(ICO)
         ABCO%IOR0 = IOR0
         ABCO%ICO  = JC
         ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABCO%VAL2(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABCO%VAL2(AB%POS(ABCO%IOR0)) = 0.0_EB
         EXIT
      ENDIF
   ENDDO

   JC = NC - L%NX
   DO IOR0 = -3, 3
      IF (JC + AB%OFFSET(IOR0) == NC) THEN
         ICO = ICO + 1
         ABCO => AB%CONDENSED(ICO)
         ABCO%IOR0 = IOR0
         ABCO%ICO  = JC
         ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABCO%VAL2(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABCO%VAL2(AB%POS(ABCO%IOR0)) = 0.0_EB
         EXIT
      ENDIF
   ENDDO

   IF (.NOT.TWO_D) THEN
      JC = NC - L%NX * L%NY
      DO IOR0 = -3, 3
         IF (JC + AB%OFFSET(IOR0) == NC) THEN
            ICO = ICO + 1
            ABCO => AB%CONDENSED(ICO)
            ABCO%IOR0 = IOR0
            ABCO%ICO  = JC
            ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
            ABCO%VAL2(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
            ABCO%VAL2(AB%POS(ABCO%IOR0)) = 0.0_EB
            EXIT
         ENDIF
      ENDDO
   ENDIF

ENDIF LAST_CELL_IN_LAST_MESH_BANDWISE_IF
 
! Cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled
! with the last cell of last mesh;
! this can be a cell on the opposite side of the own mesh or a cell on a different mesh
! if such a cell exists, store corresponding matrix entry
 
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
      ICE = GWC%ICE                               ! adjacent ghost cell number
      ICN = G%ICE_TO_ICN(ICE)                     ! get column index of neighboring offdiagonal matrix entry
      IF (ICN /= SCARC(NMESHES)%NC) CYCLE         ! if no relation to last cell in last mesh, cycle
      ICO = ICO + 1
      ABCO => AB%CONDENSED(ICO)
      ABCO%IOR0 = IOR0
      ABCO%ICO  = IC
      ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(IC, 1:AB%N_STENCIL)
      ABCO%VAL2(1:AB%N_STENCIL) = AB%VAL(IC, 1:AB%N_STENCIL)
      ABCO%VAL2(AB%POS(ABCO%IOR0)) = 0.0_EB
      EXIT
   ENDIF 
ENDDO 

AB%N_CONDENSED = ICO

END SUBROUTINE SCARC_SETUP_BMATRIX_CONDENSED


! --------------------------------------------------------------------------------------------------------------
!> \brief Setup condensed system in case of periodic or pure Neumann boundary conditions
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM_CONDENSED (NV, NL, ITYPE)
USE SCARC_POINTERS, ONLY: L, G, OG, F, OL, VC, A, ACO, AB, ABCO, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV, NL, ITYPE
INTEGER :: NM, NOM, IFACE, ICN, ICE, ICW, JC, NC, ICO, IOR0, IP, ICG, INBR

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0 .OR. &
    TYPE_PRECON == NSCARC_RELAX_FFT .OR. TYPE_PRECON == NSCARC_RELAX_FFTO) RETURN

 
! In last mesh:  subtract B*RHS(end) for internal legs of stencil
 
MESH_REAL = 0.0_EB
IF (UPPER_MESH_INDEX == NMESHES) THEN

   CALL SCARC_POINT_TO_GRID (NMESHES, NL)

   NC =  G%NC_LOCAL(NMESHES)
   VC => SCARC_POINT_TO_VECTOR(NMESHES, NL, NV)

   ! Process last column entries of all rows except of last one
   ! for those rows only one matrix entry was stored, namely that one which connects to the last cell
 
   SELECT CASE (SET_MATRIX_TYPE(NL))

      CASE (NSCARC_MATRIX_COMPACT)
         A => G%POISSON
         DO ICO = 2, A%N_CONDENSED
            ACO => A%CONDENSED(ICO)
            JC = ACO%COL(1)
            IF (JC < NC) VC(JC) = VC(JC) - ACO%VAL1(1)*VC(NC)
         ENDDO

      CASE (NSCARC_MATRIX_BANDWISE)
         AB => G%POISSONB
         DO ICO = 2, AB%N_CONDENSED
            ABCO => AB%CONDENSED(ICO)
            IP = AB%POS(ABCO%IOR0)
            JC = ABCO%ICO
            IF (JC < NC) VC(JC) = VC(JC) - ABCO%VAL1(IP)*VC(NC)
        ENDDO

   END SELECT

   MESH_REAL(NMESHES) = VC(NC)     ! store last entry of RHS
   VC(NC) = 0.0_EB                 ! set last entry of last mesh to zero

ENDIF
IF (ITYPE == 0) RETURN

! Broadcast last RHS-value of last cell in last mesh to all meshes
 
IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MESH_REAL, 1, MPI_DOUBLE_PRECISION,&
                      MPI_COMM_WORLD, IERROR)

DO NM = 1, NMESHES
   SCARC(NM)%RHS_END = MESH_REAL(NMESHES)
ENDDO

! Only in case of periodic BC's:
! Subtract B*RHS(end) for corresponding entries of all periodic communication partners
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

   SNODE = PROCESS(NM)
   RNODE = PROCESS(NMESHES)

   IF (.NOT. ARE_FACE_NEIGHBORS(NM, NMESHES)) CYCLE

   CALL SCARC_POINT_TO_OTHER_GRID (NM, NMESHES, NL)
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

 
   ! Subtract B*RHS(end) at corresponding positions
 
   DO IFACE = 1, 6                                         ! check if this face has connection to last cell

      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)

      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
         IF (NOM /= NMESHES) CYCLE                         ! only check for common matrix entries with last mesh
         CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

            ICW = OG%ICG_TO_ICW(ICG, 1)
            ICE = OG%ICG_TO_ICE(ICG, 1)
            ICN = G%ICE_TO_ICN(ICE)                        ! get column index of neighboring offdiagonal matrix entry

            IF (ICN /= SCARC(NMESHES)%NC) CYCLE            ! if no relation to last cell in last mesh, cycle

            VC(ICW) = VC(ICW) - F%INCR_FACE * SCARC(NM)%RHS_END

         ENDDO

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_SYSTEM_CONDENSED


! ----------------------------------------------------------------------------------------------------------------------
!> \brief Extract overlapping matrix parts after data exchange with neighbors and add them to main matrix
! ----------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_MATRIX_OVERLAPS (NMATRIX, NTYPE, NL)
USE SCARC_POINTERS, ONLY: G, F, OL, OG, A, OA, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN) :: NL, NMATRIX, NTYPE
INTEGER :: NM, IFACE, NOM, IOR0, ICG, ICE, IP, ICOL, INBR, ICN, ICE1, ICE2

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                 
   A => SCARC_POINT_TO_CMATRIX (NMATRIX)

   IP = A%ROW(G%NC+1)
   FACES_LOOP: DO IFACE = 1, 6               

      IOR0 = FACE_ORIENTATION(IFACE)
      F => SCARC(NM)%LEVEL(NLEVEL_MIN)%FACE(IOR0)
   
      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
         OA => SCARC_POINT_TO_OTHER_CMATRIX (NMATRIX)

         ICOL = 1
         DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
  
            ICE = OG%ICG_TO_ICE(ICG, 1)
            A%ROW(ICE) = IP 

            IF (NTYPE == 1) THEN
               ICOL = OA%ROW(ICG)
               ICN = ABS(OA%COLG(ICOL))
               A%COL(IP)  = ICE
               A%COLG(IP) = ICN
               A%VAL(IP) = OA%VAL(ICOL)
               IP = IP + 1
               DO ICOL = OA%ROW(ICG)+1, OA%ROW(ICG+1)-1
                  ICN = OA%COLG(ICOL)
                  IF (SCARC_CELL_WITHIN_MESH(G, NM, ICN)) THEN
                     A%COL(IP) = ABS(OA%COLG(ICOL)) - G%NC_OFFSET(NM)     
                  ELSE
                     A%COL(IP) = -ABS(OA%COLG(ICOL))
                     IF (ICG == OL%GHOST_FIRSTE(IOR0)) THEN
                        ICE2 = OG%ICG_TO_ICE(ICG+1, 1)
                        IF (G%LOCAL_TO_GLOBAL(ICE2) == ICN) A%COL(IP) = ICE2
                     ELSE IF (ICG == OL%GHOST_LASTW(IOR0)) THEN
                        ICE1 = OG%ICG_TO_ICE(ICG-1, 1)
                        IF (G%LOCAL_TO_GLOBAL(ICE1) == ICN) A%COL(IP) = ICE1
                     ELSE
                        ICE1 = OG%ICG_TO_ICE(ICG-1, 1)
                        ICE2 = OG%ICG_TO_ICE(ICG+1, 1)
                        IF (G%LOCAL_TO_GLOBAL(ICE1) == ICN) A%COL(IP) = ICE1
                        IF (G%LOCAL_TO_GLOBAL(ICE2) == ICN) A%COL(IP) = ICE2
                     ENDIF
                  ENDIF
                  A%COLG(IP) = ABS(OA%COLG(ICOL))      
                  A%VAL(IP)  = OA%VAL(ICOL)
                  IP = IP + 1
               ENDDO
            ELSE
               DO ICOL = OA%ROW(ICG), OA%ROW(ICG+1)-1
                  A%COL(IP) = -OA%COL(ICOL)   
                  A%COLG(IP) = ABS(OA%COLG(ICOL))      
                  A%VAL(IP) = OA%VAL(ICOL)
                  IP = IP + 1
               ENDDO
            ENDIF
         ENDDO

         A%ROW(ICE+1) = IP 
         A%N_ROW = ICE + 1
         A%N_VAL = IP - 1

      ENDDO
   ENDDO FACES_LOOP


ENDDO MESHES_LOOP

END SUBROUTINE SCARC_EXTRACT_MATRIX_OVERLAPS


! --------------------------------------------------------------------------------------------------------------------
!> \brief Extract diagonal of Poisson matrix and store it in a separate vector for further use
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_MATRIX_DIAGONAL(NL)
USE SCARC_POINTERS, ONLY: G, A, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC, JC, ICOL

CROUTINE = 'SCARC_EXTRACT_MATRIX_DIAGONAL'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)

   CALL SCARC_ALLOCATE_REAL1 (G%DIAG, 1, G%NCE2, NSCARC_INIT_ZERO, 'G%DIAG', CROUTINE)
   DO IC = 1, G%NC
      DO ICOL = A%ROW(IC), A%ROW(IC+1) - 1
         JC = A%COL(ICOL)
         IF (JC == IC) G%DIAG(IC) = A%VAL(ICOL)
      ENDDO
   ENDDO

ENDDO MESHES_LOOP

! If there are multiple meshes exchange diagonal matrix on overlapping parts

IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_DIAGS, NSCARC_NONE, NL)

END SUBROUTINE SCARC_EXTRACT_MATRIX_DIAGONAL

END MODULE SCARC_MATRICES



!=======================================================================================================================
!
! MODULE SCARC_FFT
!
!> \brief Define environment needed for the use of the FFT preconditioner based on Crayfishpak
!
!=======================================================================================================================
MODULE SCARC_FFT
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE SCARC_CONSTANTS
USE SCARC_TYPES
USE SCARC_VARIABLES
USE SCARC_MESSAGES
USE SCARC_TROUBLESHOOTING
USE SCARC_VECTORS
USE SCARC_MATRICES
USE SCARC_CONVERGENCE

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup data structures for the use of blockwise FFT methods as preconditioners
! New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFT(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SCARC_POINT_TO_GRID
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO NL = NLMIN, NLMAX
      CALL SCARC_SETUP_FFT_MESH(NM, NL)
   ENDDO 
ENDDO 

END SUBROUTINE SCARC_SETUP_FFT


! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup data structures for the use of blockwise FFT methods as preconditioners
! New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFT_MESH(NM, NL)
USE SCARC_POINTERS, ONLY: M, S, L, FFT, SCARC_POINT_TO_GRID
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NM, NL

CROUTINE = 'SCARC_SETUP_FFT'
 
! Allocate working space for FFT routine
 
CALL SCARC_POINT_TO_GRID (NM, NL)                                    
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

CALL SCARC_ALLOCATE_REAL1 (FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT%SAVE1', CROUTINE)
CALL SCARC_ALLOCATE_REAL1 (FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT%WORK', CROUTINE)

! Allocate stretching vector (set to 1)
CALL SCARC_ALLOCATE_REAL1 (FFT%HX, 0, FFT%ITRN, NSCARC_INIT_ONE, 'FFT%HX', CROUTINE)

! Allocate RHS vector for FFT routine
IF (L%NY == 1) THEN
   CALL SCARC_ALLOCATE_REAL3 (FFT%PRHS, 1, FFT%ITRN, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
ELSE
   CALL SCARC_ALLOCATE_REAL3 (FFT%PRHS, 1, FFT%ITRN, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
ENDIF

! Allocate boundary data vector for XS
IF (L%NZ>1) THEN
   IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BXS, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
   IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BXS, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
ELSE
   CALL SCARC_ALLOCATE_REAL2 (FFT%BXS, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
ENDIF

! Allocate boundary data vector for XF
IF (L%NZ>1) THEN
   IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BXF, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
   IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BXF, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
ELSE
   CALL SCARC_ALLOCATE_REAL2 (FFT%BXF, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
ENDIF

! Allocate boundary data vector for YS
IF (L%NZ > 1) THEN
   CALL SCARC_ALLOCATE_REAL2 (FFT%BYS, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYS', CROUTINE)
ELSE
   CALL SCARC_ALLOCATE_REAL2 (FFT%BYS, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYS', CROUTINE)
ENDIF

! Allocate boundary data vector for YF
IF (L%NZ > 1) THEN
   CALL SCARC_ALLOCATE_REAL2 (FFT%BYF, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYF', CROUTINE)
ELSE
   CALL SCARC_ALLOCATE_REAL2 (FFT%BYF, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYF', CROUTINE)
ENDIF

! Allocate boundary data vector for ZS
IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BZS, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZS', CROUTINE)
IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BZS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZS', CROUTINE)

! Allocate boundary data vector for ZF
IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2 (FFT%BZF, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZF', CROUTINE)
IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BZF, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZF', CROUTINE)

IF (TWO_D) THEN
   CALL H2CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
               FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
               FFT%HX,FFT%XLM,FFT%ITRN,IERROR,FFT%SAVE1)
ELSE
   CALL H3CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
               FFT%YS,FFT%YF,FFT%JBAR,FFT%MBC,&
               FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
               FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERROR,FFT%SAVE1)
ENDIF

END SUBROUTINE SCARC_SETUP_FFT_MESH


! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup data structures for the use of blockwise FFT methods as preconditioners
! New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_FFT(NM, NL)
USE SCARC_POINTERS, ONLY: M, S, L, MGM, FFT, SCARC_POINT_TO_MGM
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NM, NL

CROUTINE = 'SCARC_SETUP_MGM_FFT'
 
CALL SCARC_POINT_TO_MGM (NM, NL)                                    
FFT => MGM%FFT

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

CALL SCARC_ALLOCATE_REAL1 (FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT%SAVE1', CROUTINE)
CALL SCARC_ALLOCATE_REAL1 (FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT%WORK', CROUTINE)

! Allocate stretching vector (set to 1)
CALL SCARC_ALLOCATE_REAL1 (FFT%HX, 0, FFT%ITRN, NSCARC_INIT_ONE, 'FFT%HX', CROUTINE)

! Allocate RHS vector for FFT routine
IF (L%NY == 1) THEN
   CALL SCARC_ALLOCATE_REAL3 (FFT%PRHS, 1, FFT%ITRN, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
ELSE
   CALL SCARC_ALLOCATE_REAL3 (FFT%PRHS, 1, FFT%ITRN, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
ENDIF

IF (TWO_D) THEN
   CALL H2CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
               FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
               FFT%HX,FFT%XLM,FFT%ITRN,IERROR,FFT%SAVE1)
ELSE
   CALL H3CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
               FFT%YS,FFT%YF,FFT%JBAR,FFT%MBC,&
               FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
               FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERROR,FFT%SAVE1)
ENDIF

END SUBROUTINE SCARC_SETUP_MGM_FFT


! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for blockwise FFT methods with overlap
! New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFTO(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: M, S, L, FFT, SCARC_POINT_TO_GRID
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL

CROUTINE = 'SCARC_SEtUP_FFTO'

! Allocate working space for FFT routine
 
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    
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

      CALL SCARC_ALLOCATE_REAL1 (FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT%SAVE1', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT%WORK', CROUTINE)

      ! Allocate stretching vector (set to 1)

      CALL SCARC_ALLOCATE_REAL1 (FFT%HX, 0, FFT%ITRN, NSCARC_INIT_ONE, 'FFT', CROUTINE)

      ! Allocate RHS vector for FFT routine

      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3 (FFT%PRHS, 1, FFT%ITRN, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL3 (FFT%PRHS, 1, FFT%ITRN, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for XS

      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BXS, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BXS, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2 (FFT%BXS, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for XF

      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BXF, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BXF, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2 (FFT%BXF, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for YS

      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2 (FFT%BYS, 1, FFT%ITRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYS', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2 (FFT%BYS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BYS', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for YF

      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2 (FFT%BYF, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYF', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2 (FFT%BYF, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYF', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for ZS

      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BZS, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZS', CROUTINE)
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BZS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZS', CROUTINE)

      ! Allocate boundary data vector for ZF

      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2 (FFT%BZF, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZF', CROUTINE)
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2 (FFT%BZF, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZF', CROUTINE)

      IF (TWO_D) THEN
         CALL H2CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,IERROR,FFT%SAVE1)
      ELSE
         CALL H3CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%YS,FFT%YF,FFT%JBAR,FFT%MBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERROR,FFT%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FFTO

END MODULE SCARC_FFT



!=======================================================================================================================
!
! MODULE SCARC_GMG
!
!> \brief Define environment needed for the use of the Geometric Multigrid Method
!
!=======================================================================================================================
MODULE SCARC_GMG
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE SCARC_CONSTANTS
USE SCARC_VARIABLES

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief Perform restriction from finer to coarser grid level
!    - 'VF' corresponds to vector on fine   grid
!    - 'VC' corresponds to vector on coarse grid
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTRICTION (NVB, NVC, NLF, NLC)
USE SCARC_POINTERS, ONLY: LC, GF, GC, VF, VC, RF, &
                          SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_CMATRIX
USE SCARC_MPI, ONLY: SCARC_EXCHANGE
INTEGER, INTENT(IN) :: NVB, NVC, NLF, NLC
REAL(EB) :: DSUM
INTEGER :: NM
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-2:2,-2:2)=0
INTEGER :: IXC, IYC, IZC, ICC, IC, ICOL

IF (IS_GMG .OR. IS_CG_GMG .OR. IS_CG_ADD .OR. IS_CG_MUL1) THEN
 
! ---------- Twolevel-CG or Geometric multigrid (as main solver or preconditioner) 
 
   IF (HAS_MULTIPLE_LEVELS) THEN
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_MULTIGRID (NM, NLF, NLC)
   
         VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)
         VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
   
         IF (TWO_D) THEN
   
            SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)
   
               ! ---------- Constant Interpolation
   
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
      
               ! ---------- Bilinear Interpolation
   
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
   
         ! ---------- Constant Interpolation (Note: 3D-bilinear case is still missing)
   
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
   
! ---------- Use restriction based on smoothed aggregation method

ELSE

   IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NVB, NLF)

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      CALL SCARC_POINT_TO_MULTIGRID (NM, NLF, NLC)

      RF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_RESTRICTION)

      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)
      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
   


      DO IC = 1, GF%N_COARSE
         DSUM = 0.0_EB
         DO ICOL = RF%ROW(IC), RF%ROW(IC+1)-1                            
            DSUM =  DSUM + RF%VAL(ICOL) * VF(RF%COLG(ICOL))
         ENDDO
         VC(IC) = DSUM
      ENDDO

   ENDDO

ENDIF

END SUBROUTINE SCARC_RESTRICTION


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform prolongation from coarser to finer grid level
!    - 'VC' corresponds to coarser grid
!    - 'VF' corresponds to finer   grid
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PROLONGATION (NVC, NVB, NLC, NLF)
USE SCARC_POINTERS, ONLY: LC, GF, GC, VF, VC, PF, &
                          SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_CMATRIX
USE SCARC_MPI, ONLY: SCARC_EXCHANGE
INTEGER, INTENT(IN) :: NVC, NVB, NLC, NLF
REAL(EB) :: DSUM
INTEGER :: NM, I
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-1:1,-1:1)=0
INTEGER :: IXC, IYC, IZC, ICC, IC, ICOL

IF (IS_GMG .OR. IS_CG_GMG .OR. IS_CG_ADD .OR. IS_CG_MUL1) THEN
 
! ------------------ Twolevel CG or Geometric Multigrid 
 
   IF (HAS_MULTIPLE_LEVELS) THEN
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_MULTIGRID (NM, NLF, NLC)
   
         VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
         VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)
   
         IF (TWO_D) THEN
   
            SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)
   
               CASE (NSCARC_INTERPOL_CONSTANT)
      
                  !!$OMP PARALLEL DO PRIVATE(IXC, IZC, IXF, IZF, ICC, ICF) SCHEDULE(STATIC)
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
                  !!$OMP END PARALLEL DO 
      
               CASE (NSCARC_INTERPOL_BILINEAR)
      
                  !!$OMP PARALLEL DO PRIVATE(IXC, IZC, IXF, IZF, ICC, ICFB) SCHEDULE(STATIC)
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
                  !!$OMP END PARALLEL DO 
      
            END SELECT SELECT_INTERPOL
   
         ELSE
   
            ! ---------- Constant Interpolation (Note: 3D-bilinear case is still missing)
   
            !!$OMP PARALLEL DO PRIVATE(IXC, IYC, IZC, IXF, IYF, IZF, ICC, ICF) SCHEDULE(STATIC)
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
            !!$OMP END PARALLEL DO 
   
         ENDIF
      ENDDO
   ENDIF
   
! ---------- Use Prolongation matrix based on smoothed aggregation method
ELSE

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      CALL SCARC_POINT_TO_MULTIGRID (NM, NLF, NLC)

      PF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)

      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)

      DO IC = 1, GF%NC
         DSUM = 0.0_EB
         DO ICOL = PF%ROW(IC), PF%ROW(IC+1)-1                            
            DSUM = DSUM + PF%VAL(ICOL) * VC(PF%COL(ICOL))
         ENDDO
         VF(IC) = DSUM
      ENDDO

   ENDDO

   IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NVB, NLF)

ENDIF

END SUBROUTINE SCARC_PROLONGATION

END MODULE SCARC_GMG


#ifdef WITH_SCARC_AMG
!=======================================================================================================================
!
! MODULE SCARC_AMG
!
!> \brief Setup algebraic multigrid structures
!   Allocate needed workspace for hierarchy of system matrices, prolongation, restriction, etc.
!   Note: all used pointers end with either 'F' or 'C' where:
!       'F' corresponds to fine   level NL
!       'C' corresponds to coarse level NL+1
!   Determine mesh hierarchy based on smoothed aggregation
!   Compute QR-decomposition of nullspace vector in order to determine tentative prolongator 
!   Set nullspace for next level and perform Jacobi relaxation to get the final prolongator
!   If the maximum allowed level is not yet reached, set dimensions for next coarser level, 
!   define its nullspace and perform relaxation to define the respective Prolongation matrix
!   Define Poisson matrix on coarser level by Galerkin approach: A_coarse = R*A_fine*P
!
!=======================================================================================================================
MODULE SCARC_AMG
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE SCARC_CONSTANTS
USE SCARC_TYPES
USE SCARC_VARIABLES
USE SCARC_STORAGE
USE SCARC_UTILITIES
USE SCARC_MESSAGES
USE SCARC_CPU, ONLY: CPU
USE SCARC_MPI, ONLY: SCARC_EXCHANGE
USE SCARC_MATRICES

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! ------------------------------------------------------------------------------------------------------------------
!> \brief Setup structures needed for the use of the algebraic multigrid method
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AMG_ENVIRONMENT
INTEGER:: NL
LOGICAL:: FURTHER_COARSENING_REQUIRED

IF (.NOT.HAS_AMG_LEVELS) RETURN
FURTHER_COARSENING_REQUIRED = .TRUE.

NL = NLEVEL_MIN
!CALL  SCARC_PYTHON_MATRIX(NL, 'A')

COARSENING_LOOP: DO WHILE (FURTHER_COARSENING_REQUIRED)


   ! Determine the aggregation order among the meshes
   CALL SCARC_SETUP_AGGREGATION_ORDER                    

   ! Extract matrix diagonal from Poisson matrix A, determine strength of connection matrix and store inverted matrix diagonal 
   CALL SCARC_EXTRACT_MATRIX_DIAGONAL(NL)           
   CALL SCARC_SETUP_CONNECTION(NL)                  
   CALL SCARC_INVERT_MATRIX_DIAGONAL(NL)            

   ! Apply smoothed aggregation heuristic to specify aggregation zones and improve near null space by Jacobi relaxation step
   CALL SCARC_SETUP_AGGREGATION_ZONES(NL)           
   CALL SCARC_RELAX_NULLSPACE(NL)                   

   ! Setup final aggregation Zones matrix Z and Prolongation matrix P based on QR-decomposition
   CALL SCARC_SETUP_ZONE_OPERATOR(NL)               

      ! Setup restriction and prolongation matrices for GMG-like coarsening
   IF (TYPE_COARSENING == NSCARC_COARSENING_GMG) THEN
      CALL SCARC_SETUP_TRANSFER_GMG(NL)

   ! Setup Prolongation matrix P based on QR-decomposition, near nullspace on coarser level and corresponding Restriction matrix R
   ELSE
      CALL SCARC_SETUP_PROLONGATION_AMG(NL)
      CALL SCARC_SETUP_NULLSPACE_COARSE(NL)
      CALL SCARC_SETUP_RESTRICTION(NL)
   ENDIF

   ! First setup A*P matrix to finally build the Galerkin matrix R*A*P
   CALL SCARC_SETUP_POISSON_PROL(NL)             
   CALL SCARC_SETUP_GALERKIN(NL)             

   ! Remove workspace which is no longer used and get the next coarsening round on the wa
   CALL SCARC_CLEAN_WORKSPACE_AMG(NL)
      
   NL = NL+1
   IF (NL == NLEVEL_MAX) FURTHER_COARSENING_REQUIRED = .FALSE.

ENDDO COARSENING_LOOP

END SUBROUTINE SCARC_SETUP_AMG_ENVIRONMENT


! --------------------------------------------------------------------------------------------------------------------
!> \brief  Setup order in which aggregation is performed over mesh decomposition
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AGGREGATION_ORDER
USE SCARC_POINTERS, ONLY: SUB, S
LOGICAL, ALLOCATABLE, DIMENSION(:):: NOT_AGGREGATED
INTEGER:: NM, NOM, INBR, ICYCLE

CROUTINE = 'SCARC_SETUP_AGGREGATION_ORDER'
SUB => SUBDIVISION

CALL SCARC_ALLOCATE_INT2 (SUB%ORDER, 1, NMESHES, 1, NMESHES, NSCARC_INIT_NONE, 'SUB%ORDER', CROUTINE)
CALL SCARC_ALLOCATE_LOG1 (NOT_AGGREGATED, 1, NMESHES, NSCARC_INIT_TRUE, 'NOT_AGGREGATED', CROUTINE)

ICYCLE = 1
DO WHILE (ANY(NOT_AGGREGATED)) 
   SUB%ORDER(1:NMESHES, ICYCLE) = NSCARC_ORDER_UNASSIGNED
   DO NM = 1, NMESHES
      S => SCARC(NM)
      IF (NOT_AGGREGATED(NM) .AND. SUB%ORDER(NM, ICYCLE) /= NSCARC_ORDER_LOCKED) THEN
         SUB%ORDER(NM, ICYCLE) = NSCARC_ORDER_ACTIVE
         DO INBR = 1, SUB%N_NEIGHBORS(NM)
            NOM = SUB%NEIGHBORS(INBR, NM)
            SUB%ORDER(NOM, ICYCLE) = NSCARC_ORDER_LOCKED
         ENDDO
         NOT_AGGREGATED(NM) = .FALSE.
      ENDIF
   ENDDO
   ICYCLE = ICYCLE+1
ENDDO

SUB%N_CYCLES = ICYCLE-1

CALL SCARC_DEALLOCATE_LOG1 (NOT_AGGREGATED, 'NOT_AGGREGATED', CROUTINE)

END SUBROUTINE SCARC_SETUP_AGGREGATION_ORDER


! --------------------------------------------------------------------------------------------------------------------
!> \brief  Invert matrix diagonal which is already stored in DIAG-vector (reuse workspace)
! Scale each matrix element with inverse of diagonal and approximate spectral radius (currently disabled) 
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INVERT_MATRIX_DIAGONAL(NL)
USE SCARC_POINTERS, ONLY: G, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN):: NL
INTEGER:: NM, IC

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   DO IC = 1, G%NCE
      G%DIAG(IC) = 1.0_EB/G%DIAG(IC)
   ENDDO
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_INVERT_MATRIX_DIAGONAL


! --------------------------------------------------------------------------------------------------------------------
!> \brief  Compute a strength of connection matrix based on symmetric smoothed aggregation heuristic. 
! A nonzero connection A[i, j] is considered strong if:

!     abs(A[i, j]) >= theta*sqrt( abs(A[i, i]) * abs(A[j, j]) )

! The strength matrix S corresponds to the set of nonzero entries of A that are strong connections
! based on a strength of connection tolerance theta
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CONNECTION(NL)
USE SCARC_POINTERS, ONLY: G, A, S, C, OA, OC, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN):: NL
REAL(EB):: VAL, EPS, SCAL, CVAL_MAX, THETA
INTEGER:: NM, NOM, IC, JC, ICOL, IZONE, INBR

IF (TYPE_COARSENING == NSCARC_COARSENING_CUBIC) RETURN

CROUTINE = 'SCARC_SETUP_CONNECTION'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

   A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
   C => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_CONNECTION)


   ! Allocate workspace for strength of connection matrix (use same size as Poisson matrix)
   C%N_VAL = A%N_VAL                         
   C%N_ROW = A%N_ROW
   CALL SCARC_ALLOCATE_CMATRIX (C, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_LIGHT, 'G%CONNECTION', CROUTINE)

   IF (NL == NLEVEL_MIN) THEN
      THETA = 0.10E+0_EB
   ELSE
      THETA = SCARC_MULTIGRID_THETA
   ENDIF
   THETA = SCARC_MULTIGRID_THETA
   
   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

      OA => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_POISSON)
      OC => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_CONNECTION)

      OC%N_VAL = 2*OA%N_VAL                   ! use double layers
      OC%N_ROW = OA%N_ROW           
      CALL SCARC_ALLOCATE_CMATRIX (OC, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_LIGHT, 'OG%CONNECTION', CROUTINE)
      
   ENDDO 

   ! Check strength-of-connection criterion  
   IZONE = 1
   C%ROW(1) = 1
   DO IC = 1, G%NC
   
      EPS = THETA**2 * ABS(G%DIAG(IC))                  ! EPS = theta**2 * A_ii
      DO ICOL = A%ROW(IC), A%ROW(IC+1) - 1
         JC  = A%COL(ICOL)
         IF (JC == 0) CYCLE                                             ! omit second layer
         VAL = A%VAL(ICOL)                                              ! VAL = A_ij
   
         ! Always add the diagonal: |A_ii|  >= THETA*sqrt(|A_ii| * |A_ii|)     true!
         IF (IC == JC) THEN
            C%COL(IZONE) = JC
            C%VAL(IZONE) = VAL
            IZONE = IZONE+1

         ! Check subdiagonal entry: |A_ij|  >= THETA*sqrt(|A_ii| * |A_jj|)     ??
         ELSE IF (VAL**2 >= EPS*ABS(G%DIAG(JC))) THEN
            C%COL(IZONE) = JC
            C%VAL(IZONE) = VAL
            IZONE = IZONE+1
         ENDIF

      ENDDO
      C%ROW(IC+1) = IZONE
   ENDDO
   
   
   DO IC = 1, G%NC
      CVAL_MAX = 0.0_EB
      DO ICOL = C%ROW(IC), C%ROW(IC+1) - 1
         CVAL_MAX = MAX(ABS(C%VAL(ICOL)), CVAL_MAX)
      ENDDO
      SCAL = 1.0_EB/CVAL_MAX
      DO ICOL = C%ROW(IC), C%ROW(IC+1) - 1
         C%VAL(ICOL) = ABS(C%VAL(ICOL))*SCAL
      ENDDO
   ENDDO
   
   
ENDDO MESHES_LOOP
   
! If there are multiple meshes, exchange strength matrix on overlapping parts
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_CONNECTION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_CONNECTION, NL)
ENDIF

END SUBROUTINE SCARC_SETUP_CONNECTION
 

! --------------------------------------------------------------------------------------------------------------------
!> \brief Setup aggregation zones for Smoothed Aggregation Algebraic Multigrid Method
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AGGREGATION_ZONES(NL)
USE SCARC_POINTERS, ONLY: SUB, C, CF, G, LF, LC, GF, GC, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_MULTIGRID
INTEGER, INTENT(IN):: NL
INTEGER:: NM, NM2, ICYCLE, IC, IZL

CROUTINE = 'SCARC_SETUP_AGGREGATION_ZONES'

SUB => SUBDIVISION
MESH_INT = -1

! Allocate workspaces for coarse points, global and local aggregation zones 

MESHES_ALLOCATION_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

   CALL SCARC_ALLOCATE_INT1 (G%ZONE_CENTERS,  1, G%NCE,  NSCARC_INIT_ZERO, 'G%ZONE_CENTERS', CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (G%ZONES_GLOBAL, 1, G%NCE2, NSCARC_INIT_ZERO, 'G%ZONES_GLOBAL', CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (G%ZONES_LOCAL,  1, G%NCE2, NSCARC_INIT_ZERO, 'G%ZONES_LOCAL', CROUTINE)

ENDDO MESHES_ALLOCATION_LOOP

COARSENING_TYPE_SELECT: SELECT CASE (TYPE_COARSENING)
 
   ! ---- Default aggregation procedure for SAMG
 
   CASE (NSCARC_COARSENING_AGGREGATED)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_GRID (NM, NL)
         C => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_CONNECTION)
         CALL SCARC_SETUP_COARSENING_AGGREGATION(G, C)
         MESH_INT (NM) = G%N_ZONES
      ENDDO
      
      ! Exchange overlapping information of active meshes

      IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_NONE, NL)
   
   ! ---- Staggered aggregation procedure for SAMG
 
   CASE (NSCARC_COARSENING_AGGREGATED_S)

      CYCLES_LOOP1: DO ICYCLE = 1, SUB%N_CYCLES
   
         ! First aggregate on active meshes

         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
            IF (SUB%ORDER(NM, ICYCLE) == NSCARC_ORDER_ACTIVE) THEN
               CALL SCARC_POINT_TO_GRID (NM, NL)
               C => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_CONNECTION)
               CALL SCARC_SETUP_COARSENING_AGGREGATION(G, C)
               MESH_INT (NM) = G%N_ZONES
            ENDIF
   
         ENDDO
      
         ! Exchange overlapping information of active meshes

         IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_NONE, NL)
   
         ! Then aggregate on passive meshes (taking into account overlapping aggregate information of active meshes)

         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      
            IF (SUB%ORDER (NM, ICYCLE) /= NSCARC_ORDER_ACTIVE) THEN
               CALL SCARC_POINT_TO_GRID (NM, NL)
               C => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_CONNECTION)
               CALL SCARC_SETUP_COARSENING_AGGREGATION(G, C)
               MESH_INT (NM) = G%N_ZONES
            ENDIF
   
         ENDDO
   
         ! Exchange overlapping information of passive meshes

         IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_NONE, NL)
   
      ENDDO CYCLES_LOOP1
   
   ! ---- GMG-like aggregation procedure 
   !      In case of even cell numbers this process corresponds to the usual GMG coarsening
   !      in case of uneven cell number in a coordinate direction, on patch with 3 cells is used, the rest with patches of 2
 
   CASE (NSCARC_COARSENING_CUBIC)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)
         CALL SCARC_SETUP_COARSENING_CUBIC(LF, LC, GF, GC)
         MESH_INT(NM) = GF%N_ZONES
      ENDDO
      
      ! Exchange overlapping information 

      IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_NONE, NL)
   
   CASE (NSCARC_COARSENING_GMG)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)
         CALL SCARC_SETUP_COARSENING_GMG(LF, LC, GF, GC)
         MESH_INT(NM) = GF%N_ZONES
      ENDDO
      
END SELECT COARSENING_TYPE_SELECT

! Broadcast number of zones of all meshes

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE, 1, MPI_INTEGER, MESH_INT, COUNTS, DISPLS, MPI_INTEGER, MPI_COMM_WORLD, IERROR)
      
! Prepare grid dimensions of coarse grid level
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)         

   ! Setup grid dimensions on coarse level 
   GC%NC_LOCAL(1:NMESHES) = MESH_INT(1:NMESHES)
   GC%NC_GLOBAL = SUM(MESH_INT(1:NMESHES))
   GC%NC  = GC%NC_LOCAL(NM)
   GC%NCE = GC%NC_LOCAL(NM)
   IF (NMESHES > 1) THEN
      DO NM2 = 2, NMESHES
         GC%NC_OFFSET(NM2) = GC%NC_OFFSET(NM2-1) + GC%NC_LOCAL(NM2-1)
      ENDDO
   ENDIF                   

   ! Setup mapping from local zones to global zones

   CALL SCARC_ALLOCATE_INT1 (GC%LOCAL_TO_GLOBAL, 1, GF%NCE, NSCARC_INIT_ZERO, 'G%LOCAL_TO_GLOBAL', CROUTINE)
   DO IZL = 1, GC%NC
      GC%LOCAL_TO_GLOBAL(IZL) = IZL+GC%NC_OFFSET(NM)
   ENDDO

   DO IC = 1, GF%NC
      GF%ZONES_GLOBAL(IC) = GF%ZONES_LOCAL(IC) + GC%NC_OFFSET(NM)
   ENDDO

ENDDO

! Exchange zones information between meshes

IF (NMESHES > 1)  THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_NEIGHBORS, NSCARC_NONE, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_NEIGHBORS, NSCARC_NONE, NL)
   CALL SCARC_EXTRACT_ZONE_OVERLAPS(NL)
   CALL SCARC_EXTRACT_ZONE_POINTERS(NL)
ENDIF

! Determine final grid dimensions on coarser level and reduce zone arrays to correct length

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)         

   GC%NCE  = GF%N_ZONES
   GC%NCE2 = GF%N_ZONES

   CALL SCARC_REDUCE_INT1 (GC%LOCAL_TO_GLOBAL, GC%NCE2, 'GC%LOCAL_TO_GLOBAL', CROUTINE)
   CALL SCARC_REDUCE_INT1 (GC%ZONE_CENTERS, GC%NCE2, 'GC%ZONE_CENTERS', CROUTINE)

   GF%N_COARSE = GF%N_ZONES

   CALL SCARC_REDUCE_INT1 (GF%ZONES_LOCAL,  GF%NCE2, 'GC%LOCAL_TO_GLOBAL', CROUTINE)
   CALL SCARC_REDUCE_INT1 (GF%ZONES_GLOBAL, GF%NCE2, 'GC%LOCAL_TO_GLOBAL', CROUTINE)

   CF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_CONNECTION)
   CALL SCARC_DEALLOCATE_CMATRIX (CF, 'STRENGTH OF CONNECTION', CROUTINE)


ENDDO

END SUBROUTINE SCARC_SETUP_AGGREGATION_ZONES



! -------------------------------------------------------------------------------------------------------------
!> \brief Extract overlapping zone information (including second layers)
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_ZONE_OVERLAPS(NL)
USE SCARC_POINTERS, ONLY: GC, GF, OLF, OGF, F, &
                          SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_OTHER_MULTIGRID
INTEGER, INTENT(IN):: NL
INTEGER:: NM, INBR, IOR0, NOM, IZ, ICG, ICE, ICE2, IFOUND, IZL_CURRENT
INTEGER:: IZL1, IZL2, IZG1, IZG2, IFACE

CROUTINE = 'SCARC_EXTRACT_ZONE_OVERLAPS'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)    

   ! Clear overlapping parts and fill them with the recently exchanged data
   GF%ZONES_LOCAL (GF%NC+1: GF%NCE2) = 0
   GF%ZONES_GLOBAL(GF%NC+1: GF%NCE2) = 0

   IZL_CURRENT = GF%N_ZONES+1

   DO IFACE = 1, 6                                        

      IOR0 = FACE_ORIENTATION(IFACE)
      F => SCARC(NM)%LEVEL(NLEVEL_MIN)%FACE(IOR0)

      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_MULTIGRID (NM, NOM, NL, NL+1)

         IZ = 0
         DO ICG = OLF%GHOST_FIRSTE(IOR0), OLF%GHOST_LASTE(IOR0)

            ICE  = OGF%ICG_TO_ICE(ICG, 1)


            IZG1 = OGF%ICG_TO_GZONE(ICG)
            IFOUND = FINDLOC (GF%ZONES_GLOBAL, VALUE = IZG1, DIM = 1)
            IF (IFOUND == 0) THEN
               IZ = IZ+1
               GF%N_ZONES = GF%N_ZONES+1
               IZL1 = GF%N_ZONES
            ELSE
               IZL1 = GF%ZONES_LOCAL(IFOUND)
            ENDIF
            GF%ZONES_LOCAL(ICE)  = IZL1
            GF%ZONES_GLOBAL(ICE) = IZG1
            GC%LOCAL_TO_GLOBAL(IZL1) = IZG1

 
            IF (NL /= NLEVEL_MIN) CYCLE

            ICE2 = OGF%ICG_TO_ICE(ICG, 2)
            IZG2 = OGF%ICG_TO_GZONE(ICG+OGF%NCG)
            IFOUND = FINDLOC (GF%ZONES_GLOBAL, VALUE = IZG2, DIM = 1)
            IF (IFOUND == 0) THEN
               GF%N_ZONES = GF%N_ZONES+1
               IZL2 = GF%N_ZONES
            ELSE
               IZL2 = GF%ZONES_LOCAL(IFOUND)
            ENDIF
            GF%ZONES_LOCAL(ICE2)  = IZL2
            GF%ZONES_GLOBAL(ICE2) = IZG2
            GC%LOCAL_TO_GLOBAL(IZL2) = IZG2


         ENDDO
      ENDDO
   ENDDO 


ENDDO MESHES_LOOP

END SUBROUTINE SCARC_EXTRACT_ZONE_OVERLAPS


! -------------------------------------------------------------------------------------------------------------
!> \brief Setup pointers for overlapping zones for a pair of grid levels
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_ZONE_POINTERS(NL)
USE SCARC_POINTERS, ONLY : F, GF, GC, OLF, OLC, OGF, OGC, &
                           SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_OTHER_MULTIGRID
INTEGER, INTENT(IN):: NL
INTEGER:: INBR, IOR0, IZ, ICW1, ICW2, ICE1, ICE2, IZL1, IZL2, ICG, IZW, IZE, IFACE
INTEGER:: NM, NOM, NCGE_TOTAL = 0

CROUTINE = 'SCARC_EXTRACT_ZONE_POINTERS'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)    
   
   DO IFACE = 1, 6                                        

      IOR0 = FACE_ORIENTATION(IFACE)
      F => SCARC(NM)%LEVEL(NLEVEL_MIN)%FACE(IOR0)

      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_MULTIGRID (NM, NOM, NL, NL+1)

         CALL SCARC_ALLOCATE_INT1 (OGF%ICG_TO_IZONE, 1, 2*OGF%NCG, NSCARC_INIT_ZERO, 'OGF%ICG_TO_IZONE', CROUTINE)
         CALL SCARC_ALLOCATE_INT1 (OGF%ICG_TO_EZONE, 1, 2*OGF%NCG, NSCARC_INIT_ZERO, 'OGF%ICG_TO_EZONE', CROUTINE)

         IZ  = 0
         INTERNAL_ZONES_LOOP: DO ICG = OLF%GHOST_FIRSTW(IOR0), OLF%GHOST_LASTW(IOR0)

            ICW1 = OGF%ICG_TO_ICW(ICG, 1)
            IZL1 = GF%ZONES_LOCAL(ICW1)
            IF (FINDLOC(OGF%ICG_TO_IZONE, VALUE = IZL1, DIM = 1) == 0) THEN
               IZ = IZ+1
               OGF%ICG_TO_IZONE(IZ) = IZL1
            ENDIF

            IF (NL /= NLEVEL_MIN) CYCLE

            ICW2 = OGF%ICG_TO_ICW(ICG, 2)
            IZL2 = GF%ZONES_LOCAL(ICW2)
            IF (FINDLOC(OGF%ICG_TO_IZONE, VALUE = IZL2, DIM = 1) == 0) THEN
               IZ = IZ+1
               OGF%ICG_TO_IZONE(IZ) = IZL2
            ENDIF

         ENDDO INTERNAL_ZONES_LOOP
         OGF%NCGI = IZ
         CALL SCARC_REDUCE_INT1 (OGF%ICG_TO_IZONE, OGF%NCGI, 'OGF%ICG_TO_IZONE', CROUTINE)

         !First allocate in fine cell related length

         CALL SCARC_ALLOCATE_INT2 (OGC%ICG_TO_ICW, 1, OGF%NCGI, 1, 1, NSCARC_INIT_ZERO, 'OGF%ICG_TO_ICW', CROUTINE)
         
         IZ = 0
         DO ICG = 1, OGF%NCGI
            IZW = OGF%ICG_TO_IZONE(ICG)
            IF (FINDLOC(OGC%ICG_TO_ICW(1:OGF%NCGI, 1), VALUE = IZW, DIM = 1) == 0) THEN
               IZ = IZ+1
               OGC%ICG_TO_ICW(IZ, 1) = IZW
            ENDIF
         ENDDO
         OGC%NCG  = IZ
         OGC%NCGI = IZ

         OLC%GHOST_FIRSTW(IOR0) = 1
         OLC%GHOST_LASTW(IOR0)  = IZ

         ! Then reduce to zone related length 

         CALL SCARC_REDUCE_INT2 (OGC%ICG_TO_ICW, OGC%NCGI, 1, 'OGC%ICG_TO_ICW', CROUTINE)


         IZ  = 0 
         EXTERNAL_ZONES_LOOP: DO ICG = OLF%GHOST_FIRSTE(IOR0), OLF%GHOST_LASTE(IOR0)

            ICE1 = OGF%ICG_TO_ICE(ICG, 1)
            IZL1 = GF%ZONES_LOCAL(ICE1)
            IF (FINDLOC(OGF%ICG_TO_EZONE, VALUE = IZL1, DIM = 1) == 0) THEN
               IZ = IZ+1
               OGF%ICG_TO_EZONE(IZ) = IZL1
            ENDIF

            IF (NL /= NLEVEL_MIN) CYCLE

            ICE2 = OGF%ICG_TO_ICE(ICG, 2)
            IZL2 = GF%ZONES_LOCAL(ICE2)
            IF (FINDLOC(OGF%ICG_TO_EZONE, VALUE = IZL2, DIM = 1) == 0) THEN
               IZ = IZ+1
               OGF%ICG_TO_EZONE(IZ) = IZL2
            ENDIF

         ENDDO EXTERNAL_ZONES_LOOP
         OGF%NCGE = IZ
         CALL SCARC_REDUCE_INT1 (OGF%ICG_TO_EZONE, OGF%NCGE, 'OGF%ICG_TO_EZONE', CROUTINE)

         ! First allocate in fine cell related length

         CALL SCARC_ALLOCATE_INT2 (OGC%ICG_TO_ICE, 1, 2*OGF%NCGE, 1, 1, NSCARC_INIT_ZERO, 'OGC%ICG_TO_ICE', CROUTINE)
         
         IZ = 0
         DO ICG = 1, OGF%NCGE
            IZE = OGF%ICG_TO_EZONE(ICG)
            IF (FINDLOC(OGC%ICG_TO_ICE(1:OGF%NCGE, 1), VALUE = IZE, DIM = 1) == 0) THEN
               IZ = IZ+1
               OGC%ICG_TO_ICE(IZ, 1) = IZE
            ENDIF
         ENDDO
         OGC%NCGE = IZ
         NCGE_TOTAL = NCGE_TOTAL+OGC%NCGE

         OLC%GHOST_FIRSTE(IOR0) = 1
         OLC%GHOST_LASTE(IOR0)  = IZ

         ! Then reduce to zone related length 

         CALL SCARC_REDUCE_INT2 (OGC%ICG_TO_ICE, OGC%NCGE, 1, 'OGC%ICG_TO_ICE', CROUTINE)

         GC%N_STENCIL_MAX = 25                  ! TODO: ONLY TEMPORARILY
         OGC%NLEN_BUFFER_LAYER1  = MAX(OGC%NCGI, OGC%NCGE)
         OGC%NLEN_BUFFER_LAYER2  = OGC%NLEN_BUFFER_LAYER1*2
         OGC%NLEN_BUFFER_LAYER4  = OGC%NLEN_BUFFER_LAYER1*4
         OGC%NLEN_BUFFER_STENCIL = OGC%NLEN_BUFFER_LAYER1*GC%N_STENCIL_MAX
         OGC%NLEN_BUFFER_FULL    = OGC%NLEN_BUFFER_LAYER1*GC%N_STENCIL_MAX*2



      ENDDO 
   ENDDO 
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_EXTRACT_ZONE_POINTERS

! -------------------------------------------------------------------------------------------------------------
!> \brief Identify cells on second layer
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_IDENTIFY_LAYER2(NL)
USE SCARC_POINTERS, ONLY : S, A, OL, OG, &
                           SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN):: NL
INTEGER:: NM, NOM, ICW, ICOL, JC, JCG, INBR, IOR0, ICG, IS

CROUTINE = 'SCARC_IDENTIFY_LAYER2'
 
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)    
   A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)

   NEIGHBORS_LOOP: DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
    
      CALL SCARC_ALLOCATE_INT1 (OG%ICG_TO_ELAYER2, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'OG%ICG_TO_ELAYER2', CROUTINE)

      IS = 1
      FACE_LOOP: DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         GHOST_CELL_LOOP: DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            DO ICOL = A%ROW(ICW), A%ROW(ICW+1) - 1
               JC = A%COL(ICOL)
               IF (JC == 0) THEN
                  JCG = A%COLG(ICOL)
                  IF (FINDLOC (OG%ICG_TO_ELAYER2(1:2*OG%NCG), VALUE = JCG, DIM = 1) == 0) THEN 
                  OG%ICG_TO_ELAYER2(IS) = JCG
                  IS = IS+1
                  ENDIF
               ENDIF
            ENDDO
         ENDDO GHOST_CELL_LOOP
      ENDDO FACE_LOOP
      OL%N_LAYER2 = IS-1

   ENDDO NEIGHBORS_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_IDENTIFY_LAYER2



! --------------------------------------------------------------------------------------------------------------------
!> \brief  Standard aggregation prodecure based on strength of connection matrix
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSENING_AGGREGATION(G, C)
TYPE (SCARC_CMATRIX_TYPE), POINTER, INTENT(IN):: C
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN):: G
INTEGER:: IC, ICOL, JC, IZONE, JZONE
LOGICAL:: HAS_NEIGHBORS, HAS_AGGREGATED_NEIGHBORS

CROUTINE = 'SCARC_SETUP_COARSENING_AGGREGATION'

! 
! Pass 1 of aggregation:  Setup aggregation zones on internal cells of active mesh
! 
G%ZONES_LOCAL = 0
G%ZONES_GLOBAL = 0
G%ZONE_CENTERS = 0

PASS1_LOOP: DO IC = 1, G%NC

   IF (G%ZONES_LOCAL(IC) /= 0) CYCLE                           ! has cell already been aggregated?

   HAS_NEIGHBORS = .FALSE.
   HAS_AGGREGATED_NEIGHBORS = .FALSE.

   DO ICOL = C%ROW(IC), C%ROW(IC+1)-1                          ! are all neighbors free (not already aggregated)?
      JC = C%COL(ICOL)
      IF (JC /= 0 .AND. IC /= JC .AND. JC <= G%NC) THEN        ! only consider internal cells here
         HAS_NEIGHBORS = .TRUE.
         IF (G%ZONES_LOCAL(JC) /= 0) THEN
            HAS_AGGREGATED_NEIGHBORS = .TRUE.
            EXIT
         ENDIF
      ENDIF
   ENDDO

   IF (.NOT. HAS_NEIGHBORS) THEN                               ! do not aggregate isolated cells
      G%ZONES_LOCAL(IC) = NSCARC_HUGE_INT
   ELSE IF (.NOT. HAS_AGGREGATED_NEIGHBORS) THEN               ! build aggregate of this cell and its neighbors
      G%N_ZONES = G%N_ZONES+1
      G%ZONES_LOCAL(IC) = G%N_ZONES
      G%ZONE_CENTERS(G%N_ZONES) = IC                
      DO ICOL = C%ROW(IC), C%ROW(IC+1)-1 
         JC = C%COL(ICOL)
         IF (JC /= 0 .AND. JC <= G%NC) G%ZONES_LOCAL(C%COL(ICOL)) = G%N_ZONES
      ENDDO
   ENDIF


ENDDO PASS1_LOOP


! 
! Pass 2 of Aggregation:  Add unaggregated nodes to neighboring aggregate
! 
PASS2_LOOP: DO IC = 1, G%NC

   IF (G%ZONES_LOCAL(IC) /= 0) CYCLE            
   DO ICOL = C%ROW(IC), C%ROW(IC+1)-1
      JC = C%COL(ICOL)
      JZONE = G%ZONES_LOCAL(JC)
      IF (JZONE > 0) THEN
         IF (JC >= G%NC .OR. G%ZONE_CENTERS(JZONE)>0) THEN
            G%ZONES_LOCAL(IC) = -JZONE
            EXIT
         ENDIF
      ENDIF
   ENDDO
ENDDO PASS2_LOOP
!G%N_ZONES = G%N_ZONES-1                         !TODO: check
      


! 
! Pass 3 of Aggregation:  Process remaining nodes which have not been aggregated yet
! 
PASS3_LOOP: DO IC = 1, G%NC

   IZONE = G%ZONES_LOCAL(IC)

   ! cell IC has not been aggregated
   IF (IZONE /= 0) THEN
      IF (IZONE > 0) THEN
         G%ZONES_LOCAL(IC) = IZONE 
      ELSE IF (IZONE == NSCARC_HUGE_INT ) THEN
         G%ZONES_LOCAL(IC) = -1
      ELSE
         G%ZONES_LOCAL(IC) = -IZONE 
      ENDIF
      CYCLE PASS3_LOOP
   ENDIF

   G%ZONES_LOCAL(IC) = G%N_ZONES
   G%ZONE_CENTERS(G%N_ZONES) = IC

   DO ICOL = C%ROW(IC), C%ROW(IC+1)-1
      JC = C%COL(ICOL)
      IF (JC <= G%NC .AND. G%ZONES_LOCAL(JC) == 0) G%ZONES_LOCAL(JC) = G%N_ZONES
   ENDDO
   G%N_ZONES = G%N_ZONES+1

ENDDO PASS3_LOOP

IF (MINVAL(G%ZONES_LOCAL) < 0) THEN
   WRITE(*,*) MY_RANK+1, ':CAUTION: CELL ',MINLOC(G%ZONES_LOCAL), ' HAS NOT BEEN AGGREGATED DURING AGGREGATION'
   CALL MPI_FINALIZE(IERROR)
   STOP
ENDIF
      

END SUBROUTINE SCARC_SETUP_COARSENING_AGGREGATION


! --------------------------------------------------------------------------------------------------------------------
!> \brief Selfdefined geometric motivated aggregation procedure using cubic zones
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSENING_CUBIC(LF, LC, GF, GC)
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN):: LF, LC
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN):: GF, GC
INTEGER:: NXM, NYM, NZM, NXD, NYD, NZD, NXI, NYI, NZI
INTEGER:: IX, IY, IZ, IXZ, IYZ, IZZ, IXP, IYP, IZP, IX0, IY0, IZ0, IC
INTEGER, DIMENSION(:), ALLOCATABLE:: OFFX, OFFY, OFFZ
LOGICAL:: BFIRST

CROUTINE = 'SCARC_SETUP_COARSENING_CUBIC'

NXM = MOD(LF%NX, 2)
NXD = LF%NX/2

IF (TWO_D) THEN
   NYM = 0
   NYD = 1
ELSE
   NYM = MOD(LF%NY, 2)
   NYD = LF%NY/2
ENDIF

NZM = MOD(LF%NZ, 2)
NZD = LF%NZ/2

! Temporarily-to prevent failure of following algorithm

IF ((LF%NX < 4) .OR. (.NOT.TWO_D .AND. LF%NY < 4) .OR. (LF%NZ < 4)) THEN 
   WRITE(*,*) 'Grid dimensions too small für GMG-like aggregation'
   CALL MPI_FINALIZE(IERROR)
   STOP
ENDIF

CALL SCARC_ALLOCATE_INT1 (OFFX, 1, NXD, NSCARC_INIT_ZERO, 'OFFX', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OFFY, 1, NYD, NSCARC_INIT_ZERO, 'OFFY', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OFFZ, 1, NZD, NSCARC_INIT_ZERO, 'OFFZ', CROUTINE)

! If even number of cells in x-direction, use patch length of 2 as in default GMG
! else insert a patch of 3 after the first quarter of cells in x-direction

IF (NXM == 0) THEN
   OFFX = 2
ELSE
   NXI = MAX(NXD/4, 1)
   DO IX = 1, NXI
      OFFX(IX) = 2
   ENDDO
   OFFX(NXI+1) = 3
   DO IX = NXI+2, NXD
      OFFX(IX) = 2
   ENDDO
ENDIF

! If even number of cells in y-direction, use patch length of 2 as in default GMG
! else insert a patch of 3 after the first third of cells in x-direction

IF (TWO_D) THEN
   OFFY = 0
ELSE
   IF (NYM == 0) THEN
      OFFY = 2
   ELSE
      NYI = MAX(NYD/3, 1)
      DO IY = 1, NYI
         OFFY(IY) = 2
      ENDDO
      OFFY(NYI+1) = 3
      DO IY = NYI+2, NYD
         OFFY(IY) = 2
      ENDDO
   ENDIF
ENDIF

! If even number of cells in x-direction, use patch length of 2 as in default GMG
! else insert a patch of 3 after the first half of cells in x-direction
! the idea is to use different portions in the different coordinate direction to prevent local concentrations

IF (NZM == 0) THEN
   OFFZ = 2
ELSE
   NZI = MAX(NZD/2, 1)
   DO IZ = 1, NZI
      OFFZ(IZ) = 2
   ENDDO
   OFFZ(NZI+1) = 3
   DO IZ = NZI+2, NZD
      OFFZ(IZ) = 2
   ENDDO
ENDIF

LC%NX = NXD
LC%NY = NYD
LC%NZ = NZD

CALL SCARC_ALLOCATE_INT3 (GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER', CROUTINE)
CALL SCARC_ALLOCATE_LOG3 (LC%IS_SOLID, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_TRUE, 'LC%IS_SOLID', CROUTINE)


GF%ZONES_LOCAL = 0
GF%ZONES_GLOBAL = 0
GF%ZONE_CENTERS = 0
DIMENSION_IF: IF (TWO_D) THEN

   IZ0 = 1
   DO IZ = 1, NZD
      IX0 = 1
      DO IX = 1, NXD

         BFIRST = .TRUE.
         DO IZZ = 0, OFFZ(IZ)-1
            DO IXZ = 0, OFFX(IX)-1
               IXP = IX0+IXZ
               IZP = IZ0+IZZ
               IF (IS_UNSTRUCTURED .AND. LF%IS_SOLID(IXP, 1, IZP)) CYCLE
               IC = GF%CELL_NUMBER(IXP, 1, IZP)
               IF (BFIRST) THEN
                  GF%N_ZONES = GF%N_ZONES+1
                  BFIRST = .FALSE. 
                  GF%ZONE_CENTERS(GF%N_ZONES) = IC
                  GC%CELL_NUMBER(IX, 1, IZ) = GF%N_ZONES
                  LC%IS_SOLID(IX, 1, IZ) = .FALSE.
               ENDIF
               GF%ZONES_LOCAL(IC) = GF%N_ZONES
            ENDDO
         ENDDO
         IX0 = IX0+OFFX(IX)
      ENDDO
      IZ0 = IZ0+OFFZ(IZ)
   ENDDO

ELSE

   IZ0 = 1
   DO IZ = 1, NZD
      IY0 = 1
      DO IY = 1, NYD
         IX0 = 1
         DO IX = 1, NXD

            BFIRST = .TRUE.
            DO IZZ = 0, OFFZ(IZ)-1
               DO IYZ = 0, OFFY(IY)-1
                  DO IXZ = 0, OFFX(IX)-1
                     IXP = IX0+IXZ
                     IYP = IY0+IYZ
                     IZP = IZ0+IZZ
                     IF (IS_UNSTRUCTURED .AND. LF%IS_SOLID(IXP, IYP, IZP)) CYCLE
                     IC = GF%CELL_NUMBER(IXP, IYP, IZP)
                     IF (BFIRST) THEN
                        GF%N_ZONES = GF%N_ZONES+1
                        BFIRST = .FALSE. 
                        GF%ZONE_CENTERS(GF%N_ZONES) = IC
                        GC%CELL_NUMBER(IX, IY, IZ) = GF%N_ZONES
                        LC%IS_SOLID(IX, IY, IZ) = .FALSE.
                     ENDIF
                     IC = GF%CELL_NUMBER(IXP, IYP, IZP)
                     GF%ZONES_LOCAL(IC) = GF%N_ZONES
                  ENDDO
               ENDDO
            ENDDO
            IX0 = IX0+OFFX(IX)
         ENDDO
         IY0 = IY0+OFFY(IY)
      ENDDO
      IZ0 = IZ0+OFFZ(IZ)
   ENDDO
ENDIF DIMENSION_IF

CALL SCARC_DEALLOCATE_INT1 (OFFX, 'OFFX', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (OFFY, 'OFFY', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (OFFZ, 'OFFZ', CROUTINE)


END SUBROUTINE SCARC_SETUP_COARSENING_CUBIC

! --------------------------------------------------------------------------------------------------------------------
!> \brief Selfdefined geometric motivated aggregation procedure using cubic zones
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSENING_GMG(LF, LC, GF, GC)
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN):: LF, LC
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN):: GF, GC
INTEGER:: MODX, MODY, MODZ
INTEGER:: RELX, RELY, RELZ
INTEGER:: IX, IY, IZ, IXZ, IYZ, IZZ, IXP, IYP, IZP, IX0, IY0, IZ0, IC
INTEGER, DIMENSION(:), ALLOCATABLE:: OFFX, OFFY, OFFZ
LOGICAL:: BFIRST

CROUTINE = 'SCARC_SETUP_COARSENING_GMG'

MODX = MOD(LF%NX, 2)
LC%NX = INT(FLOOR(REAL(LF%NX/2), EB), 4) + MODX

IF (TWO_D) THEN
   MODY = 0
   LC%NY = 1
ELSE
   MODY = MOD(LF%NY, 2)
   LC%NY = INT(FLOOR(REAL(LF%NY/2), EB), 4) + MODY
ENDIF

MODZ = MOD(LF%NZ, 2)
LC%NZ = INT(FLOOR(REAL(LF%NZ/2), EB), 4) + MODZ

CALL SCARC_ALLOCATE_INT1 (OFFX, 1, LC%NX, NSCARC_INIT_ZERO, 'OFFX', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OFFY, 1, LC%NY, NSCARC_INIT_ZERO, 'OFFY', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OFFZ, 1, LC%NZ, NSCARC_INIT_ZERO, 'OFFZ', CROUTINE)

OFFX = 2;  IF (MODX /= 0) OFFX(LC%NX) = 1
IF (TWO_D) THEN
   OFFY = 0;  LC%NY = 1
ELSE
   OFFY = 2;  IF (MODY /= 0) OFFY(LC%NY) = 1
ENDIF
OFFZ = 2;  IF (MODZ /= 0) OFFZ(LC%NZ) = 1

RELX = INT(CEILING(REAL(LF%NX/LC%NX), EB), 4)
RELY = INT(CEILING(REAL(LF%NY/LC%NY), EB), 4)
RELZ = INT(CEILING(REAL(LF%NZ/LC%NZ), EB), 4)

CALL SCARC_ALLOCATE_INT3 (GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER', CROUTINE)
CALL SCARC_ALLOCATE_LOG3 (LC%IS_SOLID,    0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_TRUE,  'LC%IS_SOLID', CROUTINE)


GF%ZONES_LOCAL  = 0
GF%ZONES_GLOBAL = 0
GF%ZONE_CENTERS = 0

DIMENSION_IF: IF (TWO_D) THEN

   IZ0 = 1
   DO IZ = 1, LC%NZ
      IX0 = 1
      DO IX = 1, LC%NX

         BFIRST = .TRUE.
         DO IZZ = 0, OFFZ(IZ)-1
            DO IXZ = 0, OFFX(IX)-1
               IXP = IX0+IXZ
               IZP = IZ0+IZZ
               IF (IS_UNSTRUCTURED .AND. LF%IS_SOLID(IXP, 1, IZP)) CYCLE
               IC = GF%CELL_NUMBER(IXP, 1, IZP)
               IF (BFIRST) THEN
                  GF%N_ZONES = GF%N_ZONES+1
                  BFIRST = .FALSE. 
                  GF%ZONE_CENTERS(GF%N_ZONES) = IC
                  GC%CELL_NUMBER(IX, 1, IZ) = GF%N_ZONES
                  LC%IS_SOLID(IX, 1, IZ) = .FALSE.
               ENDIF
               GF%ZONES_LOCAL(IC) = GF%N_ZONES
            ENDDO
         ENDDO
         IX0 = IX0+OFFX(IX)
      ENDDO
      IZ0 = IZ0+OFFZ(IZ)
   ENDDO

ELSE

   IZ0 = 1
   DO IZ = 1, LC%NZ
      IY0 = 1
      DO IY = 1, LC%NY
         IX0 = 1
         DO IX = 1, LC%NX

            BFIRST = .TRUE.
            DO IZZ = 0, OFFZ(IZ)-1
               DO IYZ = 0, OFFY(IY)-1
                  DO IXZ = 0, OFFX(IX)-1
                     IXP = IX0+IXZ
                     IYP = IY0+IYZ
                     IZP = IZ0+IZZ
                     IF (IS_UNSTRUCTURED .AND. LF%IS_SOLID(IXP, IYP, IZP)) CYCLE
                     IC = GF%CELL_NUMBER(IXP, IYP, IZP)
                     IF (BFIRST) THEN
                        GF%N_ZONES = GF%N_ZONES+1
                        BFIRST = .FALSE. 
                        GF%ZONE_CENTERS(GF%N_ZONES) = IC
                        GC%CELL_NUMBER(IX, IY, IZ) = GF%N_ZONES
                        LC%IS_SOLID(IX, IY, IZ) = .FALSE.
                     ENDIF
                     IC = GF%CELL_NUMBER(IXP, IYP, IZP)
                     GF%ZONES_LOCAL(IC) = GF%N_ZONES
                  ENDDO
               ENDDO
            ENDDO
            IX0 = IX0+OFFX(IX)
         ENDDO
         IY0 = IY0+OFFY(IY)
      ENDDO
      IZ0 = IZ0+OFFZ(IZ)
   ENDDO

ENDIF DIMENSION_IF

CALL SCARC_DEALLOCATE_INT1 (OFFX, 'OFFX', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (OFFY, 'OFFY', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (OFFZ, 'OFFZ', CROUTINE)


END SUBROUTINE SCARC_SETUP_COARSENING_GMG

! --------------------------------------------------------------------------------------------------------------------
!> \brief Remove workspace on specified grid level which will no longer be needed after matrix setup
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CLEAN_WORKSPACE_SYSTEM(NL)
USE SCARC_POINTERS, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN):: NL
INTEGER:: NM

! TODO: deallocate arrays which are no longer used
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)
ENDDO

END SUBROUTINE SCARC_CLEAN_WORKSPACE_SYSTEM


! --------------------------------------------------------------------------------------------------------------------
!> \brief Remove workspace on specified grid level which will no longer be needed in SAMG method
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CLEAN_WORKSPACE_AMG(NL)
USE SCARC_POINTERS, ONLY: G, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN):: NL
INTEGER:: NM

CROUTINE = 'SCARC_CLEAN_WORKSPACE_AMG'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)
   IF (ALLOCATED(G%ZONE_CENTERS)) CALL SCARC_DEALLOCATE_INT1 (G%ZONE_CENTERS, 'G%ZONE_CENTERS', CROUTINE)
   IF (ALLOCATED(G%AUX1)) CALL SCARC_DEALLOCATE_REAL1 (G%AUX1, 'G%AUX1', CROUTINE)
   IF (ALLOCATED(G%AUX2)) CALL SCARC_DEALLOCATE_REAL1 (G%AUX2, 'G%AUX2', CROUTINE)
   IF (ALLOCATED(G%RR)) CALL SCARC_DEALLOCATE_REAL1 (G%RR, 'G%RR', CROUTINE)
   IF (ALLOCATED(G%QQ)) CALL SCARC_DEALLOCATE_REAL1 (G%QQ, 'G%QQ', CROUTINE)
ENDDO

END SUBROUTINE SCARC_CLEAN_WORKSPACE_AMG


! --------------------------------------------------------------------------------------------------------------------
!> \brief Perform relaxation of nullspac
! Perform AMG Jacobi :.. x = x-omega D^{-1} (Ax-b)
! Near-null space vector is given in vector G%NULLSPACE--> corresponds to x
! vector b is zero
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELAX_NULLSPACE(NL)
USE SCARC_POINTERS, ONLY: L, G, A, MG, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
#ifndef WITH_MKL
USE SCARC_VECTORS, ONLY: SCARC_DAXPY_CONSTANT_DOUBLE
#else
EXTERNAL:: DAXPBY
#endif
INTEGER, INTENT(IN):: NL
INTEGER:: IC, ICOL, NM, JC, JCG

CROUTINE = 'SCARC_RELAX_NULLSPACE'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL) 
   CALL SCARC_ALLOCATE_REAL1 (G%AUX1, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%AUX1', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (G%AUX2, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%AUX2', CROUTINE)
ENDDO

! For coarser levels exchange numbers and values of second layer cells with are needed for nullspace computation

IF (SCARC_MULTIGRID_RELAXING .AND. NL > NLEVEL_MIN) THEN
   CALL SCARC_IDENTIFY_LAYER2(NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_LAYER2_NUMS, NSCARC_NONE, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_LAYER2_VALS, NSCARC_NONE, NL)
ENDIF

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)         
   A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)

   MG => L%MG
   MG%OMEGA = 1.0_EB/MG%APPROX_SPECTRAL_RADIUS
   !MG%OMEGA = 4.0_EB/(3.0_EB*MG%APPROX_SPECTRAL_RADIUS)


   ! Compute defect to near-null-space vector: d = Ax-b, in this case
   !    'x' corresponds to nullspace vector consisting of only '1'-entries 
   !    'b' corresponds to zero vector 

   ! On finest level NULLSPACE vector is preset with 1, so matrix entries can simply be added
   IF (NL == NLEVEL_MIN) THEN

      CALL SCARC_ALLOCATE_REAL1 (G%NULLSPACE, 1, G%NCE2, NSCARC_INIT_ONE, 'G%NULLSPACE', CROUTINE)
      FINE_CELLS_LOOP: DO IC = 1, G%NC
         G%AUX2(IC) = 0.0_EB
         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1                          
            G%AUX2(IC) = G%AUX2(IC) + A%VAL(ICOL)
         ENDDO
      ENDDO FINE_CELLS_LOOP

   ! on coarser levels NULLSPACE vector was set in preceding coarsening loop to R-vector from QR-decomposition 
   ELSE

      G%AUX1 = G%NULLSPACE
      COARSE_CELLS_LOOP: DO IC = 1, G%NC
   
         G%AUX2(IC) = 0.0_EB
         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1                          
            JC = A%COL(ICOL)
            IF (JC /= 0) THEN
               G%AUX2(IC) =  G%AUX2(IC) + A%VAL(ICOL) * G%AUX1(JC)
            ELSE 
               JCG = A%COLG(ICOL)
               JC = FINDLOC (G%LOCAL_TO_GLOBAL(1:G%NCE2), VALUE = JCG, DIM = 1)
               IF (JC /= 0) THEN
                  G%AUX2(IC) = G%AUX2(IC) + A%VAL(ICOL) * G%AUX1(JC)
               ELSE
                  JC = FINDLOC (G%ELAYER2_NUMS(1: L%N_LAYER2_TOTAL), VALUE = JCG, DIM = 1)
                  IF (JC /= 0) G%AUX2(IC) = G%AUX2(IC) + A%VAL(ICOL) * G%ELAYER2_VALS(JC)
               ENDIF
            ENDIF
         ENDDO
   
      ENDDO COARSE_CELLS_LOOP

   ENDIF
   
   ! Scale it by parameter omega and inverse of diagonal:   d = omega D^{-1} d

   DO IC = 1, G%NC
      G%AUX2(IC) = MG%OMEGA*G%DIAG(IC) * G%AUX2(IC) 
   ENDDO
   
  ! Get new iterate:   x = x-d

#ifdef WITH_MKL
  CALL DAXPBY(G%NC, -1.0_EB, G%AUX2, 1, 1.0_EB, G%NULLSPACE, 1)
#else
  CALL SCARC_DAXPY_CONSTANT_DOUBLE(G%NC, -1.0_EB, G%AUX2, 1.0_EB, G%NULLSPACE)
#endif
   

ENDDO

IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_NULLSPACE, NSCARC_NONE, NL)

END SUBROUTINE SCARC_RELAX_NULLSPACE


! --------------------------------------------------------------------------------------------------------------
!> \brief Setup basic structure of Prolongation matrix
! This concerns the setting of the number of rows and the column pointers
! The values are still missing and are set in a later step
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_ZONE_OPERATOR(NL)
USE SCARC_POINTERS, ONLY: GF, GC, AF, ZF, SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN):: NL
INTEGER:: NM, ICC, ICCL, ICCG, IC, IP, N_ROW, N_VAL
LOGICAL:: IS_INCLUDED 

CROUTINE = 'SCARC_SETUP_ZONE_OPERATOR'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                    

   AF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
   ZF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_ZONES)

   ! First use very conservative bounds for the size of the zones operator matrix 
   ! reduce it later once the real size is known

   ZF%N_VAL = AF%N_VAL                  
   ZF%N_ROW = AF%N_VAL                  
   CALL SCARC_ALLOCATE_CMATRIX (ZF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_MINIMAL, 'GF%ZONES', CROUTINE)
  
   ! Again conservative upper bound for length-to be reduced later 

   CALL SCARC_ALLOCATE_INT1 (GF%ZONES_LOCAL, 1, GF%NCE2, NSCARC_INIT_ZERO, 'GF%ZONES_LOCAL', CROUTINE)


   ! Based on global zone numbers determine local zone numbers within mesh

   IP = 1
   ICC = 1
   ZF%ROW(ICC) = 1
   DO ICCL = 1, GC%NCE2
      ICCG = GC%LOCAL_TO_GLOBAL(ICCL)
      IS_INCLUDED = .FALSE.
      DO IC = 1, GF%NCE2
         IF (GF%ZONES_LOCAL(IC) /= ICCL) CYCLE
         IS_INCLUDED = .TRUE.
         ZF%COL(IP)  = IC
         IP = IP+1
      ENDDO
      IF (IS_INCLUDED) THEN
         ICC = ICC+1
         ZF%ROW(ICC) = IP
      ENDIF
   ENDDO

   N_ROW = ICC
   N_VAL = IP-1

   ZF%N_ROW = N_ROW
   ZF%N_VAL = N_VAL

   CALL SCARC_REDUCE_CMATRIX (ZF, 'ZF', CROUTINE)


ENDDO

END SUBROUTINE SCARC_SETUP_ZONE_OPERATOR


! --------------------------------------------------------------------------------------------------------------------
!> \brief Determine tentative prolongator for current level by computing QR-decomposition of smoothed 
! nullspace vector and set nullspace for next level
! Compute the tentative prolongator, T, which is a tentative interpolation
! matrix from the coarse-grid to the fine-grid.  T exactly interpolates  B_fine = T B_coarse.
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PROLONGATION_AMG(NL)
USE SCARC_POINTERS, ONLY: S, L, G, A, OA, P, OP, GF, GC, PF, OPF, Z, MG, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_MULTIGRID, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN):: NL
REAL(EB):: DSUM, SCAL !, TOL = 1.0E-12_EB
INTEGER:: NM, NOM, IC, JC, ICC, ICC0, ICOL, ICCOL, JCCOL, IP0, IP, JCC, IQ, INBR, NLEN

CROUTINE = 'SCARC_SETUP_PROLONGATION'

! Allocate several workspaces (with conservative bounds which will be reduced later)
!    - Prolongation matrix on internal part of mesh 
!    - for every neighbor small Prolongation matrix for corresponding overlap
!    - vectors Q and R for QR-decomposition of aggregation zones operator
! Initialize QR-decomposition
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

   A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
   P => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)
   Z => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_ZONES)

   P%N_VAL = A%N_VAL+1        
   P%N_ROW = G%NCE+1                  
   CALL SCARC_ALLOCATE_CMATRIX (P, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'G%PROLONGATION', CROUTINE)

   DO INBR = 1, S%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

      OA => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_POISSON)
      OP => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_PROLONGATION)

      OP%N_VAL = OA%N_VAL+1              ! TODO : CHECK : MUCH TOO BIG !!!
      OP%N_ROW = G%NCE+1                 ! TODO : CHECK : MUCH TOO BIG !!!
      CALL SCARC_ALLOCATE_CMATRIX (OP, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OG%PROLONGATION', CROUTINE)

   ENDDO

   CALL SCARC_ALLOCATE_REAL1 (G%QQ, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%QQ', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (G%RR, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%RR', CROUTINE)     ! TODO check length!

 
   ! Copy blocks into Q according to aggregation zones and compute norms for single ZONES
   ! In each cell corresponding to a single zone, store square-sum of entries
 
   IQ = 1
   G%AUX1 = 0.0_EB
   DO ICC = 1, Z%N_ROW-1
      DSUM = 0.0_EB
      DO ICOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
         IC = Z%COL(ICOL)
         G%QQ(IQ) = G%NULLSPACE(IC)
         DSUM = DSUM+G%QQ(IQ)**2
         IQ = IQ+1
      ENDDO
      DO ICOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
         G%AUX1(Z%COL(ICOL)) = DSUM
      ENDDO
      G%RR(ICC) = DSUM
   ENDDO

   CALL SCARC_REDUCE_REAL1 (G%QQ, IQ, 'G%QQ', CROUTINE)
   CALL SCARC_REDUCE_REAL1 (G%RR, Z%N_ROW, 'G%RR', CROUTINE)


ENDDO

 
! Exchange sums of nullspace entries within single aggregation zones
 
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_AUXILIARY, NSCARC_NONE, NL)
   

 
! Build norms over single zones and scale Q-entries by norms
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                          
   Z => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_ZONES)

   DO ICC = 1, Z%N_ROW-1
      ICOL = Z%ROW(ICC)
      IC = Z%COL(ICOL)
      G%RR(ICC) = SQRT(G%AUX1(IC))
   ENDDO

   IQ = 1
   DO ICC = 1, Z%N_ROW-1
      DO ICOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
         IC = Z%COL(ICOL)
         G%QQ(IQ) = G%QQ(IQ)/G%RR(ICC)
         IQ = IQ+1
      ENDDO
   ENDDO

ENDDO
   

! ----------- Relax prolongator:
! Smooth the tentative prolongator, so that it's accuracy is greatly improved for algebraically smooth error.
! Compute:                P =  P-A_Dinv*P   
! with:                   A_Dinv = 4/3 * 1/rho*D^{-1} A   
 
! First step: Compute P_0: = A_Dinv*Q
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                 

   A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
   P => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)
   Z => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_ZONES)

   MG => L%MG
   MG%OMEGA = 4.0_EB/3.0_EB                                         ! currently used default
   IF (SCARC_MULTIGRID_RELAXING) THEN
      SCAL = MG%OMEGA/MG%APPROX_SPECTRAL_RADIUS                     ! for testing purposes rho is set to 2 currently
   ELSE
      SCAL = 0.0_EB
   ENDIF
   
   IP = 1
   IP0 = IP
   P%ROW(1) = IP
   DO IC = 1, G%NC
      DO ICC = 1, Z%N_ROW-1
   
         DSUM = 0.0_EB
         DO ICCOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
            JC = Z%COL(ICCOL)
            ICOL = SCARC_MATCH_MATRIX_COLUMN(A, IC, JC)
            IF (ICOL /= -1) THEN
               ICC0 = ICC
               DSUM = DSUM-SCAL*G%DIAG(IC) * A%VAL(ICOL) * G%QQ(ICCOL)
            ENDIF
         ENDDO
   
         IF (ABS(DSUM) /= 0.0_EB) THEN
            P%VAL(IP) = DSUM
            P%COL(IP) = ICC
            IP = IP+1
         ENDIF
   
      ENDDO
 
      ! take care that at least one entry per fine cell is generated
      IF (IP == IP0) THEN
         P%VAL(IP) = 0.0_EB
         P%COL(IP) = ICC0
         IP = IP+1
      ENDIF
      IP0 = IP

      P%ROW(IC+1) = IP
   ENDDO
 
   P%N_VAL = IP-1

ENDDO
   
 
! Second step: Compute P: = P-P_0
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   P => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)
   Z => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_ZONES)


   DO ICC = 1, Z%N_ROW-1
      DO ICCOL = Z%ROW(ICC), Z%ROW(ICC+1) - 1
         IC = Z%COL(ICCOL)

         IF (IC > G%NC) CYCLE
         DO JCCOL = P%ROW(IC), P%ROW(IC+1) - 1
            JCC = P%COL(JCCOL)
            IF (JCC == ICC) THEN
               P%VAL(JCCOL) = P%VAL(JCCOL) + G%QQ(ICCOL)
               CYCLE
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   
ENDDO

! Determine global columns array for Prolongation matrix
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                
   PF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)

   DO IC = 1, GF%NC
      DO JCCOL = PF%ROW(IC), PF%ROW(IC+1) - 1
         JCC = PF%COL(JCCOL)
         PF%COLG(JCCOL) = GC%LOCAL_TO_GLOBAL(JCC)
      ENDDO
   ENDDO
   

ENDDO

 
! Exchange resulting columns and values of Prolongation matrix and extract exchanged data from 
! overlapping parts with single neighbors and attach them to main matrix
 
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLSG, NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_PROLONGATION, 0, NL)
ENDIF
   
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1) 
   PF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)

   NLEN = MAX(GC%NCE2, 4 * (GC%NCE2-GC%NC+2) + 50)           ! TODO check length
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_LOCAL,  1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_LOCAL',  CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_GLOBAL, 1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_GET_CELL_DEPENDENCIES_GALERKIN(GC, GF, PF, NLEN)

   CALL SCARC_REDUCE_INT1 (GC%CELLS_LOCAL, GC%NC_GALERKIN, 'GC%CELLS_LOCAL', CROUTINE)
   CALL SCARC_REDUCE_INT1 (GC%CELLS_GLOBAL, GC%NC_GALERKIN, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_REDUCE_CMATRIX (PF, 'P%PROLONGATION', CROUTINE)
   DO INBR = 1, S%N_NEIGHBORS
      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
      OPF => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_PROLONGATION)
      CALL SCARC_REDUCE_CMATRIX (OPF, 'OP%PROLONGATION', CROUTINE)
   ENDDO


ENDDO

END SUBROUTINE SCARC_SETUP_PROLONGATION_AMG


! -------------------------------------------------------------------------------------------------------------
!> \brief Setup restriction and prolongation matrices in case of GMG-like coarsening
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TRANSFER_GMG(NL)
USE SCARC_POINTERS, ONLY: S, LC, LF, GF, GC, AF, RF, ZF, PF, OPF, &
                          SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_MULTIGRID, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN):: NL
INTEGER:: NM, IC, JC, ICC, ICF, ICCOL, IRCOL, IP, JCC, JCCOL, NOM, INBR, NLEN
INTEGER:: IS, IXC, IZC, IXF, IZF, IOFFX, IOFFZ, ICC0(4)
INTEGER:: STENCIL(16) = 0
INTEGER:: RN, RE, RS, RW

CROUTINE = 'SCARC_SETUP_PROLONGATION_GMG'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                                    

   AF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
   ZF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_ZONES)
   RF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_RESTRICTION)
   PF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)


   SELECT CASE (TYPE_INTERPOL)

      CASE (NSCARC_INTERPOL_CONSTANT)

         RF%N_VAL = GF%NCE2
         RF%N_ROW = GC%NCE2+1
         CALL SCARC_ALLOCATE_CMATRIX (RF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%RESTRICTION', CROUTINE)

         PF%N_VAL = AF%N_VAL+1        
         PF%N_ROW = GF%NCE2+1
         CALL SCARC_ALLOCATE_CMATRIX (PF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%PROLONGATION', CROUTINE)

         IP = 1
         RF%ROW(1) = IP
         DO ICC = 1, ZF%N_ROW-1
            DO ICCOL = ZF%ROW(ICC), ZF%ROW(ICC+1)-1
               JC = ZF%COL(ICCOL)
               RF%COL(IP) = ZF%COL(ICCOL)
               RF%COLG(IP) = GF%LOCAL_TO_GLOBAL(ZF%COL(ICCOL))
               RF%VAL(IP) = 1.00_EB
               IP = IP+1
            ENDDO
            RF%ROW(ICC+1) = IP
            !WRITE(*,*) ICC, ':SUM(RF)=',SUM(RF%VAL(RF%ROW(ICC):RF%ROW(ICC+1)-1))
         ENDDO
      
         IP = 1
         PF%ROW(1) = IP
         DO IC = 1, GF%NC
            DO ICC = 1, ZF%N_ROW-1
               COLUMN_LOOP: DO ICCOL = ZF%ROW(ICC), ZF%ROW(ICC+1)-1
                  IF (ZF%COL(ICCOL) == IC) THEN
                     PF%COL(IP) = ICC
                     PF%COLG(IP) = GC%LOCAL_TO_GLOBAL(ICC)
                     PF%VAL(IP) = 0.25_EB
                     IP = IP+1
                     EXIT COLUMN_LOOP
                  ENDIF
               ENDDO COLUMN_LOOP
            ENDDO
            PF%ROW(IC+1) = IP
            !WRITE(*,*) ICC, ':SUM(PF)=',SUM(PF%VAL(PF%ROW(IC):PF%ROW(IC+1)-1))
         ENDDO

      CASE (NSCARC_INTERPOL_BILINEAR)

         RF%N_VAL = 16*GF%NC
         RF%N_ROW = GC%NC+1
         CALL SCARC_ALLOCATE_CMATRIX (RF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%RESTRICTION', CROUTINE)

         PF%N_VAL = AF%N_VAL+1        
         PF%N_ROW = GF%NC+1
         CALL SCARC_ALLOCATE_CMATRIX (PF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%PROLONGATION', CROUTINE)

         IP = 1
         RF%ROW(1) = IP
         DO IZC = 1, LC%NZ
            DO IXC = 1, LC%NX

               ICC = (IZC-1) * LC%NX+IXC

               IXF = 2*IXC
               IZF = 2*IZC
               
               IF (IXC == 1 .AND. IZC == 1) THEN
                  RN = 1; RE = 1; RS = 0; RW = 0
               ELSE IF (IXC == 1 .AND. IZC == LC%NZ) THEN
                  RN = 0; RE = 1; RS = 1; RW = 0
               ELSE IF (IXC == LC%NX .AND. IZC == 1) THEN
                  RN = 1; RE = 0; RS = 0; RW = 1
               ELSE IF (IXC == LC%NX .AND. IZC == LC%NZ) THEN
                  RN = 0; RE = 0; RS = 1; RW = 1
               ELSE IF (IXC == 1) THEN
                  RN = 1; RE = 1; RS = 1; RW = 0
               ELSE IF (IXC == LC%NX) THEN
                  RN = 1; RE = 0; RS = 1; RW = 1
               ELSE IF (IZC == 1) THEN
                  RN = 1; RE = 1; RS = 0; RW = 1
               ELSE IF (IZC == LC%NZ) THEN
                  RN = 0; RE = 1; RS = 1; RW = 1
               ELSE
                  RN = 1; RE = 1; RS = 1; RW = 1
               ENDIF

               STENCIL(13:16) =  (/    RN*RW,    RN *(2+RW),    RN *(2+RE),    RN*RE /)
               STENCIL( 9:12) =  (/ (2+RN)*RW, (2+RN)*(2+RW), (2+RN)*(2+RE), (2+RN)*RE /)
               STENCIL( 5: 8) =  (/ (2+RS)*RW, (2+RS)*(2+RW), (2+RS)*(2+RE), (2+RS)*RE /)
               STENCIL( 1: 4) =  (/    RS*RW,    RS *(2+RW),    RS *(2+RE),    RS*RE /)

               IS = 1
               DO IOFFZ = -2, 1
                  DO IOFFX = -2, 1
                     CALL PROCESS_FINE_CELL (IXF+IOFFX, 1, IZF+IOFFZ, IP, STENCIL(IS))
                     IS = IS+1
                  ENDDO
               ENDDO

               RF%ROW(ICC+1) = IP
               !WRITE(*,*) ICC, ':SUM(RF)=',SUM(RF%VAL(RF%ROW(ICC):RF%ROW(ICC+1)-1))


            ENDDO
         ENDDO

         IP = 1
         PF%ROW(1) = IP
         DO IC = 1, GF%NC
            DO ICC = 1, GC%NC 
               COLUMN_LOOP2: DO IRCOL = RF%ROW(ICC), RF%ROW(ICC+1)-1
                  IF (RF%COL(IRCOL) == IC) THEN
                     PF%COL(IP) = ICC
                     PF%COLG(IP) = GC%LOCAL_TO_GLOBAL(ICC)
                     PF%VAL(IP) = RF%VAL(IRCOL)
                     IP = IP+1
                     EXIT COLUMN_LOOP2
                  ENDIF
               ENDDO COLUMN_LOOP2
            ENDDO
            PF%ROW(IC+1) = IP
            !WRITE(*,*) ICC, ':SUM(PF)=',SUM(PF%VAL(PF%ROW(IC):PF%ROW(IC+1)-1))
         ENDDO

         PF%VAL = PF%VAL/16.0_EB
         IF (TWO_D) THEN
            RF%VAL = RF%VAL/4.0_EB
         ELSE
            RF%VAL = RF%VAL/2.0_EB
         ENDIF

      CASE (NSCARC_INTERPOL_BILINEAR2)

         RF%N_VAL = 16*GF%NC
         RF%N_ROW = GC%NC+1
         CALL SCARC_ALLOCATE_CMATRIX (RF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%RESTRICTION', CROUTINE)

         PF%N_VAL = 4*GF%NC
         PF%N_ROW = GF%NC+1
         CALL SCARC_ALLOCATE_CMATRIX (PF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%PROLONGATION', CROUTINE)

         IP = 1
         PF%ROW(1) = IP
         DO IZF = 1, LF%NZ
            DO IXF = 1, LF%NX

               ICF = GF%CELL_NUMBER(IXF, 1, IZF)

               !IXC = CEILING(REAL((IXF+1)/2), EB) 
               !IZC = CEILING(REAL((IZF+1)/2), EB)
               IXC = MOD(IXF, 2) + 1
               IZC = MOD(IZF, 2) + 1

               ICC0(1) = (IZC-1) * LC%NX+IXC
               ICC0(2) = (IZC-1) * LC%NX+IXC+1
               ICC0(3) = IZC*LC%NX+IXC
               ICC0(4) = IZC*LC%NX+IXC+1

 WRITE(*,*) 'IXF, IZF, ICF, IXC, IZC, ICC0(1:4):', IXF, IZF, ICF, IXC, IZC, ICC0(1:4)
               
 IF (IXC /= 0 .AND. IZC /= 0) CALL PROCESS_COARSE_CELL (ICC0(1), IP, 9)
 IF (IXC /= 0 .AND. IZC /= 0) CALL PROCESS_COARSE_CELL (ICC0(2), IP, 3)
 IF (IXC /= 0 .AND. IZC /= 0) CALL PROCESS_COARSE_CELL (ICC0(3), IP, 3)
 IF (IXC /= 0 .AND. IZC /= 0) CALL PROCESS_COARSE_CELL (ICC0(4), IP, 1)

               PF%ROW(ICF+1) = IP

            ENDDO
         ENDDO

         IP = 1
         PF%ROW(1) = IP
         DO IC = 1, GF%NC
            DO ICC = 1, GC%NC 
               COLUMN_LOOP3: DO IRCOL = RF%ROW(ICC), RF%ROW(ICC+1)-1
                  IF (RF%COL(IRCOL) == IC) THEN
                     PF%COL(IP) = ICC
                     PF%COLG(IP) = GC%LOCAL_TO_GLOBAL(ICC)
                     PF%VAL(IP) = RF%VAL(IRCOL)
                     IP = IP+1
                     EXIT COLUMN_LOOP3
                  ENDIF
               ENDDO COLUMN_LOOP3
            ENDDO
            PF%ROW(IC+1) = IP
            WRITE(*,*) ICC, ':SUM(PF)=',SUM(PF%VAL(PF%ROW(IC):PF%ROW(IC+1)-1))
         ENDDO

         PF%VAL = PF%VAL/16.0_EB
         IF (TWO_D) THEN
            RF%VAL = RF%VAL/4.0_EB
         ELSE
            RF%VAL = RF%VAL/2.0_EB
         ENDIF

   END SELECT


ENDDO

! Determine global columns array for Prolongation matrix
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                
   PF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)

   DO IC = 1, GF%NC
      DO JCCOL = PF%ROW(IC), PF%ROW(IC+1) - 1
         JCC = PF%COL(JCCOL)
         PF%COLG(JCCOL) = GC%LOCAL_TO_GLOBAL(JCC)
      ENDDO
   ENDDO
   

ENDDO

 
! Exchange resulting columns and values of Prolongation matrix and extract exchanged data from 
! overlapping parts with single neighbors and attach them to main matrix
 
IF (NMESHES > 1 .AND. TYPE_COARSENING /= NSCARC_COARSENING_GMG) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLSG, NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_PROLONGATION, 0, NL)
ENDIF
   
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1) 
   PF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)

   NLEN = MAX(GC%NCE2, 4 * (GC%NCE2-GC%NC+2) + 50)           ! TODO check length
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_LOCAL,  1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_LOCAL',  CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_GLOBAL, 1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_GET_CELL_DEPENDENCIES_GALERKIN(GC, GF, PF, NLEN)

   CALL SCARC_REDUCE_INT1 (GC%CELLS_LOCAL, GC%NC_GALERKIN, 'GC%CELLS_LOCAL', CROUTINE)
   CALL SCARC_REDUCE_INT1 (GC%CELLS_GLOBAL, GC%NC_GALERKIN, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_REDUCE_CMATRIX (PF, 'P%PROLONGATION', CROUTINE)
   IF (TYPE_COARSENING /= NSCARC_COARSENING_GMG) THEN
      DO INBR = 1, S%N_NEIGHBORS
         NOM = S%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
         OPF => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_PROLONGATION)
         CALL SCARC_REDUCE_CMATRIX (OPF, 'OP%PROLONGATION', CROUTINE)
      ENDDO
   ENDIF


ENDDO

END SUBROUTINE SCARC_SETUP_TRANSFER_GMG

SUBROUTINE PROCESS_FINE_CELL(IXF, IYF, IZF, IP, VAL)
USE SCARC_POINTERS, ONLY: RF, GF
INTEGER, INTENT(IN):: IXF, IYF, IZF, VAL
INTEGER, INTENT(INOUT):: IP
IF (VAL == 0) RETURN
RF%COL(IP) = ICF (IXF, IYF, IZF)
RF%COLG(IP) = GF%LOCAL_TO_GLOBAL(RF%COL(IP))
RF%VAL(IP) = VAL
IP = IP+1
END SUBROUTINE PROCESS_FINE_CELL

SUBROUTINE PROCESS_COARSE_CELL(ICC, IP, VAL)
USE SCARC_POINTERS, ONLY: PF, GC
INTEGER, INTENT(IN):: ICC, VAL
INTEGER, INTENT(INOUT):: IP
IF (VAL == 0) RETURN
PF%COL(IP) = ICC
PF%COLG(IP) = GC%LOCAL_TO_GLOBAL(ICC)
PF%VAL(IP) = VAL
IP = IP+1
END SUBROUTINE PROCESS_COARSE_CELL

INTEGER FUNCTION ICF (IXF, IYF, IZF)
USE SCARC_POINTERS, ONLY : LF
INTEGER, INTENT(IN):: IXF, IYF, IZF
ICF = (IZF-1) * LF%NX*LF%NY + (IYF-1) * LF%NX+IXF
RETURN
END FUNCTION

INTEGER FUNCTION ICC (IXC, IYC, IZC)
USE SCARC_POINTERS, ONLY : LC
INTEGER, INTENT(IN):: IXC, IYC, IZC
ICC = (IZC-1) * LC%NX*LC%NY + (IYC-1) * LC%NX+IXC
RETURN
END FUNCTION

! -------------------------------------------------------------------------------------------------------------
!> \brief Determine on which overlapping global coarse cells are given mesh depends (also considering diagonal connections)
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_CELL_DEPENDENCIES_GALERKIN(GC, GF, PF, NLEN)
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN):: GC, GF
TYPE (SCARC_CMATRIX_TYPE), POINTER, INTENT(IN):: PF
INTEGER, INTENT(IN):: NLEN
INTEGER:: IZL, IZG, IP, ICOL, IC, IFOUND1, IFOUND2, IFOUND3

IP = 1
PROLONGATION_CELLS_LOOP: DO IC = 1, GF%NCE

   ! Check if zone number used in given row of Prolongation matrix is already accounted for
   DO ICOL = PF%ROW(IC), PF%ROW(IC+1) - 1

      IZL = PF%COL(ICOL)
      IZG = PF%COLG(ICOL)

      IFOUND1 = FINDLOC (GC%CELLS_GLOBAL(1:NLEN),  VALUE = IZG, DIM = 1)
      IFOUND2 = FINDLOC (GF%ZONES_GLOBAL(1:GF%NC), VALUE = IZG, DIM = 1)
      IFOUND3 = FINDLOC (GC%LOCAL_TO_GLOBAL(1:GC%NCE2), VALUE = IZG, DIM = 1)
      IF (IFOUND1 <= 0 .AND. IFOUND2 <= 0) THEN  
         GC%CELLS_LOCAL(IP)  = IFOUND3
         GC%CELLS_GLOBAL(IP) = IZG
         IP = IP+1
      ENDIF
   ENDDO

ENDDO PROLONGATION_CELLS_LOOP
GC%NC_GALERKIN = IP-1


END SUBROUTINE SCARC_GET_CELL_DEPENDENCIES_GALERKIN


! --------------------------------------------------------------------------------------------------------------------
!> \brief Define nullspace for next coarser level, if coarsest level isn't reached yet
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_NULLSPACE_COARSE(NL)
USE SCARC_POINTERS, ONLY: GC, GF, ZF, SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN):: NL
INTEGER:: NM

CROUTINE = 'SCARC_SETUP_NULLSPACE_COARSE'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                   ! Sets pointers GC and GF
   ZF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_ZONES)

   IF (NL < NLEVEL_MAX) THEN
      GC%N_FINE = GC%NC_LOCAL(NM)
      CALL SCARC_ALLOCATE_REAL1 (GC%NULLSPACE, 1, GC%NCE, NSCARC_INIT_ZERO, 'GC%NULLSPACE', CROUTINE)
      GC%NULLSPACE(1:GC%NCE) = GF%RR(1:GC%NCE)
      CALL SCARC_REDUCE_INT1 (GC%LOCAL_TO_GLOBAL, GC%NCE, 'GC%LOCAL_TO_GLOBAL', CROUTINE)
   ENDIF


ENDDO
   
END SUBROUTINE SCARC_SETUP_NULLSPACE_COARSE


! --------------------------------------------------------------------------------------------------------------------
!> \brief Determine which columns of system matrix are involved in multiplication with tentative prolongator
! --------------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_MATCH_MATRIX_COLUMN(A, IC, JC)
TYPE (SCARC_CMATRIX_TYPE), POINTER:: A
INTEGER, INTENT(IN):: IC, JC
INTEGER:: ICOL
SCARC_MATCH_MATRIX_COLUMN = -1
DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
   IF (A%COL(ICOL) == JC) THEN
      SCARC_MATCH_MATRIX_COLUMN = ICOL
      RETURN 
   ENDIF
ENDDO
END FUNCTION SCARC_MATCH_MATRIX_COLUMN


! --------------------------------------------------------------------------------------------------------------------
!> \brief Setup Restriction matrix: Build transpose of Prolongation matrix
! Compute the Restriction matrix, R, which interpolates from the fine-grid to the coarse-grid.
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_RESTRICTION(NL)
USE SCARC_POINTERS, ONLY: GC, GF, RF, PF, SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN):: NL
INTEGER:: NM, IC, IRROW, IRCOL, IPCOL, IPC, ICCL, ICCG, IFOUND 
LOGICAL:: IS_INCLUDED

CROUTINE = 'SCARC_SETUP_RESTRICTION'

! Allocate Restriction matrix R on internal mesh part
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)            ! Sets pointers GF and GC to fine and coarse level

   PF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)
   RF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_RESTRICTION)

   RF%N_VAL = PF%N_VAL+100
   RF%N_ROW = PF%N_ROW+100

   CALL SCARC_ALLOCATE_CMATRIX (RF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%RESTRICTION', CROUTINE)

   IRROW = 1                             ! counter of current row of Restriction matrix-corresponds to coarse cells
   IRCOL = 1                             ! counter of current column of Restriction matrix
   RF%ROW(IRROW) = IRCOL

   LOCAL_COARSE_CELLS_LOOP: DO ICCL = 1, GC%NCE

      ICCG = GC%LOCAL_TO_GLOBAL(ICCL)                     ! corresponding global coarse cell

      IFOUND = -1
      IFOUND = FINDLOC(PF%COLG, VALUE = ICCG, DIM = 1)
      IF (IFOUND == -1) CYCLE

      IS_INCLUDED = .FALSE.

      FINE_CELLS_LOOP: DO IC = 1, GF%NCE                  ! counter of fine cell (including overlaps)

         ROW_LOOP: DO IPCOL = PF%ROW(IC), PF%ROW(IC+1)-1
            IPC = PF%COLG(IPCOL)
            IF (IPC == ICCG) THEN
               IS_INCLUDED = .TRUE.
               RF%VAL(IRCOL) = PF%VAL(IPCOL)
               RF%COLG(IRCOL) = IC
               IRCOL = IRCOL+1
               EXIT ROW_LOOP
            ENDIF
         ENDDO ROW_LOOP
      ENDDO FINE_CELLS_LOOP

      IF (IS_INCLUDED) THEN 
         RF%ROW(IRROW+1) = IRCOL
         IRROW = IRROW+1
      ENDIF

      RF%N_ROW = IRROW 
      RF%N_VAL = IRCOL-1

   ENDDO LOCAL_COARSE_CELLS_LOOP

   CALL SCARC_REDUCE_CMATRIX (RF, 'GF%RESTRICTION', CROUTINE)


ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_RESTRICTION


! --------------------------------------------------------------------------------------------------------------------
!> \brief Find matching column index during matrix-matrix multiplication of compact matrices
! --------------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_FIND_MATCHING_COLUMN(P, JC, ICCG)
TYPE (SCARC_CMATRIX_TYPE), POINTER:: P
INTEGER, INTENT(IN):: JC, ICCG
INTEGER:: IPCOL

SCARC_FIND_MATCHING_COLUMN = -1
IF (JC == 0) RETURN
DO IPCOL = P%ROW(JC), P%ROW(JC+1)-1
   IF (P%COLG(IPCOL) == ICCG) THEN
      SCARC_FIND_MATCHING_COLUMN = IPCOL
      RETURN 
   ENDIF
ENDDO

END FUNCTION SCARC_FIND_MATCHING_COLUMN


! --------------------------------------------------------------------------------------------------------------------
!> \brief Find matching components to multiply row of Poisson matrix with column of Prolongation matrix
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MULTIPLY_POISSON_PROL(A, P, PP, ICC, ICC0, IP, IC)
TYPE (SCARC_CMATRIX_TYPE), POINTER, INTENT(IN):: A, P, PP
INTEGER, INTENT(IN):: ICC, IC
INTEGER, INTENT(INOUT):: IP, ICC0
REAL(EB):: DSUM, TOL = 1E-12_EB
INTEGER:: IACOL, IPCOL, JC

DSUM = 0.0_EB
DO IACOL = A%ROW(IC), A%ROW(IC+1)-1
   JC = A%COL(IACOL)
   IF (JC == 0) CYCLE
   IPCOL = SCARC_FIND_MATCHING_COLUMN(P, JC, ICC)
   IF (JC < 0 .OR. IPCOL <= 0) CYCLE
   ICC0 = ICC
   DSUM = DSUM+A%VAL(IACOL) * P%VAL(IPCOL)
ENDDO

IF (ABS(DSUM) > TOL) THEN
   PP%COL(IP)  = ICC
   PP%COLG(IP) = ICC
   PP%VAL(IP)  = DSUM
   IP = IP+1
ENDIF

END SUBROUTINE SCARC_MULTIPLY_POISSON_PROL


! --------------------------------------------------------------------------------------------------------------------
!> \brief Perform matrix multiplication between fine Poisson matrix and Prolongation matrix 
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSON_PROL(NL)
USE SCARC_POINTERS, ONLY: GC, GF, AF, PF, PPF, OA, OPP, &
                          SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_MULTIGRID, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN):: NL
INTEGER  :: NM, NOM, IC, IP, IP0, INBR, ICC, ICC0 = -1
REAL(EB):: TNOW, TSUM = 0.0_EB

CROUTINE = 'SCARC_SETUP_POISSON_PROL'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)

   AF  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)          
   PF  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_PROLONGATION)     
   PPF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_PROL)         

   !PPF%N_ROW = AF%N_ROW
   PPF%N_ROW = GF%NCE+1
   PPF%N_VAL = PPF%N_ROW*30            ! TODO: only temporarily
   CALL SCARC_ALLOCATE_CMATRIX (PPF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%PPF', CROUTINE)

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = SCARC(NM)%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

      OA  => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_POISSON)
      OPP => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_POISSON_PROL)

      OPP%N_VAL = AF%N_VAL              ! TODO : CHECK : MUCH TOO BIG !!!
      OPP%N_ROW = GF%NCE+1            
      CALL SCARC_ALLOCATE_CMATRIX (OPP, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OGF%PP', CROUTINE)

   ENDDO

   IP = 1
   IP0 = IP
   PPF%ROW(1) = IP


   FINE_CELLS_LOOP: DO IC = 1, GF%NC

      ! TODO: Better time measurement!
      TNOW = CURRENT_TIME()
      IF (MY_RANK == 0 .AND. MOD(IC, 1000) == 0) WRITE(*,*) 'ScaRC-AMG-Setup: Processing cell ',IC, ' of ',GF%NC

      INTERNAL_COARSE_CELLS_LOOP: DO ICC = 1, GC%NC
         CALL SCARC_MULTIPLY_POISSON_PROL(AF, PF, PPF, GC%LOCAL_TO_GLOBAL(ICC), ICC0, IP, IC)
      ENDDO INTERNAL_COARSE_CELLS_LOOP

      EXTERNAL_COARSE_CELLS_LOOP: DO ICC = 1, GC%NC_GALERKIN
         CALL SCARC_MULTIPLY_POISSON_PROL(AF, PF, PPF, GC%CELLS_GLOBAL(ICC), ICC0, IP, IC)
      ENDDO EXTERNAL_COARSE_CELLS_LOOP

      ! take care that at least one entry per fine cell is generated
      IF (IP == IP0) THEN
         PPF%VAL(IP)  = 0.0_EB
         PPF%COL(IP)  = ICC0
         PPF%COLG(IP) = ICC0
         IP = IP+1
      ENDIF
      IP0 = IP

      PPF%ROW(IC+1) = IP

CPU(MY_RANK)%AMG = CPU(MY_RANK)%AMG+CURRENT_TIME()-TNOW
TSUM = TSUM+CPU(MY_RANK)%AMG

   ENDDO FINE_CELLS_LOOP


ENDDO

! Exchange overlapping parts of Prolongation matrix and extract exchanged data
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_POISSON_PROL, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLSG, NSCARC_MATRIX_POISSON_PROL, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_POISSON_PROL, NL)
   CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_POISSON_PROL, 0, NL)
ENDIF

! Reduce workspace for Prolongation matrix to really needed size
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)
   PPF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_PROL)         


   CALL SCARC_REDUCE_CMATRIX (PPF, 'GF%POISSON-PROL', CROUTINE)
   DO INBR = 1, SCARC(NM)%N_NEIGHBORS
      NOM = SCARC(NM)%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)
      OPP => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_POISSON_PROL)
      CALL SCARC_REDUCE_CMATRIX (OPP, 'OGF%POISSON_PROL', CROUTINE)
   ENDDO

ENDDO

END SUBROUTINE SCARC_SETUP_POISSON_PROL


! --------------------------------------------------------------------------------------------------------------------
!> \brief Find matching components to multiply row of Restriction matrix with column of Poisson-Prol matrix
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MULTIPLY_GALERKIN(PPF, RF, AC, ICCL, JCCL, JCCG, JCCG0, IP)
USE SCARC_POINTERS, ONLY: GF
TYPE (SCARC_CMATRIX_TYPE), POINTER, INTENT(IN):: PPF, RF, AC
INTEGER, INTENT(IN):: ICCL, JCCL, JCCG
INTEGER, INTENT(INOUT):: IP, JCCG0
INTEGER:: IAPCOL, IRCOL, JC
REAL(EB):: DSUM, TOL = 1E-12_EB

DSUM = 0.0_EB
DO IRCOL = RF%ROW(ICCL), RF%ROW(ICCL+1)-1
   JC = RF%COLG(IRCOL)
   IF (TYPE_COARSENING == NSCARC_COARSENING_GMG .AND. JC > GF%NC) RETURN
   IAPCOL = SCARC_FIND_MATCHING_COLUMN(PPF, JC, JCCG) 
   IF (IAPCOL > 0) THEN
      JCCG0 = JCCG
      DSUM = DSUM+RF%VAL(IRCOL) * PPF%VAL(IAPCOL)
   ENDIF
ENDDO

IF (ABS(DSUM) > TOL) THEN
   AC%COL(IP)  = JCCL
   AC%COLG(IP) = JCCG
   AC%VAL(IP) = DSUM
   IP = IP+1
ENDIF

END SUBROUTINE SCARC_MULTIPLY_GALERKIN


! --------------------------------------------------------------------------------------------------------------------
!> \brief Setup Galerkin matrix on coarser grid level (AMG only)
! Note: Matrix POISSON_PROL corresponds to POISSON x PROLONGATION  ~ AP
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GALERKIN(NL)
USE SCARC_POINTERS, ONLY: GF, GC, PPF, RF, AC, OAC, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                          SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_OTHER_MULTIGRID, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN):: NL
INTEGER  :: NM, IP, IP0, INBR, NOM, NLEN, ICCL, ICCG, JCC, JCCL, JCCG, JCCG0 = -1
REAL(EB):: TNOW, TSUM = 0.0_EB

CROUTINE = 'SCARC_SETUP_GALERKIN'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)

   PPF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_PROL)         
   RF  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_RESTRICTION)         
   AC  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)         


   IF (.NOT.ALLOCATED (AC%VAL)) THEN
      AC%N_ROW = GC%NCE+1
      AC%N_VAL = AC%N_ROW**2             ! only temporarily TODO TOO BIG
      CALL SCARC_ALLOCATE_CMATRIX (AC, NL+1, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GC%POISSON', CROUTINE)
   ENDIF

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = SCARC(NM)%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL+1)

      OAC => SCARC_POINT_TO_OTHER_CMATRIX (NSCARC_MATRIX_POISSON)

      OAC%N_VAL = AC%N_VAL              ! TODO : CHECK : MUCH TOO BIG !!!
      OAC%N_ROW = GC%NCE2+1           ! TODO : CHECK : MUCH TOO BIG !!!
      CALL SCARC_ALLOCATE_CMATRIX (OAC, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OAC%POISSON', CROUTINE)

   ENDDO

   CALL SCARC_DEALLOCATE_INT1 (GC%CELLS_LOCAL,  'GC%CELLS_LOCAL',  CROUTINE)
   CALL SCARC_DEALLOCATE_INT1 (GC%CELLS_GLOBAL, 'GC%CELLS_GLOBAL', CROUTINE)

   NLEN = 4 * (GC%NCE2-GC%NC+2)

   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_LOCAL, 1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_LOCAL', CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_GLOBAL, 1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_GET_CELL_DEPENDENCIES_GALERKIN(GC, GF, PPF, NLEN)

   CALL SCARC_REDUCE_INT1 (GC%CELLS_LOCAL, GC%NC_GALERKIN, 'GC%CELLS_LOCAL', CROUTINE)
   CALL SCARC_REDUCE_INT1 (GC%CELLS_GLOBAL, GC%NC_GALERKIN, 'GC%CELLS_GLOBAL', CROUTINE)

ENDDO

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)

   PPF => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_PROL)         
   RF  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_RESTRICTION)         
   AC  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)         

   IP = 1
   IP0 = IP
   AC%ROW(1) = IP

   LOCAL_COARSE_CELLS_LOOP: DO ICCL = 1, GC%NC

TNOW = CURRENT_TIME()
      ICCG = GC%LOCAL_TO_GLOBAL(ICCL)               ! corresponding global coarse cell

      INTERNAL_COARSE_CELLS_LOOP: DO JCCL = 1, GC%NC
         JCCG = GC%LOCAL_TO_GLOBAL(JCCL)
         CALL SCARC_MULTIPLY_GALERKIN(PPF, RF, AC, ICCL, JCCL, JCCG, JCCG0, IP)
      ENDDO INTERNAL_COARSE_CELLS_LOOP

      !EXTERNAL_COARSE_CELLS_LOOP: DO JCCL = GC%NC+1, GC%NCE2
      EXTERNAL_COARSE_CELLS_LOOP: DO JCC = 1, GC%NC_GALERKIN
         JCCL = GC%CELLS_LOCAL(JCC)
         JCCG = GC%CELLS_GLOBAL(JCC)
         CALL SCARC_MULTIPLY_GALERKIN(PPF, RF, AC, ICCL, JCCL, JCCG, JCCG0, IP)
      ENDDO EXTERNAL_COARSE_CELLS_LOOP

      ! take care that at least one entry per fine cell is generated
      IF (IP == IP0) THEN
         AC%VAL(IP)  = 0.0_EB
         AC%COL(IP)  = JCCG0
         AC%COLG(IP) = JCCG0
         IP = IP+1
      ENDIF
      IP0 = IP

      AC%ROW(ICCL+1) = IP

! TODO: better time measurement
CPU(MY_RANK)%AMG = CPU(MY_RANK)%AMG+CURRENT_TIME()-TNOW
TSUM = TSUM+CPU(MY_RANK)%AMG

   ENDDO LOCAL_COARSE_CELLS_LOOP
   AC%N_ROW = ICCL 


   CALL SCARC_REDUCE_CMATRIX (AC, 'POISSON-COARSE', CROUTINE)
   CALL SCARC_GET_MATRIX_STENCIL_MAX(AC, GC%NC)

ENDDO

CALL SCARC_RESORT_MATRIX_ROWS(NL+1)             

MESH_INT = 0                            
RANK_INT = 0

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL+1)                       ! TODO: Check if correct coarse level is used
   AC => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)         
   MESH_INT(NM) = AC%N_STENCIL_MAX
   RANK_INT = MAX(RANK_INT, MESH_INT(NM))
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_INT, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERROR)

#ifdef WITH_MKL
IF (TYPE_MKL(NL+1) == NSCARC_MKL_LOCAL .OR. TYPE_MKL(NL+1) == NSCARC_MKL_GLOBAL) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      CALL SCARC_SETUP_MATRIX_MKL(NSCARC_MATRIX_POISSON, NM, NL+1)
   ENDDO
ENDIF
#endif

END SUBROUTINE SCARC_SETUP_GALERKIN


! --------------------------------------------------------------------------------------------------------------------
!> \brief Compute entry of Poisson times Prolongation matrix at specified position
! This consists of a summation over the entries:    P(:,ICC)*A(IC, :) 
! Thus, it must be checked, if-for a given entry of A in row IC-the Prolongation matrix
! has a corresponding non-zero value
! --------------------------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION SCARC_VALUE_RAP(IC, ICC)
USE SCARC_POINTERS, ONLY : GC, AF, PF
INTEGER, INTENT(IN):: IC, ICC
INTEGER:: JC, IA, IP, JCC
REAL(EB):: DSUM

DSUM = 0.0_EB
DO IA = AF%ROW(IC), AF%ROW(IC+1) - 1
   JC = AF%COL(IA)
   IF (JC < 0) CYCLE
   DO IP = PF%ROW(JC), PF%ROW(JC+1) -1
      JCC = PF%COL(IP) 
      IF (JCC == GC%LOCAL_TO_GLOBAL(ICC)) THEN
         DSUM = DSUM+AF%VAL(IA)*PF%VAL(IP)
         CYCLE
      ENDIF
   ENDDO
ENDDO
SCARC_VALUE_RAP = DSUM


END FUNCTION SCARC_VALUE_RAP


! --------------------------------------------------------------------------------------------------------------------
!> \brief Resort matrix entries such that diagonal entry comes first (compact storage technique only)
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESORT_MATRIX_ROWS(NL)
USE SCARC_POINTERS, ONLY: G, A, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN):: NL
INTEGER, ALLOCATABLE, DIMENSION(:):: COL_AUX, COLG_AUX
REAL(EB), ALLOCATABLE, DIMENSION(:):: VAL_AUX
INTEGER:: NM, NCOL, ICOL, JCOL, KCOL, IC
LOGICAL:: COLG_IS_DEFINED = .FALSE.

CROUTINE = 'SCARC_RESORT_MATRIX_ROWS'

! TODO: use correct length of COL
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)
   A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)

   CALL SCARC_ALLOCATE_INT1  (COL_AUX,  1, A%N_STENCIL_MAX, NSCARC_INIT_ZERO, 'COL_AUX',  CROUTINE)
   CALL SCARC_ALLOCATE_INT1  (COLG_AUX, 1, A%N_STENCIL_MAX, NSCARC_INIT_ZERO, 'COLG_AUX', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (VAL_AUX,  1, A%N_STENCIL_MAX, NSCARC_INIT_ZERO, 'VAL_AUX',  CROUTINE)

   IF (ALLOCATED(A%COLG)) COLG_IS_DEFINED = .TRUE.

   !DO IC = 1, G%NCE2
   DO IC = 1, G%NC
      COL_AUX = 0
      JCOL = 1
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         COL_AUX(JCOL) = A%COL(ICOL)
         IF (COLG_IS_DEFINED)  COLG_AUX(JCOL) = A%COLG(ICOL)
         VAL_AUX(JCOL) = A%VAL(ICOL)
         JCOL = JCOL+1
      ENDDO
      NCOL = JCOL-1

      ! Find column index of diagonal element
      JCOL = 0
      DO WHILE (JCOL <= NCOL)
        JCOL = JCOL+1
        IF (COL_AUX(JCOL) == IC) EXIT
      ENDDO

      ! Store corresponding index and value in first matrix element of that row
      ICOL = A%ROW(IC)
      A%COL(ICOL)  = COL_AUX(JCOL)
      A%COLG(ICOL) = COLG_AUX(JCOL)
      A%VAL(ICOL)  = VAL_AUX(JCOL)

      COL_AUX(JCOL) = 99999999
      COLG_AUX(JCOL) = 99999999

      IF (COLG_IS_DEFINED) THEN
         JCOL = MINLOC(COLG_AUX(1:NCOL), DIM = 1)
      ELSE
         JCOL = MINLOC(COL_AUX(1:NCOL), DIM = 1)
      ENDIF
      KCOL = 1
      ICOL = ICOL+1
      DO WHILE (KCOL < NCOL)
         A%COLG(ICOL) = COLG_AUX(JCOL)
         A%COL(ICOL)  = COL_AUX(JCOL)
         A%VAL(ICOL)  = VAL_AUX(JCOL)
         IF (COLG_IS_DEFINED) THEN
            COLG_AUX(JCOL) = 99999999
            JCOL = MINLOC(COLG_AUX(1:NCOL), DIM = 1)
         ELSE
            COL_AUX(JCOL) = 99999999
            JCOL = MINLOC(COL_AUX(1:NCOL), DIM = 1)
         ENDIF
         KCOL = KCOL+1
         ICOL = ICOL+1
      ENDDO

      IF (ICOL /= A%ROW(IC+1)) WRITE(*,*) 'ERROR IN RESORT_MATRIX_ROWS'

   ENDDO


   CALL SCARC_DEALLOCATE_INT1  (COL_AUX,  'COL_AUX',  CROUTINE)
   CALL SCARC_DEALLOCATE_INT1  (COLG_AUX, 'COLG_AUX', CROUTINE)
   CALL SCARC_DEALLOCATE_REAL1 (VAL_AUX,  'VAL_AUX',  CROUTINE)

ENDDO

END SUBROUTINE SCARC_RESORT_MATRIX_ROWS

END MODULE SCARC_AMG



#endif
!=======================================================================================================================
!
! MODULE SCARC_MGM
! 
!> \brief Setup environment for call of McKeeney-Greengard-Mayo method (still experimental )
!
!=======================================================================================================================
MODULE SCARC_MGM

USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS
USE MESH_VARIABLES
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME, GET_FILE_NUMBER, SHUTDOWN
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_STORAGE
USE SCARC_CONVERGENCE

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate vectors and define variables needed for McKeeney-Greengard-Mayo method
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM (NL)
USE SCARC_POINTERS, ONLY: L, G, MGM, LO, UP, SCARC_POINT_TO_MGM, SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
USE SCARC_CONVERGENCE, ONLY: VELOCITY_ERROR_MGM, NIT_MGM
INTEGER, INTENT(IN):: NL
INTEGER:: NM

CROUTINE = 'SCARC_SETUP_MGM'
IS_MGM = .TRUE.

CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)

   ! Initialize MGM related variables
    
   MGM%NCS = L%STRUCTURED%NC
   MGM%NCU = L%UNSTRUCTURED%NC

   MGM%NWE = L%N_WALL_CELLS_EXT
   MGM%NWI = L%N_WALL_CELLS_INT
   MGM%NW1 = L%N_WALL_CELLS_EXT+1
   MGM%NW2 = L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT

   VELOCITY_ERROR_MGM = SCARC_MGM_ACCURACY
   NIT_MGM = SCARC_MGM_ITERATIONS

   ! Allocate workspace for the storage of the different vectors in the MGM methods

   ! SIP   : structured inhomogeneous Poisson solution (pass 1)
   ! UIP   : unstructured inhomogeneous Poisson solution (merge)
   ! UHL   : unstructured homogeneous Laplace solution (pass 2)
   ! UHL2  : unstructured homogeneous Laplace solution of previous time step (extrapolation BCs only)

   CALL SCARC_ALLOCATE_REAL3 (MGM%SIP, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%SIP', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%UIP, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%UIP', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%UHL, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%UHL', CROUTINE)

   ! OUIP  : other unstructured inhomogeneous Poisson solution on boundary
   ! OUHL  : other unstructured homogeneous Laplace solution on boundary
   ! OUHL2 : other unstructured homogeneous Laplace solution of previous time step on boundary

   CALL SCARC_ALLOCATE_REAL1 (MGM%OUIP, 1, MGM%NWE, NSCARC_INIT_ZERO, 'MGM%OUIP', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (MGM%OUHL, 1, MGM%NWE, NSCARC_INIT_ZERO, 'MGM%OUHL', CROUTINE)

   IF (TYPE_MGM_BC == NSCARC_MGM_BC_EXPOL) THEN
      CALL SCARC_ALLOCATE_REAL3 (MGM%UHL2, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%UHL2', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (MGM%OUHL2, 1, MGM%NWE, NSCARC_INIT_ZERO, 'MGM%OUHL2', CROUTINE)
   ENDIF

   ! BXS, BXF, BYS, BYF, BZS, BZF: boundary value vectors for the different faces

   CALL SCARC_ALLOCATE_REAL2 (MGM%BXS, 1, L%NY, 1, L%NZ, NSCARC_INIT_ZERO, 'MGM%BXS', CROUTINE)
   CALL SCARC_ALLOCATE_REAL2 (MGM%BXF, 1, L%NY, 1, L%NZ, NSCARC_INIT_ZERO, 'MGM%BXF', CROUTINE)
   CALL SCARC_ALLOCATE_REAL2 (MGM%BYS, 1, L%NX, 1, L%NZ, NSCARC_INIT_ZERO, 'MGM%BYS', CROUTINE)
   CALL SCARC_ALLOCATE_REAL2 (MGM%BYF, 1, L%NX, 1, L%NZ, NSCARC_INIT_ZERO, 'MGM%BYF', CROUTINE)
   CALL SCARC_ALLOCATE_REAL2 (MGM%BZS, 1, L%NX, 1, L%NY, NSCARC_INIT_ZERO, 'MGM%BZS', CROUTINE)
   CALL SCARC_ALLOCATE_REAL2 (MGM%BZF, 1, L%NX, 1, L%NY, NSCARC_INIT_ZERO, 'MGM%BZF', CROUTINE)

   ! UVEL,   VVEL,  WVEL : u-, v- and w-velocity components 
   ! OUVEL, OVVEL, OWVEL : u-, v- and w-velocity components of other mesh

   CALL SCARC_ALLOCATE_REAL3 (MGM%UVEL, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%UVEL', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%VVEL, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%VVEL', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%WVEL, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%WVEL', CROUTINE)

   CALL SCARC_ALLOCATE_REAL1 (MGM%OUVEL, 1, MGM%NWE, NSCARC_INIT_ZERO, 'MGM%OUVEL', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (MGM%OVVEL, 1, MGM%NWE, NSCARC_INIT_ZERO, 'MGM%OVVEL', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (MGM%OWVEL, 1, MGM%NWE, NSCARC_INIT_ZERO, 'MGM%OWVEL', CROUTINE)

   ! SCARC  : ScaRC solution (structured inhomogeneous by construction)
   ! USCARC : UScaRC solution (unstructured inhomogeneous by construction)
   ! DSCARC : difference of UScaRC and ScaRC solution 

   CALL SCARC_ALLOCATE_REAL3 (MGM%SCARC,  0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%SCARC',  CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%USCARC, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%USCARC', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%DSCARC, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%DSCARC', CROUTINE)

   ! UIP_VS_USCARC : difference vector of unstructured inhomogeneous Poisson versus UScaRC
   ! UHL_VS_DSCARC : difference vector of unstructured homogeneous Laplace versus difference UScaRC-ScaRC

   CALL SCARC_ALLOCATE_REAL3 (MGM%UIP_VS_USCARC, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%UIP_VS_US', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%UHL_VS_DSCARC, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%UHL_VS_DSCARC', CROUTINE)

   ! U1, V1, W1: u-, v- and w-velocity components in first MGM pass
   ! U2, V2, W2: u-, v- and w-velocity components in second MGM pass

   CALL SCARC_ALLOCATE_REAL3 (MGM%U1, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%U1', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%V1, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%V1', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%W1, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%W1', CROUTINE)
 
   CALL SCARC_ALLOCATE_REAL3 (MGM%U2, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%U2', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%V2, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%V2', CROUTINE)
   CALL SCARC_ALLOCATE_REAL3 (MGM%W2, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'MGM%W2', CROUTINE)

   ! Configure boundary cell counters and weights for 'True Approximate' boundary setting
   ! BTYPE  : Type of boundary condition in single boundary cells (Dirichlet/Neumann/Internal)
   ! WEIGHT : Weight for true approximate setting in single boundary cells

   IF (TYPE_MGM_BC == NSCARC_MGM_BC_TRUE) THEN

      CALL SCARC_ALLOCATE_INT2 (MGM%BTYPE, 1, MGM%NWE, -3, 3, NSCARC_INIT_NONE, 'MGM%BTYPE', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (MGM%WEIGHT, 1, MGM%NWE, NSCARC_INIT_ZERO, 'MGM%WEIGHT', CROUTINE)

      CALL SCARC_SETUP_MGM_TRUE_APPROXIMATE 

   ENDIF

   ! Allocate workspace for MGM solution, RHS and auxiliary vectors if not solver by CG

   IF (TYPE_MGM_LAPLACE /= NSCARC_MGM_LAPLACE_CG) THEN

      CALL SCARC_ALLOCATE_REAL1 (MGM%X, 1, G%NC, NSCARC_INIT_ZERO, 'X', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (MGM%B, 1, G%NC, NSCARC_INIT_ZERO, 'B', CROUTINE)

#ifdef WITH_MKL
      IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
         CALL SCARC_ALLOCATE_REAL1_FB (MGM%X_FB, 1, G%NC, NSCARC_INIT_ZERO, 'X', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1_FB (MGM%B_FB, 1, G%NC, NSCARC_INIT_ZERO, 'B', CROUTINE)
      ENDIF
#endif

   ENDIF

   ! The following code is still experimental and addresses the solution of the LU method compactly stored matrices

   IF (TYPE_MGM_LAPLACE == NSCARC_MGM_LAPLACE_LU .OR. TYPE_MGM_LAPLACE == NSCARC_MGM_LAPLACE_LUPERM) THEN

      LO => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LOWER)
      UP => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_UPPER)

      CALL SCARC_ALLOCATE_REAL1 (MGM%Y, 1, G%NC, NSCARC_INIT_ZERO, 'Y', CROUTINE)

      CALL SCARC_SETUP_MGM_PASS2_SIZES(NM, NLEVEL_MIN)             ! TODO
      CALL SCARC_ALLOCATE_CMATRIX (LO, NLEVEL_MIN, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_LIGHT, 'LO', CROUTINE)
      CALL SCARC_ALLOCATE_CMATRIX (UP, NLEVEL_MIN, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_LIGHT, 'UP', CROUTINE)
   
      CALL SCARC_SETUP_MGM_PASS2(NM, NLEVEL_MIN)

   ENDIF

ENDDO

END SUBROUTINE SCARC_SETUP_MGM


! -------------------------------------------------------------------------------------------------------------
!> \brief Setup structures for the true approximate boundary setting in MGM
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_TRUE_APPROXIMATE 
USE SCARC_POINTERS, ONLY: L, G, GWC, MGM
INTEGER:: IW, I, J, K, IOR0, IC 
INTEGER, ALLOCATABLE, DIMENSION(:,:):: CNT
REAL(EB):: SX, SY, SZ

! temporarily allocate counter

CALL SCARC_ALLOCATE_INT2 (CNT, 1, G%NC, -3, 3, NSCARC_INIT_NONE, 'CNT', CROUTINE)
CNT = 2;  CNT(:,0) = 0

DO IW = 1, MGM%NWE

   GWC => G%WALL(IW)

   I = GWC%IXW
   J = GWC%IYW
   K = GWC%IZW

   IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

   IOR0 = GWC%IOR
   IC = G%CELL_NUMBER(I, J, K)
 
   ! Counter for the main diagonal entries to be considered

   IF (GWC%BTYPE == DIRICHLET) THEN
      CNT(IC, IOR0) = CNT(IC, IOR0) + 1 
   ELSE IF (GWC%BTYPE == NEUMANN) THEN
      CNT(IC, IOR0) = CNT(IC, IOR0) - 1 
   ENDIF

ENDDO

! Configure boundary type information for 'True Approximate' boundary setting

DO IW = 1, MGM%NWE

   GWC => G%WALL(IW)
   IF (GWC%BTYPE /= INTERNAL) CYCLE

   I = GWC%IXW
   J = GWC%IYW
   K = GWC%IZW

   IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
   IC = G%CELL_NUMBER(I, J, K)

   DO IOR0 = -3, 3
      IF (CNT(IC, IOR0) > 2) THEN
         MGM%BTYPE(IW, IOR0) = DIRICHLET
      ELSE IF (CNT(IC, IOR0) < 2) THEN
         MGM%BTYPE(IW, IOR0) = NEUMANN
      ELSE
         MGM%BTYPE(IW, IOR0) = INTERNAL
      ENDIF
   ENDDO

   IF (TWO_D) THEN
      SX = REAL(CNT(IC, 1) + CNT(IC, -1) - 2, EB)
      SZ = REAL(CNT(IC, 3) + CNT(IC, -3) - 2, EB)
      MGM%WEIGHT(IW) = 1.0_EB/(SX*L%DXI2+SZ*L%DZI2)
   ELSE
      SX = REAL(CNT(IC, 1) + CNT(IC, -1) - 2, EB)
      SY = REAL(CNT(IC, 2) + CNT(IC, -2) - 2, EB)
      SZ = REAL(CNT(IC, 3) + CNT(IC, -3) - 2, EB)
      MGM%WEIGHT(IW) = 1.0_EB/(SX*L%DXI2+SY*L%DYI2+SZ*L%DZI2)
   ENDIF
   
ENDDO

CALL SCARC_DEALLOCATE_INT2 (CNT, 'CNT', CROUTINE)

END SUBROUTINE SCARC_SETUP_MGM_TRUE_APPROXIMATE


! -------------------------------------------------------------------------------------------------------------
!> \brief Setup sizes for LU-decomposition of McKeeney-Greengard-Mayo method
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_PASS2_SIZES(NM, NL)
USE SCARC_POINTERS, ONLY: G, LO, UP, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN):: NM, NL
INTEGER:: IC, JC, NMAX_U, NMAX_L

CROUTINE = 'SCARC_SETUP_MGM_PASS2_SIZES'

CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
CALL SCARC_POINT_TO_GRID (NM, NL)

! Preset pointers for LO and UP with one-value rows (corresponding to initialization with diagonal element)
 
NMAX_U = G%NC
NMAX_L = G%NC

UP%N_ROW = G%NC+1
LO%N_ROW = G%NC+1
UP%N_VAL = G%NC
LO%N_VAL = G%NC

ROW_LOOP: DO IC = 1, G%NC  
   COL_LOOP: DO JC = IC, G%NC
      UP%N_VAL = UP%N_VAL + 1
      LO%N_VAL = LO%N_VAL + 1
   ENDDO COL_LOOP
ENDDO ROW_LOOP

END SUBROUTINE SCARC_SETUP_MGM_PASS2_SIZES


! -------------------------------------------------------------------------------------------------------------
!> \brief Setup LU-decomposition for McKeeney-Greengard-Mayo method
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_PASS2(NM, NL)
USE SCARC_POINTERS, ONLY: G, A, LO, UP, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN):: NM, NL
INTEGER:: IC0, IC, JC, KC, NMAX_U, NMAX_L
REAL (EB):: SCAL, VL = 0.0_EB, VU = 0.0_EB, VAL
INTEGER:: TYPE_SCOPE_SAVE

CROUTINE = 'SCARC_SETUP_MGM_PASS2'

TYPE_SCOPE_SAVE = TYPE_SCOPE(0)
TYPE_SCOPE(0) = NSCARC_SCOPE_LOCAL

CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
CALL SCARC_POINT_TO_GRID (NM, NL)

A  => G%LAPLACE

! Preset pointers for LM and UM with one-value rows (corresponding to initialization with diagonal element)
!
DO IC = 1, G%NC
   UP%ROW(IC) = IC ;  UP%COL(IC) = IC
   LO%ROW(IC) = IC ;  LO%COL(IC) = IC
ENDDO
UP%ROW(G%NC+1) = G%NC+1
LO%ROW(G%NC+1) = G%NC+1

NMAX_U = G%NC
NMAX_L = G%NC

ROW_LOOP: DO IC0 = 1, G%NC  

   IC = IC0

   ! Set main diagonal element of L to 1.0
   VAL = 1.0_EB
   CALL SCARC_INSERT_TO_CMATRIX (LO, VAL, IC, IC, G%NC, NMAX_L, 'LO')

   COL_LOOP: DO JC = IC, G%NC

      SCAL = 0.0_EB
      DO KC = 1, IC-1
         VL = SCARC_EVALUATE_CMATRIX (LO, IC, KC)
         VU = SCARC_EVALUATE_CMATRIX (UP, KC, JC)
         SCAL = SCAL+VL*VU
      ENDDO

      VAL = SCARC_EVALUATE_CMATRIX(A, IC, JC)  - SCAL
      IF (ABS(VAL) > TWO_EPSILON_EB) CALL SCARC_INSERT_TO_CMATRIX (UP, VAL, IC, JC, G%NC, NMAX_U, 'UM')

      SCAL = 0.0_EB
      DO KC = 1, IC-1
         VL = SCARC_EVALUATE_CMATRIX (LO, JC, KC)
         VU = SCARC_EVALUATE_CMATRIX (UP, KC, IC)
         SCAL = SCAL+VL*VU
      ENDDO
      VAL = (SCARC_EVALUATE_CMATRIX(A, JC, IC) - SCAL)/SCARC_EVALUATE_CMATRIX(UP, IC, IC)
      IF (ABS(VAL) > TWO_EPSILON_EB) CALL SCARC_INSERT_TO_CMATRIX (LO, VAL, JC, IC, G%NC, NMAX_L, 'LM')

   ENDDO COL_LOOP

ENDDO ROW_LOOP

CALL SCARC_REDUCE_CMATRIX (LO, 'LO', CROUTINE)
CALL SCARC_REDUCE_CMATRIX (UP, 'UP', CROUTINE)

TYPE_SCOPE(0) = TYPE_SCOPE_SAVE

END SUBROUTINE SCARC_SETUP_MGM_PASS2


! --------------------------------------------------------------------------------------------------------------------
!> \brief Convergence state of MGM method
! --------------------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_MGM_CONVERGENCE_STATE(ITE_MGM, NTYPE)
USE SCARC_POINTERS, ONLY: MGM, SCARC_POINT_TO_MGM
INTEGER, INTENT(IN):: ITE_MGM, NTYPE
INTEGER:: NM

SCARC_MGM_CONVERGENCE_STATE = NSCARC_MGM_FAILURE
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM(NM, NLEVEL_MIN)

   SCARC_MGM_ACCURACY   = VELOCITY_ERROR_GLOBAL  ! Store achieved MGM accuracy for statistics in chid.out
   SCARC_MGM_ITERATIONS = ITE_MGM                ! Store required MGM iterations for statistics in chid.out

   SELECT CASE (NTYPE)

      ! Initialization - after first structured inhomogeneous Poisson solution

      CASE (0)

         MGM%ITE = 0
         MGM%ITE_LAPLACE = 0
         MGM%ITE_POISSON = ITE                   ! ITE, CAPPA contain statistics of preceding structured CG-solution
         MGM%CAPPA_POISSON = CAPPA

      ! MGM iteration - after each unstructured homogeneous Laplace solution

      CASE (1)

         MGM%ITE = ITE_MGM
         IF (TYPE_MGM_LAPLACE == NSCARC_MGM_LAPLACE_CG .AND. ITE > MGM%ITE_LAPLACE) THEN
            MGM%ITE_LAPLACE = MAX(ITE, MGM%ITE_LAPLACE)            ! Store worst Laplace-CG statistics
            MGM%CAPPA_LAPLACE = CAPPA
         ENDIF                         

      ! Termination - after whole MGM solution

      CASE (-1)

         CAPPA = MGM%CAPPA_POISSON   ! Reset to Krylov statistics of Poisson solution for statistics in chid.out
         ITE   = MGM%ITE_POISSON


   END SELECT
   IF (VELOCITY_ERROR_GLOBAL <= VELOCITY_ERROR_MGM) SCARC_MGM_CONVERGENCE_STATE = NSCARC_MGM_SUCCESS

ENDDO

END FUNCTION SCARC_MGM_CONVERGENCE_STATE


! --------------------------------------------------------------------------------------------------------------
!> \brief Set correct boundary values at external and internal boundaries
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_UPDATE_GHOSTCELLS(NTYPE)
USE SCARC_POINTERS, ONLY: M, L, G, GWC, HP, MGM, BXS, BXF, BYS, BYF, BZS, BZF, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN):: NTYPE
INTEGER:: NM, IW, IOR0, IXG, IYG, IZG, IXW, IYW, IZW 

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)                                   

   SELECT CASE(NTYPE)

      ! Update ghostcells for local Laplace problems
      ! Along external boundaries use zero Dirichlet or Neumann BC's
      ! Along mesh interfaces use Dirichlet BC's corresponding to MGM interface settings 
        
      CASE (NSCARC_MGM_LAPLACE)

         HP => MGM%UHL
         BXS => MGM%BXS ; BXF => MGM%BXF
         BYS => MGM%BYS ; BYF => MGM%BYF
         BZS => MGM%BZS ; BZF => MGM%BZF
    
         WALL_CELLS_LOOP_LAPLACE: DO IW = 1, L%N_WALL_CELLS_EXT
      
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
                  IF (GWC%BTYPE == INTERNAL) THEN
                     HP(IXG, IYW, IZW) = -HP(IXW, IYW, IZW) + 2.0_EB*BXS(IYW, IZW)
                  ELSE IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXG, IYW, IZW) = -HP(IXW, IYW, IZW) 
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXG, IYW, IZW) =  HP(IXW, IYW, IZW) 
                  ENDIF
               CASE (-1)
                  IF (GWC%BTYPE == INTERNAL) THEN
                     HP(IXG, IYW, IZW) = -HP(IXW, IYW, IZW) + 2.0_EB*BXF(IYW, IZW)
                  ELSE IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXG, IYW, IZW) = -HP(IXW, IYW, IZW) 
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXG, IYW, IZW) =  HP(IXW, IYW, IZW) 
                  ENDIF
               CASE ( 2)
                  IF (GWC%BTYPE == INTERNAL) THEN
                     HP(IXW, IYG, IZW) = -HP(IXW, IYW, IZW) + 2.0_EB*BYS(IXW, IZW)
                  ELSE IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXW, IYG, IZW) = -HP(IXW, IYW, IZW) 
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXW, IYG, IZW) =  HP(IXW, IYW, IZW) 
                  ENDIF
               CASE (-2)
                  IF (GWC%BTYPE == INTERNAL) THEN
                     HP(IXW, IYG, IZW) = -HP(IXW, IYW, IZW) + 2.0_EB*BYF(IXW, IZW)
                  ELSE IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXW, IYG, IZW) = -HP(IXW, IYW, IZW) 
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXW, IYG, IZW) =  HP(IXW, IYW, IZW) 
                  ENDIF
               CASE ( 3)
                  IF (GWC%BTYPE == INTERNAL) THEN
                     HP(IXW, IYW, IZG) = -HP(IXW, IYW, IZW) + 2.0_EB*BZS(IXW, IYW)
                  ELSE IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXW, IYW, IZG) = -HP(IXW, IYW, IZW) 
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXW, IYW, IZG) =  HP(IXW, IYW, IZW) 
                  ENDIF
               CASE (-3)
                  IF (GWC%BTYPE == INTERNAL) THEN
                     HP(IXW, IYW, IZG) = -HP(IXW, IYW, IZW) + 2.0_EB*BZF(IXW, IYW)
                  ELSE IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXW, IYW, IZG) = -HP(IXW, IYW, IZW) 
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXW, IYW, IZG) =  HP(IXW, IYW, IZW) 
                  ENDIF
            END SELECT
      
         ENDDO WALL_CELLS_LOOP_LAPLACE

      ! All other cases: 
      ! Updating the ghost cells for the inhomogeneous structured Poisson as well as 
      ! the (optional) ScaRC and UScaRC solutions in MGM and when terminating the current MGM run

      CASE  (NSCARC_MGM_POISSON, NSCARC_MGM_SCARC, NSCARC_MGM_USCARC, NSCARC_MGM_TERMINATE) 

         SELECT CASE (NTYPE)
            CASE (NSCARC_MGM_POISSON)
               HP => MGM%SIP
            CASE (NSCARC_MGM_SCARC) 
               HP => MGM%SCARC
            CASE (NSCARC_MGM_USCARC)
               HP => MGM%USCARC
            CASE (NSCARC_MGM_TERMINATE)
               HP => MGM%UIP
         END SELECT

         BXS => M%BXS ; BXF => M%BXF
         BYS => M%BYS ; BYF => M%BYF
         BZS => M%BZS ; BZF => M%BZF
       
         WALL_CELLS_LOOP: DO IW = 1, L%N_WALL_CELLS_EXT
      
            GWC => G%WALL(IW)
      
            IF (GWC%BTYPE == INTERNAL) CYCLE

            IXG = GWC%IXG
            IYG = GWC%IYG
            IZG = GWC%IZG
      
            IXW = GWC%IXW
            IYW = GWC%IYW
            IZW = GWC%IZW
      
            IOR0 = GWC%IOR
      
            SELECT CASE (IOR0)
               CASE ( 1)
                  IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXG, IYW, IZW) = -HP(IXW, IYW, IZW) + 2.0_EB*BXS(IYW, IZW)
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXG, IYW, IZW) =  HP(IXW, IYW, IZW) - L%DX*BXS(IYW, IZW)
                  ENDIF
               CASE (-1)
                  IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXG, IYW, IZW) = -HP(IXW, IYW, IZW) + 2.0_EB*BXF(IYW, IZW)
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXG, IYW, IZW) =  HP(IXW, IYW, IZW) + L%DX*BXF(IYW, IZW)
                  ENDIF
               CASE ( 2)
                  IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXW, IYG, IZW) = -HP(IXW, IYW, IZW) + 2.0_EB*BYS(IXW, IZW)
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXW, IYG, IZW) =  HP(IXW, IYW, IZW) - L%DY*BYS(IXW, IZW)
                  ENDIF
               CASE (-2)
                  IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXW, IYG, IZW) = -HP(IXW, IYW, IZW) + 2.0_EB*BYF(IXW, IZW)
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXW, IYG, IZW) =  HP(IXW, IYW, IZW) + L%DY*BYF(IXW, IZW)
                  ENDIF
               CASE ( 3)
                  IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXW, IYW, IZG) = -HP(IXW, IYW, IZW) + 2.0_EB*BZS(IXW, IYW)
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXW, IYW, IZG) =  HP(IXW, IYW, IZW) - L%DZ*BZS(IXW, IYW)
                  ENDIF
               CASE (-3)
                  IF (GWC%BTYPE == DIRICHLET) THEN
                     HP(IXW, IYW, IZG) = -HP(IXW, IYW, IZW) + 2.0_EB*BZF(IXW, IYW)
                  ELSE IF (GWC%BTYPE == NEUMANN) THEN
                     HP(IXW, IYW, IZG) =  HP(IXW, IYW, IZW) + L%DZ*BZF(IXW, IYW)
                  ENDIF
            END SELECT
      
         ENDDO WALL_CELLS_LOOP

   END SELECT

ENDDO

END SUBROUTINE SCARC_MGM_UPDATE_GHOSTCELLS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Copy specified vectors in McKeeney-Greengard-Mayo method 
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_COPY(NTYPE)
USE SCARC_POINTERS, ONLY: L, MGM, SCARC_POINT_TO_MGM
INTEGER, INTENT(IN):: NTYPE
INTEGER:: NM
INTEGER:: IX, IY, IZ

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM(NM, NLEVEL_MIN)

   SELECT CASE(NTYPE)

      ! Copy structured ihomogeneous Poisson MGM solution to unstructured ihomogeneous Poisson MGM solution (initialization) 

      CASE (NSCARC_MGM_SIP_TO_UIP) 

         MGM%UIP = MGM%SIP

      ! Copy structured inhomogeneous ScaRC solution to structured inhomogeneous Poisson MGM solution

      CASE (NSCARC_MGM_SCARC_TO_SIP)                    

         MGM%SIP = MGM%SCARC

      ! Copy unstructured inhomogeneous ScaRC solution to unstructured inhomogeneous Poisson MGM solution

      CASE (NSCARC_MGM_USCARC_TO_UIP)                     

         MGM%UIP = MGM%USCARC
         IF (TWO_D) THEN
            DO IZ = 1, L%NZ
               DO IX = 1, L%NX
                  IF (L%IS_SOLID(IX, 1, IZ)) MGM%UIP(IX, 0:2, IZ) = 0.0_EB
               ENDDO
            ENDDO
         ELSE
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  DO IX = 1, L%NX
                     IF (L%IS_SOLID(IX, IY, IZ)) MGM%UIP(IX, IY, IZ) = 0.0_EB
                  ENDDO
               ENDDO
            ENDDO
         ENDIF


      ! Copy difference between ScaRC and UScaRC to unstructured homogeneous Laplace MGM solution

      CASE (NSCARC_MGM_DSCARC_TO_UHL)    

         MGM%UHL = MGM%DSCARC

         IF (TWO_D) THEN
            DO IZ = 1, L%NZ
               DO IX = 1, L%NX
                  IF (L%IS_SOLID(IX, 1, IZ)) MGM%UHL(IX, 0:2, IZ) = 0.0_EB
               ENDDO
            ENDDO
         ELSE
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  DO IX = 1, L%NX
                     IF (L%IS_SOLID(IX, IY, IZ)) MGM%UHL(IX, IY, IZ) = 0.0_EB
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      ! Copy difference between ScaRC and UScaRC to previous unstructured homogeneous Laplace solution UHL2

      CASE (NSCARC_MGM_DSCARC_TO_UHL2)  

         MGM%UHL2 = MGM%DSCARC

      ! Copy unstructured homogeneous Laplace MGM solution UHL to previous UHL2

      CASE (NSCARC_MGM_UHL_TO_UHL2)  

         MGM%UHL2 = MGM%UHL

      ! Copy other unstructured homogeneous Laplace MGM solution to previous 

      CASE (NSCARC_MGM_OUHL_TO_OUHL2) 

         MGM%OUHL2 = MGM%OUHL


   END SELECT

ENDDO

END SUBROUTINE SCARC_MGM_COPY


! ------------------------------------------------------------------------------------------------------------------
!> \brief Build difference of specified vectors in McKeeney-Greengard-Mayo method
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_DIFF(NTYPE)
USE SCARC_POINTERS, ONLY: L, G, ST, MGM, GWC, SCARC_POINT_TO_MGM
INTEGER, INTENT(IN):: NTYPE
INTEGER:: NM, IX, IY, IZ, IOR0, IW

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM(NM, NLEVEL_MIN)
   ST  => L%STAGE(NSCARC_STAGE_ONE)

   SELECT CASE(NTYPE)

      ! Build difference between ScaRC and UScaRC solution

      CASE (NSCARC_MGM_USCARC_VS_SCARC)

         G  => L%STRUCTURED
         MGM%DSCARC = MGM%USCARC-MGM%SCARC

         IF (TWO_D) THEN
            !DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT
            DO IW = 1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT
               GWC => G%WALL(IW)
               IOR0 = GWC%IOR
               IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE
               IX = GWC%IXG;  IZ = GWC%IZG
               IF (L%IS_SOLID(IX, 1, IZ)) MGM%DSCARC(IX, 0:2, IZ) = 0.0_EB
            ENDDO
         ELSE
            DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT
               GWC => G%WALL(IW)
               IOR0 = GWC%IOR
               IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE
               IX = GWC%IXG;  IY = GWC%IYG;  IZ = GWC%IZG
               IF (L%IS_SOLID(IX, IY, IZ)) MGM%DSCARC(IX, IY, IZ) = 0.0_EB
            ENDDO
         ENDIF

      ! Build difference of unstructured Laplace MGM solution and difference of ScaRC and UScaRC

      CASE (NSCARC_MGM_UHL_VS_DSCARC)            

         MGM%UHL_VS_DSCARC = MGM%UHL-MGM%DSCARC

      ! Build difference of unstructured Poisson MGM solution and UScaRC solution

      CASE (NSCARC_MGM_UIP_VS_USCARC)                     

         MGM%UIP_VS_USCARC = MGM%UIP-MGM%USCARC


   END SELECT

ENDDO

END SUBROUTINE SCARC_MGM_DIFF


! ------------------------------------------------------------------------------------------------------------------
!> \brief Store specified type of vector in McKeeney-Greengard-Mayo method
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_STORE(NTYPE)
USE SCARC_POINTERS, ONLY: L, G, ST, MGM, GWC, M, SCARC_POINT_TO_MGM
INTEGER, INTENT(IN):: NTYPE
INTEGER:: NM, IX, IY, IZ, ICS, ICU, ICE, IOR0, IW

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM(NM, NLEVEL_MIN)
   ST  => L%STAGE(NSCARC_STAGE_ONE)

   SELECT CASE(NTYPE)

      ! ---------- Store ScaRC solution MGM%SCARC

      CASE (NSCARC_MGM_SCARC)

         G  => L%STRUCTURED

         DO IZ = 1, L%NZ
            DO IY = 1, L%NY
               DO IX = 1, L%NX
                  ICS = G%CELL_NUMBER(IX, IY, IZ)  
                  MGM%SCARC(IX, IY, IZ) = ST%X(ICS) 
               ENDDO
            ENDDO
         ENDDO

         DO IW = 1, L%N_WALL_CELLS_EXT
            GWC => G%WALL(IW)
            IF (GWC%BTYPE /= INTERNAL) CYCLE
            IOR0 = GWC%IOR
            IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE
            IX = GWC%IXG;  IY = GWC%IYG;  IZ = GWC%IZG
            ICE = L%STRUCTURED%CELL_NUMBER(IX, IY, IZ)
            MGM%SCARC(IX, IY, IZ) = ST%X(ICE) 
         ENDDO

      ! ---------- Store UScaRC solution MGM%USCARC

      CASE (NSCARC_MGM_USCARC)

         G  => L%UNSTRUCTURED

         MGM%USCARC = 0.0_EB
         DO IZ = 1, L%NZ
            DO IY = 1, L%NY
               DO IX = 1, L%NX
                  IF (L%IS_SOLID(IX, IY, IZ)) CYCLE
                  ICS = G%CELL_NUMBER(IX, IY, IZ)           ! unstructured cell number
                  MGM%USCARC(IX, IY, IZ) = ST%X(ICS) 
               ENDDO
            ENDDO
         ENDDO

         DO IW = 1, L%N_WALL_CELLS_EXT
            GWC => G%WALL(IW)
            IF (GWC%BTYPE /= INTERNAL) CYCLE
            IOR0 = GWC%IOR
            IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE
            IX = GWC%IXG;  IY = GWC%IYG;  IZ = GWC%IZG
            ICE = G%CELL_NUMBER(IX, IY, IZ)
            MGM%USCARC(IX, IY, IZ) = ST%X(ICE) 
         ENDDO


      ! ---------- Store structured inhomogeneous Poisson solution MGM%SIP

      CASE (NSCARC_MGM_POISSON)

         G  => L%STRUCTURED

         !MGM%SIP = 0.0_EB
         DO IZ = 1, L%NZ
            DO IY = 1, L%NY
               DO IX = 1, L%NX
                  ICS = G%CELL_NUMBER(IX, IY, IZ)              ! structured cell number
                  MGM%SIP(IX, IY, IZ) = ST%X(ICS) 
               ENDDO
            ENDDO
         ENDDO

         DO IW = 1, L%N_WALL_CELLS_EXT
            GWC => G%WALL(IW)
            IF (GWC%BTYPE /= INTERNAL) CYCLE
            IOR0 = GWC%IOR
            IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE
            IX = GWC%IXG;  IY = GWC%IYG;  IZ = GWC%IZG
            ICE = L%STRUCTURED%CELL_NUMBER(IX, IY, IZ)
            MGM%SIP(IX, IY, IZ) = ST%X(ICE) 
         ENDDO

      ! ---------- Store homogeneous unstructured Laplace solution MGM%UHL

      CASE (NSCARC_MGM_LAPLACE)

         G => L%UNSTRUCTURED

         IF (TYPE_MGM_BC == NSCARC_MGM_BC_EXPOL) THEN
            DO IZ = 0, L%NZ+1
               DO IY = 0, L%NY+1
                  DO IX = 0, L%NX+1
                     IF (L%IS_SOLID(IX, IY, IZ)) CYCLE
                     ICU = G%CELL_NUMBER(IX, IY, IZ)                  ! unstructured cell number
                     MGM%UHL2(IX, IY, IZ) = MGM%UHL(IX, IY, IZ)       ! also store second level
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         MGM%UHL = 0.0_EB
         SELECT CASE (TYPE_MGM_LAPLACE)
            CASE (NSCARC_MGM_LAPLACE_CG)
               DO IZ = 1, L%NZ
                  DO IY = 1, L%NY
                     DO IX = 1, L%NX
                        IF (L%IS_SOLID(IX, IY, IZ)) CYCLE
                        ICU = G%CELL_NUMBER(IX, IY, IZ)               ! unstructured cell number
                        MGM%UHL(IX, IY, IZ) = ST%X(ICU)               ! solution contained in ST%X
                     ENDDO
                  ENDDO
               ENDDO
            CASE (NSCARC_MGM_LAPLACE_LUPERM)
               DO IZ = 1, L%NZ
                  DO IY = 1, L%NY
                     DO IX = 1, L%NX
                        IF (L%IS_SOLID(IX, IY, IZ)) CYCLE
                        ICU = G%PERM_FW(G%CELL_NUMBER(IX, IY, IZ))    ! unstructured permuted cell number
                        MGM%UHL(IX, IY, IZ) = MGM%X(ICU)              ! solution contained in MGM%X
                     ENDDO
                  ENDDO
               ENDDO
            CASE DEFAULT
               DO IZ = 1, L%NZ
                  DO IY = 1, L%NY
                     DO IX = 1, L%NX
                        IF (L%IS_SOLID(IX, IY, IZ)) CYCLE
                        ICU = G%CELL_NUMBER(IX, IY, IZ)                ! unstructured cell number
                        MGM%UHL(IX, IY, IZ) = MGM%X(ICU)               ! solution contained in MGM%X
                     ENDDO
                  ENDDO
               ENDDO
          END SELECT


      ! ---------- Merge structured inhomogeneous Poisson and unstructured homogeneous Laplace solutions

      CASE (NSCARC_MGM_MERGE)


         DO IZ = 0, L%NZ+1
            DO IY = 0, L%NY+1
               DO IX = 0, L%NX+1
                  MGM%UIP(IX, IY, IZ) = MGM%SIP(IX, IY, IZ) + MGM%UHL(IX, IY, IZ)              ! Variant A
               ENDDO
            ENDDO
         ENDDO




      ! ---------- Terminate MGM method and extract predictor/corrector solution for FDS code

      CASE (NSCARC_MGM_TERMINATE)

         IF (PREDICTOR) THEN

            DO IZ = 0, L%NZ+1
               DO IY = 0, L%NY+1
                  DO IX = 0, L%NX+1
                     M%H(IX, IY, IZ) = MGM%UIP(IX, IY, IZ) 
                  ENDDO
               ENDDO
            ENDDO
            IF (TWO_D) THEN
               DO IZ = 1, L%NZ
                  DO IX = 1, L%NX
                     IF (L%IS_SOLID(IX, 1, IZ)) M%H(IX, 0:2, IZ) = 0.0_EB
                  ENDDO
               ENDDO
            ELSE
               DO IZ = 1, L%NZ
                  DO IY = 1, L%NY
                     DO IX = 1, L%NX
                        IF (L%IS_SOLID(IX, IY, IZ)) M%H(IX, IY, IZ) = 0.0_EB
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

         ELSE

            DO IZ = 0, L%NZ+1
               DO IY = 0, L%NY+1
                  DO IX = 0, L%NX+1
                     M%HS(IX, IY, IZ) = MGM%UIP(IX, IY, IZ) 
                  ENDDO
               ENDDO
            ENDDO
            IF (TWO_D) THEN
               DO IZ = 1, L%NZ
                  DO IX = 1, L%NX
                     IF (L%IS_SOLID(IX, 1, IZ)) M%HS(IX, 0:2, IZ) = 0.0_EB
                  ENDDO
               ENDDO
            ELSE
               DO IZ = 1, L%NZ
                  DO IY = 1, L%NY
                     DO IX = 1, L%NX
                        IF (L%IS_SOLID(IX, IY, IZ)) M%HS(IX, IY, IZ) = 0.0_EB
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF

   END SELECT

ENDDO

!#ifdef WITH_SCARC_DEBUG
!1000 FORMAT (A, ': IX, IY, IZ =', 3I6, ': ICU =', I6, ': HP =', E14.6)
!#endif
END SUBROUTINE SCARC_MGM_STORE


! --------------------------------------------------------------------------------------------------------------
!> \brief Setup workspace for McKeeney-Greengard-Mayo method
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_WORKSPACE(NL)
USE SCARC_POINTERS, ONLY: M, MGM, SCARC_POINT_TO_MGM
INTEGER, INTENT(IN):: NL
INTEGER  :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM(NM, NL)

   IF (PREDICTOR) THEN
      MGM%U1 = M%U
      MGM%V1 = M%V
      MGM%W1 = M%W
   ELSE
      MGM%U1 = M%US
      MGM%V1 = M%VS
      MGM%W1 = M%WS
   ENDIF

   MGM%SIP = 0.0_EB
   !MGM%UHL = 0.0_EB
   MGM%UIP = 0.0_EB
   !IF (TYPE_MGM_BC == NSCARC_MGM_BC_EXPOL) MGM%UHL2 = 0.0_EB
   MGM%UIP_VS_USCARC = 0.0_EB
   MGM%UHL_VS_DSCARC = 0.0_EB

   MGM%SCARC = 0.0_EB
   MGM%USCARC = 0.0_EB
   MGM%DSCARC = 0.0_EB

ENDDO

END SUBROUTINE SCARC_SETUP_MGM_WORKSPACE


! --------------------------------------------------------------------------------------------------------------
!> \brief Set interface boundary conditions for unstructured, homogeneous part of McKeeney-Greengard-Mayo method
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_SET_INTERFACES(VB, NM)
USE SCARC_POINTERS, ONLY: L, F, G, OL, OG, MGM, UHL, UHL2, OUHL, OUHL2, BTYPE, &
                          SCARC_POINT_TO_MGM, SCARC_POINT_TO_OTHER_GRID
INTEGER, INTENT(IN):: NM
REAL(EB), DIMENSION(:),   INTENT(IN), POINTER :: VB
INTEGER:: I, J, K, IOR0, IFACE, INBR, NOM, ICG, ICW, IWG, ITYPE
REAL(EB):: VAL, HB(-3:3) = 0.0_EB

ITYPE = TYPE_MGM_BC
IF (TYPE_MGM_BC == NSCARC_MGM_BC_EXPOL .AND. TOTAL_PRESSURE_ITERATIONS <= 2) ITYPE = NSCARC_MGM_BC_MEAN

!MGM_MESH_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

!   CALL SCARC_POINT_TO_MGM(NM, NL)
!   G  => L%UNSTRUCTURED
   BTYPE => MGM%BTYPE


   MGM_FACE_LOOP: DO IFACE = 1, 6

      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)

      MGM_NBR_LOOP: DO INBR = 1, F%N_NEIGHBORS
         
         NOM = F%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NLEVEL_MIN)
         
         MGM_CELL_LOOP: DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

            IWG = OG%ICG_TO_IWG(ICG)
            ICW = OG%ICG_TO_ICW(ICG, 1)

            IF (ICW == -1) CYCLE

            I = G%ICX(ICW) 
            J = G%ICY(ICW) 
            K = G%ICZ(ICW) 


            SELECT CASE (ITYPE)

               ! Boundary setting along interfaces by simple mean values

               CASE (NSCARC_MGM_BC_MEAN)

                  SELECT CASE (IOR0)
                     CASE ( 1)
                        VAL = 0.5_EB * (UHL(I, J, K) + OUHL(IWG)) 
                     CASE (-1)
                        VAL = 0.5_EB * (UHL(I, J, K) + OUHL(IWG)) 
                     CASE ( 2)
                        VAL = 0.5_EB * (UHL(I, J, K) + OUHL(IWG)) 
                     CASE (-2)
                        VAL = 0.5_EB * (UHL(I, J, K) + OUHL(IWG)) 
                     CASE ( 3)
                        VAL = 0.5_EB * (UHL(I, J, K) + OUHL(IWG)) 
                     CASE (-3)
                        VAL = 0.5_EB * (UHL(I, J, K) + OUHL(IWG)) 
                  END SELECT

               ! Boundary setting along interfaces by extrapolation

               CASE (NSCARC_MGM_BC_EXPOL)

                  SELECT CASE (IOR0)
                     CASE ( 1)
                        VAL = UHL(I, J, K) + OUHL(IWG) - 0.5_EB*(UHL2(I, J, K) + OUHL2(IWG))  
                     CASE (-1)
                        VAL = UHL(I, J, K) + OUHL(IWG) - 0.5_EB*(UHL2(I, J, K) + OUHL2(IWG))  
                     CASE ( 2)
                        VAL = UHL(I, J, K) + OUHL(IWG) - 0.5_EB*(UHL2(I, J, K) + OUHL2(IWG))  
                     CASE (-2)
                        VAL = UHL(I, J, K) + OUHL(IWG) - 0.5_EB*(UHL2(I, J, K) + OUHL2(IWG))  
                     CASE ( 3)
                        VAL = UHL(I, J, K) + OUHL(IWG) - 0.5_EB*(UHL2(I, J, K) + OUHL2(IWG))  
                     CASE (-3)
                        VAL = UHL(I, J, K) + OUHL(IWG) - 0.5_EB*(UHL2(I, J, K) + OUHL2(IWG))  
                  END SELECT

               ! Boundary setting along interfaces by simple mean values

               CASE (NSCARC_MGM_BC_TRUE)

                  HB = 0.0_EB
                  MGM_TRUE_IOR_SELECT: SELECT CASE (IOR0)

                     ! ---------------------------------------
                     CASE( 1)
                        HB( 1) = 0.5_EB * (OUHL(IWG)    + OUHL2(IWG) )
                        HB(-1) = 0.5_EB * (UHL(1, J, K) + UHL(2, J, K))
                        IF (.NOT.TWO_D) THEN
                           IF (BTYPE(IWG, 2) == INTERNAL) THEN
                              HB( 2) = 0.5_EB * (UHL(1, J-1, K) + OUHL(IWG+F%INCR_STENCIL(2)))
                           ENDIF
                           IF (BTYPE(IWG, -2) == INTERNAL) THEN
                              HB(-2) = 0.5_EB * (UHL(1, J+1, K) + OUHL(IWG+F%INCR_STENCIL(-2)))
                           ENDIF
                        ENDIF
                        IF (BTYPE(IWG, 3) == INTERNAL) THEN
                           HB( 3) = 0.5_EB * (UHL(1, J, K-1) + OUHL(IWG+F%INCR_STENCIL(3)))
                        ENDIF
                        IF (BTYPE(IWG, -3) == INTERNAL) THEN
                           HB(-3) = 0.5_EB *(UHL(1, J, K+1) + OUHL(IWG+F%INCR_STENCIL(-3)))
                        ENDIF


                     ! ---------------------------------------
                     CASE(-1)
                        HB( 1) = 0.5_EB * (UHL(I-1, J, K) + UHL(I, J, K))
                        HB(-1) = 0.5_EB * (OUHL(IWG)      + OUHL2(IWG) )
                        IF (.NOT.TWO_D) THEN
                           IF (BTYPE(IWG, 2) == INTERNAL) THEN
                              HB( 2) = 0.5_EB * (UHL(I, J-1, K) + OUHL(IWG+F%INCR_STENCIL(2)))
                           ENDIF
                           IF (BTYPE(IWG, -2) == INTERNAL) THEN
                              HB(-2) = 0.5_EB * (UHL(I, J+1, K) + OUHL(IWG+F%INCR_STENCIL(-2)))
                           ENDIF
                        ENDIF
                        IF (BTYPE(IWG, 3) == INTERNAL) THEN   
                           HB( 3) = 0.5_EB * (UHL(I, J, K-1) + OUHL(IWG+F%INCR_STENCIL(3)))
                        ENDIF
                        IF (BTYPE(IWG, -3) == INTERNAL) THEN  
                           HB(-3) = 0.5_EB * (UHL(I, J, K+1) + OUHL(IWG+F%INCR_STENCIL(-3)))
                        ENDIF

                     ! ---------------------------------------
                     CASE( 3)
                        WRITE(*,*) 'Not yet done'
                     ! ---------------------------------------
                     CASE(-3)
                        WRITE(*,*) 'Not yet done'
                  END SELECT MGM_TRUE_IOR_SELECT

                  IF (TWO_D) THEN
                     VAL = (L%DXI2*(HB(1)+HB(-1)) + L%DZI2*(HB(3)+HB(-3))) * MGM%WEIGHT(IWG)
                  ELSE
                     VAL = (L%DXI2*(HB(1)+HB(-1)) + L%DYI2*(HB(2)+HB(-2)) + L%DZI2*(HB(3)+HB(-3))) * MGM%WEIGHT(IWG)
                  ENDIF


            END SELECT

            SELECT CASE (IOR0)
               CASE ( 1)
                  MGM%BXS(J,K) = VAL
               CASE (-1)
                  MGM%BXF(J,K) = VAL
               CASE ( 2)
                  MGM%BYS(I,K) = VAL
               CASE (-2)
                  MGM%BYF(I,K) = VAL
               CASE ( 3)
                  MGM%BZS(I,J) = VAL
               CASE (-3)
                  MGM%BZF(I,J) = VAL
            END SELECT

            VB(ICW) = VB(ICW) + F%SCAL_DIRICHLET * VAL    

         ENDDO MGM_CELL_LOOP
      ENDDO MGM_NBR_LOOP
   ENDDO MGM_FACE_LOOP
!ENDDO MGM_MESH_LOOP

END SUBROUTINE SCARC_MGM_SET_INTERFACES


! --------------------------------------------------------------------------------------------------------------
!> \brief Set BC's along internal obstructions for MGM method
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_SET_OBSTRUCTIONS(VB)
USE SCARC_POINTERS, ONLY: L, G, MGM, UU, VV, WW, GWC, SCARC_POINT_TO_MGM
REAL(EB), DIMENSION(:), POINTER, INTENT(IN) :: VB
INTEGER:: IW, I, J, K, IOR0, IC
REAL(EB):: VAL

!MGM_MESH_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

!   CALL SCARC_POINT_TO_MGM(NM, NL)
!   G => L%UNSTRUCTURED


   MGM_OBST_LOOP: DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT
      
      GWC => G%WALL(IW)

      UU  => MGM%U1
      VV  => MGM%V1
      WW  => MGM%W1
      
      I = GWC%IXW
      J = GWC%IYW
      K = GWC%IZW
      
      IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
      
      IOR0 = GWC%IOR
      IC   = G%CELL_NUMBER(I, J, K)
      
      SELECT CASE (IOR0)
         CASE(1)
            VAL =  L%DXI*DTI*UU(I-1, J, K)
         CASE(-1)
            VAL = -L%DXI*DTI*UU(I, J, K)
         CASE(2)
            VAL =  L%DYI*DTI*VV(I, J-1, K)

         CASE(-2)
            VAL = -L%DYI*DTI*VV(I, J, K)
         CASE(3)
            VAL =  L%DZI*DTI*WW(I, J, K-1)
         CASE(-3)
            VAL = -L%DZI*DTI*WW(I, J, K)
      END SELECT

      !IF (BFIRST_WORKSPACE) VB(IC) = VB(IC) + VAL             ! Variant A
      VB(IC) = VB(IC) + VAL                                    ! Variant B
      
   ENDDO MGM_OBST_LOOP


!ENDDO MGM_MESH_LOOP

END SUBROUTINE SCARC_MGM_SET_OBSTRUCTIONS


! --------------------------------------------------------------------------------------------------------------
!> \brief Update velocities after either the first or second pass of the MGM method
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_UPDATE_VELOCITY(NTYPE)
USE SCARC_POINTERS, ONLY: M, L, GWC, MGM, UU, VV, WW, HP, SCARC_POINT_TO_MGM
INTEGER, INTENT(IN):: NTYPE
INTEGER  :: NM, I, J, K, IW, IOR0

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM(NM, NLEVEL_MIN)

   MGM_PART_SELECT: SELECT CASE (NTYPE)
      
      ! Update velocity with new information of previous structured inhomogeneous Poisson solution (SIP)
      
      CASE (NSCARC_MGM_POISSON)

         HP => MGM%SIP
         IF (PREDICTOR) THEN
            UU => M%U
            VV => M%V
            WW => M%W
         ELSE
            UU => M%US
            VV => M%VS
            WW => M%WS
         ENDIF

         DO K = 1, M%KBAR
            DO J = 1, M%JBAR
               DO I = 0, M%IBAR
                  MGM%U1(I, J, K) = UU(I, J, K) - DT*( M%FVX(I, J, K) + M%RDXN(I)*(HP(I+1, J, K)-HP(I, J, K)) )
               ENDDO
            ENDDO
         ENDDO
            
         DO K = 1, M%KBAR
            DO J = 0, M%JBAR
               DO I = 1, M%IBAR
                  MGM%V1(I, J, K) = VV(I, J, K) - DT*( M%FVY(I, J, K) + M%RDYN(J)*(HP(I, J+1, K)-HP(I, J, K)) )
               ENDDO
            ENDDO
         ENDDO
            
         DO K = 0, M%KBAR
            DO J = 1, M%JBAR
               DO I = 1, M%IBAR
                  MGM%W1(I, J, K) = WW(I, J, K) - DT*( M%FVZ(I, J, K) + M%RDZN(K)*(HP(I, J, K+1)-HP(I, J, K)) )
               ENDDO
            ENDDO
         ENDDO
            
         MGM%UVEL = MGM%U1
         MGM%VVEL = MGM%V1
         MGM%WVEL = MGM%W1 

      ! Update velocity with new information of previous unstructured homogeneous Laplace solution (UHL)
      
      CASE (NSCARC_MGM_LAPLACE)

         HP => MGM%UHL
            
         DO K = 1, M%KBAR
            DO J = 1, M%JBAR
               DO I = 0, M%IBAR
                  MGM%U2(I, J, K) = - DT*M%RDXN(I)*(HP(I+1, J, K)-HP(I, J, K))
               ENDDO
            ENDDO
         ENDDO

            
         DO K = 1, M%KBAR
            DO J = 0, M%JBAR
               DO I = 1, M%IBAR
                  MGM%V2(I, J, K) = - DT*M%RDYN(J)*(HP(I, J+1, K)-HP(I, J, K))
               ENDDO
            ENDDO
         ENDDO
            

         DO K = 0, M%KBAR
            DO J = 1, M%JBAR
               DO I = 1, M%IBAR
                  MGM%W2(I, J, K) = - DT*M%RDZN(K)*(HP(I, J, K+1)-HP(I, J, K))
               ENDDO
            ENDDO
         ENDDO

            
         ! Recompute velocities on obstruction cells, such that correct normal derivative of Laplace solution is used 
         DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT

            GWC => L%UNSTRUCTURED%WALL(IW)
            IF (GWC%BOUNDARY_TYPE /= SOLID_BOUNDARY) CYCLE

            IOR0 = GWC%IOR
            I = GWC%IXW
            J = GWC%IYW
            K = GWC%IZW
            
            SELECT CASE(IOR0)
               CASE( 1)
                  MGM%U2(I-1, J, K) = - MGM%U1(I-1, J, K)
               CASE(-1)
                  MGM%U2(I, J, K)   = - MGM%U1(I, J, K) 
               CASE( 2)
                  MGM%V2(I, J-1, K) = - MGM%V1(I, J-1, K) 
               CASE(-2)
                  MGM%V2(I, J, K)   = - MGM%V1(I, J, K) 
               CASE( 3)
                  MGM%W2(I, J, K-1) = - MGM%W1(I, J, K-1) 
               CASE(-3)
                  MGM%W2(I, J, K)   = - MGM%W1(I, J, K) 
            END SELECT

         ENDDO

         MGM%UVEL = MGM%U1+MGM%U2
         MGM%VVEL = MGM%V1+MGM%V2
         MGM%WVEL = MGM%W1+MGM%W2

   END SELECT MGM_PART_SELECT
ENDDO
 
END SUBROUTINE SCARC_MGM_UPDATE_VELOCITY


! --------------------------------------------------------------------------------------------------------------
!> \brief Set internal boundary conditions for unstructured, homogeneous part of McKeeney-Greengard-Mayo method
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_COMPUTE_VELOCITY_ERROR(NTYPE)
USE SCARC_POINTERS, ONLY: M, L, MGM, GWC, EWC, HP, SCARC_POINT_TO_MGM
INTEGER, INTENT(IN) ::  NTYPE
INTEGER:: NM, I, J, K, IW, IOR0, IIO1, IIO2, JJO1, JJO2, KKO1, KKO2, IIO, JJO, KKO, ITYPE
REAL(EB):: UN_NEW_OTHER, UN_NEW, DUDT, DVDT, DWDT
TYPE(MESH_TYPE), POINTER:: M2
TYPE(OMESH_TYPE), POINTER:: OM

MESH_REAL = 0.0_EB                            
RANK_REAL = 0.0_EB
UN_NEW_OTHER = 0.0_EB

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM(NM, NLEVEL_MIN)

   IF (NTYPE == NSCARC_MGM_POISSON) THEN
      HP  => MGM%SIP
   ELSE
      HP  => MGM%UIP
   ENDIF

   MGM%VELOCITY_ERROR = 0.0_EB

   SELECT CASE (NTYPE)

      ! ------------------------- Poisson case
       
      CASE (NSCARC_MGM_POISSON) 

         WALLCELLS_POISSON_LOOP: DO IW = 1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT

            GWC => L%STRUCTURED%WALL(IW)                        ! point to structured grid

            IF (GWC%BOUNDARY_TYPE /= SOLID_BOUNDARY         .AND. &
                GWC%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY) CYCLE

            IOR0 = GWC%IOR

            I = GWC%IXG
            J = GWC%IYG
            K = GWC%IZG

            ! Update normal component of velocity at the mesh boundary

            SELECT CASE(IOR0)
               CASE( 1)
                  UN_NEW = M%U(I, J, K)   - DT*(M%FVX(I, J, K)   + M%RDXN(I)  *(HP(I+1, J, K)-HP(I, J, K)))
               CASE(-1)
                  UN_NEW = M%U(I-1, J, K) - DT*(M%FVX(I-1, J, K) + M%RDXN(I-1)*(HP(I, J, K)-HP(I-1, J, K)))
               CASE( 2)
                  UN_NEW = M%V(I, J, K)   - DT*(M%FVY(I, J, K)   + M%RDYN(J)  *(HP(I, J+1, K)-HP(I, J, K)))
               CASE(-2)
                  UN_NEW = M%V(I, J-1, K) - DT*(M%FVY(I, J-1, K) + M%RDYN(J-1)*(HP(I, J, K)-HP(I, J-1, K)))
               CASE( 3)
                  UN_NEW = M%W(I, J, K)   - DT*(M%FVZ(I, J, K)   + M%RDZN(K)  *(HP(I, J, K+1)-HP(I, J, K)))
               CASE(-3)
                  UN_NEW = M%W(I, J, K-1) - DT*(M%FVZ(I, J, K-1) + M%RDZN(K-1)*(HP(I, J, K)-HP(I, J, K-1)))
            END SELECT

            IF (M%WALL(IW)%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY) THEN
         
               UN_NEW_OTHER = 0._EB
         
               EWC => M%EXTERNAL_WALL(IW)

               OM => M%OMESH(EWC%NOM)
               M2 => MESHES(EWC%NOM)

               IIO1 = EWC%IIO_MIN
               JJO1 = EWC%JJO_MIN
               KKO1 = EWC%KKO_MIN
               IIO2 = EWC%IIO_MAX
               JJO2 = EWC%JJO_MAX
               KKO2 = EWC%KKO_MAX
         
               IOR_SELECT_1: SELECT CASE(IOR0)
                  CASE( 1)
                     DO KKO = KKO1, KKO2
                        DO JJO = JJO1, JJO2
                           DO IIO = IIO1, IIO2
                              DUDT = -OM%FVX(IIO, JJO, KKO)   - M2%RDXN(IIO)  *(HP(I+1, J, K)-HP(I, J, K))
                              UN_NEW_OTHER = UN_NEW_OTHER+OM%U(IIO, JJO, KKO)   + DT*DUDT
                           ENDDO
                        ENDDO
                     ENDDO
                  CASE(-1)
                     DO KKO = KKO1, KKO2
                        DO JJO = JJO1, JJO2
                           DO IIO = IIO1, IIO2
                              DUDT = -OM%FVX(IIO-1, JJO, KKO) - M2%RDXN(IIO-1)*(HP(I, J, K)-HP(I-1, J, K))
                              UN_NEW_OTHER = UN_NEW_OTHER+OM%U(IIO-1, JJO, KKO) + DT*DUDT
                           ENDDO
                        ENDDO
                     ENDDO
                  CASE( 2)
                     DO KKO = KKO1, KKO2
                        DO JJO = JJO1, JJO2
                           DO IIO = IIO1, IIO2
                              DVDT = -OM%FVY(IIO, JJO, KKO)   - M2%RDYN(JJO)  *(HP(I, J+1, K)-HP(I, J, K))
                              UN_NEW_OTHER = UN_NEW_OTHER+OM%V(IIO, JJO, KKO)   + DT*DVDT
                           ENDDO
                        ENDDO
                     ENDDO
                  CASE(-2)
                     DO KKO = KKO1, KKO2
                        DO JJO = JJO1, JJO2
                           DO IIO = IIO1, IIO2
                              DVDT = -OM%FVY(IIO, JJO-1, KKO) - M2%RDYN(JJO-1)*(HP(I, J, K)-HP(I, J-1, K))
                              UN_NEW_OTHER = UN_NEW_OTHER+OM%V(IIO, JJO-1, KKO) + DT*DVDT
                           ENDDO
                        ENDDO
                     ENDDO
                  CASE( 3)
                     DO KKO = KKO1, KKO2
                        DO JJO = JJO1, JJO2
                           DO IIO = IIO1, IIO2
                              DWDT = -OM%FVZ(IIO, JJO, KKO)   - M2%RDZN(KKO)  *(HP(I, J, K+1)-HP(I, J, K))
                              UN_NEW_OTHER = UN_NEW_OTHER+OM%W(IIO, JJO, KKO)   + DT*DWDT
                           ENDDO
                        ENDDO
                     ENDDO
                  CASE(-3)
                     DO KKO = KKO1, KKO2
                        DO JJO = JJO1, JJO2
                           DO IIO = IIO1, IIO2
                              DWDT = -OM%FVZ(IIO, JJO, KKO-1) - M2%RDZN(KKO-1)*(HP(IIO, JJO, KKO)-HP(IIO, JJO, KKO-1))
                              UN_NEW_OTHER = UN_NEW_OTHER+OM%W(IIO, JJO, KKO-1) + DT*DWDT
                           ENDDO
                        ENDDO
                     ENDDO
               END SELECT IOR_SELECT_1
            ENDIF

            IF (M%WALL(IW)%BOUNDARY_TYPE == SOLID_BOUNDARY) THEN
               UN_NEW_OTHER = -SIGN(1._EB, REAL(IOR0, EB))*MESHES(NM)%WALL(IW)%ONE_D%U_NORMAL_S
            ENDIF

            ! Compute velocity difference

            MGM%VELOCITY_ERROR = MAX(MGM%VELOCITY_ERROR, ABS(UN_NEW-UN_NEW_OTHER))


         ENDDO WALLCELLS_POISSON_LOOP

      ! ------------------------- Laplace case
       
      CASE (NSCARC_MGM_LAPLACE) 

         WALLCELLS_LAPLACE_LOOP: DO IW = 1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT

            GWC => L%UNSTRUCTURED%WALL(IW)                     ! point to unstructured grid

            ITYPE = GWC%BOUNDARY_TYPE
            IF (.NOT. ((IW > L%N_WALL_CELLS_EXT .AND. ITYPE == SOLID_BOUNDARY) .OR. ITYPE == INTERPOLATED_BOUNDARY)) CYCLE

            IOR0 = GWC%IOR

            I = GWC%IXG
            J = GWC%IYG
            K = GWC%IZG

            ! Update normal component of velocity at the mesh boundary

            UN_NEW = 0.0_EB
            UN_NEW_OTHER = 0.0_EB
            SELECT CASE(IOR0)
               CASE( 1)
                  !UN_NEW = M%U(I, J, K)   - DT*(M%FVX(I, J, K)   + M%RDXN(I)  *(HP(I+1, J, K)-HP(I, J, K))*DHFCT)
                  UN_NEW = MGM%UVEL(I, J, K)
               CASE(-1)
                  !UN_NEW = M%U(I-1, J, K) - DT*(M%FVX(I-1, J, K) + M%RDXN(I-1)*(HP(I, J, K)-HP(I-1, J, K))*DHFCT)
                  UN_NEW = MGM%UVEL(I-1, J, K)
               CASE( 2)
                  !UN_NEW = M%V(I, J, K)   - DT*(M%FVY(I, J, K)   + M%RDYN(J)  *(HP(I, J+1, K)-HP(I, J, K))*DHFCT)
                  UN_NEW = MGM%VVEL(I, J, K)
               CASE(-2)
                  !UN_NEW = M%V(I, J-1, K) - DT*(M%FVY(I, J-1, K) + M%RDYN(J-1)*(HP(I, J, K)-HP(I, J-1, K))*DHFCT)
                  UN_NEW = MGM%VVEL(I, J-1, K)
               CASE( 3)
                  !UN_NEW = M%W(I, J, K)   - DT*(M%FVZ(I, J, K)   + M%RDZN(K)  *(HP(I, J, K+1)-HP(I, J, K))*DHFCT)
                  UN_NEW = MGM%WVEL(I, J, K) 
               CASE(-3)
                  !UN_NEW = M%W(I, J, K-1) - DT*(M%FVZ(I, J, K-1) + M%RDZN(K-1)*(HP(I, J, K)-HP(I, J, K-1))*DHFCT)
                  UN_NEW = MGM%WVEL(I, J, K-1) 
            END SELECT

            IF (GWC%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY) THEN
               SELECT CASE(ABS(IOR0))
                  CASE( 1)
                     UN_NEW_OTHER = MGM%OUVEL(IW)
                  CASE( 2)
                     UN_NEW_OTHER = MGM%OVVEL(IW)
                  CASE( 3)
                     UN_NEW_OTHER = MGM%OWVEL(IW)
               END SELECT
            ENDIF

            IF (GWC%BOUNDARY_TYPE == SOLID_BOUNDARY) THEN
               UN_NEW_OTHER = -SIGN(1._EB, REAL(IOR0, EB))*MESHES(NM)%WALL(IW)%ONE_D%U_NORMAL_S
            ENDIF

            ! Compute velocity difference

            MGM%VELOCITY_ERROR = MAX(MGM%VELOCITY_ERROR, ABS(UN_NEW-UN_NEW_OTHER))


         ENDDO WALLCELLS_LAPLACE_LOOP

   END SELECT

   MESH_REAL(NM) = MGM%VELOCITY_ERROR
   RANK_REAL = MAX(RANK_REAL, MESH_REAL(NM))

ENDDO MESHES_LOOP

IF (N_MPI_PROCESSES > 1) & 
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_REAL, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
VELOCITY_ERROR_GLOBAL = RANK_REAL

END SUBROUTINE SCARC_MGM_COMPUTE_VELOCITY_ERROR

END MODULE SCARC_MGM



!=======================================================================================================================
! 
! MODULE SCARC_METHODS
!
!> \brief Collection of available ScaRC/UScaRC solvers : 
!
!  - Krylov method (without and with coarse grid correction)
!  - McKeeney-Greengard-Mayo method (still experimental)
!  - Geometric Multigrid method
!  - Algebraic Multigrid method
!  - FFT method (Crayfish Pak)
!  - IntelMKL methods (Pardiso/Cluster_Sparse_Solver)
!
!=======================================================================================================================
MODULE SCARC_METHODS
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE SCARC_CONSTANTS
USE SCARC_TYPES
USE SCARC_VARIABLES
USE SCARC_STACK
USE SCARC_VECTORS
USE SCARC_CONVERGENCE

IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief  Setup environment for global Krylov method to solve the overall Poisson problem
! This includes 
!  - environment for the global Krylov method 
!  - environment for the (mostly) local or global preconditioners 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_KRYLOV_ENVIRONMENT
USE SCARC_FFT, ONLY: SCARC_SETUP_FFT, SCARC_SETUP_FFTO
#ifdef WITH_MKL
USE SCARC_MKL, ONLY: SCARC_SETUP_PARDISO, SCARC_SETUP_CLUSTER
#endif
INTEGER :: NSTACK

NSTACK = NSCARC_STACK_ROOT
STACK(NSTACK)%SOLVER => POISSON_SOLVER
CALL SCARC_SETUP_STACK_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

! Setup preconditioner for Krylov solver, all acting locally by default
 
NSTACK = NSTACK + 1
STACK(NSTACK)%SOLVER => POISSON_PRECON
SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

   CASE (NSCARC_RELAX_JAC)                                          ! Jacobi preconditioning
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

   CASE (NSCARC_RELAX_SSOR)                                         ! SSOR preconditioning
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

   CASE (NSCARC_RELAX_MJAC)                                         ! Jacobi preconditioning in matrix form
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MJAC(NLEVEL_MIN, NLEVEL_MAX)

   CASE (NSCARC_RELAX_MGS)                                          ! GS preconditioning in matrix form
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MGS(NLEVEL_MIN, NLEVEL_MAX)

   CASE (NSCARC_RELAX_MSGS)                                         ! SGS preconditioning in matrix form
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSGS(NLEVEL_MIN, NLEVEL_MAX)

   CASE (NSCARC_RELAX_MSOR)                                         ! SOR preconditioning in matrix form
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

   CASE (NSCARC_RELAX_MSSOR)                                        ! SSOR preconditioning in matrix form
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

   CASE (NSCARC_RELAX_LU)                                           ! LU preconditioning in matrix form
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_LU(NLEVEL_MIN, NLEVEL_MAX)

   CASE (NSCARC_RELAX_ILU)                                          ! ILU(0) preconditioning in matrix form
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_ILU(NLEVEL_MIN, NLEVEL_MAX)

   CASE (NSCARC_RELAX_FFT)                                          ! FFT preconditioning
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MIN)

   CASE (NSCARC_RELAX_FFTO)                                         ! FFT preconditioning with overlap
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_FFTO(NLEVEL_MIN, NLEVEL_MIN)

#ifdef WITH_MKL
   CASE (NSCARC_RELAX_MKL)                                          ! IntelMKL preconditioning
      SELECT CASE(TYPE_SCOPE(1))
         CASE (NSCARC_SCOPE_GLOBAL)                                 ! by global CLUSTER_SPARSE_SOLVER
            CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_GLOBAL)
            CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)
         CASE (NSCARC_SCOPE_LOCAL)                                  ! by local PARDISO solvers
            CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
      END SELECT

   CASE (NSCARC_RELAX_OPTIMIZED)                                    ! mixture of IntelMKL and FFT preconditioning
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_OPTIMIZED(NLEVEL_MIN, NLEVEL_MIN)
#endif
 
   CASE (NSCARC_RELAX_GMG)                                          ! Preconditioning by complete GMG method

      CALL SCARC_SETUP_STACK_PRECON(NSTACK, TYPE_SCOPE(1))
      CALL SCARC_SETUP_STACK_MULTIGRID(NSCARC_SOLVER_PRECON, TYPE_SCOPE(1), NSCARC_STAGE_TWO, NSTACK, &
                                       NLEVEL_MIN, NLEVEL_MAX)

      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => POISSON_SMOOTH
      SELECT CASE (TYPE_SMOOTH)

         CASE (NSCARC_RELAX_JAC)                                    ! Jacobi smoothing
            CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

         CASE (NSCARC_RELAX_SSOR)                                   ! SSOR smoothing
            CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

         CASE (NSCARC_RELAX_FFT)                                    ! FFT smoothing
            CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MAX-1)

#ifdef WITH_MKL
         CASE (NSCARC_RELAX_MKL)                                    ! IntelMKL smoothing (local PARDISO's)
            STACK(NSTACK)%SOLVER => POISSON_SMOOTH
            CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
#endif
      END SELECT

      ! Coarse grid solver (same scope of action as calling GMG)

      NSTACK = NSTACK + 1
      CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_TWO, TYPE_SCOPE(1), NSTACK, NLEVEL_MAX, NLEVEL_MAX)

END SELECT SELECT_KRYLOV_PRECON

 
! If two-level Krylov, allocate intermediate structures for interpolation and workspace for global coarse solver
 
IF (HAS_XMEAN_LEVELS) THEN
   CALL SCARC_SETUP_PRECON_XMEAN(NLEVEL_MIN)
ELSE IF (HAS_TWO_LEVELS) THEN
   IF (.NOT.IS_CG_AMG) CALL SCARC_SETUP_INTERPOLATION(NSCARC_STAGE_ONE, NLEVEL_MIN+1, NLEVEL_MAX)
   NSTACK = NSTACK + 1
   CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)
ENDIF

! Store final number of stacks

N_STACK_TOTAL = NSTACK                                 

END SUBROUTINE SCARC_SETUP_KRYLOV_ENVIRONMENT


#ifdef WITH_MKL
! --------------------------------------------------------------------------------------------------------------
!> \brief  Setup optimized local preconditioner, either local FFT or PARDISO, depending on structure of mesh
! If mesh does not contain obstructions, the faster FFT preconditioner is used, otherwise PARDISO from IntelMKL
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OPTIMIZED(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: L, SCARC_POINT_TO_LEVEL
USE SCARC_FFT, ONLY: SCARC_SETUP_FFT_MESH
USE SCARC_MKL, ONLY: SCARC_SETUP_PARDISO_MESH
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL

DO NL = NLMIN, NLMAX
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      CALL SCARC_POINT_TO_LEVEL(NM, NL)
      IF (L%HAS_OBSTRUCTIONS) THEN
         CALL SCARC_SETUP_PARDISO_MESH(NM, NL)
      ELSE
         CALL SCARC_SETUP_FFT_MESH(NM, NL)
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_OPTIMIZED
#endif


! --------------------------------------------------------------------------------------------------------------
!> \brief  Setup environment for global multigrid method to solve the overall Poisson problem
! This includes 
!  - environment for the global multigrid method 
!  - environment for the (mostly) local or global smoothers 
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MULTIGRID_ENVIRONMENT
USE SCARC_FFT, ONLY: SCARC_SETUP_FFT, SCARC_SETUP_FFTO
#ifdef WITH_MKL
USE SCARC_MKL, ONLY: SCARC_SETUP_CLUSTER, SCARC_SETUP_PARDISO
#endif
INTEGER :: NSTACK

NSTACK = NSCARC_STACK_ROOT
STACK(NSTACK)%SOLVER => POISSON_SOLVER
CALL SCARC_SETUP_STACK_MULTIGRID(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MAX)

NSTACK = NSTACK + 1
STACK(NSTACK)%SOLVER => POISSON_SMOOTH
SELECT CASE(TYPE_SMOOTH)

   CASE (NSCARC_RELAX_JAC)                                       ! Jacobi smoothing
      CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

   CASE (NSCARC_RELAX_SSOR)                                      ! SSOR smoothing
      CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

   CASE (NSCARC_RELAX_MJAC)                                      ! Jacobi smoothing in matrix form
      CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MJAC(NLEVEL_MIN, NLEVEL_MAX)

   CASE (NSCARC_RELAX_MGS)                                       ! GS smoothing in matrix form
      CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MGS(NLEVEL_MIN, NLEVEL_MAX)

   CASE (NSCARC_RELAX_MSGS)                                      ! SGS smoothing in matrix form
      CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSGS(NLEVEL_MIN, NLEVEL_MAX)

   CASE (NSCARC_RELAX_MSOR)                                      ! SOR smoothing in matrix form
      CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

   CASE (NSCARC_RELAX_MSSOR)                                     ! SSOR smoothing in matrix form
      CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

   CASE (NSCARC_RELAX_FFT)                                       ! FFT smoothing
      CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MAX-1)

#ifdef WITH_MKL
   CASE (NSCARC_RELAX_MKL)                                       ! IntelMKL smoothing
      CALL SCARC_SETUP_STACK_SMOOTH(NSTACK, TYPE_SCOPE(2))
      SELECT CASE(TYPE_SCOPE(2))
         CASE (NSCARC_SCOPE_GLOBAL)                              ! Globally acting:  Use global CLUSTER_SPARSE_SOLVER
            CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)
         CASE (NSCARC_SCOPE_LOCAL)                               ! Locally acting: Use local PARDISO solvers
            CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
      END SELECT
#endif

END SELECT

! Globally acting coarse grid solver

NSTACK = NSTACK + 1
CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)

! Store final number of stacks

N_STACK_TOTAL = NSTACK

END SUBROUTINE SCARC_SETUP_MULTIGRID_ENVIRONMENT


! --------------------------------------------------------------------------------------------------------------
!> \brief Setup environment needed for the McKenney-Greengard-Mayo method to solve the overall Poisson problem
! This includes
!  - environment for structured global solver in pass 1 of MGM (global structured Krylov & preconditioners)
!  - if required, environment for unstructured global solver in pass 1 of MGM (global unstructured Krylov & preconditioners)
!  - environment for unstructured local solvers in pass 2 of MGM
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_ENVIRONMENT
USE SCARC_POINTERS, ONLY: SCARC_POINT_TO_MGM
#ifdef WITH_MKL
USE SCARC_POINTERS, ONLY: L
USE SCARC_MKL, ONLY: SCARC_SETUP_PARDISO, SCARC_SETUP_MGM_PARDISO
USE SCARC_MATRICES, ONLY: SCARC_SETUP_MATRIX_MKL
#endif
USE SCARC_MGM, ONLY: SCARC_SETUP_MGM
USE SCARC_FFT, ONLY: SCARC_SETUP_FFT, SCARC_SETUP_MGM_FFT
INTEGER :: NSTACK, TYPE_MATRIX_SAVE
#ifdef WITH_MKL
INTEGER :: NM
#endif

! Allocate workspace and define variables for the different boundary settings in MGM-method

CALL SCARC_SETUP_MGM (NLEVEL_MIN)

! ------- First pass: Setup solver for inhomogeneous Poisson global problem on structured discretization
!         By default, a global CG-method with FFT-preconditioning is used
!         Goal is, to use PFFT or CG with PFFT-preconditioning later

TYPE_PRECON = NSCARC_RELAX_FFT
CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_STRUCTURED)

NSTACK = NSCARC_STACK_ROOT
STACK(NSTACK)%SOLVER => POISSON_SOLVER_STRUCTURED
CALL SCARC_SETUP_STACK_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

NSTACK = NSTACK + 1
STACK(NSTACK)%SOLVER => POISSON_PRECON_STRUCTURED
CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MIN)

! If exact initialization or comparison versus exact solution (that is, the difference USCARC-SCARC) is required,
! the global unstructured Poisson matrix was already assembled in SETUP_SYSTEMS, 
! thus also initialize the related preconditioners 

TYPE_MATRIX_SAVE = TYPE_MATRIX             
TYPE_MATRIX = NSCARC_MATRIX_COMPACT        
CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)

IF (SCARC_MGM_CHECK_LAPLACE .OR. SCARC_MGM_EXACT_INITIAL) THEN

   NSTACK = NSTACK + 1
   STACK(NSTACK)%SOLVER => POISSON_SOLVER_UNSTRUCTURED
   CALL SCARC_SETUP_STACK_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

   ! If IntelMKL-library is available, use local PARDISO preconditioners, otherwise local SSOR preconditioners 
   ! for the global unstructured Poisson problem 

   NSTACK = NSTACK + 1
   STACK(NSTACK)%SOLVER => POISSON_PRECON_UNSTRUCTURED
#ifdef WITH_MKL
   TYPE_PRECON = NSCARC_RELAX_MKL
   CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
   CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)        
#else
   TYPE_PRECON = NSCARC_RELAX_SSOR
   CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
   CALL SCARC_WARNING (NSCARC_WARNING_NO_MKL_PRECON, SCARC_NONE, NSCARC_NONE)
#endif
      
ENDIF

! ------- Second pass: Setup solvers for local homogeneous Laplace problems on unstructured discretization
!         Different local solvers are available (CG, own LU, own permuted LU, PARDISO, even FFT if mesh happens to be structured)

TYPE_MATRIX = NSCARC_MATRIX_COMPACT
CALL SCARC_SET_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)

NSTACK = NSTACK + 1
STACK(NSTACK)%SOLVER => LAPLACE_SOLVER_UNSTRUCTURED
N_STACK_LAPLACE = NSTACK

#ifdef WITH_MKL
TYPE_MKL (NLEVEL_MIN) = NSCARC_MKL_LOCAL
#endif

SELECT_LAPLACE_SOLVER: SELECT CASE (TYPE_MGM_LAPLACE)

   ! Setup CG-solvers for the solution of the local Laplace problems 
   ! If IntelMKL is available use PARDISO preconditioners by default, otherwise SSOR

   CASE (NSCARC_MGM_LAPLACE_CG) 

      CALL SCARC_SETUP_STACK_KRYLOV(NSCARC_SOLVER_MGM, NSCARC_SCOPE_LOCAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)
   
      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => LAPLACE_PRECON_UNSTRUCTURED

#ifdef WITH_MKL
      TYPE_PRECON = NSCARC_RELAX_MKL
      CALL SCARC_SETUP_STACK_MGM (NSTACK, NSCARC_SCOPE_LOCAL)
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_MGM (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_MGM_PARDISO (NM, NLEVEL_MIN)
      ENDDO
#else
      TYPE_PRECON = NSCARC_RELAX_SSOR
      CALL SCARC_SETUP_STACK_MGM(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_WARNING (NSCARC_WARNING_NO_MKL_PRECON, SCARC_NONE, NSCARC_NONE)
#endif

   ! Setup local LU or permuted LU solvers

   CASE (NSCARC_MGM_LAPLACE_LU, NSCARC_MGM_LAPLACE_LUPERM)
      CALL SCARC_SETUP_STACK_MGM(NSTACK, NSCARC_SCOPE_LOCAL)

#ifdef WITH_MKL

   ! Setup local IntelMKL-PARDISO solvers 

   CASE (NSCARC_MGM_LAPLACE_PARDISO)

      CALL SCARC_SETUP_STACK_MGM(NSTACK, NSCARC_SCOPE_LOCAL)
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_MGM (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_MGM_PARDISO (NM, NLEVEL_MIN)
      ENDDO

   ! Setup mixture of local Crayfishpak-FFT or IntelMKL-PARDISO solvers, depending on wether the mesh is structured or not

   CASE (NSCARC_MGM_LAPLACE_OPTIMIZED)

      CALL SCARC_SETUP_STACK_MGM(NSTACK, NSCARC_SCOPE_LOCAL)
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_MGM (NM, NLEVEL_MIN)
         IF (L%HAS_OBSTRUCTIONS) THEN
            CALL SCARC_SETUP_MGM_PARDISO (NM, NLEVEL_MIN)
         ELSE
            CALL SCARC_SETUP_MGM_FFT (NM, NLEVEL_MIN)
         ENDIF

      ENDDO
#endif

END SELECT SELECT_LAPLACE_SOLVER

TYPE_MATRIX = TYPE_MATRIX_SAVE

! Store final number of stacks

N_STACK_TOTAL = NSTACK       

END SUBROUTINE SCARC_SETUP_MGM_ENVIRONMENT

! ------------------------------------------------------------------------------------------------------------------
!> \brief Store Jacobi preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
! the MJAC-preconditioner in matrix form is defined
!           M_MJAC = D 
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MJAC(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, A, AB, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC

CROUTINE = 'SCARC_SETUP_MJAC'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    

      SELECT CASE(TYPE_MATRIX)
         CASE (NSCARC_MATRIX_COMPACT)
            A => G%POISSON
            CALL SCARC_ALLOCATE_REAL1 (A%RELAX, 1, G%NC, NSCARC_INIT_ZERO, 'G%POISSON%RELAX', CROUTINE)
            DO IC = 1, G%NC
               A%RELAX(IC) = 1.0_EB/A%VAL(A%ROW(IC))
            ENDDO 
         CASE (NSCARC_MATRIX_BANDWISE)
            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL1 (AB%RELAXD, 1, AB%N_DIAG,  NSCARC_INIT_ZERO, 'G%POISSON%RELAXD', CROUTINE)
            DO IC = 1, G%NC
               AB%RELAXD(IC) = 1.0_EB/AB%VAL(IC, AB%POS(0))
            ENDDO 
       END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MJAC


! ------------------------------------------------------------------------------------------------------------------
!> \brief Store GS preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
!          -E is the strictly lower part
!          -F is the strictly upper part
! the SGS-preconditioner in matrix form is defined
!           M_MGS = (D - E) = (I - E D^{-1}) D
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGS(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, A, AB, SCARC_POINT_TO_GRID
USE SCARC_UTILITIES, ONLY: SET_MATRIX_TYPE
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, IPTR

CROUTINE = 'SCARC_SETUP_MGS'
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    

      SELECT CASE (SET_MATRIX_TYPE(NL))

         ! ---------- Matrix in compact storage technique
 
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSON
            CALL SCARC_ALLOCATE_REAL1 (A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSON%RELAX', CROUTINE)

            DO IC = 1, G%NC
               DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                  JC = A%COL(IPTR)
                  IF (JC <  IC) A%RELAX(IPTR) = A%VAL(IPTR) / A%VAL(A%ROW(JC))
                  IF (JC == IC) A%RELAX(IPTR) = A%VAL(IPTR)
               ENDDO
            ENDDO 

 
         ! ---------- Matrix in bandwise storage technique
 
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2 (AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%POISSONB%RELAX', CROUTINE)
            WRITE(*,*) 'SCARC_SETUP_MGS: BANDWISE: NOT FINISHED YET'

      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MGS


! ------------------------------------------------------------------------------------------------------------------
!> \brief Store symmetric Gauss-Seidel preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
!          -E is the strictly lower part
!          -F is the strictly upper part
! the SGS-preconditioner in matrix form is defined
!           M_MSGS = (D - E) D^{-1} (D - F)  =  (I - E D^{-1}) (D - F)
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MSGS(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, A, AB, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, IPTR, I, IS, IL, IOR0

CROUTINE = 'SCARC_SETUP_MSGS'
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    

      SELECT CASE(TYPE_MATRIX)
 
         ! ---------- Matrix in compact storage technique
 
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSON
            CALL SCARC_ALLOCATE_REAL1 (A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSON%RELAX', CROUTINE)
            A%RELAX = A%VAL

            DO IC = 1, G%NC
               !  l(i,j) = a(i,j)/a(j,j)
               DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                  JC = A%COL(IPTR)
                  IF (JC < IC)  A%RELAX(IPTR) = A%VAL(IPTR) / A%VAL(A%ROW(JC))
               ENDDO
            ENDDO 
 
         ! ---------- Matrix in bandwise storage technique
 
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2 (AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%AB%RELAX', CROUTINE)
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


! ------------------------------------------------------------------------------------------------------------------
!> \brief Store SOR preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
!          -E is the strictly lower part
!          -F is the strictly upper part
! the SOR-preconditioner in matrix form is defined
!           M_MSOR = (D−ωE) = (I−ωE D^{-1}) D
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MSOR(NLMIN, NLMAX, NSTACK)
USE SCARC_POINTERS, ONLY: G, A, AB, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NLMIN, NLMAX, NSTACK
REAL (EB) :: OMEGA
INTEGER :: NM, NL, IC, JC, IPTR

CROUTINE = 'SCARC_SETUP_MSOR'
OMEGA = STACK(NSTACK)%SOLVER%OMEGA

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    

      SELECT CASE(TYPE_MATRIX)

         ! ---------- Matrix in compact storage technique
 
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSON
            CALL SCARC_ALLOCATE_REAL1 (A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSON%RELAX', CROUTINE)
      
            DO IC = 1, G%NC
      
               !DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
               !   JC = A%COL(IPTR)
               !   IF (JC <  IC) A%RELAX(IPTR) = OMEGA * A%VAL(IPTR) / A%VAL(A%ROW(JC))
               !   IF (JC == IC) A%RELAX(IPTR) = A%VAL(IPTR)
               !ENDDO
               DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                  JC = A%COL(IPTR)
                  IF (JC < IC) THEN
                     A%RELAX(IPTR) = OMEGA * A%VAL(IPTR) / A%VAL(A%ROW(JC))
                  ELSE IF (JC == IC) THEN
                     A%RELAX(IPTR) = A%VAL(IPTR)
                  ELSE IF (JC > IC) THEN
                     A%RELAX(IPTR) = OMEGA * A%VAL(IPTR) 
                  ENDIF
               ENDDO
      
            ENDDO 

         ! ---------- Matrix in bandwise storage technique
 
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2 (AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%POISSONB%RELAX', CROUTINE)
            WRITE(*,*) 'SCARC_SETUP_MSOR: BANDWISE: NOT FINISHED YET'
      
      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MSOR


! ------------------------------------------------------------------------------------------------------------------
!> \brief Store SSOR preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
!          -E is the strictly lower part
!          -F is the strictly upper part
! the SSOR-preconditioner in matrix form is defined
!           B_SSOR = 1/(omega * (2-omega)) * (D - omega * E) D^{-1} (D - omega * F)
!                  = (I - omega E D^{-1}) * [1/(omega * (2-omega)) * D -  1/(2-omega) * F]
! Defining the triangular matrices
!               L  = I - omega E D^{-1}
!               U  = 1/(omega * (2-omega)) * D -  1/(2-omega) * F
! the SSOR-preconditioning can be thought as the solution of two triangular systems
! Both matrices can be stored as a single matrix that occupies the same amount of storage as A
! where the same row and column pointers can be used as for A (identical pattern)
! Note that the diagonal elements of L are 1 (and are omitted, only the diagonal of U is stored there)
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MSSOR(NLMIN, NLMAX, NSTACK)
USE SCARC_POINTERS, ONLY: G, A, AB, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NLMIN, NLMAX, NSTACK
INTEGER :: NM, NL, IC, JC, IPTR, INCR, IS, IL, IOR0, I
REAL(EB) :: OMEGA, SCAL1, SCAL2

CROUTINE = 'SCARC_SETUP_MSSOR'

OMEGA = STACK(NSTACK)%SOLVER%OMEGA
SCAL1  = 1.0_EB / (OMEGA * (2.0_EB - OMEGA))
SCAL2  = 1.0_EB / (2.0_EB - OMEGA)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    

      SELECT CASE(TYPE_MATRIX)

         ! ---------- Matrix in compact storage technique
 
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSON
            CALL SCARC_ALLOCATE_REAL1 (A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSON%RELAX', CROUTINE)
      
            DO IC = 1, G%NC
      
               DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                  JC = A%COL(IPTR)
                  IF (JC < IC) THEN
                     A%RELAX(IPTR) = OMEGA * A%VAL(IPTR) / A%VAL(A%ROW(JC))
                  ELSE IF (JC == IC) THEN
                     A%RELAX(IPTR) = SCAL1 * A%VAL(IPTR)
                  ELSE IF (JC > IC) THEN
                     A%RELAX(IPTR) = SCAL2 * A%VAL(IPTR)
                  ENDIF
               ENDDO

            ENDDO 

         ! ---------- Matrix in bandwise storage technique
 
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2 (AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%POISSONB%RELAX', CROUTINE)

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
!               L  = I - omega E D^{-1}
!               U  = 1/(omega * (2-omega)) * D -  1/(2-omega) * F
      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

!STACK(NSTACK)%SOLVER%OMEGA = 1.0_EB

END SUBROUTINE SCARC_SETUP_MSSOR


! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize ILU(0) decomposition of Poisson matrix
! L- and U-parts are stored in the same array, diagonal elements of L are supposed to be 1
! Based on Saad-algorithm 10.4 from 'Iterative Methods for Sparse Linear Systems':
!   for i = 2 , ... , n do
!      for k = 1 , ... , i-1 and for (i,k) in NZ(A) do
!         compute a_ik = a_ik / a_kk
!         for j = k+1 , ... , n and for (i,j) in NZ(A) do
!            compute a_ij = a_ij - a_ik a_kj
!         enddo
!      enddo
!   enddo
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LU(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, A, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, KC, IPTR, JPTR, KPTR, KPTR0

CROUTINE = 'SCARC_SETUP_LU'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    

      A => G%POISSON
      CALL SCARC_ALLOCATE_REAL1 (A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSON%RELAX', CROUTINE)
      A%RELAX = A%VAL
    
      CELL_LOOP: DO IC = 2, G%NC
   
         COLUMN_LOOP: DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
  
            KC = A%COL(IPTR)                         ! get number of neighboring cell
            IF (KC >= IC) CYCLE                      ! only consider neighbors with lower cell numbers than IC
            IF (A%RELAX(IPTR) == 0) CYCLE
 
            KPTR = A%ROW(KC)                         ! get diagonal entry of neighbor
            A%RELAX(IPTR) = A%RELAX(IPTR)/A%RELAX(KPTR)

            DO JPTR = A%ROW(IC), A%ROW(IC+1)-1

               JC = A%COL(JPTR)
               IF (JC<=KC) CYCLE                     ! only consider neighbors with higher cell numbers than IC
               IF (A%RELAX(JPTR) == 0) CYCLE

               KPTR = -1
               DO KPTR0 = A%ROW(KC), A%ROW(KC+1)-1
                  IF (A%COL(KPTR0) == JC) THEN
                    KPTR = KPTR0
                  ENDIF
               ENDDO
               IF (KPTR>0) A%RELAX(JPTR) = A%RELAX(JPTR) - A%RELAX(IPTR) * A%RELAX(KPTR)

            ENDDO

         ENDDO COLUMN_LOOP
      ENDDO CELL_LOOP
   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_LU


! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize LU decomposition of Poisson matrix
! L- and U-parts are stored in the same array, diagonal elements of L are supposed to be 1
!   for i = 2 , ... , n do
!      for k = 1 , ... , i-1 and for (i,k) do
!         compute a_ik = a_ik / a_kk
!         for j = k+1 , ... , n and for (i,j) do
!            compute a_ij = a_ij - a_ik a_kj
!         enddo
!      enddo
!   enddo
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_ILU(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, A, AB, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, KC, IPTR, JPTR, KPTR, KPTR0, IOR0, JOR0, KOR0
LOGICAL :: BFOUND

CROUTINE = 'SCARC_SETUP_ILU'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    

      SELECT CASE(TYPE_MATRIX)

         ! ---------- Matrix in compact storage technique
 
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSON
            CALL SCARC_ALLOCATE_REAL1 (A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSON%RELAX', CROUTINE)
            A%RELAX = A%VAL
      
            CELL_LOOP: DO IC = 2, G%NC
      
               COLUMN_LOOP: DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
      
                  KC = A%COL(IPTR)                        ! get number of neighboring cell
                  IF (KC >= IC) CYCLE                      ! only consider neighbors with lower cell numbers than IC
                  IF (A%RELAX(IPTR) == 0) CYCLE
      
                  KPTR = A%ROW(KC)                        ! get diagonal entry of neighbor
                  A%RELAX(IPTR) = A%RELAX(IPTR)/A%RELAX(KPTR)
      
                  DO JPTR = A%ROW(IC), A%ROW(IC+1)-1
      
                     JC = A%COL(JPTR)
                     IF (JC<=KC) CYCLE                     ! only consider neighbors with higher cell numbers than IC
                     IF (A%RELAX(JPTR) == 0) CYCLE
      
                     KPTR = -1
                     DO KPTR0 = A%ROW(KC), A%ROW(KC+1)-1
                        IF (A%COL(KPTR0) == JC) THEN
                          KPTR = KPTR0
                        ENDIF
                     ENDDO
                     IF (KPTR>0) A%RELAX(JPTR) = A%RELAX(JPTR) - A%RELAX(IPTR) * A%RELAX(KPTR)
      
                  ENDDO
      
               ENDDO COLUMN_LOOP
            ENDDO CELL_LOOP
 
         ! ---------- Matrix in bandwise storage technique
 
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2 (AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%POISSONB%RELAX', CROUTINE)
            AB%RELAX = AB%VAL
      
            CELL_BANDWISE_LOOP: DO IC = 2, G%NC
      
               COLUMN_BANDWISE_LOOP: DO IOR0 = 3, -3, -1 
      
                  IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE

                  KC = IC + AB%OFFSET(IOR0)                           ! get number of neighboring cell
                  IF (KC<=0  .OR. KC >= IC) CYCLE                     ! only consider neighbors with lower cell numbers than IC
                  IF (AB%RELAX(IC, AB%POS(IOR0)) == 0) CYCLE
      
                  AB%RELAX(IC,AB%POS(IOR0)) = AB%RELAX(IC, AB%POS(IOR0))/AB%RELAX(KC,AB%POS(0))
      
                  DO JOR0 = 3, -3, -1

                     IF (TWO_D .AND. ABS(JOR0) == 2) CYCLE
      
                     JC = IC + AB%OFFSET(JOR0)                  ! get number of neighboring cell

                     IF (JC<=KC .OR. JC >G%NC) CYCLE            ! only consider neighbors with higher cell numbers than IC
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
      
               ENDDO COLUMN_BANDWISE_LOOP
            ENDDO CELL_BANDWISE_LOOP
      
      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_ILU


! --------------------------------------------------------------------------------------------------------------------
!> \brief Setup onedirectional Poisson-preconditioner
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PRECON_XMEAN(NL)
USE SCARC_POINTERS, ONLY: M, L, PRE, SCARC_POINT_TO_LEVEL
INTEGER, INTENT(IN) :: NL
INTEGER :: II, IX, NM

CROUTINE = 'SCARC_SETUP_PRECON_XMEAN'

NX_GLOBAL = 0
DO NM = 1, NMESHES
   NX_GLOBAL = NX_GLOBAL + MESHES(NM)%IBAR
ENDDO

CALL SCARC_ALLOCATE_INT1 (NX_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'NX_OFFSET', CROUTINE)
NX_OFFSET(1) = 0
DO NM = 2, NMESHES
   NX_OFFSET(NM) = NX_OFFSET(NM-1) + MESHES(NM-1)%IBAR
ENDDO

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_LEVEL (NM, NL)                                    
   PRE => L%PRECON

   CALL SCARC_ALLOCATE_REAL1 (PRE%DIAG, 1, NX_GLOBAL, NSCARC_INIT_ZERO, 'PRE%DIAG', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (PRE%LOW , 1, NX_GLOBAL, NSCARC_INIT_ZERO, 'PRE%LOW',  CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (PRE%UP  , 1, NX_GLOBAL, NSCARC_INIT_ZERO, 'PRE%UP',   CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (PRE%RHS , 1, NX_GLOBAL, NSCARC_INIT_ZERO, 'PRE%UP',   CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (PRE%AUX , 1, NX_GLOBAL, NSCARC_INIT_ZERO, 'PRE%AUX',  CROUTINE)

   DO IX = 1, L%NX
      II = NX_OFFSET(NM) + IX  
      PRE%DIAG(II) = -M%RDX(IX)*(M%RDXN(IX)+M%RDXN(IX-1))  ! Main diagonal
      PRE%UP(II)   =  M%RDX(IX)*M%RDXN(IX)                 ! Upper diagonal 
      PRE%LOW(II)  =  M%RDX(IX)*M%RDXN(IX-1)               ! Lower diagonal 
   ENDDO
      
   ! Apply boundary conditions at end of tunnel to the matrix components

   IF (NM==1) THEN
      IF (M%LBC==FISHPAK_BC_NEUMANN_NEUMANN .OR. M%LBC==FISHPAK_BC_NEUMANN_DIRICHLET) THEN     ! Neumann 
         PRE%DIAG(1) = PRE%DIAG(1) + PRE%LOW(1)
      ELSE                                                                                     ! Dirichlet
         PRE%DIAG(1) = PRE%DIAG(1) - PRE%LOW(1)
      ENDIF
   ENDIF

   IF (NM==NMESHES) THEN
      IF (M%LBC==FISHPAK_BC_NEUMANN_NEUMANN .OR. M%LBC==FISHPAK_BC_DIRICHLET_NEUMANN) THEN     ! Neumann BC
         PRE%DIAG(NX_GLOBAL) = PRE%DIAG(NX_GLOBAL) + PRE%UP(NX_GLOBAL)
      ELSE                                                                                     ! Dirichet BC
         PRE%DIAG(NX_GLOBAL) = PRE%DIAG(NX_GLOBAL) - PRE%UP(NX_GLOBAL)
      ENDIF
   ENDIF

   ! If not MPI process 0 send matrix components to process 0

   IF (MY_RANK>0) THEN  
   
      CALL MPI_GATHERV(PRE%UP(DISPLS_TP(MY_RANK)+1),COUNTS_TP(MY_RANK),MPI_DOUBLE_PRECISION,PRE%UP,COUNTS_TP,DISPLS_TP,&
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_GATHERV(PRE%LOW(DISPLS_TP(MY_RANK)+1),COUNTS_TP(MY_RANK),MPI_DOUBLE_PRECISION,PRE%LOW,COUNTS_TP,DISPLS_TP,&
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_GATHERV(PRE%DIAG(DISPLS_TP(MY_RANK)+1),COUNTS_TP(MY_RANK),MPI_DOUBLE_PRECISION,PRE%DIAG,COUNTS_TP,DISPLS_TP,&
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
      
   ! If MPI process 0 receive matrix components from slaves and compose global matrix

   ELSE  
   
      CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,PRE%UP,COUNTS_TP,DISPLS_TP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,PRE%LOW,COUNTS_TP,DISPLS_TP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,PRE%DIAG,COUNTS_TP,DISPLS_TP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
      
   ENDIF
      
ENDDO MESHES_LOOP
   
END SUBROUTINE SCARC_SETUP_PRECON_XMEAN


! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for MKL-methods
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSE_SOLVER(NSTAGE, NSCOPE, NSTACK, NLMIN, NLMAX)
#ifdef WITH_MKL
USE SCARC_MKL, ONLY: SCARC_SETUP_MKL, SCARC_SETUP_CLUSTER, SCARC_SETUP_PARDISO
#endif
INTEGER, INTENT(IN)    :: NSCOPE, NSTAGE, NLMIN, NLMAX
INTEGER, INTENT(INOUT) :: NSTACK

STACK(NSTACK)%SOLVER => COARSE_SOLVER         

SELECT_COARSE: SELECT CASE (TYPE_COARSE)

   ! -------------- CG-method is used as iterative coarse grid solver

   CASE (NSCARC_COARSE_ITERATIVE)

      CALL SCARC_SETUP_STACK_KRYLOV(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)

      NSTACK = NSTACK + 1
      TYPE_PRECON = NSCARC_RELAX_SSOR
      STACK(NSTACK)%SOLVER => COARSE_PRECON
      CALL SCARC_SETUP_STACK_PRECON(NSTACK, NSCOPE)

   ! -------------- LU-decomposition (from MKL) for coarse Poisson matrix is used as direct coarse grid solver

#ifdef WITH_MKL 
   CASE (NSCARC_COARSE_DIRECT)

      CALL SCARC_SETUP_MKL(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
      IF (NMESHES > 1) THEN
         CALL SCARC_SETUP_CLUSTER(NLMIN, NLMAX)                ! Use CLUSTER_SPARSE_SOLVER from IntelMKL
      ELSE
         CALL SCARC_SETUP_PARDISO(NLMIN, NLMAX)                ! Use PARDISO solver from IntelMKL
      ENDIF
#endif

   ! -------------- Otherwise: print error message
   CASE DEFAULT
      CALL SCARC_ERROR(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_COARSE)

END SELECT SELECT_COARSE
END SUBROUTINE SCARC_SETUP_COARSE_SOLVER


! ------------------------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for additive or multiplicative coarse grid
! (corresponding to Schwarz domain decomposition method)
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERPOLATION(NSTAGE, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, ST, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NSTAGE, NLMIN, NLMAX
INTEGER :: NM, NL

CROUTINE = 'SCARC_SETUP_INTERPOLATION'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                    
      ST => SCARC(NM)%LEVEL(NL)%STAGE(NSTAGE)

      CALL SCARC_ALLOCATE_REAL1 (ST%X, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%X', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (ST%B, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%B', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (ST%V, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Q', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (ST%R, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%W', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (ST%Y, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Y', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1 (ST%Z, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Z', CROUTINE)

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_INTERPOLATION


! ------------------------------------------------------------------------------------------------------------------
!> \brief Set initial solution corresponding to boundary data in BXS, BXF, ...
! ------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WORKSPACE(NS, NL)
USE SCARC_POINTERS, ONLY: M, L, F, G, SV, ST, STP, GWC, PRHS, HP, MGM, VB, VX, &
                          BXS, BXF, BYS, BYF, BZS, BZF, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_MGM
USE SCARC_MGM, ONLY: SCARC_MGM_SET_OBSTRUCTIONS, SCARC_MGM_SET_INTERFACES
INTEGER, INTENT(IN) :: NS, NL
REAL(EB) :: VAL
INTEGER  :: NM, IW, IW1, IW2, IOR0, I, J, K, IC


SV => STACK(NS)%SOLVER

SELECT_SOLVER_TYPE: SELECT CASE (SV%TYPE_SOLVER)

   ! ---------- If used as main solver use values from pressure-routine as initialization
 
   CASE (NSCARC_SOLVER_MAIN)
   
      MAIN_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_GRID (NM, NL)                                    
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
   
         PRHS => M%PRHS
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF

         ! Get right hand side (PRHS from pres.f90) and initial vector (H or HS from last time step)

         BXS => M%BXS   ;  BXF => M%BXF
         BYS => M%BYS   ;  BYF => M%BYF
         BZS => M%BZS   ;  BZF => M%BZF

         !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
         DO IC = 1, G%NC
            ST%X(IC) = HP(G%ICX(IC), G%ICY(IC), G%ICZ(IC))        ! use last iterate as initial solution
            ST%B(IC) = PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC))      ! get new RHS from surrounding code
         ENDDO                         
         !$OMP END PARALLEL DO

         !!$OMP PARALLEL 
         MAIN_INHOMOGENEOUS_LOOP: DO IOR0 = -3, 3, 1 
      
            IF (IOR0 == 0) CYCLE
            F => SCARC(NM)%LEVEL(NL)%FACE(IOR0)
            
            IW1 = F%NCW0
            IW2 = F%NCW0 + F%NCW - 1
      
            !!$OMP DO PRIVATE(IW, GWC, I, J, K, IC, VAL) SCHEDULE(STATIC)
            MAIN_FACE_INHOMOGENEOUS_LOOP: DO IW = IW1, IW2
      
               GWC => G%WALL(IW)
   
               I = GWC%IXW
               J = GWC%IYW
               K = GWC%IZW
   
               IF (TWO_D .AND. J /= 1) CALL SCARC_ERROR(NSCARC_ERROR_GRID_INDEX, SCARC_NONE, J)
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
   
               IC = G%CELL_NUMBER(I,J,K)
   
               ! ---------- Dirichlet BC's:
               ! these are based on the SETTING in BTYPE
               ! in the structured case this corresponds to the face-wise SETTING according to the FFT
               ! (this allows to use local FFT's as preconditioners)
               ! in the unstructured case only open boundary cells lead to Dirichlet BC's

               MAIN_IF_DIRICHLET: IF (GWC%BTYPE == DIRICHLET) THEN
   
                  SELECT CASE (IOR0)
                     CASE (1)
                        VAL =  BXS(J,K)
                     CASE (-1)
                        VAL =  BXF(J,K)
                     CASE (2)
                        VAL =  BYS(I,K)
                     CASE (-2)
                        VAL =  BYF(I,K)
                     CASE (3)
                        VAL =  BZS(I,J)
                     CASE (-3)
                        VAL =  BZF(I,J)
                  END SELECT
   
                  ST%B(IC) = ST%B(IC) + F%SCAL_DIRICHLET * VAL
               ENDIF MAIN_IF_DIRICHLET
   
               ! ---------- Neumann BC's:
               ! Note for the unstructured case only:
               ! Here, the matrix also contains Neumann BC's for those cells which have a
               ! PRESSURE_BC_INDEX == DIRICHLET but are NOT open; these cells must be excluded below,
               ! because BXS, BXF, ... contain the Dirichlet information from pres.f90 there;
               ! excluding them corresponds to a homogeneous Neumann condition for these cells

               MAIN_IF_NEUMANN: IF (GWC%BTYPE == NEUMANN) THEN
   
                  IF (IS_UNSTRUCTURED .AND. M%WALL(IW)%PRESSURE_BC_INDEX /= NEUMANN) CYCLE
   
                  SELECT CASE (IOR0)
                     CASE (1)
                        VAL =  BXS(J,K)
                     CASE (-1)
                        VAL =  BXF(J,K)
                     CASE (2)
                        VAL =  BYS(I,K)
                     CASE (-2)
                        VAL =  BYF(I,K)
                     CASE (3)
                        VAL =  BZS(I,J)
                     CASE (-3)
                        VAL =  BZF(I,J)
                  END SELECT
   
                  ST%B(IC) = ST%B(IC) + F%SCAL_NEUMANN * VAL

               ENDIF MAIN_IF_NEUMANN
      
            ENDDO MAIN_FACE_INHOMOGENEOUS_LOOP
            !!$OMP END DO

         ENDDO MAIN_INHOMOGENEOUS_LOOP
         !!$OMP END PARALLEL 
   
   
      ENDDO MAIN_MESHES_LOOP
      
      ! In case of a Krylov method clear overlapping parts of auxiliary vectors

      IF (IS_CG.OR.IS_MGM.OR.HAS_TWO_LEVELS) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            ST%D(1:G%NCE) = 0.0_EB
            ST%R(1:G%NCE) = 0.0_EB
            ST%V(1:G%NCE) = 0.0_EB
            ST%Y(1:G%NCE) = 0.0_EB
            ST%Z(1:G%NCE) = 0.0_EB
         ENDDO
      ENDIF
   
      ! In case of a multigrid method as main solver clear
      ! overlapping parts of auxiliary vectors and coarse grid solver vectors

      IF (IS_GMG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            ST%V(1:G%NCE) = 0.0_EB
            ST%Z(1:G%NCE) = 0.0_EB
         ENDDO
      ENDIF
   
      ! In case of pure Neumann or periodic BCs, broadcast RHS(end) from last mesh
      ! to all and store it on all meshes

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
   
 
   ! ---------- If used as second pass in a MGM method, use basically zero RHS with special boundary values along
   !            mesh interfaces and internal obstructions
 
   CASE (NSCARC_SOLVER_MGM)

      MGM_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_MGM (NM, NL)                                    
         G => L%UNSTRUCTURED
   
         BXS => MGM%BXS ;  BXS = 0.0_EB
         BXF => MGM%BXF ;  BXF = 0.0_EB
         BYS => MGM%BYS ;  BYS = 0.0_EB
         BYF => MGM%BYF ;  BYF = 0.0_EB
         BZS => MGM%BZS ;  BZS = 0.0_EB
         BZF => MGM%BZF ;  BZF = 0.0_EB

         IF (TYPE_MGM_LAPLACE == NSCARC_MGM_LAPLACE_CG) THEN
            VB => ST%B ;  VB = 0.0_EB
            VX => ST%X ;  VX = 0.0_EB
         ELSE
            VB => MGM%B ;  VB = 0.0_EB
            VX => MGM%X ;  VX = 0.0_EB
         ENDIF

         CALL SCARC_MGM_SET_INTERFACES (VB, NM)
         CALL SCARC_MGM_SET_OBSTRUCTIONS (VB)

         !!$OMP PARALLEL 
         MGM_INHOMOGENEOUS_LOOP: DO IOR0 = -3, 3, 1 
      
            IF (IOR0 == 0) CYCLE
            F => SCARC(NM)%LEVEL(NL)%FACE(IOR0)
            
            IW1 = F%NCW0
            IW2 = F%NCW0 + F%NCW - 1
      
            !!$OMP DO PRIVATE(IW, GWC, I, J, K, IC, VAL) SCHEDULE(STATIC)
            MGM_FACE_INHOMOGENEOUS_LOOP: DO IW = IW1, IW2
      
               GWC => G%WALL(IW)
   
               I = GWC%IXW
               J = GWC%IYW
               K = GWC%IZW
   
               IF (TWO_D .AND. J /= 1) CALL SCARC_ERROR(NSCARC_ERROR_GRID_INDEX, SCARC_NONE, J)
               IF (L%IS_SOLID(I, J, K)) CYCLE
   
               IC = G%CELL_NUMBER(I,J,K)
   
               MGM_IF_DIRICHLET: IF (GWC%BTYPE == DIRICHLET) THEN
   
                  SELECT CASE (IOR0)
                     CASE (1)
                        VAL =  BXS(J,K)
                     CASE (-1)
                        VAL =  BXF(J,K)
                     CASE (2)
                        VAL =  BYS(I,K)
                     CASE (-2)
                        VAL =  BYF(I,K)
                     CASE (3)
                        VAL =  BZS(I,J)
                     CASE (-3)
                        VAL =  BZF(I,J)
                  END SELECT
   
                  VB(IC) = VB(IC) + F%SCAL_DIRICHLET * VAL
   
               ENDIF MGM_IF_DIRICHLET
   
               MGM_IF_NEUMANN: IF (GWC%BTYPE == NEUMANN) THEN
   
                  IF (IS_UNSTRUCTURED .AND. M%WALL(IW)%PRESSURE_BC_INDEX /= NEUMANN) CYCLE
   
                  SELECT CASE (IOR0)
                     CASE (1)
                        VAL =  BXS(J,K)
                     CASE (-1)
                        VAL =  BXF(J,K)
                     CASE (2)
                        VAL =  BYS(I,K)
                     CASE (-2)
                        VAL =  BYF(I,K)
                     CASE (3)
                        VAL =  BZS(I,J)
                     CASE (-3)
                        VAL =  BZF(I,J)
                  END SELECT
   
                  VB(IC) = VB(IC) + F%SCAL_NEUMANN * VAL

   
               ENDIF MGM_IF_NEUMANN
      
            ENDDO MGM_FACE_INHOMOGENEOUS_LOOP
            !!$OMP END DO

         ENDDO MGM_INHOMOGENEOUS_LOOP
         !!$OMP END PARALLEL 
   
      ENDDO MGM_MESHES_LOOP
      
   ! ---------- If MG is used as Krylov preconditioner, vector R of main Krylov is the new RHS for MG
 
   CASE (NSCARC_SOLVER_PRECON)
   
      IF (IS_CG_MG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L   => SCARC(NM)%LEVEL(NL)
            ST  => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)              ! current stage
            STP => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)           ! parent stage
            ST%X(1:G%NCE) = 0.0_EB
            ST%B(1:G%NCE) = STP%R(1:G%NCE)                                                
            ST%V(1:G%NCE) = 0.0_EB
            ST%Z(1:G%NCE) = 0.0_EB
         ENDDO
      ENDIF
   
   ! ---------- If used as coarse grid solver start with zero initialization
 
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


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform requested coarse grid solver (iterative/direct)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_COARSE(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL

SELECT CASE (TYPE_COARSE)

   CASE (NSCARC_COARSE_ITERATIVE)
      CALL SCARC_METHOD_KRYLOV (NSTACK, NPARENT, NLEVEL)

   CASE (NSCARC_COARSE_DIRECT)
#ifdef WITH_MKL
      !IF (STACK(NPARENT)%SOLVER%TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL .AND. NMESHES > 1) THEN
      IF (NMESHES > 1) THEN
         CALL SCARC_METHOD_CLUSTER (NSCARC_MATRIX_POISSON_SYM, NSTACK, NPARENT, NLEVEL)
      ELSE
         CALL SCARC_METHOD_PARDISO (NSCARC_MATRIX_POISSON_SYM, NSTACK, NPARENT, NLEVEL)
      ENDIF
#else
      CALL SCARC_ERROR(NSCARC_ERROR_DIRECT_NOMKL, SCARC_NONE, NLEVEL)
#endif

END SELECT

END SUBROUTINE SCARC_METHOD_COARSE


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform global conjugate gradient method based on global Possion-matrix
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_KRYLOV(NSTACK, NPARENT, NLEVEL)
USE SCARC_MATRICES, ONLY: SCARC_SETUP_SYSTEM_CONDENSED
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER :: NSTATE, NS, NP, NL, NG
REAL (EB) :: ALPHA0, BETA0, SIGMA0, SIGMA1=0.0_EB
REAL (EB) :: TNOW, TNOWI

TNOW = CURRENT_TIME()
ITE_CG = 0

TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS = NSTACK
NP = NPARENT
NL = NLEVEL
NG = TYPE_GRID                              ! TODO: Why? Forgot the reason ...

! ---------- Initialization
!   - Get parameters for current scope (note: NL denotes the finest level)
!   - Get right hand side vector and clear solution vectors

CALL SCARC_SETUP_SCOPE(NS, NP)
TYPE_GRID = NG                              ! TODO: Why? Forgot the reason ...

CALL SCARC_SETUP_WORKSPACE(NS, NL)

!CALL SCARC_PRESET_VECTOR(B, NL)

! In case of pure Neumann boundary conditions setup condensed system

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_VECTOR_INIT (X, 0.0_EB, NL)                    
   CALL FILTER_MEANVALUE(B, NL)                       
   CALL SCARC_SETUP_SYSTEM_CONDENSED (B, NL, 1)            
ENDIF

! Compute initial residual 

IF (IS_MGM .AND. NSTACK == 3) THEN
   TYPE_MATVEC = NSCARC_MATVEC_LOCAL
ELSE
   TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
ENDIF
CALL SCARC_MATVEC_PRODUCT (X, R, NL)                         !  r^0 := A*x^0
CALL SCARC_VECTOR_SUM     (B, R, -1.0_EB, 1.0_EB, NL)        !  r^0 := r^0 - b     corresponds to  A*x^0 - b

RES    = SCARC_L2NORM (R, NL)                                !  res   := ||r^0||
RESIN  = RES                                                 !  resin := res
NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                 !  res < tolerance ?

! Perform initial preconditioning

IF (NSTATE /= NSCARC_STATE_CONV_INITIAL) THEN                !  if no convergence yet, call intial preconditioner
   CALL SCARC_PRECONDITIONER(NS, NS, NL)                     !  v^0 := Precon(r^0)
   SIGMA1 = SCARC_SCALAR_PRODUCT(R, V, NL)                   !  SIGMA1 := (r^0,v^0)
   CALL SCARC_VECTOR_COPY (V, D, -1.0_EB, NL)                !  d^0 := -v^0
ENDIF

! ---------- Perform conjugate gradient looping

CG_LOOP: DO ITE = 1, NIT

   TNOWI = CURRENT_TIME()
   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   !TYPE_MATVEC = NSCARC_MATVEC_LOCAL
   CALL SCARC_MATVEC_PRODUCT (D, Y, NL)                        !  y^k := A*d^k

   ALPHA0 = SCARC_SCALAR_PRODUCT (D, Y, NL)                    !  ALPHA0 := (d^k,y^k)     corresponds to   (d^k,A*d^k)
   ALPHA0 = SIGMA1/ALPHA0                                      !  ALPHA0 := (r^k,v^k)/(d^k,A*d^k)

   CALL SCARC_VECTOR_SUM (D, X, ALPHA0, 1.0_EB, NL)            !  x^{k+1} := x^k + ALPHA0 * d^k
   CALL SCARC_VECTOR_SUM (Y, R, ALPHA0, 1.0_EB, NL)            !  r^{k+1} := r^k + ALPHA0 * y^k   ~  r^k + ALPHA0 * A * d^k

   RES = SCARC_L2NORM (R, NL)                                  !  res := ||r^{k+1}||
   NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                !  res < tolerance ??
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP

   CALL SCARC_PRECONDITIONER(NS, NS, NL)                       !  v^{k+1} := Precon(r^{k+1})

   SIGMA0 = SCARC_SCALAR_PRODUCT (R, V, NL)                    !  SIGMA0 := (r^{k+1},v^{k+1})
   BETA0  = SIGMA0/SIGMA1                                      !  BETA0  := (r^{k+1},v^{k+1})/(r^k,v^k)
   SIGMA1 = SIGMA0                                             !  save last SIGMA0

   CALL SCARC_VECTOR_SUM (V, D, -1.0_EB, BETA0, NL)            !  d^{k+1} := -v^{k+1} + BETA0 * d^{k+1}

   CPU(MY_RANK)%ITERATION=MAX(CPU(MY_RANK)%ITERATION,CURRENT_TIME()-TNOWI)

ENDDO CG_LOOP

! ---------- Determine convergence rate and print corresponding information
! In case of CG as main solver:
!   - Transfer ScaRC solution vector X to FDS pressure vector
!   - Set ghost cell values along external boundaries
!   - Exchange values along internal boundaries

CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL RESTORE_LAST_CELL(X, NL)
   CALL FILTER_MEANVALUE(X, NL)
ENDIF


IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN .AND. .NOT.IS_MGM) THEN
   CALL SCARC_UPDATE_MAINCELLS(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF


CALL SCARC_RELEASE_SCOPE(NS, NP)

END SUBROUTINE SCARC_METHOD_KRYLOV


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform McKeeney-Greengard-Mayo (MGM) method
! Note that the MGM method only works on finest grid level NLEVEL_MIN by default
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MGM(NSTACK)
USE SCARC_MGM
USE SCARC_CONVERGENCE, ONLY: NIT_MGM
INTEGER, INTENT(IN) :: NSTACK
INTEGER :: ITE_MGM = 0, STATE_MGM
LOGICAL :: USE_CORRECT_INITIALIZATION

CALL SCARC_SETUP_MGM_WORKSPACE(NLEVEL_MIN)
! ------------- Pass 1: Compute global structured inhomogeneous Poisson solution SIP

CALL SCARC_SET_SYSTEM_TYPE (NSCARC_GRID_STRUCTURED, NSCARC_MATRIX_POISSON)
CALL SCARC_METHOD_KRYLOV (NSTACK, NSCARC_STACK_ZERO, NLEVEL_MIN)

CALL SCARC_MGM_STORE (NSCARC_MGM_POISSON)                   ! store this structured inhomogeneous Poisson solution in MGM%SIP
CALL SCARC_MGM_UPDATE_GHOSTCELLS (NSCARC_MGM_POISSON)       ! update ghost cell values correspondingly (global solution!)

CALL SCARC_MGM_COPY (NSCARC_MGM_SIP_TO_UIP)                 ! Initialize unstructured inhomogeneous Poisson UIP with SIP

CALL SCARC_MGM_UPDATE_VELOCITY (NSCARC_MGM_POISSON)         ! update velocity based on SIP
CALL SCARC_MGM_COMPUTE_VELOCITY_ERROR (NSCARC_MGM_POISSON)  ! compute related velocity error

! Determine if correct initialization of interface boundary values is required
! This is only needed in the multi-mesh case and only in the first or first two (in case of extrapolated BCs) pressure iterations 

USE_CORRECT_INITIALIZATION = NMESHES > 1 .AND. SCARC_MGM_EXACT_INITIAL .AND. &
                             ((TOTAL_PRESSURE_ITERATIONS <= 1) .OR. &
                              (TOTAL_PRESSURE_ITERATIONS <= 2  .AND.TYPE_MGM_BC == NSCARC_MGM_BC_EXPOL))

! The upper computation of the structured inhomogeneous Poisson (SIP) solution corresponds to the usual ScaRC solution
! (To this end the structured Poisson matrix was assembled during the setup phase)
! If comparison with exact solution is required, also compute UScaRC solution
! (In this case the unstructured Poisson matrix was also assembled during the setup phase)
! In both cases inhomogeneous external BC's are used and the ghost cells are set correspondingly (both are global solutions)
! Compute the difference DScaRC of UScaRC and ScaRC 

IF (SCARC_MGM_CHECK_LAPLACE .OR. USE_CORRECT_INITIALIZATION) THEN

   CALL SCARC_MGM_STORE (NSCARC_MGM_SCARC)                                 
   CALL SCARC_MGM_UPDATE_GHOSTCELLS (NSCARC_MGM_SCARC)    

   CALL SCARC_SET_SYSTEM_TYPE (NSCARC_GRID_UNSTRUCTURED, NSCARC_MATRIX_POISSON)
   CALL SCARC_METHOD_KRYLOV (NSTACK, NSCARC_STACK_ZERO, NLEVEL_MIN)             ! compute UScaRC with unstructured CG-method 

   CALL SCARC_MGM_STORE (NSCARC_MGM_USCARC)                                
   CALL SCARC_MGM_UPDATE_GHOSTCELLS (NSCARC_MGM_USCARC)                   

   CALL SCARC_MGM_DIFF (NSCARC_MGM_USCARC_VS_SCARC)                             ! build difference DSCARC of USCARC and SCARC

ENDIF

! Only if correct initialization is required:
! In the very first pressure iteration ever, or - in case of extrapolated MGM BCs - in the two very first pressure iterations 
!    - use UScaRC solution as unstructured inhomogeneous Poisson UIP solution of this MGM pass 
!    - use difference DScaRC of UScaRC and ScaRC as unstructured homogeneous Laplace UHL solution of this MGM pass
! Exchange the interface values of the local Laplace problems for the later BC setting in the next MGM call
! In this case the requested velocity tolerance has been reached by default here and this MGM call can be left
! Otherwise: Check if the requested velocity tolerance has already been reached by pass 1

IF (USE_CORRECT_INITIALIZATION) THEN

   CALL SCARC_MGM_COPY (NSCARC_MGM_USCARC_TO_UIP)        
   CALL SCARC_MGM_COPY (NSCARC_MGM_DSCARC_TO_UHL)        

   SELECT CASE (TYPE_MGM_BC)
      CASE (NSCARC_MGM_BC_MEAN)
          CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MGM_SINGLE, NSCARC_NONE, NLEVEL_MIN)
      CASE (NSCARC_MGM_BC_EXPOL)
          CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MGM_SINGLE, NSCARC_NONE, NLEVEL_MIN)
          IF (TOTAL_PRESSURE_ITERATIONS == 1) THEN
             CALL SCARC_MGM_COPY (NSCARC_MGM_UHL_TO_UHL2)        ! store also second last values for UHL
             CALL SCARC_MGM_COPY (NSCARC_MGM_OUHL_TO_OUHL2)      ! store also second last values for other UHL
          ENDIF
      CASE (NSCARC_MGM_BC_TRUE)
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MGM_DOUBLE, NSCARC_NONE, NLEVEL_MIN)
   END SELECT
 
   STATE_MGM = NSCARC_MGM_SUCCESS

ELSE

   STATE_MGM = SCARC_MGM_CONVERGENCE_STATE(0, 0)

ENDIF

! ------------- Pass 2: Solve local unstructured homogeneous Laplace solutions UHL
! This is only necessary if the requested accuracy has not already been achieved by pass 1

IF (STATE_MGM /= NSCARC_MGM_SUCCESS) THEN

   CALL SCARC_SET_SYSTEM_TYPE (NSCARC_GRID_UNSTRUCTURED, NSCARC_MATRIX_LAPLACE)

   MGM_CORRECTION_LOOP: DO ITE_MGM = 1, NIT_MGM

      ! Compute local Laplace problems either by (permuted) LU- or CG-method
      ! In both cases the following is done within the solver:
      ! - definition of  BC's along obstructions according to MGM-algorithm 
      ! - definition of  BC's along interfaces by 'MEAN', 'EXTRAPOLATION' or 'TRUE' based on previous Laplace solutions

      SELECT CASE (TYPE_MGM_LAPLACE)
         CASE (NSCARC_MGM_LAPLACE_CG)
            CALL SCARC_METHOD_KRYLOV (N_STACK_LAPLACE, NSCARC_STACK_ZERO, NLEVEL_MIN)
         CASE (NSCARC_MGM_LAPLACE_LU, NSCARC_MGM_LAPLACE_LUPERM)
            CALL SCARC_METHOD_MGM_LU(N_STACK_LAPLACE, NSCARC_STACK_ZERO, NLEVEL_MIN)
#ifdef WITH_MKL
         CASE (NSCARC_MGM_LAPLACE_PARDISO)
            CALL SCARC_METHOD_MGM_PARDISO(N_STACK_LAPLACE, NSCARC_STACK_ZERO, NLEVEL_MIN)
         CASE (NSCARC_MGM_LAPLACE_OPTIMIZED)
            CALL SCARC_METHOD_MGM_OPTIMIZED(N_STACK_LAPLACE, NSCARC_STACK_ZERO, NLEVEL_MIN)
#endif
      END SELECT
   
      CALL SCARC_MGM_STORE (NSCARC_MGM_LAPLACE)            

      ! Exchange interface data between neighboring meshes according to chosen boundary method
       
      SELECT CASE (TYPE_MGM_BC)
         CASE (NSCARC_MGM_BC_MEAN)
            CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MGM_SINGLE, NSCARC_NONE, NLEVEL_MIN)
         CASE (NSCARC_MGM_BC_TRUE)
            CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MGM_DOUBLE, NSCARC_NONE, NLEVEL_MIN)
         CASE (NSCARC_MGM_BC_EXPOL)
            CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MGM_SINGLE, NSCARC_NONE, NLEVEL_MIN)
            CALL SCARC_MGM_COPY (NSCARC_MGM_UHL_TO_UHL2)                            ! store second time level
            CALL SCARC_MGM_COPY (NSCARC_MGM_OUHL_TO_OUHL2)
      END SELECT

      CALL SCARC_MGM_UPDATE_GHOSTCELLS (NSCARC_MGM_LAPLACE)
      CALL SCARC_MGM_STORE (NSCARC_MGM_MERGE)
   
      ! Get new velocities based on local Laplace solutions and compute corresponding velocity error

      CALL SCARC_MGM_UPDATE_VELOCITY (NSCARC_MGM_LAPLACE)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MGM_VELO, NSCARC_NONE, NLEVEL_MIN)
      CALL SCARC_MGM_COMPUTE_VELOCITY_ERROR (NSCARC_MGM_LAPLACE)

      IF (SCARC_MGM_CHECK_LAPLACE) THEN
         CALL SCARC_MGM_DIFF (NSCARC_MGM_UHL_VS_DSCARC)       ! unstructured homogeneous Laplace vs difference UScaRC-ScaRC
         CALL SCARC_MGM_DIFF (NSCARC_MGM_UIP_VS_USCARC)       ! unstructured inhomogeneous Poisson vs UScaRC
      ENDIF

      STATE_MGM = SCARC_MGM_CONVERGENCE_STATE(ITE_MGM, 1)
      IF (STATE_MGM == NSCARC_MGM_SUCCESS) EXIT MGM_CORRECTION_LOOP

      IF (TYPE_MGM_BC == NSCARC_MGM_BC_TRUE) THEN
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MGM_DOUBLE, NSCARC_NONE, NLEVEL_MIN)   ! also 2. layer of interface adjacent cells
      ELSE 
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MGM_SINGLE, NSCARC_NONE, NLEVEL_MIN)   ! only 1. layer of interface adjacent cells
      ENDIF

   ENDDO MGM_CORRECTION_LOOP
      
   STATE_MGM = SCARC_MGM_CONVERGENCE_STATE(ITE_MGM-1, -1)
   
ENDIF

! Reset method type (which has been changed during Krylov method) to MGM
TYPE_METHOD = NSCARC_METHOD_MGM
CALL SCARC_MGM_STORE (NSCARC_MGM_TERMINATE)

END SUBROUTINE SCARC_METHOD_MGM


! -------------------------------------------------------------------------------------------------------------
!> \brief Perform LU-decompositions for local unstructured Laplace matrices 
! Two different variants are available:
! - 'LU'     : self-coded LU decomposition 
! - 'LUPERM' : self-coded LU decomposition with additional grid permutation (non-zero RHS entries at the end)
! Note: in both cases, L and U are stored as compact matrices
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MGM_LU(NS, NP, NL)
USE SCARC_POINTERS, ONLY: L, G, MGM, A, LO, UP, ST, SCARC_POINT_TO_MGM, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN):: NS, NP, NL
INTEGER:: IC, JC, NM

CALL SCARC_SETUP_SCOPE(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL)

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM (NM, NL)   
   G => L%UNSTRUCTURED

   A   => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE)
   LO  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LOWER)
   UP  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_UPPER)

   ST  => L%STAGE(STACK(NS)%SOLVER%TYPE_STAGE)


   DO IC = G%NONZERO, G%NC
      MGM%Y(IC) = MGM%B(G%PERM_BW(IC))
      DO JC = 1, IC-1
         MGM%Y(IC) = MGM%Y(IC) - SCARC_EVALUATE_CMATRIX(LO, IC, JC) * MGM%Y(JC)
      ENDDO
   ENDDO

   DO IC = G%NC, 1, -1
      MGM%X(IC) = MGM%Y(IC)
      DO JC = IC+1, G%NC
         MGM%X(IC) = MGM%X(IC) - SCARC_EVALUATE_CMATRIX(UP, IC, JC) * MGM%X(JC)
      ENDDO
      MGM%X(IC) = MGM%X(IC)/SCARC_EVALUATE_CMATRIX(UP, IC, IC)
   ENDDO

ENDDO

CALL SCARC_RELEASE_SCOPE(NS, NP)
END SUBROUTINE SCARC_METHOD_MGM_LU


#ifdef WITH_MKL
! --------------------------------------------------------------------------------------------------------------
!> \brief Perform solution of local Laplace problems by IntelMKL Pardiso methods on each mesh
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MGM_PARDISO(NS, NP, NL)
USE SCARC_POINTERS, ONLY: L, G, MGM, MKL, AS, ST, SCARC_POINT_TO_MGM, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NS, NP, NL
INTEGER :: NM
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

CALL SCARC_SETUP_SCOPE(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM (NM, NL)                                    
   G   => L%UNSTRUCTURED
   AS  => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE_SYM)
   MKL => L%MKL
   ST  => L%STAGE(STACK(NS)%SOLVER%TYPE_STAGE)

   MKL%PHASE  = 33                    ! only solving

   IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_DOUBLE) THEN
      MGM%X = 0.0_EB
      CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                     AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, MGM%B, MGM%X, MKL%ERROR)
   ELSE
      MGM%B_FB = REAL(MGM%B, FB)
      MGM%X_FB = 0.0_EB
      CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                     AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, MGM%B_FB, MGM%X_FB, MKL%ERROR)
      MGM%X = REAL(MGM%X_FB,EB)
   ENDIF
   
   IF (MKL%ERROR /= 0) CALL SCARC_ERROR(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

ENDDO MESHES_LOOP

CALL SCARC_RELEASE_SCOPE(NS, NP)
END SUBROUTINE SCARC_METHOD_MGM_PARDISO


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform solution of local Laplace problems with Pardiso or FFT, depending on grid type:
!- if mesh happens to be unstructured : Use IntelMKL Pardiso
!- if mesh happens to be structured   : Use Crayfishpak FFT
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MGM_OPTIMIZED (NS, NP, NL)
USE SCARC_POINTERS, ONLY: L, G, MGM, MKL, FFT, AS, ST, SCARC_POINT_TO_MGM, SCARC_POINT_TO_CMATRIX
USE POIS, ONLY: H2CZSS, H3CZSS
INTEGER, INTENT(IN) :: NS, NP, NL
INTEGER :: NM, IC
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

CALL SCARC_SETUP_SCOPE(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MGM (NM, NL)                                    
   G => L%UNSTRUCTURED
   ST  => L%STAGE(STACK(NS)%SOLVER%TYPE_STAGE)

   SELECT CASE (L%HAS_OBSTRUCTIONS)

      ! If mesh contains obstructions, then the grid is really unstructured and PARDISO from IntelMKL is used
       
      CASE (.TRUE.)

         AS => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE_SYM)
      
         MKL => L%MKL
         MKL%PHASE  = 33         ! only solving

         IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_DOUBLE) THEN
            MGM%X = 0.0_EB
            CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                           AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                           MKL%MSGLVL, MGM%B, MGM%X, MKL%ERROR)
         ELSE
            MGM%B_FB = REAL(MGM%B, FB)
            MGM%X_FB = 0.0_EB
            CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                           AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                           MKL%MSGLVL, MGM%B_FB, MGM%X_FB, MKL%ERROR)
            MGM%X = REAL(MGM%X_FB, EB)
         ENDIF
      
         IF (MKL%ERROR /= 0) CALL SCARC_ERROR(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)


      ! If mesh contains obstructions, then the grid is structured and the FFT from Crayfishpak is used 
       
      CASE (.FALSE.)

         FFT => L%FFT

         FFT%PRHS = 0.0_EB
         !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
         DO IC = 1, G%NC
            FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = MGM%B(IC)
         ENDDO
         !$OMP END PARALLEL DO

         IF (TWO_D) THEN
            CALL H2CZSS (MGM%BXS,  MGM%BXF, MGM%BZS, MGM%BZF, FFT%ITRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ELSE
            CALL H3CZSS (MGM%BXS,  MGM%BXF, MGM%BYS, MGM%BYF, MGM%BZS, MGM%BZF, FFT%ITRN, FFT%JTRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ENDIF

         !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
         DO IC = 1, G%NC
            MGM%X(IC) = FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) 
         ENDDO
         !$OMP END PARALLEL DO 

   END SELECT

ENDDO MESHES_LOOP

CALL SCARC_RELEASE_SCOPE(NS, NP)
END SUBROUTINE SCARC_METHOD_MGM_OPTIMIZED
#endif


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform geometric multigrid method based on global possion-matrix
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MULTIGRID(NSTACK, NPARENT, NLEVEL)
USE SCARC_UTILITIES, ONLY: CYCLING_CONTROL
USE SCARC_GMG, ONLY: SCARC_RESTRICTION, SCARC_PROLONGATION
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER :: NS, NP, NL
INTEGER :: NSTATE, ICYCLE
REAL (EB) :: TNOW, TNOW_COARSE

TNOW = CURRENT_TIME()
ITE_MG = 0

! Store current and parent stack position and current level

TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS = NSTACK
NP = NPARENT
NL = NLEVEL

 
! ---------- Initialization:
!   - Save SETTING (in case that subsequent solvers with different SETTING are called)
!   - Define parameters for current scope (note: NL denotes the finest level)
!   - Initialize solution, right hand side vector
  
CALL SCARC_SETUP_SCOPE(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL)

  
! ---------- Compute initial defect:  
!            RESIN := || B - A*X ||
!   - Initialize cycle counts for MG-iteration
!   - Perform initial matrix-vector product on finest level
!   - calculate norm of initial residual on finest level
  

CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                  !  V := A*X
CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                     !  V := B - V


RES    = SCARC_L2NORM (V, NL)                                         !  RESIN := ||V||
RESIN  = RES
NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                          !  RES < TOL already ??

IF (TYPE_SOLVER == NSCARC_SOLVER_PRECON .AND. RESIN <= 1E-6_EB) THEN
   CALL SCARC_VECTOR_SUM (V, X, 1.0_EB, 1.0_EB, NL)                    !  x := omega * v + x
   CALL SCARC_UPDATE_PRECONDITIONER(NLEVEL_MIN)
   CALL SCARC_RELEASE_SCOPE(NS, NP)
   RETURN
ENDIF

ICYCLE = CYCLING_CONTROL(NSCARC_CYCLING_SETUP, NL)

  
! ---------- Perform multigrid-looping (start each iteration on finest level)
  
MULTIGRID_LOOP: DO ITE = 1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   NL = NLEVEL_MIN
   ICYCLE = CYCLING_CONTROL(NSCARC_CYCLING_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLING_EXIT)

      ! Presmoothing  (smoothing/restriction till coarsest level is reached)
      ! initial and final residual are passed via vector V by default
 
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)
         !IF (ITE /= 1) CALL SCARC_SMOOTHER (NSCARC_CYCLING_PRESMOOTH, NS+1, NS, NL)         ! D_fine   := Smooth(defect)
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_PRESMOOTH, NS+1, NS, NL)         ! D_fine   := Smooth(defect)
         CALL SCARC_RESTRICTION (V, B, NL, NL+1)                              ! B_coarse := Rest(D_fine)
         CALL SCARC_VECTOR_CLEAR (X, NL+1)                                    ! use zero initial guess on coarse level
         NL = NL + 1                                                          ! set coarser level
      ENDDO PRESMOOTHING_LOOP

 
      ! Coarse grid solver
 
      TNOW_COARSE = CURRENT_TIME()
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)                          ! X_coarse := exact_sol(.)
      CPU(MY_RANK)%COARSE =CPU(MY_RANK)%COARSE+CURRENT_TIME()-TNOW_COARSE

 
      ! Postsmoothing (smoothing/restriction till finest level is reached again)
 
      POSTSMOOTHING_LOOP: DO WHILE (NL > NLEVEL_MIN)
         NL=NL-1
         CALL SCARC_PROLONGATION (X, V, NL+1, NL)                             ! V_fine := Prol(X_coarse)
         CALL SCARC_VECTOR_SUM (V, X, 1.0_EB, 1.0_EB, NL)                     ! X_fine := V_fine + X_fine

         CALL SCARC_SMOOTHER (NSCARC_CYCLING_POSTSMOOTH, NS+1, NS, NL)        ! V_fine := Smooth(defect)
         ICYCLE = CYCLING_CONTROL(NSCARC_CYCLING_PROCEED, NL)           ! perform requested cycle
         IF (ICYCLE /= NSCARC_CYCLING_POSTSMOOTH) CYCLE CYCLE_LOOP
      ENDDO POSTSMOOTHING_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) CALL SCARC_ERROR(NSCARC_ERROR_MULTIGRID_LEVEL, SCARC_NONE, NL)

 
   ! Compute norm of new residual on finest level and  leave loop correspondingly
 
   CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                       ! V := A*X
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                          ! V := F - V

   RES = SCARC_L2NORM (V, NL)                                                 ! RES := ||V||
   NSTATE = SCARC_CONVERGENCE_STATE(0, NS, NL)                                ! convergence ?
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT MULTIGRID_LOOP

ENDDO MULTIGRID_LOOP

  
! ---------- Determine convergence rate and print corresponding information:
! In case of MG as main solver:
!   - Transfer ScaRC solution vector X to FDS pressure vector
!   - Set ghost cell values along external boundaries
!   - Exchange values along internal boundaries (consistency!)
  
CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)


SELECT CASE (TYPE_SOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      CALL SCARC_UPDATE_MAINCELLS(NLEVEL_MIN)
      CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
   CASE (NSCARC_SOLVER_PRECON)
      CALL SCARC_UPDATE_PRECONDITIONER(NLEVEL_MIN)
END SELECT

CALL SCARC_RELEASE_SCOPE(NS, NP)

END SUBROUTINE SCARC_METHOD_MULTIGRID


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform requested MKL solver (global/local)
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MKL(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
#ifndef WITH_MKL
INTEGER :: NDUMMY
#endif

#ifdef WITH_MKL
SELECT_MKL: SELECT CASE (TYPE_MKL(0))
   CASE (NSCARC_MKL_GLOBAL)
      CALL SCARC_METHOD_CLUSTER(NSCARC_MATRIX_POISSON_SYM, NSTACK, NPARENT, NLEVEL)
   CASE (NSCARC_MKL_LOCAL)
      CALL SCARC_METHOD_PARDISO(NSCARC_MATRIX_POISSON_SYM, NSTACK, NPARENT, NLEVEL)
END SELECT SELECT_MKL
#else
NDUMMY = NPARENT * NLEVEL
CALL SCARC_ERROR(NSCARC_ERROR_MKL_STACK, SCARC_NONE, NSTACK)
#endif

END SUBROUTINE SCARC_METHOD_MKL


! --------------------------------------------------------------------------------------------------------------------
!> \brief Perform preceding FFT method to improve start solution for ScaRC
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_FFT
USE MESH_POINTERS
USE POIS, ONLY: H2CZSS, H3CZSS
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
INTEGER :: NM, I, J, K
LOGICAL :: WITH_BDRY = .FALSE.

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   !CALL POINT_TO_MESH(NM)
   
   IF (PREDICTOR) THEN
      HP => H
   ELSE
      HP => HS
   ENDIF
   
   ! Call the Poisson solver
 
   IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
   IF (TWO_D .AND. .NOT.CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
   
   DO K=1,KBAR
     DO J=1,JBAR
        DO I=1,IBAR
            HP(I,J,K) = PRHS(I,J,K)
        ENDDO
      ENDDO
   ENDDO
   
   ! Apply boundary conditions to H
 
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


#ifdef WITH_MKL
! --------------------------------------------------------------------------------------------------------------
!> \brief Perform global Pardiso-method based on MKL
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CLUSTER(NMATRIX, NSTACK, NPARENT, NLEVEL)
USE SCARC_POINTERS, ONLY: L, G, MKL, V1, V2, AS, V1_FB, V2_FB, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_VECTOR_FB, &
                          SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NMATRIX, NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NS, NP, NL
REAL (EB) :: TNOW

NS = NSTACK
NP = NPARENT
NL = NLEVEL

TNOW = CURRENT_TIME()

CALL SCARC_SETUP_SCOPE(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   AS => SCARC_POINT_TO_CMATRIX (NMATRIX)

   V1 => SCARC_POINT_TO_VECTOR (NM, NL, B)
   V2 => SCARC_POINT_TO_VECTOR (NM, NL, X)

   MKL => L%MKL
   MKL%PHASE  = 33                                ! only solving


   IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

      V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, B)
      V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, X)

      V1_FB = REAL(V1, FB)
      V2_FB = 0.0_FB
      CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                   AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)
      V2 = REAL(V2_FB, EB)

   ELSE

      CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                   AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)
   ENDIF

   IF (MKL%ERROR /= 0) CALL SCARC_ERROR(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

ENDDO MESHES_LOOP

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, X, NL)

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SCOPE(NS, NP)

END SUBROUTINE SCARC_METHOD_CLUSTER


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform global Pardiso-method based on MKL
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_PARDISO(NMATRIX, NSTACK, NPARENT, NLEVEL)
USE SCARC_POINTERS, ONLY: L, G, MKL, AS, V1, V2, V1_FB, V2_FB, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_VECTOR_FB, &
                          SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NMATRIX, NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NS, NP, NL
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

NS = NSTACK
NP = NPARENT
NL = NLEVEL

CALL SCARC_SETUP_SCOPE(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   AS => SCARC_POINT_TO_CMATRIX (NMATRIX)

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
                     AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

      V2 = REAL(V2_FB, EB)

   ELSE


      V2 = 0.0_EB
      CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                     AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, V1, V2, MKL%ERROR)
   ENDIF

   IF (MKL%ERROR /= 0) CALL SCARC_ERROR(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

ENDDO MESHES_LOOP

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SCOPE(NSTACK, NPARENT)

END SUBROUTINE SCARC_METHOD_PARDISO
#endif


! ----------------------------------------------------------------------------------------------------------------------
!> \brief Copy final solution from GMG (as preconditioner) to corresponding vector of CG (as main solver)
! ----------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRECONDITIONER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%V = SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%X
ENDDO
END SUBROUTINE SCARC_UPDATE_PRECONDITIONER


! ----------------------------------------------------------------------------------------------------------------------
!> \brief Finalize data for pressure vector (predictor/corrector) when local ScaRC solver has finished
! ----------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_MAINCELLS(NL)
USE SCARC_POINTERS, ONLY: M, G, L, ST, HP, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC 

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    
   ST  => L%STAGE(NSCARC_STAGE_ONE)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF


   HP = 0.0_EB
   !!$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
   DO IC = 1, G%NC
      HP (G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = ST%X(IC)
   ENDDO
   !!$OMP END PARALLEL DO 


ENDDO

END SUBROUTINE SCARC_UPDATE_MAINCELLS


! --------------------------------------------------------------------------------------------------------------
!> \brief Set correct boundary values at external and internal boundaries
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_GHOSTCELLS(NL)
USE SCARC_POINTERS, ONLY: M, L, G, GWC, HP, SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IXG, IYG, IZG, IXW, IYW, IZW 

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

   ! Compute ghost cell values
 
   !!$OMP PARALLEL DO SHARED(HP, M, L, G) PRIVATE(IW, IXG, IYG, IZG, IXW, IYW, IZW, IOR0, GWC) SCHEDULE(STATIC)
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
   !!$OMP END PARALLEL DO


ENDDO

! Perform data exchange to achieve consistency of ghost values along internal boundaries
! Note: this is most probably no longer necessary because MESH_EXCHANGE(5) is used after the call of ScaRC

CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_PRESSURE, NSCARC_NONE, NL)

   
END SUBROUTINE SCARC_UPDATE_GHOSTCELLS
   
! -------------------------------------------------------------------------------------------------------------
!> \brief Preconditioning method which is based on the following input and output convention:
!  - the residual which has to be preconditioned is passed in via vector R
!  - the result of preconditioning is passed out via vector V
!  - for several variants Y and Z are used as auxiliary vectors
!  - in the comments: call is based on current grid level l (mostly the finest one)
!  -                  l=1 denotes the finest  grid level NLEVEL_MIN
!  -                  l=L denotes the coarset grid level NLEVEL_MAX
! -------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRECONDITIONER(NS, NP, NL)
USE SCARC_GMG, ONLY: SCARC_RESTRICTION, SCARC_PROLONGATION
INTEGER, INTENT(IN) :: NS, NP, NL     
INTEGER :: IL

SELECT_PRECON_TYPE: SELECT CASE (TYPE_TWOLEVEL)

   ! ---------- Classical one-level preconditioning
   ! Solve local fine Poisson problems by selected relaxation method (SSOR, FFT, PARDISO, etc.)
 
   CASE (NSCARC_TWOLEVEL_NONE)

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !  v := r
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !  v := Relax(r)
 
   ! ---------- Additive two-level preconditioning : COARSE_ADD
   ! Local fine Poisson problems and global coarse Poisson problem are based on same defect
 
   CASE (NSCARC_TWOLEVEL_COARSE_ADD)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !  Use r^l as right hand side for preconditioner
      DO IL = NL, NLEVEL_MAX-1                                    !  successively restrict to coarser levels up to coarsest
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !  b^{l+1} := Restriction(r^l)
      ENDDO
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Z, 1.0_EB, NLEVEL_MAX)           !  z^L := x^L
      DO IL = NLEVEL_MAX-1, NL, -1                                !  successively interpolate to finer levels up to finest
         CALL SCARC_PROLONGATION(Z, Z, IL+1, IL)                  !  z^l := Prolongation(z^{l+1})
      ENDDO
      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !  v^l := r^l
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !  v^l := Relax(r^l)
      CALL SCARC_VECTOR_SUM (Z, V, SCARC_COARSE_OMEGA, 1.0_EB, NL)    !  v^l := z^l + v^l

   ! ---------- Multiplicative two-level preconditioning : COARSE_MUL1
   ! First solve global coarse Poisson problem, build new defect, then solve local fine Poisson problems
 
   CASE (NSCARC_TWOLEVEL_COARSE_MUL1)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !  Use r^l as right hand side for preconditioner
      DO IL = NL, NLEVEL_MAX-1
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !  b^{l+1} := Restriction(r^l)
      ENDDO
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !  y^L := x^L
      DO IL = NLEVEL_MAX-1, NLEVEL_MIN, -1
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !  y^l := Prolongation(y^{l+1})
      ENDDO
      CALL SCARC_MATVEC_PRODUCT (Y, Z, NL)                        !  z^l := A^l * y^l
      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !  z^l := r^l - z^l
      CALL SCARC_VECTOR_COPY (Z, V, 1.0_EB, NL)                   !  v^l := z^l
      CALL SCARC_RELAXATION (Z, V, NS+1, NP, NL)                  !  v^l := Relax(z^l)
      CALL SCARC_VECTOR_SUM (Y, V, SCARC_COARSE_OMEGA, 1.0_EB, NL)   !  v^l := y^l - z^l

   ! ---------- Multiplicative two-level preconditioning : COARSE_MUL2
   ! First solve local fine Poisson problems, build new defect, then solve global coarse Poisson problem

   CASE (NSCARC_TWOLEVEL_COARSE_MUL2)

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !  v^l := r^l
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !  v^l := Relax(r^l)
      CALL SCARC_MATVEC_PRODUCT (V, Z, NL)                        !  z^l := A^{l} * v^l
      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !  z^l := r^l - z^l
      CALL SCARC_RESTRICTION (Z, B, NL, NL+1)                     !  b^{l+1} := rest(R^{l})
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !  x^{l+1} := A^{l+1}^{-1}(b^{l+1})
      CALL SCARC_PROLONGATION (X, Z, NL+1, NL)                    !  v^l := Prolongation(x^{l+1})
      CALL SCARC_VECTOR_SUM (Z, V, SCARC_COARSE_OMEGA, 1.0_EB, NL)       !  z^l := r^l - z^l
 
   ! ---------- Only coarse grid preconditioner : COARSE_ONLY
   ! Only solve global coarse Poisson problem
 
   CASE (NSCARC_TWOLEVEL_COARSE_ONLY)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !  Use r^l as right hand side for preconditioner
      DO IL = NL, NLEVEL_MAX-1                                    !  successively restrict to coarser levels up to coarsest
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !  b^{l+1} := Restriction(b^l)
      ENDDO
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !  y^L := x^L
      DO IL = NLEVEL_MAX-1, NL, -1                                !  successively interpolate to finer levels up to finest
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !  y^l := Prolongation(y^{l+1})
      ENDDO
      CALL SCARC_VECTOR_COPY (Y, V, 1.0_EB, NL)                   !  v^l := y^l

   ! ---------- Twolevel preconditioning using meanvalues in x-direction : XMEAN_ADD
   ! Global 1D-meanvalue problem and local fine Poisson problems are based on same defect
 
   CASE (NSCARC_TWOLEVEL_XMEAN_ADD)

      CALL SCARC_VECTOR_COPY (R, Z, 1.0_EB, NL)                   !  z := r                   (copy defect)
      CALL SCARC_PRECON_XMEAN (Z, NL)                             !  z := A1d^{-1} * z_mean   (solve global 1D mean problem)
      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !  v := r 
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !  v := Relax(v)            (solve local Poisson problems)
      CALL SCARC_VECTOR_SUM (Z, V, SCARC_COARSE_OMEGA, 1.0_EB, NL) !                          (relax influence of 1D-meanvalues)

   ! ---------- Twolevel preconditioning using meanvalues in x-direction : XMEAN_MUL1
   ! First solve global 1D-meanvalue problem, build new defect, then solve local fine Poisson problems
 
   CASE (NSCARC_TWOLEVEL_XMEAN_MUL1)

      CALL SCARC_VECTOR_COPY (R, Y, 1.0_EB, NL)                   !  y := r                   (copy residual)
      CALL SCARC_PRECON_XMEAN(Y, NL)                              !  y := A1d^{-1} * y_mean   (solve global 1D mean problem)
      CALL SCARC_MATVEC_PRODUCT (Y, Z, NL)                        !  z = A * y
      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !  z = r - z  = r - A * y   (new defect)
      CALL SCARC_VECTOR_COPY (Z, V, 1.0_EB, NL)                   !  v := z
      CALL SCARC_RELAXATION (Z, V, NS+1, NP, NL)                  !  v := Relax(z)            (solve local Poisson problems)
      CALL SCARC_VECTOR_SUM (Y, V, SCARC_COARSE_OMEGA, 1.0_EB, NL) !                          (relax influence of 1D-meanvalues)

   ! ---------- Twolevel preconditioning using meanvalues in x-direction : XMEAN_MUL2
   ! First solve local Poisson problems, build new defect, then solve global 1D-meanvalue problem
 
   CASE (NSCARC_TWOLEVEL_XMEAN_MUL2)

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !  v := r                   (copy residual)
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !  v := Relax(r)            (solve local Poisson problems)
      CALL SCARC_MATVEC_PRODUCT (V, Z, NL)                        !  z := A * v
      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !  z := r - z  = r - A * v  (new defect)
      CALL SCARC_PRECON_XMEAN(Z, NL)                              !  z := A1d^{-1} * z_mean   (solve global 1D mean problem)
      CALL SCARC_VECTOR_SUM (Z, V, SCARC_COARSE_OMEGA, 1.0_EB, NL) !                          (relax influence of 1D-meanvalues)

END SELECT SELECT_PRECON_TYPE

END SUBROUTINE SCARC_PRECONDITIONER

! --------------------------------------------------------------------------------------------------------------
!> \brief Perform smoothing based on specified relaxation method
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SMOOTHER(NTYPE, NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NTYPE, NSTACK, NPARENT, NLEVEL
INTEGER :: NSTATE=0, NS, NP, NL, ITE0
REAL(EB) :: TNOW
LOGICAL :: BMATVEC, BL2NORM, BVERBOSE
 
! ---------- Initialization
 
TNOW = CURRENT_TIME()
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS = NSTACK
NP = NPARENT
NL = NLEVEL

CALL SCARC_SETUP_SCOPE(NS, NP)
 
! Calculate initial defect on l2-norm on level NL (only if BMATVEC and Bl2NORM are set to .TRUE.)
! Because initial vector in MG is set to zero, this defect corresponds to F
 
ITE0 = 0
BVERBOSE = .TRUE.
IF (BVERBOSE) THEN
   BL2NORM  = .TRUE.
   BMATVEC  = .TRUE.
ELSE
   BL2NORM  = .FALSE.
   BMATVEC  = .FALSE.
ENDIF

IF (BMATVEC) THEN
   CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                 !  v := A*x
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                    !  v := b - v   
ENDIF

IF (BL2NORM) THEN
   RESIN = SCARC_L2NORM (V, NL)                                         !  resin := ||v||
ELSE
   RESIN = SCARC_RESIDUAL
ENDIF
IF (BVERBOSE) NSTATE = SCARC_CONVERGENCE_STATE(NTYPE, NS, NL)

 
! ---------- Smoothing loop - only temporarily
 
!IF (NTYPE == NSCARC_CYCLING_PRESMOOTH) THEN
!   NIT = SCARC_MULTIGRID_PRESMOOTH
!ELSE
!   NIT = SCARC_MULTIGRID_POSTSMOOTH
!ENDIF

SMOOTH_LOOP: DO ITE0=1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

#ifdef WITH_MKL
   IF (TYPE_SMOOTH == NSCARC_RELAX_MKL) THEN
      CALL SCARC_VECTOR_COPY(V, Z, 1.0_EB, NL)                          !  use additional auxiliary vector Z
      CALL SCARC_RELAXATION (Z, V, NS, NP, NL)                          !  v := Relax(z)
   ELSE
      CALL SCARC_RELAXATION (V, V, NS, NP, NL)                          !  v := Relax(v)
   ENDIF
#else
   CALL SCARC_RELAXATION (V, V, NS, NP, NL)                             !  v := Relax(v)
#endif


   CALL SCARC_VECTOR_SUM (V, X, OMEGA, 1.0_EB, NL)                      !  x := omega * v + x
   CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                 !  v := A*x
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                    !  v := b - v


   IF (BL2NORM) THEN
      RES = SCARC_L2NORM (V, NL)                                        !  res := ||v||
      IF (BVERBOSE) THEN
         NSTATE = SCARC_CONVERGENCE_STATE(NTYPE, NS, NL)
         IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP
      ENDIF
   ENDIF

ENDDO SMOOTH_LOOP

CALL SCARC_RELEASE_SCOPE(NS, NP)

CPU(MY_RANK)%SMOOTHER = CPU(MY_RANK)%SMOOTHER + CURRENT_TIME() - TNOW
END SUBROUTINE SCARC_SMOOTHER


! --------------------------------------------------------------------------------------------------------------
!> \brief Perform preconditioning based on requested local solvers
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELAXATION (NV1, NV2, NS, NP, NL)
USE SCARC_POINTERS, ONLY: L, G, A, AB, FFT, V1, V2, &
                          SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_VECTOR_FB, &
                          SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_BMATRIX
#ifdef WITH_MKL
USE SCARC_POINTERS, ONLY: AS, MKL, V1_FB, V2_FB
#endif
USE SCARC_UTILITIES, ONLY: SET_MATRIX_TYPE
USE POIS, ONLY: H2CZSS, H3CZSS
REAL(EB) :: AUX, OMEGA_SSOR = 1.5_EB
REAL (EB) :: TNOW
INTEGER, INTENT(IN) :: NV1, NV2, NS, NP, NL
INTEGER :: NM, NP0, IC, JC, ICOL, ITYPE, IDIAG, IPTR, INCR, IOR0, IC0, IY, IZ

TNOW = CURRENT_TIME()
ITYPE = STACK(NS-1)%SOLVER%TYPE_RELAX


IF ((IS_AMG .OR. IS_CG_AMG) .AND. NL > NLEVEL_MIN .AND. ITYPE == NSCARC_RELAX_FFT) ITYPE = NSCARC_RELAX_SSOR

SELECT CASE (ITYPE)

   ! --------- Preconditioning by blockwise Jacobi
 
   CASE (NSCARC_RELAX_JAC)

      JACOBI_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                    
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         SELECT CASE (SET_MATRIX_TYPE(NL))
            
            ! ---------- Matrix in compact storage technique
 
            CASE (NSCARC_MATRIX_COMPACT)

               IF (IS_LAPLACE) THEN
                  A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE)
               ELSE
                  A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
               ENDIF
               !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  V2(IC) = V2(IC) / A%VAL(A%ROW(IC))
               ENDDO
               !$OMP END PARALLEL DO

            ! ---------- Matrix in bandwise storage technique
 
            CASE (NSCARC_MATRIX_BANDWISE)

               AB => SCARC_POINT_TO_BMATRIX (NSCARC_MATRIX_POISSON)
               !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  V2(IC) = V2(IC) / AB%VAL(IC, AB%POS(0))
               ENDDO
               !$OMP END PARALLEL DO

         END SELECT 

      ENDDO JACOBI_MESHES_LOOP

 
   ! --------- Preconditioning by blockwise SSOR
 
   CASE (NSCARC_RELAX_SSOR)

      SSOR_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                    
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         SELECT CASE (SET_MATRIX_TYPE(NL))

            ! ---------- Matrix in compact storage technique
 
            CASE (NSCARC_MATRIX_COMPACT)

               IF (IS_LAPLACE) THEN
                  A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_LAPLACE)
               ELSE
                  A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
               ENDIF

               IF (NL == NLEVEL_MIN) THEN
               SSOR_FORWARD_COMPACT_LOOP: DO IC = 1, G%NC                                  ! forward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     IF (A%COL(ICOL) >= IC) CYCLE                                          ! only process lower diags
                     IF (A%COL(ICOL) <= G%NC) THEN
                        AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))  ! ignore overlaps
                     ENDIF
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / A%VAL(A%ROW(IC))
               ENDDO SSOR_FORWARD_COMPACT_LOOP

               SSOR_BACKWARD_COMPACT_LOOP: DO IC = G%NC-1, 1, -1                           ! backward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     IF (A%COL(ICOL) <= IC) CYCLE                                          ! only process upper diags
                     IF (A%COL(ICOL) <= G%NC) THEN                                         ! ignore overlaps
                        AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))
                     ENDIF
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / A%VAL(A%ROW(IC))
               ENDDO SSOR_BACKWARD_COMPACT_LOOP

               ELSE

               SSOR_FORWARD_COMPACT_LOOP_COARSE: DO IC = 1, G%NC                           ! forward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     IF (A%COL(ICOL) >= IC .OR. A%COL(ICOL) == 0) CYCLE                    ! only process lower diags
                     IF (A%COL(ICOL) <= G%NC) THEN
                        AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))  ! ignore overlaps
                     ENDIF
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / A%VAL(A%ROW(IC))
               ENDDO SSOR_FORWARD_COMPACT_LOOP_COARSE

               SSOR_BACKWARD_COMPACT_LOOP_COARSE: DO IC = G%NC-1, 1, -1                     ! backward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     IF (A%COL(ICOL) <= IC .OR. A%COL(ICOL) == 0) CYCLE                     ! only process upper diags
                     IF (A%COL(ICOL) <= G%NC) THEN                                          ! ignore overlaps
                        AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))
                     ENDIF
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / A%VAL(A%ROW(IC))
               ENDDO SSOR_BACKWARD_COMPACT_LOOP_COARSE

               ENDIF

 
         ! ---------- Matrix in bandwise storage technique
 
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => SCARC_POINT_TO_BMATRIX (NSCARC_MATRIX_POISSON)

 
            ! 2D version
 
            IF (TWO_D) THEN

               SSOR_FORWARD_BANDWISE_2D_LOOP: DO IC = 1, G%NC                 ! forward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = 1, 3, 2                                           ! only process lower x- and z-diag
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= 1 .AND. JC <= G%NC) AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_FORWARD_BANDWISE_2D_LOOP

               SSOR_BACKWARD_BANDWISE_2D_LOOP: DO IC = G%NC-1, 1, -1          ! backward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = -1, -3, -2                                        ! only process upper x- and z-diag
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC <= G%NC) THEN
                        AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                     ENDIF
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_BACKWARD_BANDWISE_2D_LOOP

 
            ! 3D version
 
            ELSE

               SSOR_FORWARD_BANDWISE_3D_LOOP: DO IC = 1, G%NC                 ! forward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = 1, 3                                              ! only process lower diags
                     IF (AB%POS(IOR0) == 0) CYCLE                             ! no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= 1 .AND. JC <= G%NC) AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_FORWARD_BANDWISE_3D_LOOP

               SSOR_BACKWARD_BANDWISE_3D_LOOP: DO IC = G%NC-1, 1, -1          ! backward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = -1, -3, -1                                        ! only process upper diags
                     IF (AB%POS(IOR0) == 0) CYCLE                             ! no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= IC .AND. JC <= G%NC) AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_BACKWARD_BANDWISE_3D_LOOP

            ENDIF

         END SELECT 

      ENDDO SSOR_MESHES_LOOP

 
   ! --------- Preconditioning by Jacobi in matrix form
 
   CASE (NSCARC_RELAX_MJAC)

      MJAC_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                    

         V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         SELECT CASE(TYPE_MATRIX)

            ! ------------ Matrix in compact storage technique

            CASE (NSCARC_MATRIX_COMPACT)
               A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
               CALL SCARC_SCALING_VARIABLE(G%NC, A%RELAX, V1, V2)

            ! ------------ Matrix in bandwise storage technique

            CASE (NSCARC_MATRIX_BANDWISE)
               AB => SCARC_POINT_TO_BMATRIX (NSCARC_MATRIX_POISSON)
               CALL SCARC_SCALING_VARIABLE(G%NC, AB%RELAXD, V1, V2)

         END SELECT

      ENDDO MJAC_MESHES_LOOP

 
   ! --------- Preconditioning by different matrix-form preconditioners
   ! in all cases the preconditioner is given as separate matrix which is based
   ! on the same storage technique as the matrix AC itself;
   ! two tridiagonal systems have to be solved
   ! V1 contains the RHS to be solved for, V2 will contain the solution
 
   CASE (NSCARC_RELAX_MGS, NSCARC_RELAX_MSGS, NSCARC_RELAX_MSOR, NSCARC_RELAX_MSSOR, NSCARC_RELAX_ILU)

      LU_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                    

         V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)
      
         SELECT CASE(TYPE_MATRIX)

            ! ------------ Matrix in compact storage technique

            CASE (NSCARC_MATRIX_COMPACT)

               A => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON)
      
               ! Forward solve:   Solve V2 = L^-1 V1
               ! Compute sol(i) = rhs(i) - sum L(i,j) x sol(j)

               DO IC = 1, G%NC
                  V2(IC) = V1(IC)
                  DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                     JC = A%COL(IPTR)
                     IF (JC >= IC) CYCLE
                     V2(IC) = V2(IC) - A%RELAX(IPTR) * V2(JC)
                  ENDDO
               ENDDO
      
               ! If preconditioner is not symmetric, upper matrix U is zero and nothing has to be solved

               IF (ITYPE == NSCARC_RELAX_MGS .OR. ITYPE == NSCARC_RELAX_MSOR) CYCLE
      
               ! Backward solve : Compute sol: inv(U) sol

               DO IC = G%NC, 1, -1
      
                  DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                     JC = A%COL(IPTR)
                     IF (JC <= IC) CYCLE
                     V2(IC) = V2(IC) - A%RELAX(IPTR) * V2(JC)
                  ENDDO
      
                  ! Compute sol(i) = sol(i)/U(i,i)

                  IDIAG = A%ROW(IC)
                  V2(IC) = V2(IC)/A%RELAX(IDIAG)
      
               ENDDO
      

            ! ---------- Matrix in bandwise storage technique
 
            CASE (NSCARC_MATRIX_BANDWISE)

               AB => SCARC_POINT_TO_BMATRIX (NSCARC_MATRIX_POISSON)
      
               IF (TWO_D) THEN
                  INCR = -2
               ELSE 
                  INCR = -1
               ENDIF
               
               ! Forward solve:   V2 = L^-1 V1
               ! Compute sol(i) = rhs(i) - sum L(i,j) x sol(j)

               !!$OMP PARALLEL DO PRIVATE(IC, JC, IOR0) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  V2(IC) = V1(IC)
                  DO IOR0 = 3, 1, INCR
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC <= 0) CYCLE
                     V2(IC) = V2(IC) - AB%RELAX(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
               ENDDO
               !!$OMP END PARALLEL DO 
      
               ! If preconditioner is not symmetric, upper matrix U is zero and nothing has to be solved

               IF (ITYPE == NSCARC_RELAX_MGS .OR. ITYPE == NSCARC_RELAX_MSOR) CYCLE
      
               ! Backward solve
               ! Compute sol: inv(U) sol

               !!$OMP PARALLEL DO PRIVATE(IC, JC, IOR0) SCHEDULE(STATIC)
               DO IC = G%NC, 1, -1
      
                  DO IOR0 = -1, -3, INCR
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC > G%NC) CYCLE
                     V2(IC) = V2(IC) - AB%RELAX(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
      
                  ! Compute sol(i) = sol(i)/U(i,i)
                  V2(IC) = V2(IC)/AB%RELAX(IC, AB%POS(0))
               ENDDO
               !!$OMP END PARALLEL DO 
      
         END SELECT
      
      ENDDO LU_MESHES_LOOP
      
   ! --------- Preconditioning by blockwise Geometric Multigrid
 
   CASE (NSCARC_RELAX_GMG)

      NP0 = NP               ! Only dummy command until multigrid is used again
      !CALL SCARC_METHOD_MULTIGRID (NS, NP, NLEVEL_MIN)

 
   ! --------- Preconditioning by blockwise FFT based on Crayfishpak
 
   CASE (NSCARC_RELAX_FFT)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                    
         FFT => L%FFT

         V1  => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2  => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         ! Feed corresponding right hand sides for FFT

         !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
         DO IC = 1, G%NC
            FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = V1(IC)
         ENDDO
         !$OMP END PARALLEL DO

         ! Call corresponding FFT solver
 
         IF (TWO_D) THEN
            CALL H2CZSS (FFT%BXS,  FFT%BXF, FFT%BZS, FFT%BZF, FFT%ITRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ELSE
            CALL H3CZSS (FFT%BXS,  FFT%BXF, FFT%BYS, FFT%BYF, FFT%BZS, FFT%BZF, FFT%ITRN, FFT%JTRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ENDIF

         ! Extract computed solution which is contained in FFT%PRHS
 
         !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
         DO IC = 1, G%NC
            V2(IC) = FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) 
         ENDDO
         !$OMP END PARALLEL DO 

      ENDDO

 
   ! --------- Preconditioning by blockwise overlapping FFT based on Crayfishpak 
   !           still test-version for tunnel-shaped geometries of type Mx1
 
   CASE (NSCARC_RELAX_FFTO)

      ! Exchange overlapping parts

      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NV1, NL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                    
         FFT => L%FFT

         V1  => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2  => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         ! Feed corresponding right hand sides for FFT
 
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

         ! Call corresponding FFT solver
 
         IF (TWO_D) THEN
            CALL H2CZSS (FFT%BXS,  FFT%BXF, FFT%BZS, FFT%BZF, FFT%ITRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ELSE
            CALL H3CZSS (FFT%BXS,  FFT%BXF, FFT%BYS, FFT%BYF, FFT%BZS, FFT%BZF, FFT%ITRN, FFT%JTRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ENDIF

         ! Extract computed solution which is contained in FFT%PRHS
 
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

      ! Exchange overlapping parts

      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_MEAN, NV2, NL)


#ifdef WITH_MKL
 

   ! --------- Preconditioning by LU-decomposition
 
   CASE (NSCARC_RELAX_MKL)

      ! Preconditioning by Cluster Sparse Solver from MKL
 
      MKL_SCOPE_IF: IF (STACK(NS)%SOLVER%TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL) THEN

         MKL_SCOPE_GLOBAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

            CALL SCARC_POINT_TO_GRID (NM, NL)                                    
            MKL => L%MKL
            AS => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_SYM)

            MKL%PHASE  = 33                            ! only solving

            V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
            V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)


            IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

               V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV1)
               V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV2)

               V1_FB(1:G%NC) = REAL(V1(1:G%NC), FB)
               V2_FB(1:G%NC) = 0.0_FB

               CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                            AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                            MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)

               V2(1:G%NC) = REAL(V2_FB(1:G%NC), EB)

            ELSE
               CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                            AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                            MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)
            ENDIF
            IF (MKL%ERROR /= 0) CALL SCARC_ERROR(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

         ENDDO MKL_SCOPE_GLOBAL_LOOP

 
      ! Preconditioning by Pardiso Solver from MKL
 
      ELSE MKL_SCOPE_IF

         MKL_SCOPE_LOCAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

            CALL SCARC_POINT_TO_GRID (NM, NL)                                    
            MKL => L%MKL
            AS => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_SYM)

            MKL%PHASE  = 33                            ! only solving

            V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
            V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

            IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

               V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV1)
               V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV2)

               V1_FB(1:G%NC) = REAL(V1(1:G%NC), FB)
               V2_FB(1:G%NC) = 0.0_FB
   
               CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                              AS%VAL_FB, AS%ROW, AS%COL, &
                              MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

               V2(1:G%NC) = REAL(V2_FB(1:G%NC), EB)
            ELSE

               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

               CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                              AS%VAL, AS%ROW, AS%COL, &
                              MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1, V2, MKL%ERROR)

            ENDIF
            IF (MKL%ERROR /= 0) CALL SCARC_ERROR(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

         ENDDO MKL_SCOPE_LOCAL_LOOP

      ENDIF MKL_SCOPE_IF

      
   ! --------- Preconditioning by optimized use of FFT or PARDISO, depending on structure of mesh
 
   CASE (NSCARC_RELAX_OPTIMIZED)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                    

         ! If mesh happens to contain obstructions, PARDISO preconditioner is used
          
         IF (L%HAS_OBSTRUCTIONS) THEN

            MKL => L%MKL
            AS => SCARC_POINT_TO_CMATRIX (NSCARC_MATRIX_POISSON_SYM)

            MKL%PHASE  = 33                            ! only solving

            V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
            V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

            IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

               V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV1)
               V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV2)

               V1_FB(1:G%NC) = REAL(V1(1:G%NC), FB)
               V2_FB(1:G%NC) = 0.0_FB
   
               CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                              AS%VAL_FB, AS%ROW, AS%COL, &
                              MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

               V2(1:G%NC) = REAL(V2_FB(1:G%NC), EB)
            ELSE

               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

               CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                              AS%VAL, AS%ROW, AS%COL, &
                              MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1, V2, MKL%ERROR)

            ENDIF
            IF (MKL%ERROR /= 0) CALL SCARC_ERROR(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

         ! If mesh happens to be purely structured, FFT preconditioner is used

         ELSE

            FFT => L%FFT
   
            V1  => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
            V2  => SCARC_POINT_TO_VECTOR(NM, NL, NV2)
   
            ! Feed corresponding right hand sides for FFT
   
            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = V1(IC)
            ENDDO
            !$OMP END PARALLEL DO
   
            ! Call corresponding FFT solver
    
            IF (TWO_D) THEN
               CALL H2CZSS (FFT%BXS,  FFT%BXF, FFT%BZS, FFT%BZF, FFT%ITRN, &
                            FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
            ELSE
               CALL H3CZSS (FFT%BXS,  FFT%BXF, FFT%BYS, FFT%BYF, FFT%BZS, FFT%BZF, FFT%ITRN, FFT%JTRN, &
                            FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
            ENDIF
   
            ! Extract computed solution which is contained in FFT%PRHS
    
            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               V2(IC) = FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) 
            ENDDO
            !$OMP END PARALLEL DO 
   
         ENDIF
      ENDDO
#endif

END SELECT

CPU(MY_RANK)%RELAXATION =CPU(MY_RANK)%RELAXATION+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_RELAXATION


! --------------------------------------------------------------------------------------------------------------
!> \brief Twolevel relaxation by meanvalues in x-direction
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRECON_XMEAN (NV2, NL)
USE SCARC_POINTERS, ONLY: M, L, G, V2, PRE, SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV2, NL
INTEGER :: II, IX, IY, IZ, IC, I, NM
REAL(EB) :: XMEAN, VAL
 
MEAN1D_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                    

   PRE => L%PRECON
   V2  => SCARC_POINT_TO_VECTOR(NM, NL, NV2)
      
   DO IX = 1, L%NX
      II = NX_OFFSET(NM) + IX  
      PRE%RHS(II) = 0._EB
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            IF (L%IS_SOLID(IX, IY, IZ)) CYCLE
            IC = G%CELL_NUMBER(IX, IY, IZ)
            XMEAN = V2(IC)*M%DY(IY)*M%DZ(IZ)
            PRE%RHS(II) = PRE%RHS(II) + V2(IC)*M%DY(IY)*M%DZ(IZ)
         ENDDO
      ENDDO
      PRE%RHS(II) = PRE%RHS(II)/((M%YF-M%YS)*(M%ZF-M%ZS)) 
   ENDDO

   ! If not MPI process 0 send RHS components to process 0

   IF (MY_RANK>0) THEN  
      
      CALL MPI_GATHERV(PRE%RHS(DISPLS_TP(MY_RANK)+1),COUNTS_TP(MY_RANK),MPI_DOUBLE_PRECISION,PRE%RHS,COUNTS_TP,DISPLS_TP,&
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
      
   ! If MPI process 0 receive RHS components form slaves and solve tridiagonal Poisson system

   ELSE  
      
      CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,PRE%RHS,COUNTS_TP,DISPLS_TP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

      PRE%AUX = PRE%DIAG
      TRIDIAGONAL_SOLVER_1: DO I=2,NX_GLOBAL
         VAL = PRE%LOW(I)/PRE%AUX(I-1)
         PRE%AUX(I) = PRE%AUX(I) - VAL*PRE%UP(I-1)
         PRE%RHS(I) = PRE%RHS(I) - VAL*PRE%RHS(I-1)
      ENDDO TRIDIAGONAL_SOLVER_1
      PRE%RHS(NX_GLOBAL)  = PRE%RHS(NX_GLOBAL)/PRE%AUX(NX_GLOBAL)
      TRIDIAGONAL_SOLVER_2: DO I=NX_GLOBAL-1,1,-1
         PRE%RHS(I) = (PRE%RHS(I) - PRE%UP(I)*PRE%RHS(I+1))/PRE%AUX(I)
      ENDDO TRIDIAGONAL_SOLVER_2

   ENDIF
      
   ! The solution to the tri-diagonal linear system is PRE%RHS. Broadcast this to all the MPI processes.

   CALL MPI_BCAST(PRE%RHS,NX_GLOBAL,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

   DO IX = 1, L%NX
      DO IY = 1, L%NY
         DO IZ = 1, L%NZ
            IC = G%CELL_NUMBER(IX, IY, IZ)
            V2(IC) = PRE%RHS(NX_OFFSET(NM)+IX)
         ENDDO
      ENDDO
   ENDDO

ENDDO MEAN1D_MESHES_LOOP
 
END SUBROUTINE SCARC_PRECON_XMEAN

END MODULE SCARC_METHODS



!=======================================================================================================================
!
! MODULE SCRC
!
!  Basic setup and call of different variants of ScaRC/UScaRC 
!
!=======================================================================================================================
MODULE SCRC

USE PRECISION_PARAMETERS, ONLY: EB
USE GLOBAL_CONSTANTS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE SCARC_CONSTANTS
USE SCARC_METHODS

IMPLICIT NONE (TYPE,EXTERNAL)

PUBLIC :: SCARC_SETUP, SCARC_SOLVER

CONTAINS

! --------------------------------------------------------------------------------------------------------------
!> \brief Initialize ScaRC structures based on SCARC-input parameters from &PRES namelist
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP
USE SCARC_PARSER
USE SCARC_GRIDS
USE SCARC_STORAGE, ONLY: SCARC_SETUP_STORAGE
#ifdef WITH_SCARC_DEBUG
USE SCARC_MESSAGES, ONLY: SCARC_SETUP_MESSAGES
#endif
USE SCARC_MPI, ONLY: SCARC_SETUP_EXCHANGES, SCARC_SETUP_GLOBALS
USE SCARC_CPU, ONLY: SCARC_SETUP_CPU
USE SCARC_STACK, ONLY: SCARC_SETUP_VECTORS
USE SCARC_MATRICES, ONLY: SCARC_SETUP_SYSTEMS
#ifdef WITH_MKL
USE SCARC_MKL, ONLY: SCARC_SETUP_MKL_ENVIRONMENT
#endif
#ifdef WITH_SCARC_AMG
USE SCARC_AMG, ONLY: SCARC_SETUP_AMG_ENVIRONMENT
#endif
#ifdef WITH_SCARC_POSTPROCESSING
USE SCARC_POSTPROCESSING, ONLY: SCARC_SETUP_PRESSURE
#endif
USE SCARC_METHODS, ONLY: SCARC_SETUP_KRYLOV_ENVIRONMENT, SCARC_SETUP_MULTIGRID_ENVIRONMENT, &
                         SCARC_SETUP_MGM_ENVIRONMENT
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

! Setup mechanisms for own memory management, different messaging services and CPU-time measurements
 
CALL SCARC_SETUP_STORAGE
CALL SCARC_SETUP_MESSAGES
CALL SCARC_SETUP_CPU

! Parse ScaRC related input parameters in &PRES namelist

CALL SCARC_PARSE_INPUT                      ; IF (STOP_STATUS==SETUP_STOP) RETURN

! Setup different components of ScaRC solver
 
CALL SCARC_SETUP_LEVELS                               ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_BASICS                               ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_GRIDS                                ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_GLOBALS                              ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_NEIGHBORS                            ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_FACES                                ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_SUBDIVISION                          ; IF (STOP_STATUS==SETUP_STOP) RETURN

! Setup wall information according to specified discretization type/method
 
IF (HAS_MULTIPLE_GRIDS) THEN
   CALL SCARC_SETUP_WALLS (NSCARC_GRID_STRUCTURED)    ; IF (STOP_STATUS==SETUP_STOP) RETURN
   CALL SCARC_SETUP_WALLS (NSCARC_GRID_UNSTRUCTURED)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
ELSE
   CALL SCARC_SETUP_WALLS (TYPE_GRID)                 ; IF (STOP_STATUS==SETUP_STOP) RETURN
ENDIF

! Setup information for data exchanges and matrix systems
 
CALL SCARC_SETUP_EXCHANGES                            ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_SYSTEMS                              ; IF (STOP_STATUS==SETUP_STOP) RETURN

! Setup information for algebraic multigrid if needed as preconditioner or main solver

#ifdef WITH_SCARC_AMG
IF (HAS_AMG_LEVELS) CALL SCARC_SETUP_AMG_ENVIRONMENT          
#endif

! Setup environment for requested solver

SELECT_METHOD: SELECT CASE(TYPE_METHOD)

   ! Global Krylov method
   CASE (NSCARC_METHOD_KRYLOV)
      CALL SCARC_SETUP_KRYLOV_ENVIRONMENT

   ! Global multigrid method
    CASE (NSCARC_METHOD_MULTIGRID)
       CALL SCARC_SETUP_MULTIGRID_ENVIRONMENT

   ! Global McKeeney-Greengard-Mayo method
   CASE (NSCARC_METHOD_MGM)
       CALL SCARC_SETUP_MGM_ENVIRONMENT

#ifdef WITH_MKL
   ! Global Intel-MKL related method
   CASE (NSCARC_METHOD_LU)
       CALL SCARC_SETUP_MKL_ENVIRONMENT
#endif
END SELECT SELECT_METHOD

! Setup vector structures for requested solver

CALL SCARC_SETUP_VECTORS                              ; IF (STOP_STATUS==SETUP_STOP) RETURN
 
CPU(MY_RANK)%SETUP   = CPU(MY_RANK)%SETUP   + CURRENT_TIME() - TNOW
CPU(MY_RANK)%OVERALL = CPU(MY_RANK)%OVERALL + CURRENT_TIME() - TNOW

END SUBROUTINE SCARC_SETUP


! --------------------------------------------------------------------------------------------------------------------
!> \brief Call of requested ScaRC solver 
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SOLVER(DT_CURRENT)
USE SCARC_CONVERGENCE
USE SCARC_METHODS, ONLY: SCARC_METHOD_KRYLOV, SCARC_METHOD_MULTIGRID, SCARC_METHOD_MGM
#ifdef WITH_MKL
USE SCARC_METHODS, ONLY: SCARC_METHOD_MKL
#endif
REAL (EB), INTENT(IN) :: DT_CURRENT
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

CALL SCARC_SET_ITERATION_STATE (DT_CURRENT)


SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   CASE (NSCARC_METHOD_KRYLOV)
      CALL SCARC_METHOD_KRYLOV (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
   
   CASE (NSCARC_METHOD_MULTIGRID)
      CALL SCARC_METHOD_MULTIGRID(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
   
   CASE (NSCARC_METHOD_MGM)
      CALL SCARC_METHOD_MGM(NSCARC_STACK_ROOT)

#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LU)
      CALL SCARC_METHOD_MKL(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
#endif
   
END SELECT SELECT_METHOD


IF (STOP_STATUS==SETUP_STOP) RETURN

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
CPU(MY_RANK)%SOLVER =CPU(MY_RANK)%SOLVER+CURRENT_TIME()-TNOW
CPU(MY_RANK)%OVERALL=CPU(MY_RANK)%OVERALL+CURRENT_TIME()-TNOW

END SUBROUTINE SCARC_SOLVER

END MODULE SCRC


