!> \brief Intel Math Kernel Library linear solver interface

! MKL Solver, defined in makefile:
#ifdef WITH_MKL
!!! #include "mkl_pardiso.f90"
#define __MKL_PARDISO_F90
MODULE MKL_PARDISO_PRIVATE
   TYPE MKL_PARDISO_HANDLE
      INTEGER(KIND=8) DUMMY
   END TYPE
   INTEGER, PARAMETER :: PARDISO_OOC_FILE_NAME = 1
END MODULE MKL_PARDISO_PRIVATE
MODULE MKL_PARDISO
   USE MKL_PARDISO_PRIVATE
   INTERFACE PARDISO
      SUBROUTINE PARDISO_D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
         USE MKL_PARDISO_PRIVATE
         TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(INOUT) :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(KIND=8),     INTENT(IN)    :: A(*)
         REAL(KIND=8),     INTENT(INOUT) :: B(*)
         REAL(KIND=8),     INTENT(OUT)   :: X(*)
      END SUBROUTINE PARDISO_D

      SUBROUTINE PARDISO_S( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
         USE MKL_PARDISO_PRIVATE
         TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(INOUT) :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(KIND=4),     INTENT(IN)    :: A(*)
         REAL(KIND=4),     INTENT(INOUT) :: B(*)
         REAL(KIND=4),     INTENT(OUT)   :: X(*)
      END SUBROUTINE PARDISO_S

      SUBROUTINE PARDISO_D_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
         USE MKL_PARDISO_PRIVATE
         TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(INOUT) :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(KIND=8),     INTENT(IN)    :: A(*)
         REAL(KIND=8),     INTENT(INOUT) :: B(N,*)
         REAL(KIND=8),     INTENT(OUT)   :: X(N,*)
      END SUBROUTINE PARDISO_D_2D

      SUBROUTINE PARDISO_S_2D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
         USE MKL_PARDISO_PRIVATE
         TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(INOUT) :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(KIND=4),     INTENT(IN)    :: A(*)
         REAL(KIND=4),     INTENT(INOUT) :: B(N,*)
         REAL(KIND=4),     INTENT(OUT)   :: X(N,*)
      END SUBROUTINE PARDISO_S_2D
   END INTERFACE
END MODULE MKL_PARDISO

!!! #include "mkl_cluster_sparse_solver.f90"
#define __MKL_CLUSTER_SPARSE_SOLVER_F90
MODULE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
   TYPE MKL_CLUSTER_SPARSE_SOLVER_HANDLE
      INTEGER(KIND=8) DUMMY
   END TYPE
END MODULE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
MODULE MKL_CLUSTER_SPARSE_SOLVER
   USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
   INTERFACE CLUSTER_SPARSE_SOLVER
      SUBROUTINE CLUSTER_SPARSE_SOLVER_D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
         USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
         USE MPI_F08
         TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(IN)    :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(KIND=8),     INTENT(IN)    :: A(*)
         REAL(KIND=8),     INTENT(INOUT) :: B(*)
         REAL(KIND=8),     INTENT(OUT)   :: X(*)
         TYPE(MPI_COMM),   INTENT(IN)    :: COMM
      END SUBROUTINE CLUSTER_SPARSE_SOLVER_D

      SUBROUTINE CLUSTER_SPARSE_SOLVER_S(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
         USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
         USE MPI_F08
         TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(IN)    :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(KIND=4),     INTENT(IN)    :: A(*)
         REAL(KIND=4),     INTENT(INOUT) :: B(*)
         REAL(KIND=4),     INTENT(OUT)   :: X(*)
         TYPE(MPI_COMM),   INTENT(IN)    :: COMM
      END SUBROUTINE CLUSTER_SPARSE_SOLVER_S

      SUBROUTINE CLUSTER_SPARSE_SOLVER_D_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
         USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
         USE MPI_F08
         TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(IN)    :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(KIND=8),     INTENT(IN)    :: A(*)
         REAL(KIND=8),     INTENT(INOUT) :: B(N,*)
         REAL(KIND=8),     INTENT(OUT)   :: X(N,*)
         TYPE(MPI_COMM),   INTENT(IN)    :: COMM
      END SUBROUTINE CLUSTER_SPARSE_SOLVER_D_2D

      SUBROUTINE CLUSTER_SPARSE_SOLVER_S_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
         USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
         USE MPI_F08
         TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(IN)    :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(KIND=4),     INTENT(IN)    :: A(*)
         REAL(KIND=4),     INTENT(INOUT) :: B(N,*)
         REAL(KIND=4),     INTENT(OUT)   :: X(N,*)
         TYPE(MPI_COMM),   INTENT(IN)    :: COMM
      END SUBROUTINE CLUSTER_SPARSE_SOLVER_S_2D
   END INTERFACE
END MODULE MKL_CLUSTER_SPARSE_SOLVER
#endif /* WITH_MKL */

#ifdef WITH_HYPRE
MODULE HYPRE_INTERFACE
   IMPLICIT NONE(TYPE,EXTERNAL)
   INCLUDE 'HYPREf.h' ! defines the integer parameters
   INTEGER :: HYPRE_IERR = 0
   INTEGER, PARAMETER :: HYPRE_SOLVER_ID = 1                 ! Preconditioned Conjugate Gradient (PCG) solver
   INTEGER, PARAMETER :: HYPRE_PRECOND_ID = 2                ! Algebraic Multi-Grid (AMG) preconditioner
   INTEGER, PARAMETER :: HYPRE_SOLVER_MAXIT = 1000           ! Max iteratations of PCG solver
   REAL(KIND=8), PARAMETER :: HYPRE_SOLVER_TOL = 1.D-12      ! Solver tolerance
   INTEGER, PARAMETER :: HYPRE_SOLVER_SETTWONORM = 1         ! 0=use L_Infty norm (max error) for convergence, 1=use L2 norm
   INTEGER, PARAMETER :: HYPRE_SOLVER_SETPRINTLEVEL = 0      ! 0=no output, 1=minimal, 2=verbose
   INTEGER, PARAMETER :: HYPRE_SOLVER_SETLOGGING = 0         ! 0=no logging, 1=solver stores intermediate info, norms, etc.
   INTEGER, PARAMETER :: HYPRE_PRECOND_SETPRINTLEVEL = 0     ! 0=no output, 1=minimal, 2=verbose
   INTEGER, PARAMETER :: HYPRE_PRECOND_COARSENINGTYPE = 6    ! 0   CLJP (Cleary-Luby-Jones-Plassmann) parallel coarsening
                                                             ! 1   Classical Ruge-Stüben (RS) coarsening
                                                             ! 3   Modified RS coarsening
                                                             ! 6   Falgout coarsening (a combination of CLJP and RS)
                                                             ! 8   PMIS (Parallel Modified Independent Set) coarsening
                                                             ! 10  HMIS (Hybrid Modified Independent Set) coarsening
                                                             ! 21  Falgout coarsening with aggressive coarsening
   INTEGER, PARAMETER :: HYPRE_PRECOND_SETRELAXTYPE = 6      ! 0   Jacobi (default)
                                                             ! 1   Gauss-Seidel, sequential (very slow in parallel)
                                                             ! 2   Gauss-Seidel, interior points first (parallel variant)
                                                             ! 3   Hybrid Gauss-Seidel or SOR (symmetric in parallel)
                                                             ! 6   L1-scaled Jacobi
                                                             ! 8   L1-scaled hybrid symmetric Gauss-Seidel/SOR
                                                             ! 13  Two-step Jacobi
                                                             ! 16  Chebyshev smoothing (useful for difficult problems)
   INTEGER, PARAMETER :: HYPRE_PRECOND_NUMSWEEPS = 1         ! Number of sweeps on each level of preconditioner
   REAL(KIND=8), PARAMETER :: HYPRE_PRECOND_TOL = 0.D0       ! Preconditioner convergence tolerance
   INTEGER, PARAMETER :: HYPRE_PRECOND_MAXITER = 1           ! Max number of iterations for preconditioner


   TYPE HYPRE_ZM_TYPE
      INTEGER(KIND=8) :: A_H                                 ! Matrix handle
      INTEGER(KIND=8) :: F_H                                 ! RHS handle
      INTEGER(KIND=8) :: X_H                                 ! Solution handle
      INTEGER(KIND=8) :: PARCSR_A_H                          ! Parallel Compressed Sparse Row matrix object
      INTEGER(KIND=8) :: PAR_F_H                             ! RHS object
      INTEGER(KIND=8) :: PAR_X_H                             ! Solution object
      INTEGER(KIND=8) :: SOLVER                              ! Solver handle
      INTEGER(KIND=8) :: PRECOND                             ! Preconditioner handle
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INDICES          ! Row indices of rhs and solution vectors
      INTEGER :: NUM_ITERATIONS                              ! Output number of iterations
      REAL(KIND=8) :: FINAL_RES_NORM                         ! Final residual norm
   END TYPE HYPRE_ZM_TYPE

   INTERFACE

      SUBROUTINE HYPRE_INITIALIZE(IERR)
         INTEGER, INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_INITIALIZE

      SUBROUTINE HYPRE_FINALIZE(IERR)
         INTEGER, INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_FINALIZE

      SUBROUTINE HYPRE_IJMATRIXCREATE(COMM,I1,I2,I3,I4,A,IERR)
         USE MPI_F08
         TYPE(MPI_COMM),  INTENT(IN)  :: COMM
         INTEGER,         INTENT(IN)  :: I1,I2,I3,I4
         INTEGER(KIND=8), INTENT(OUT) :: A
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJMATRIXCREATE

      SUBROUTINE HYPRE_IJMATRIXSETOBJECTTYPE(A,HYPRE_PARCSR_INT,IERR)
         INTEGER(KIND=8), INTENT(IN)  :: A
         INTEGER,         INTENT(IN)  :: HYPRE_PARCSR_INT
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJMATRIXSETOBJECTTYPE

      SUBROUTINE HYPRE_IJMATRIXINITIALIZE(A,IERR)
         INTEGER(KIND=8), INTENT(IN)  :: A
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJMATRIXINITIALIZE

      SUBROUTINE HYPRE_IJMATRIXSETVALUES(A, NROWS, NCOLS, ROWS, COLS, VALUES, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: A
         INTEGER,         INTENT(IN)  :: NROWS
         INTEGER,         INTENT(IN)  :: NCOLS(*)
         INTEGER,         INTENT(IN)  :: ROWS(*)
         INTEGER,         INTENT(IN)  :: COLS(*)
         REAL(KIND=8),    INTENT(IN)  :: VALUES(*)
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJMATRIXSETVALUES

      SUBROUTINE HYPRE_IJMATRIXASSEMBLE(A,IERR)
         INTEGER(KIND=8), INTENT(IN)  :: A
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJMATRIXASSEMBLE

      SUBROUTINE HYPRE_IJMATRIXGETOBJECT(A,PARCSR_A,IERR)
         INTEGER(KIND=8), INTENT(IN)  :: A
         INTEGER(KIND=8), INTENT(OUT) :: PARCSR_A
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJMATRIXGETOBJECT

      SUBROUTINE HYPRE_IJVECTORCREATE(COMM, ILOWER, IUPPER, X, IERR)
         USE MPI_F08
         TYPE(MPI_COMM),  INTENT(IN)  :: COMM
         INTEGER,         INTENT(IN)  :: ILOWER,IUPPER
         INTEGER(KIND=8), INTENT(OUT) :: X
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJVECTORCREATE

      SUBROUTINE HYPRE_IJVECTORSETOBJECTTYPE(X, HYPRE_PARCSR_INT, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: X
         INTEGER,         INTENT(IN)  :: HYPRE_PARCSR_INT
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJVECTORSETOBJECTTYPE

      SUBROUTINE HYPRE_IJVECTORINITIALIZE(X, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: X
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJVECTORINITIALIZE

      SUBROUTINE HYPRE_IJVECTORSETVALUES(X, LOCAL_SIZE, ROWS, VALUES, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: X
         INTEGER,         INTENT(IN)  :: LOCAL_SIZE
         INTEGER,         INTENT(IN)  :: ROWS(*)
         REAL(KIND=8),    INTENT(IN)  :: VALUES(*)
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJVECTORSETVALUES

      SUBROUTINE HYPRE_IJVECTORGETVALUES(X, LOCAL_SIZE, ROWS, VALUES, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: X
         INTEGER,         INTENT(IN)  :: LOCAL_SIZE
         INTEGER,         INTENT(IN)  :: ROWS(*)
         REAL(KIND=8),    INTENT(OUT) :: VALUES(*)
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJVECTORGETVALUES

      SUBROUTINE HYPRE_IJVECTORASSEMBLE(X, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: X
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJVECTORASSEMBLE

      SUBROUTINE HYPRE_IJVECTORGETOBJECT(X, PAR_X, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: X
         INTEGER(KIND=8), INTENT(OUT) :: PAR_X
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_IJVECTORGETOBJECT

      SUBROUTINE HYPRE_PARCSRPCGCREATE(COMM, SOLVER, IERR)
         USE MPI_F08
         TYPE(MPI_COMM),  INTENT(IN)  :: COMM
         INTEGER(KIND=8), INTENT(OUT) :: SOLVER
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGCREATE

      SUBROUTINE HYPRE_PARCSRPCGSETMAXITER(SOLVER, N, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: SOLVER
         INTEGER,         INTENT(IN)  :: N
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGSETMAXITER

      SUBROUTINE HYPRE_PARCSRPCGSETTOL(SOLVER, TOL, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: SOLVER
         REAL(KIND=8),    INTENT(IN)  :: TOL
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGSETTOL

      SUBROUTINE HYPRE_PARCSRPCGSETATOL(SOLVER, TOL, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: SOLVER
         REAL(KIND=8),    INTENT(IN)  :: TOL
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGSETATOL

      SUBROUTINE HYPRE_PARCSRPCGSETTWONORM(SOLVER, IFLAG, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: SOLVER
         INTEGER,         INTENT(IN)  :: IFLAG
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGSETTWONORM

      SUBROUTINE HYPRE_PARCSRPCGSETPRINTLEVEL(SOLVER, IFLAG, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: SOLVER
         INTEGER,         INTENT(IN)  :: IFLAG
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGSETPRINTLEVEL

      SUBROUTINE HYPRE_PARCSRPCGSETLOGGING(SOLVER, IFLAG, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: SOLVER
         INTEGER,         INTENT(IN)  :: IFLAG
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGSETLOGGING

      SUBROUTINE HYPRE_BOOMERAMGCREATE(PRECOND, IERR)
         INTEGER(KIND=8), INTENT(OUT) :: PRECOND
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_BOOMERAMGCREATE

      SUBROUTINE HYPRE_BOOMERAMGSETPRINTLEVEL(PRECOND, IFLAG, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: PRECOND
         INTEGER,         INTENT(IN)  :: IFLAG
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_BOOMERAMGSETPRINTLEVEL

      SUBROUTINE HYPRE_BOOMERAMGSETCOARSENTYPE(PRECOND, IFLAG, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: PRECOND
         INTEGER,         INTENT(IN)  :: IFLAG
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_BOOMERAMGSETCOARSENTYPE

      SUBROUTINE HYPRE_BOOMERAMGSETOLDDEFAULT(PRECOND, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: PRECOND
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_BOOMERAMGSETOLDDEFAULT

      SUBROUTINE HYPRE_BOOMERAMGSETRELAXTYPE(PRECOND, IFLAG, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: PRECOND
         INTEGER,         INTENT(IN)  :: IFLAG
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_BOOMERAMGSETRELAXTYPE

      SUBROUTINE HYPRE_BOOMERAMGSETNUMSWEEPS(PRECOND, N, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: PRECOND
         INTEGER,         INTENT(IN)  :: N
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_BOOMERAMGSETNUMSWEEPS

      SUBROUTINE HYPRE_BOOMERAMGSETTOL(PRECOND, TOL, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: PRECOND
         REAL(KIND=8),    INTENT(IN)  :: TOL
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_BOOMERAMGSETTOL

      SUBROUTINE HYPRE_BOOMERAMGSETMAXITER(PRECOND, N, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: PRECOND
         INTEGER,         INTENT(IN)  :: N
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_BOOMERAMGSETMAXITER

      SUBROUTINE HYPRE_PARCSRPCGSETPRECOND(SOLVER, PRECOND_ID, PRECOND, IERR)
         INTEGER(KIND=8), INTENT(IN)  :: SOLVER
         INTEGER,         INTENT(IN)  :: PRECOND_ID
         INTEGER(KIND=8), INTENT(IN)  :: PRECOND
         INTEGER,         INTENT(OUT) :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGSETPRECOND

      SUBROUTINE HYPRE_PARCSRPCGSETUP(SOLVER, PARCSR_A, PAR_B, PAR_X, IERR)
         INTEGER(KIND=8), INTENT(IN)     :: SOLVER
         INTEGER(KIND=8), INTENT(IN)     :: PARCSR_A
         INTEGER(KIND=8), INTENT(IN)     :: PAR_B
         INTEGER(KIND=8), INTENT(IN)     :: PAR_X
         INTEGER,         INTENT(OUT)    :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGSETUP

      SUBROUTINE HYPRE_PARCSRPCGSOLVE(SOLVER, PARCSR_A, PAR_B, PAR_X, IERR)
         INTEGER(KIND=8), INTENT(IN)     :: SOLVER
         INTEGER(KIND=8), INTENT(IN)     :: PARCSR_A
         INTEGER(KIND=8), INTENT(IN)     :: PAR_B
         INTEGER(KIND=8), INTENT(IN)     :: PAR_X
         INTEGER,         INTENT(OUT)    :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGSOLVE

      SUBROUTINE HYPRE_PARCSRPCGGETNUMITERATIONS(SOLVER, NUM_ITERATIONS, IERR)
         INTEGER(KIND=8), INTENT(IN)     :: SOLVER
         INTEGER,         INTENT(OUT)    :: NUM_ITERATIONS
         INTEGER,         INTENT(OUT)    :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGGETNUMITERATIONS

      SUBROUTINE HYPRE_PARCSRPCGGETFINALRELATIVE(SOLVER, FINAL_RES_NORM, IERR)
         INTEGER(KIND=8), INTENT(IN)     :: SOLVER
         REAL(KIND=8),    INTENT(OUT)    :: FINAL_RES_NORM
         INTEGER,         INTENT(OUT)    :: IERR
      END SUBROUTINE HYPRE_PARCSRPCGGETFINALRELATIVE

   END INTERFACE

   ! HYPRE info for ZONE_MESH defined in HYPRE_ZM, ZONE_MESH_TYPE, type.f90.
   PRIVATE
   PUBLIC :: HYPRE_IERR,                          &  ! introduced by FDS interface
             HYPRE_SOLVER_ID,                     &  ! introduced by FDS interface
             HYPRE_PRECOND_ID,                    &  ! introduced by FDS interface
             HYPRE_SOLVER_MAXIT,                  &  ! introduced by FDS interface
             HYPRE_SOLVER_TOL,                    &  ! introduced by FDS interface
             HYPRE_SOLVER_SETTWONORM,             &  ! introduced by FDS interface
             HYPRE_SOLVER_SETPRINTLEVEL,          &  ! introduced by FDS interface
             HYPRE_SOLVER_SETLOGGING,             &  ! introduced by FDS interface
             HYPRE_PRECOND_SETPRINTLEVEL,         &  ! introduced by FDS interface
             HYPRE_PRECOND_COARSENINGTYPE,        &  ! introduced by FDS interface
             HYPRE_PRECOND_SETRELAXTYPE,          &  ! introduced by FDS interface
             HYPRE_PRECOND_NUMSWEEPS,             &  ! introduced by FDS interface
             HYPRE_PRECOND_TOL,                   &  ! introduced by FDS interface
             HYPRE_PRECOND_MAXITER,               &  ! introduced by FDS interface
             HYPRE_ZM_TYPE,                       &  ! introduced by FDS interface
             HYPRE_UNITIALIZED,                   &  ! from HYPREf.h
             HYPRE_PETSC_MAT_PARILUT_SOLVER,      &  ! from HYPREf.h
             HYPRE_PARILUT,                       &  ! from HYPREf.h
             HYPRE_STRUCT,                        &  ! from HYPREf.h
             HYPRE_SSTRUCT,                       &  ! from HYPREf.h
             HYPRE_PARCSR,                        &  ! from HYPREf.h
             HYPRE_ISIS,                          &  ! from HYPREf.h
             HYPRE_PETSC,                         &  ! from HYPREf.h
             HYPRE_PFMG,                          &  ! from HYPREf.h
             HYPRE_SMG,                           &  ! from HYPREf.h
             HYPRE_MEMORY_HOST,                   &  ! from HYPREf.h
             HYPRE_MEMORY_DEVICE,                 &  ! from HYPREf.h
             HYPRE_EXEC_HOST,                     &  ! from HYPREf.h
             HYPRE_EXEC_DEVICE,                   &  ! from HYPREf.h
             HYPRE_INITIALIZE,                    &  ! subroutine in HYPRE library
             HYPRE_FINALIZE,                      &  ! subroutine in HYPRE library
             HYPRE_IJMATRIXCREATE,                &  ! subroutine in HYPRE library
             HYPRE_IJMATRIXSETOBJECTTYPE,         &  ! subroutine in HYPRE library
             HYPRE_IJMATRIXINITIALIZE,            &  ! subroutine in HYPRE library
             HYPRE_IJMATRIXSETVALUES,             &  ! subroutine in HYPRE library
             HYPRE_IJMATRIXASSEMBLE,              &  ! subroutine in HYPRE library
             HYPRE_IJMATRIXGETOBJECT,             &  ! subroutine in HYPRE library
             HYPRE_IJVECTORCREATE,                &  ! subroutine in HYPRE library
             HYPRE_IJVECTORSETOBJECTTYPE,         &  ! subroutine in HYPRE library
             HYPRE_IJVECTORINITIALIZE,            &  ! subroutine in HYPRE library
             HYPRE_IJVECTORSETVALUES,             &  ! subroutine in HYPRE library
             HYPRE_IJVECTORGETVALUES,             &  ! subroutine in HYPRE library
             HYPRE_IJVECTORASSEMBLE,              &  ! subroutine in HYPRE library
             HYPRE_IJVECTORGETOBJECT,             &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGCREATE,               &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGSETMAXITER,           &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGSETTOL,               &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGSETATOL,              &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGSETTWONORM,           &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGSETPRINTLEVEL,        &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGSETLOGGING,           &  ! subroutine in HYPRE library
             HYPRE_BOOMERAMGCREATE,               &  ! subroutine in HYPRE library
             HYPRE_BOOMERAMGSETPRINTLEVEL,        &  ! subroutine in HYPRE library
             HYPRE_BOOMERAMGSETCOARSENTYPE,       &  ! subroutine in HYPRE library
             HYPRE_BOOMERAMGSETOLDDEFAULT,        &  ! subroutine in HYPRE library
             HYPRE_BOOMERAMGSETRELAXTYPE,         &  ! subroutine in HYPRE library
             HYPRE_BOOMERAMGSETNUMSWEEPS,         &  ! subroutine in HYPRE library
             HYPRE_BOOMERAMGSETTOL,               &  ! subroutine in HYPRE library
             HYPRE_BOOMERAMGSETMAXITER,           &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGSETPRECOND,           &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGSETUP,                &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGSOLVE,                &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGGETNUMITERATIONS,     &  ! subroutine in HYPRE library
             HYPRE_PARCSRPCGGETFINALRELATIVE         ! subroutine in HYPRE library

END MODULE HYPRE_INTERFACE
#endif /* WITH_HYPRE */
