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
         TYPE (MPI_COMM),  INTENT(IN)    :: COMM
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
         TYPE (MPI_COMM),  INTENT(IN)    :: COMM
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
         TYPE (MPI_COMM),  INTENT(IN)    :: COMM
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
         TYPE (MPI_COMM),  INTENT(IN)    :: COMM
      END SUBROUTINE CLUSTER_SPARSE_SOLVER_S_2D
   END INTERFACE
END MODULE MKL_CLUSTER_SPARSE_SOLVER
#endif /* WITH_MKL */
