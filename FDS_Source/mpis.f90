MODULE MPI

USE PRECISION_PARAMETERS, ONLY : DPC, EB

       integer LAM_MAJOR_VERSION, LAM_MINOR_VERSION
       integer LAM_RELEASE_VERSION
       integer LAM_ALPHA_VERSION, LAM_BETA_VERSION
       integer LAM_SVN_VERSION
       integer MPI_IN_PLACE
       parameter (LAM_MAJOR_VERSION=7)
       parameter (LAM_MINOR_VERSION=1)
       parameter (LAM_RELEASE_VERSION=2)
       parameter (LAM_ALPHA_VERSION=0)
       parameter (LAM_BETA_VERSION=0)
       parameter (LAM_SVN_VERSION=0)
       integer MPI_VERSION, MPI_SUBVERSION
       parameter (MPI_VERSION=1)
       parameter (MPI_SUBVERSION=2)
       integer MPI_SUCCESS, MPI_ANY_SOURCE, MPI_ANY_TAG
       integer MPI_PROC_NULL, MPI_MAX_PROCESSOR_NAME
       integer MPI_MAX_ERROR_STRING, MPI_UNDEFINED
       integer MPI_CART, MPI_GRAPH, MPI_KEYVAL_INVALID
       integer MPI_STATUS_SIZE, MPI_SOURCE, MPI_TAG, MPI_ERROR
       integer MPI_TAG_UB, MPI_HOST, MPI_IO, MPI_WTIME_IS_GLOBAL
       integer MPI_UNIVERSE_SIZE, MPI_APPNUM, MPI_WIN_BASE
       integer MPI_WIN_SIZE, MPI_WIN_DISP_UNIT, MPI_BSEND_OVERHEAD
       integer MPI_MAX_INFO_KEY, MPI_MAX_INFO_VAL
       integer MPI_MAX_PORT_NAME, MPI_MAX_OBJECT_NAME
       integer MPI_ORDER_C, MPI_ORDER_FORTRAN
       integer MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC
       integer MPI_DISTRIBUTE_NONE, MPI_DISTRIBUTE_DFLT_DARG
       integer MPI_ROOT, MPI_CANCEL_SOURCE

       parameter (MPI_SUCCESS=0)
       parameter (MPI_ANY_SOURCE=-1)
       parameter (MPI_ANY_TAG=-1)
       parameter (MPI_PROC_NULL=-2)
       parameter (MPI_CANCEL_SOURCE=-3)
       parameter (MPI_ROOT=-4)
       parameter (MPI_MAX_PROCESSOR_NAME=255)
       parameter (MPI_MAX_ERROR_STRING=255)
       parameter (MPI_UNDEFINED=-32766)
       parameter (MPI_CART=1)
       parameter (MPI_GRAPH=2)
       parameter (MPI_KEYVAL_INVALID=-1)
       parameter (MPI_STATUS_SIZE=4)
       parameter (MPI_SOURCE=1)
       parameter (MPI_TAG=2)
       parameter (MPI_ERROR=3)
       parameter (MPI_TAG_UB=0)
       parameter (MPI_HOST=1)
       parameter (MPI_IO=2)
       parameter (MPI_WTIME_IS_GLOBAL=3)
       parameter (MPI_UNIVERSE_SIZE=4)
       parameter (MPI_APPNUM=5)
       parameter (MPI_WIN_BASE=6)
       parameter (MPI_WIN_SIZE=7)
       parameter (MPI_WIN_DISP_UNIT=8)
       parameter (MPI_BSEND_OVERHEAD=40)
       parameter (MPI_MAX_INFO_KEY=35)
       parameter (MPI_MAX_INFO_VAL=255)
       parameter (MPI_MAX_PORT_NAME=35)
       parameter (MPI_MAX_OBJECT_NAME=63)
       parameter (MPI_ORDER_C=0)
       parameter (MPI_ORDER_FORTRAN=1)
       parameter (MPI_DISTRIBUTE_BLOCK=0)
       parameter (MPI_DISTRIBUTE_CYCLIC=1)
       parameter (MPI_DISTRIBUTE_NONE=2)
       parameter (MPI_DISTRIBUTE_DFLT_DARG=-1)
       complex(DPC) :: MPI_BOTTOM, MPI_ARGV_NULL
       complex(DPC) :: MPI_ARGVS_NULL, MPI_ERRCODES_IGNORE
       integer :: MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
     ! complex(DPC) :: MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
       common/mpi_bottom/MPI_BOTTOM
       common/mpi_argv_null/MPI_ARGV_NULL
       common/mpi_argvs_null/MPI_ARGVS_NULL
       common/mpi_errcodes_ignore/MPI_ERRCODES_IGNORE
       common/mpi_status_ignore/MPI_STATUS_IGNORE
       common/mpi_statuses_ignore/MPI_STATUSES_IGNORE
       integer MPI_GROUP_NULL, MPI_COMM_NULL, MPI_DATATYPE_NULL
       integer MPI_REQUEST_NULL, MPI_OP_NULL, MPI_ERRHANDLER_NULL
       integer MPI_INFO_NULL

       parameter (MPI_GROUP_NULL=-1)
       parameter (MPI_COMM_NULL=-1)
       parameter (MPI_DATATYPE_NULL=-1)
       parameter (MPI_REQUEST_NULL=-1)
       parameter (MPI_OP_NULL=-1)
       parameter (MPI_ERRHANDLER_NULL=-1)
       parameter (MPI_INFO_NULL=-1)
       integer MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED
       integer MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE

       parameter (MPI_THREAD_SINGLE=0)
       parameter (MPI_THREAD_FUNNELED=1)
       parameter (MPI_THREAD_SERIALIZED=2)
       parameter (MPI_THREAD_MULTIPLE=3)
       integer MPI_ERR_BUFFER, MPI_ERR_COUNT, MPI_ERR_TYPE
       integer MPI_ERR_TAG, MPI_ERR_COMM, MPI_ERR_RANK
       integer MPI_ERR_REQUEST, MPI_ERR_ROOT, MPI_ERR_GROUP
       integer MPI_ERR_OP, MPI_ERR_TOPOLOGY, MPI_ERR_DIMS
       integer MPI_ERR_ARG, MPI_ERR_UNKNOWN, MPI_ERR_TRUNCATE
       integer MPI_ERR_OTHER, MPI_ERR_INTERN, MPI_ERR_IN_STATUS
       integer MPI_ERR_PENDING, MPI_ERR_SYSRESOURCE
       integer MPI_ERR_LOCALDEAD, MPI_ERR_REMOTEDEAD
       integer MPI_ERR_VALUE, MPI_ERR_FLAGS, MPI_ERR_SERVICE
       integer MPI_ERR_NAME, MPI_ERR_SPAWN, MPI_ERR_KEYVAL
       integer MPI_ERR_INFO_NOKEY, MPI_ERR_WIN
       integer MPI_ERR_EPOCH, MPI_ERR_TYPENOTSUP
       integer MPI_ERR_INFO_KEY, MPI_ERR_INFO_VALUE
       integer MPI_ERR_NO_MEM, MPI_ERR_BASE
       integer MPI_ERR_LASTCODE

       parameter (MPI_ERR_BUFFER=1)
       parameter (MPI_ERR_COUNT=2)
       parameter (MPI_ERR_TYPE=3)
       parameter (MPI_ERR_TAG=4)
       parameter (MPI_ERR_COMM=5)
       parameter (MPI_ERR_RANK=6)
       parameter (MPI_ERR_REQUEST=7)
       parameter (MPI_ERR_ROOT=8)
       parameter (MPI_ERR_GROUP=9)
       parameter (MPI_ERR_OP=10)
       parameter (MPI_ERR_TOPOLOGY=11)
       parameter (MPI_ERR_DIMS=12)
       parameter (MPI_ERR_ARG=13)
       parameter (MPI_ERR_UNKNOWN=14)
       parameter (MPI_ERR_TRUNCATE=15)
       parameter (MPI_ERR_OTHER=16)
       parameter (MPI_ERR_INTERN=17)
       parameter (MPI_ERR_IN_STATUS=18)
       parameter (MPI_ERR_PENDING=19)
       parameter (MPI_ERR_SYSRESOURCE=20)
       parameter (MPI_ERR_LOCALDEAD=21)
       parameter (MPI_ERR_REMOTEDEAD=22)
       parameter (MPI_ERR_VALUE=23)
       parameter (MPI_ERR_FLAGS=24)
       parameter (MPI_ERR_SERVICE=25)
       parameter (MPI_ERR_NAME=26)
       parameter (MPI_ERR_SPAWN=27)
       parameter (MPI_ERR_KEYVAL=28)
       parameter (MPI_ERR_INFO_NOKEY=29)
       parameter (MPI_ERR_WIN=30)
       parameter (MPI_ERR_EPOCH=31)
       parameter (MPI_ERR_TYPENOTSUP=32)
       parameter (MPI_ERR_INFO_KEY=33)
       parameter (MPI_ERR_INFO_VALUE=34)
       parameter (MPI_ERR_NO_MEM=35)
       parameter (MPI_ERR_BASE=36)
       parameter (MPI_ERR_LASTCODE=37)
       integer MPI_IDENT, MPI_CONGRUENT, MPI_SIMILAR, MPI_UNEQUAL

       parameter (MPI_IDENT=1)
       parameter (MPI_CONGRUENT=2)
       parameter (MPI_SIMILAR=3)
       parameter (MPI_UNEQUAL=4)
       integer MPI_COMM_WORLD, MPI_COMM_SELF
       integer MPI_GROUP_EMPTY
       integer MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN

       parameter (MPI_COMM_WORLD=0)
       parameter (MPI_COMM_SELF=1)
       parameter (MPI_GROUP_EMPTY=2)
       parameter (MPI_ERRORS_ARE_FATAL=3)
       parameter (MPI_ERRORS_RETURN=4)

       integer MPI_INTEGER, MPI_INTEGER1 
       integer MPI_INTEGER2, MPI_INTEGER4, MPI_INTEGER8
       integer MPI_REAL, MPI_REAL4, MPI_REAL8
       integer MPI_REAL16, MPI_DOUBLE_PRECISION
       integer MPI_COMPLEX, MPI_LOGICAL, MPI_CHARACTER
       integer MPI_BYTE, MPI_PACKED, MPI_UB, MPI_LB, MPI_2REAL
       integer MPI_2DOUBLE_PRECISION, MPI_2INTEGER
       integer MPI_DOUBLE_COMPLEX

       parameter (MPI_BYTE=5)
       parameter (MPI_PACKED=6)
       parameter (MPI_UB=7)
       parameter (MPI_LB=8)
       parameter (MPI_CHARACTER=9)
       parameter (MPI_LOGICAL=10)
       parameter (MPI_INTEGER=11)
       parameter (MPI_INTEGER1=12) 
       parameter (MPI_INTEGER2=13)
       parameter (MPI_INTEGER4=14) 
       parameter (MPI_INTEGER8=15)  
       parameter (MPI_REAL=16)
       parameter (MPI_REAL4=17)
       parameter (MPI_REAL8=18)
       parameter (MPI_REAL16=-1)
       parameter (MPI_DOUBLE_PRECISION=20)
       parameter (MPI_COMPLEX=21)
       parameter (MPI_DOUBLE_COMPLEX=22)
       parameter (MPI_2REAL=23)
       parameter (MPI_2DOUBLE_PRECISION=24)
       parameter (MPI_2INTEGER=25)

       integer MPI_MAX, MPI_MIN, MPI_SUM, MPI_PROD, MPI_LAND
       integer MPI_BAND, MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR
       integer MPI_MAXLOC, MPI_MINLOC, MPI_REPLACE

       parameter (MPI_MAX=26)
       parameter (MPI_MIN=27)
       parameter (MPI_SUM=28)
       parameter (MPI_PROD=29)
       parameter (MPI_LAND=30)
       parameter (MPI_BAND=31)
       parameter (MPI_LOR=32)
       parameter (MPI_BOR=33)
       parameter (MPI_LXOR=34)
       parameter (MPI_BXOR=35)
       parameter (MPI_MAXLOC=36)
       parameter (MPI_MINLOC=37)
       parameter (MPI_REPLACE=38)

  integer mpi_failure
  parameter ( mpi_failure = 1 )
  integer mpi_product
  parameter ( mpi_product = 4 )
  
  integer dummy
  logical dummyl
  character dummyc

  interface mpi_scatter
    module procedure mpi_scatter_real2   , mpi_scatter_real3
  end interface mpi_scatter


  interface mpi_gatherv
    module procedure mpi_gatherv_int1    , mpi_gatherv_int2,    mpi_gatherv_int3, &
                     mpi_gatherv_real1   , mpi_gatherv_real2,   mpi_gatherv_real3, mpi_gatherv_real4, &
                     mpi_gatherv_logical1, mpi_gatherv_logical2,mpi_gatherv_logical3
  end interface mpi_gatherv

  interface mpi_allgather
    module procedure mpi_allgather_int1    , mpi_allgather_int2,  &
                     mpi_allgather_real1   , mpi_allgather_real2,   mpi_allgather_real3, &
                     mpi_allgather_logical1, mpi_allgather_logical2
  end interface mpi_allgather

  interface mpi_allgatherv
    module procedure mpi_allgatherv_int1    , mpi_allgatherv_int2, &
                     mpi_allgatherv_real1, mpi_allgatherv_real1a   , mpi_allgatherv_real2, mpi_allgatherv_real3, &
                     mpi_allgatherv_logical1, mpi_allgatherv_logical2
  end interface mpi_allgatherv

  interface mpi_reduce
    module procedure mpi_reduce_int0    , mpi_reduce_int1,   &
                     mpi_reduce_real0   , mpi_reduce_real1,  &
                     mpi_reduce_logical0, mpi_reduce_logical1
  end interface mpi_reduce

  interface mpi_allreduce
    module procedure mpi_allreduce_int0    , mpi_allreduce_int1,   &
                     mpi_allreduce_real0   , mpi_allreduce_real1,  &
                     mpi_allreduce_logical0, mpi_allreduce_logical1, mpi_allreduce_logical2
  end interface mpi_allreduce

  interface mpi_recv_init
    module procedure mpi_recv_init_int0    , mpi_recv_init_int1,    &
                     mpi_recv_init_real0   , mpi_recv_init_real1,   &
                     mpi_recv_init_logical0, mpi_recv_init_logical1,&
                     mpi_recv_init_char0   , mpi_recv_init_char1   
  end interface mpi_recv_init

  interface mpi_send_init
    module procedure mpi_send_init_int0    , mpi_send_init_int1,    &
                     mpi_send_init_real0   , mpi_send_init_real1,   &
                     mpi_send_init_logical0, mpi_send_init_logical1,&
                     mpi_send_init_char0   , mpi_send_init_char1   
  end interface mpi_send_init

  interface mpi_recv
    module procedure mpi_recv_int0    , mpi_recv_int1,    &
                     mpi_recv_real0   , mpi_recv_real1,   &
                     mpi_recv_logical0, mpi_recv_logical1,&
                     mpi_recv_char0   , mpi_recv_char1   
  end interface mpi_recv

  interface mpi_send
    module procedure mpi_send_int0    , mpi_send_int1,     &
                     mpi_send_real0   , mpi_send_real1,    &
                     mpi_send_logical0, mpi_send_logical1, &
                     mpi_send_char0   , mpi_send_char1
  end interface mpi_send

  interface mpi_bsend
    module procedure mpi_bsend_int0    , mpi_bsend_int1,     &
                     mpi_bsend_real0   , mpi_bsend_real1,    &
                     mpi_bsend_logical0, mpi_bsend_logical1, &
                     mpi_bsend_char0   , mpi_bsend_char1
  end interface mpi_bsend

  interface mpi_isend
    module procedure mpi_isend_int0    , mpi_isend_int1,     &
                     mpi_isend_real0   , mpi_isend_real1,    &
                     mpi_isend_logical0, mpi_isend_logical1, &
                     mpi_isend_char0   , mpi_isend_char1   
  end interface mpi_isend

  interface mpi_irecv
    module procedure mpi_irecv_int0    , mpi_irecv_int1,     &
                     mpi_irecv_real0   , mpi_irecv_real1,    &
                     mpi_irecv_logical0, mpi_irecv_logical1, &
                     mpi_irecv_char0   , mpi_irecv_char1
  end interface mpi_irecv

  interface mpi_bcast
    module procedure mpi_bcast_int0   , mpi_bcast_int1,      &
                     mpi_bcast_real0   , mpi_bcast_real1,    &
                     mpi_bcast_logical0, mpi_bcast_logical1
  end interface mpi_bcast

  interface mpi_startall
     module procedure mpi_startall_int0, mpi_startall_int1
  end interface mpi_startall

  interface mpi_waitall
     module procedure mpi_waitall_int0, mpi_waitall_int1, mpi_waitall_int1b
  end interface mpi_waitall

  interface mpi_testall
     module procedure mpi_testall
  end interface mpi_testall

  contains

subroutine mpi_scatter_real2 ( data1, nsend, sendtype, data2, nrecv, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:,:) :: data2
  integer:: ierror
  integer:: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1) + nrecv + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_scatter_real3 ( data1, nsend, sendtype, data2, nrecv, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:,:,:) :: data2
  integer:: ierror
  integer :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1,1) + nrecv + recvtype + recvnode + comm + ierror
end subroutine

subroutine mpi_gatherv_int1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  integer:: data1
  integer, dimension(:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1) + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_gatherv_int2 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  integer:: data1
  integer, dimension(:,:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1) + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_gatherv_int3 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  integer:: data1
  integer, dimension(:,:,:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1,1) + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_gatherv_real1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1) + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_gatherv_real2 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:,:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1) + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_gatherv_real3 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:,:,:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1,1) + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_gatherv_real4 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:,:,:,:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1,1,1) + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_gatherv_logical1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  logical:: data1
  logical, dimension(:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummyl = data1
  dummyl = data2(1) 
  dummy = nsend + sendtype + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_gatherv_logical2 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  logical:: data1
  logical, dimension(:,:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummyl = data1
  dummyl = data2(1,1) 
  dummy = nsend + sendtype + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine
subroutine mpi_gatherv_logical3 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  logical:: data1
  logical, dimension(:,:,:) :: data2
  integer:: ierror
  integer, dimension(:) :: ndispls
  integer, dimension(:) :: nrecv
  integer:: recvtype
  integer:: recvnode
  integer:: sendtype
  dummyl = data1
  dummyl = data2(1,1,1)  
  dummy = nsend + sendtype + nrecv(1) + ndispls(1) + recvtype + recvnode + comm + ierror
end subroutine


subroutine mpi_allgather_int1 ( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  integer:: data1
  integer, dimension(:) :: data2
  integer:: ierror
  integer:: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1) + nrecv + recvtype + comm + ierror
end subroutine
subroutine mpi_allgather_int2 ( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  integer:: data1
  integer, dimension(:,:) :: data2
  integer:: ierror
  integer:: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1) + nrecv + recvtype + comm + ierror
end subroutine
subroutine mpi_allgather_real1 ( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:) :: data2
  integer:: ierror
  integer:: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1) + nrecv + recvtype + comm + ierror
end subroutine
subroutine mpi_allgather_real2 ( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:,:) :: data2
  integer:: ierror
  integer:: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1) + nrecv + recvtype + comm + ierror
end subroutine
subroutine mpi_allgather_real3 ( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:,:,:) :: data2
  integer:: ierror
  integer:: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1,1) + nrecv + recvtype + comm + ierror
end subroutine
subroutine mpi_allgather_logical1 ( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  logical:: data1
  logical, dimension(:) :: data2
  integer:: ierror
  integer:: nrecv
  integer:: recvtype
  integer:: sendtype
  dummyl = data1
  dummyl = data2(1) 
  dummy = nsend + sendtype + nrecv + recvtype + comm + ierror
end subroutine
subroutine mpi_allgather_logical2 ( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )
  implicit none
  integer:: nsend
  integer:: comm
  logical:: data1
  logical, dimension(:,:) :: data2
  integer:: ierror
  integer:: nrecv
  integer:: recvtype
  integer:: sendtype
  dummyl = data1
  dummyl = data2(1,1) 
  dummy = nsend + sendtype + nrecv + recvtype + comm + ierror
end subroutine


subroutine mpi_allgatherv_int1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
  integer:: nsend
  integer:: comm
  integer:: data1
  integer, dimension(:):: data2
  integer:: ierror
  integer, dimension(:):: ndispls
  integer, dimension(:):: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1) + nrecv(1) + ndispls(1) + recvtype + comm + ierror
end subroutine
subroutine mpi_allgatherv_int2 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
  integer:: nsend
  integer:: comm
  integer:: data1
  integer, dimension(:,:):: data2
  integer:: ierror
  integer, dimension(:):: ndispls
  integer, dimension(:):: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1) + nrecv(1) + ndispls(1) + recvtype + comm + ierror
end subroutine
subroutine mpi_allgatherv_real1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
  integer:: nsend
  integer:: comm
! real(eb):: data1
  integer:: data1
  real(eb), dimension(:):: data2
  integer:: ierror
  integer, dimension(:):: ndispls
  integer, dimension(:):: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1) + nrecv(1) + ndispls(1) + recvtype + comm + ierror
end subroutine
subroutine mpi_allgatherv_real1a ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
  integer:: nsend
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:):: data2
  integer:: ierror
  integer, dimension(:):: ndispls
  integer, dimension(:):: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1) + nrecv(1) + ndispls(1) + recvtype + comm + ierror
end subroutine
subroutine mpi_allgatherv_real2 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
  integer:: nsend
  integer:: comm
! real(eb):: data1
  integer:: data1
  real(eb), dimension(:,:):: data2
  integer:: ierror
  integer, dimension(:):: ndispls
  integer, dimension(:):: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1) + nrecv(1) + ndispls(1) + recvtype + comm + ierror
end subroutine
subroutine mpi_allgatherv_real3 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
  integer:: nsend
  integer:: comm
! real(eb):: data1
  integer:: data1
  real(eb), dimension(:,:,:):: data2
  integer:: ierror
  integer, dimension(:):: ndispls
  integer, dimension(:):: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1 + nsend + sendtype + data2(1,1,1) + nrecv(1) + ndispls(1) + recvtype + comm + ierror
end subroutine
subroutine mpi_allgatherv_logical1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
  integer:: nsend
  integer:: comm
! logical:: data1
  integer:: data1
  logical, dimension(:):: data2
  integer:: ierror
  integer, dimension(:):: ndispls
  integer, dimension(:):: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1
  dummyl = data2(1) 
  dummy = nsend + sendtype + nrecv(1) + ndispls(1) + recvtype + comm + ierror
end subroutine
subroutine mpi_allgatherv_logical2 ( data1, nsend, sendtype, data2, nrecv, ndispls, &
  recvtype, comm, ierror )
  integer:: nsend
  integer:: comm
! logical:: data1
  integer:: data1
  logical, dimension(:,:):: data2
  integer:: ierror
  integer, dimension(:):: ndispls
  integer, dimension(:):: nrecv
  integer:: recvtype
  integer:: sendtype
  dummy = data1
  dummyl = data2(1,1) 
  dummy = nsend + sendtype + nrecv(1) + ndispls(1) + recvtype + comm + ierror
end subroutine


subroutine mpi_reduce_int0( data1, data2, n, datatype, operation, receiver, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer:: data1
  integer:: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  integer:: receiver
  dummy = data1 + data2 + n + datatype + operation + receiver + comm + ierror
end subroutine
subroutine mpi_reduce_int1 ( data1, data2, n, datatype, operation, receiver, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:) :: data1
  integer:: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  integer:: receiver
  dummy = data1(1) + data2 + n + datatype + operation + receiver + comm + ierror
end subroutine
subroutine mpi_reduce_real0( data1, data2, n, datatype, operation, receiver, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data1
  real(eb):: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  integer:: receiver
  dummy = data1 + data2 + n + datatype + operation + receiver + comm + ierror
end subroutine
subroutine mpi_reduce_real1( data1, data2, n, datatype, operation, receiver, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data1
  real(eb), dimension(:) :: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  integer:: receiver
  dummy = data1 + data2(1) + n + datatype + operation + receiver + comm + ierror
end subroutine
subroutine mpi_reduce_logical0( data1, data2, n, datatype, operation, receiver, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data1
  logical:: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  integer:: receiver
  dummyl = data1
  dummyl = data2 
  dummy = n + datatype + operation + receiver + comm + ierror
end subroutine
subroutine mpi_reduce_logical1( data1, data2, n, datatype, operation, receiver, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data1
  logical, dimension(:) :: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  integer:: receiver
  dummyl = data1
  dummyl = data2(1) 
  dummy = n + datatype + operation + receiver + comm + ierror
end subroutine


subroutine mpi_allreduce_int0 ( data1, data2, n, datatype, operation, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer:: data1
  integer:: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  dummy = data1 + data2 + n + datatype + operation + comm + ierror
end subroutine
subroutine mpi_allreduce_int1 ( data1, data2, n, datatype, operation, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:) :: data1
  integer:: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  dummy = data1(1) + data2 + n + datatype + operation + comm + ierror
end subroutine
subroutine mpi_allreduce_int2 ( data1, data2, n, datatype, operation, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:,:) :: data1
  integer:: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  dummy = data1(1,1) + data2 + n + datatype + operation + comm + ierror
end subroutine
subroutine mpi_allreduce_real0( data1, data2, n, datatype, operation, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data1
  real(eb):: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  dummy = data1 + data2 + n + datatype + operation + comm + ierror
end subroutine
subroutine mpi_allreduce_real1( data1, data2, n, datatype, operation, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb), dimension(:) :: data1
  real(eb):: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  dummy = data1(1) + data2 + n + datatype + operation + comm + ierror
end subroutine
subroutine mpi_allreduce_logical0 ( data1, data2, n, datatype, operation, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data1
  logical:: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  dummyl = data1
  dummyl = data2
  dummy = n + datatype + operation + comm + ierror
end subroutine
subroutine mpi_allreduce_logical1( data1, data2, n, datatype, operation, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical, dimension(:) :: data1
  logical:: data2
  integer:: datatype
  integer:: ierror
  integer:: operation
  dummyl = data1(1)
  dummyl = data2
  dummy = n + datatype + operation + comm + ierror
end subroutine
subroutine mpi_allreduce_logical2( data1, data2, n, datatype, operation, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical :: data2
  integer :: data1
  integer:: datatype
  integer:: ierror
  integer:: operation
  dummy  = data1
  dummyl = data2
  dummy = n + datatype + operation + comm + ierror
end subroutine

subroutine mpi_recv_init_int0( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_recv_init_int1( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_recv_init_real0( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_recv_init_real1( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb), dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_recv_init_logical0( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummyl = data
  dummy = n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_recv_init_logical1( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummyl = data(1)
  dummy = n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_recv_init_char0( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  character:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummyc = data
  dummy = n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_recv_init_char1( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  character, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummyc = data(1)
  dummy = n + datatype + iproc + itag + comm + irequest + ierror
end subroutine


subroutine mpi_send_init_int0 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummy = data + n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_send_init_int1 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummy = data(1) + n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_send_init_real0 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummy = data + n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_send_init_real1 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb), dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummy = data(1) + n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_send_init_logical0 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummyl = data
  dummy = n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_send_init_logical1 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummyl = data(1)
  dummy = n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_send_init_char0 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  character:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummyc = data
  dummy = n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_send_init_char1 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  character, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummyc = data(1)
  dummy = n + datatype + iproc + itag + comm + request + ierror
end subroutine


subroutine mpi_recv_int0 ( data, n, datatype, iproc, itag, comm, istatus, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer, dimension(:) :: istatus
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + istatus(1) + ierror
end subroutine
subroutine mpi_recv_int1 ( data, n, datatype, iproc, itag, comm, istatus, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer, dimension(:) :: istatus
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + istatus(1) + ierror
end subroutine
subroutine mpi_recv_real0 ( data, n, datatype, iproc, itag, comm, istatus, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer, dimension(:) :: istatus
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + istatus(1) + ierror
end subroutine
subroutine mpi_recv_real1 ( data, n, datatype, iproc, itag, comm, istatus, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb), dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer, dimension(:) :: istatus
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + istatus(1) + ierror
end subroutine
subroutine mpi_recv_logical0 ( data, n, datatype, iproc, itag, comm, istatus, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer, dimension(:) :: istatus
  integer:: itag
  dummyl = data
  dummy = n + datatype + iproc + itag + comm + istatus(1) + ierror
end subroutine
subroutine mpi_recv_logical1 ( data, n, datatype, iproc, itag, comm, istatus, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer, dimension(:) :: istatus
  integer:: itag
  dummyl = data(1)
  dummy = n + datatype + iproc + itag + comm + istatus(1) + ierror
end subroutine
subroutine mpi_recv_char0 ( data, n, datatype, iproc, itag, comm, istatus, ierror )
  implicit none
  integer:: n
  integer:: comm
  character:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer, dimension(:) :: istatus
  integer:: itag
  dummyc = data
  dummy = n + datatype + iproc + itag + comm + istatus(1) + ierror
end subroutine
subroutine mpi_recv_char1 ( data, n, datatype, iproc, itag, comm, istatus, ierror )
  implicit none
  integer:: n
  integer:: comm
  character, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer, dimension(:) :: istatus
  integer:: itag
  dummyc = data(1)
  dummy = n + datatype + iproc + itag + comm + istatus(1) + ierror
end subroutine


subroutine mpi_send_int0( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_send_int1( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_send_real0( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_send_real1( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb), dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_send_logical0( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummyl = data
  dummy = n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_send_logical1( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical, dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummyl = data(1)
  dummy = n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_send_char0( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  character:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummyc = data
  dummy = n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_send_char1( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  character, dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummyc = data(1)
  dummy = n + datatype + iproc + itag + comm + ierror
end subroutine

subroutine mpi_bsend_int0( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_bsend_int1( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_bsend_real0( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_bsend_real1( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb), dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_bsend_logical0( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummyl = data
  dummy = n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_bsend_logical1( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical, dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummyl = data(1)
  dummy = n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_bsend_char0( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  character:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummyc = data
  dummy = n + datatype + iproc + itag + comm + ierror
end subroutine
subroutine mpi_bsend_char1( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer:: n
  integer:: comm
  character, dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  dummyc = data(1)
  dummy = n + datatype + iproc + itag + comm + ierror
end subroutine

subroutine mpi_irecv_int0 ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_irecv_int1 ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_irecv_real0 ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummy = data + n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_irecv_real1 ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb), dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummy = data(1) + n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_irecv_logical0 ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummyl = data
  dummy = n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_irecv_logical1 ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummyl = data(1)
  dummy = n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_irecv_char0 ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  character:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummyc = data
  dummy = n + datatype + iproc + itag + comm + irequest + ierror
end subroutine
subroutine mpi_irecv_char1 ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer:: n
  integer:: comm
  character, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: irequest
  integer:: itag
  dummyc = data(1)
  dummy = n + datatype + iproc + itag + comm + irequest + ierror
end subroutine


subroutine mpi_isend_int0 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummy = data + n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_isend_int1 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  integer, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummy = data(1) + n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_isend_real0 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb):: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummy = data + n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_isend_real1 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  real(eb), dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummy = data(1) + n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_isend_logical0 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummyl = data
  dummy = n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_isend_logical1 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  logical, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummyl = data(1)
  dummy = n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_isend_char0 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  character:: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummyc = data
  dummy = n + datatype + iproc + itag + comm + request + ierror
end subroutine
subroutine mpi_isend_char1 ( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer:: n
  integer:: comm
  character, dimension(:) :: data
  integer:: datatype
  integer:: ierror
  integer:: iproc
  integer:: itag
  integer:: request
  dummyc = data(1)
  dummy = n + datatype + iproc + itag + comm + request + ierror
end subroutine


subroutine mpi_bcast_int0 ( data, n, datatype, node, comm, ierror )
  integer:: n
  integer:: comm
  integer:: data
  integer:: datatype
  integer:: ierror
  integer:: node
  dummy = data + n + datatype + node + comm + ierror
end subroutine
subroutine mpi_bcast_int1 ( data, n, datatype, node, comm, ierror )
  integer:: n
  integer:: comm
  integer, dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: node
  dummy = data(1) + n + datatype + node + comm + ierror
end subroutine
subroutine mpi_bcast_real0 ( data, n, datatype, node, comm, ierror )
  integer:: n
  integer:: comm
  real(eb):: data
  integer:: datatype
  integer:: ierror
  integer:: node
  dummy = data + n + datatype + node + comm + ierror
end subroutine
subroutine mpi_bcast_real1 ( data, n, datatype, node, comm, ierror )
  integer:: n
  integer:: comm
  real(eb), dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: node
  dummy = data(1) + n + datatype + node + comm + ierror
end subroutine
subroutine mpi_bcast_logical0 ( data, n, datatype, node, comm, ierror )
  integer:: n
  integer:: comm
  logical:: data
  integer:: datatype
  integer:: ierror
  integer:: node
  dummyl = data
  dummy =  n + datatype + node + comm + ierror
end subroutine
subroutine mpi_bcast_logical1 ( data, n, datatype, node, comm, ierror )
  integer:: n
  integer:: comm
  logical, dimension(:):: data
  integer:: datatype
  integer:: ierror
  integer:: node
  dummyl = data(1)
  dummy = n + datatype + node + comm + ierror
end subroutine


subroutine mpi_startall_int0 ( icount, irequest, ierror )
  implicit none
  integer:: icount
  integer:: ierror
  integer:: irequest
  dummy = icount + irequest + ierror
end subroutine
subroutine mpi_startall_int1 ( icount, irequest, ierror )
  implicit none
  integer:: icount
  integer:: ierror
  integer, dimension(:):: irequest
  dummy = icount + irequest(1) + ierror
end subroutine


subroutine mpi_waitall_int0 ( icount, irequest, istatus, ierror )
  integer:: icount
  integer:: ierror
  integer:: irequest
  integer:: istatus
  dummy = icount + irequest + istatus + ierror
end subroutine
subroutine mpi_waitall_int1 ( icount, irequest, istatus, ierror )
  integer:: icount
  integer:: ierror
  integer, dimension(:)   :: irequest
  integer, dimension(:,:) :: istatus
  dummy = icount + irequest(1) + istatus(1,1) + ierror
end subroutine
subroutine mpi_waitall_int1b ( icount, irequest, istatus, ierror )
  integer:: icount
  integer:: ierror
  integer, dimension(:)   :: irequest
  integer :: istatus
  dummy = icount + irequest(1) + istatus + ierror
end subroutine

subroutine mpi_testall ( icount, irequest, flag, istatus, ierror )
  integer:: icount
  integer:: ierror
  logical:: flag
  integer, dimension(:)   :: irequest
  integer :: istatus
  flag=.true.
  dummy = icount + irequest(1) + istatus + ierror
end subroutine



subroutine mpi_get_processor_name(pname, pnamelen, ierror)
  character(*) :: pname
  integer pnamelen
  integer ierror
  pname = 'null'
  dummy = pnamelen + ierror
end subroutine


subroutine mpi_barrier ( comm, ierror )
  integer comm
  integer ierror
  dummy = comm + ierror
end subroutine


subroutine mpi_comm_size ( comm, nprocs, ierror )
  integer comm
  integer ierror
  integer nprocs
  nprocs = 1
  dummy = comm + ierror
end subroutine


subroutine mpi_comm_rank ( comm, me, ierror )
  integer comm
  integer ierror
  integer me
  me = 0
  dummy = comm + ierror
end subroutine


subroutine mpi_finalize ( ierror )
  integer ierror
  dummy = ierror
end subroutine


subroutine mpi_abort ( comm, errorcode, ierror )
  integer ierror,comm,errorcode
  dummy = ierror
  errorcode = ierror
  comm = ierror
end subroutine


subroutine mpi_init ( ierror )
  integer ierror
  dummy = ierror
end subroutine


subroutine mpi_buffer_attach( buff, icount, ierror)
  integer buff(*),icount,ierror,errorcode
  dummy = buff(1)
  errorcode = ierror
  dummy = icount
end subroutine


subroutine mpi_buffer_detach( buff, icount, ierror)
  integer buff(*),icount,ierror,errorcode
  dummy = buff(1)
  errorcode = ierror
  dummy = icount
end subroutine


subroutine mpi_init_thread (required,provided,ierror)
  integer, intent(in) :: required
  integer ierror,provided
  provided = required
  dummy = ierror
end subroutine


subroutine mpi_wait ( irequest, istatus, ierror )
  integer ierror
  integer irequest
  integer istatus
  dummy = irequest + istatus + ierror 
 end subroutine


function mpi_wtime ( )
  integer count
  integer count_max
  integer count_rate
  real ( kind = 8 ) mpi_wtime
  call system_clock ( count, count_rate, count_max )
  mpi_wtime = real ( count, kind = 8 ) / real ( count_rate, kind = 8 )  
end function


END MODULE MPI
