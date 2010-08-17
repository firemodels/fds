MODULE MPI

USE PRECISION_PARAMETERS, ONLY : DPC

       integer LAM_MAJOR_VERSION, LAM_MINOR_VERSION
       integer LAM_RELEASE_VERSION
       integer LAM_ALPHA_VERSION, LAM_BETA_VERSION
       integer LAM_SVN_VERSION
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
       complex(DPC) :: MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
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
      double precision MPI_WTIME, MPI_WTICK !, PMPI_WTIME, PMPI_WTICK
      external MPI_WTIME, MPI_WTICK !, PMPI_WTIME, PMPI_WTICK
  integer mpi_failure
  parameter ( mpi_failure = 1 )
  integer mpi_product
  parameter ( mpi_product = 4 )

END MODULE MPI


subroutine mpi_get_processor_name(pname, pnamelen, ierror)
  implicit none
  character(*) :: pname
  integer pnamelen
  integer ierror
  pname = 'null'
end subroutine


subroutine mpi_allreduce ( data1, data2, n, datatype, operation, comm, ierror )
  implicit none
  integer n
  integer comm
  integer data1(n)
  integer data2(n)
  integer datatype
  integer ierror
  integer operation
end subroutine



subroutine mpi_allgather ( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )

  implicit none

  integer nsend
  integer comm
  integer data1(nsend)
  integer data2(nsend)
  integer ierror
  integer nrecv
  integer recvtype
  integer sendtype

end subroutine



subroutine mpi_allgatherv ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )

  implicit none

  integer nsend
  integer comm
  integer data1(nsend)
  integer data2(nsend)
  integer ierror
  integer ndispls
  integer nrecv
  integer recvtype
  integer sendtype

end subroutine



subroutine mpi_gatherv ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, recvnode, comm, ierror )

  implicit none

  integer nsend
  integer comm
  integer data1(nsend)
  integer data2(nsend)
  integer ierror
  integer ndispls
  integer nrecv
  integer recvtype
  integer recvnode
  integer sendtype

end subroutine



subroutine mpi_barrier ( comm, ierror )

  implicit none

  integer comm
  integer ierror

end subroutine



subroutine mpi_bcast ( data, n, datatype, node, comm, ierror )

  implicit none

  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer node

end subroutine


subroutine mpi_comm_size ( comm, nprocs, ierror )
  implicit none

  integer comm
  integer ierror
  integer nprocs

  nprocs = 1

end subroutine



subroutine mpi_comm_rank ( comm, me, ierror )

  implicit none

  integer comm
  integer ierror
  integer me

  me = 0

end subroutine



subroutine mpi_finalize ( ierror )

  implicit none

  integer ierror

end subroutine



subroutine mpi_init ( ierror )

  implicit none

  integer ierror

end subroutine


subroutine mpi_recv_init ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none
  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer irequest
  integer itag
end subroutine mpi_recv_init


subroutine mpi_irecv ( data, n, datatype, iproc, itag, comm, irequest, ierror )

  implicit none

  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer irequest
  integer itag


end subroutine


subroutine mpi_send_init( data, n, datatype, iproc, itag, comm, request, ierror )
  implicit none
  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer request
end subroutine mpi_send_init


subroutine mpi_isend ( data, n, datatype, iproc, itag, comm, request, ierror )

  implicit none

  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer request

  request = 0

end subroutine



subroutine mpi_recv ( data, n, datatype, iproc, itag, comm, istatus, ierror )

  implicit none

  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer istatus
  integer itag

end subroutine



subroutine mpi_reduce ( data1, data2, n, datatype, operation, receiver, comm, ierror )

  implicit none

  integer n
  integer comm
  integer data1(n)
  integer data2
  integer datatype
  integer ierror
  integer operation
  integer receiver

end subroutine



subroutine mpi_rsend ( data, n, datatype, iproc, itag, comm, ierror )

  implicit none

  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag

end subroutine



subroutine mpi_send ( data, n, datatype, iproc, itag, comm, ierror )

  implicit none

  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag

end subroutine



subroutine mpi_wait ( irequest, istatus, ierror )

  implicit none

  integer ierror
  integer irequest
  integer istatus

 end subroutine


subroutine mpi_startall ( icount, irequest, ierror )
  implicit none
  integer icount
  integer ierror
  integer irequest
end subroutine mpi_startall



subroutine mpi_waitall ( icount, irequest, istatus, ierror )

  implicit none

  integer icount
  integer ierror
  integer irequest
  integer istatus

end subroutine



function mpi_wtime ( )

  implicit none

  integer count
  integer count_max
  integer count_rate
  real ( kind = 8 ) mpi_wtime

  call system_clock ( count, count_rate, count_max )
  mpi_wtime = real ( count, kind = 8 ) / real ( count_rate, kind = 8 )
  
end function


