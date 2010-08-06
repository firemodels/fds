MODULE MPI

! -*- fortran -*-
!
! Copyright (c) 2001-2004 The Trustees of Indiana University.  
!                         All rights reserved.
! Copyright (c) 1998-2001 University of Notre Dame. 
!                         All rights reserved.
! Copyright (c) 1994-1998 The Ohio State University.  
!                         All rights reserved.
! 
! This file is part of the LAM/MPI software package.  For license
! information, see the LICENSE file in the top level directory of the
! LAM/MPI source distribution.
! 
! $HEADER$
!
! $Id: mpif.h.in,v 1.5 2003/04/03 17:56:16 jsquyres Exp $
!
!	Function:	- LAM/MPI F77 header file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do ***not*** copy this file to the directory where your Fortran 
! fortran application is compiled unless it is absolutely necessary!  Most
! modern Fortran compilers now support the -I command line flag, which
! tells the compiler where to find .h files (specifically, this one).  For
! example:
!
!      unix% mpif77 foo.f -o foo -I$LAMHOME/include
!
! will probably do the trick (assuming that you have set LAMHOME 
! properly).
!
! That being said, LAMs mpif77 wrapper compiler should
! automatically include the -I option for you.  The following command
! should be equivalent to the command listed above:
!
!      unix% mpif77 foo.f -o foo
!
! You should not copy this file to your local directory because it is
! possible that this file will be changed between versions of LAM/MPI.
! Indeed, this mpif.h is incompatible with the mpif.f of other 
! implementations of MPI.  Using this mpif.h with other implementations 
! of MPI, or with other versions of LAM/MPI will result in undefined
! behavior (to include incorrect results, segmentation faults, 
! unexplainable hanging in your application, etc.).  Always use the
! -I command line option instead (or let mpif77 do it for you).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE PRECISION_PARAMETERS, ONLY : DPC

!
! LAM version
! This file is generated from configure; do not edit it manually.
!
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
!
! MPI version
!
       integer MPI_VERSION, MPI_SUBVERSION

       parameter (MPI_VERSION=1)
       parameter (MPI_SUBVERSION=2)
!
! misc. constants
!
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
!
! global variables
!
       complex(DPC) :: MPI_BOTTOM, MPI_ARGV_NULL
       complex(DPC) :: MPI_ARGVS_NULL, MPI_ERRCODES_IGNORE
       complex(DPC) :: MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
!       double complex MPI_BOTTOM, MPI_ARGV_NULL
!       double complex MPI_ARGVS_NULL, MPI_ERRCODES_IGNORE
!       double complex MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
       common/mpi_bottom/MPI_BOTTOM
       common/mpi_argv_null/MPI_ARGV_NULL
       common/mpi_argvs_null/MPI_ARGVS_NULL
       common/mpi_errcodes_ignore/MPI_ERRCODES_IGNORE
       common/mpi_status_ignore/MPI_STATUS_IGNORE
       common/mpi_statuses_ignore/MPI_STATUSES_IGNORE
!
! NULL handles (indices)
!
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
!
! MPI_Init_thread constants
!
       integer MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED
       integer MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE

       parameter (MPI_THREAD_SINGLE=0)
       parameter (MPI_THREAD_FUNNELED=1)
       parameter (MPI_THREAD_SERIALIZED=2)
       parameter (MPI_THREAD_MULTIPLE=3)
!
! error classes
!
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
!
! comparison results
!
       integer MPI_IDENT, MPI_CONGRUENT, MPI_SIMILAR, MPI_UNEQUAL

       parameter (MPI_IDENT=1)
       parameter (MPI_CONGRUENT=2)
       parameter (MPI_SIMILAR=3)
       parameter (MPI_UNEQUAL=4)
!
! lookup table indices
!
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
!
! attribute functions
!
     ! external MPI_NULL_COPY_FN, MPI_NULL_DELETE_FN
     ! external MPI_COMM_NULL_COPY_FN, MPI_COMM_NULL_DELETE_FN
     ! external MPI_TYPE_NULL_COPY_FN, MPI_TYPE_NULL_DELETE_FN
     ! external MPI_WIN_NULL_COPY_FN, MPI_WIN_NULL_DELETE_FN
     ! external MPI_DUP_FN, MPI_COMM_DUP_FN
     ! external MPI_TYPE_DUP_FN, MPI_WIN_DUP_FN
!
! double precision functions
!
      double precision MPI_WTIME, MPI_WTICK !, PMPI_WTIME, PMPI_WTICK
      external MPI_WTIME, MPI_WTICK !, PMPI_WTIME, PMPI_WTICK
!
! extra items added by McGrattan
!
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
subroutine mpi_abort ( comm, errorcode, ierror )

!*****************************************************************************80
!
!! MPI_ABORT shuts down the processes in a given communicator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Input, integer ERRORCODE, the error code to be returned.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer comm
  integer errorcode
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_ABORT:'
  write ( *, '(a,i12)' ) '  Shut down with error code = ', errorcode

  stop
end subroutine
subroutine mpi_allgather ( data1, nsend, sendtype,data2, nrecv, recvtype, &
  comm, ierror )

!*****************************************************************************80
!
!! MPI_ALLGATHER gathers data from all the processes in a communicator.
!
!  Discussion:
!
!    Copy values from DATA1 to DATA2.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  use mpi
  implicit none

! include "mpi_stubs_f90.h"

  integer nsend

  integer comm
  integer data1(nsend)
  integer data2(nsend)
  integer ierror
  integer nrecv
  integer recvtype
  integer sendtype

  ierror = MPI_SUCCESS

  if ( sendtype == mpi_double_precision ) then
  ! call mpi_copy_double_precision ( data1, data2, nsend, ierror )
  else if ( sendtype == mpi_integer ) then
  ! call mpi_copy_integer ( data1, data2, nsend, ierror )
  else if ( sendtype == mpi_real ) then
  ! call mpi_copy_real ( data1, data2, nsend, ierror )
  else
    ierror = MPI_FAILURE
  end if

  return
end subroutine
subroutine mpi_allgatherv ( data1, nsend, sendtype, data2, nrecv, ndispls, &
  recvtype, comm, ierror )

!*****************************************************************************80
!
!! MPI_ALLGATHERV gathers data from all the processes in a communicator.
!
!  Discussion:
!
!    Copy values from DATA1 to DATA2.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  use mpi
  implicit none

! include "mpi_stubs_f90.h"

  integer nsend

  integer comm
  integer data1(nsend)
  integer data2(nsend)
  integer ierror
  integer ndispls
  integer nrecv
  integer recvtype
  integer sendtype

  ierror = MPI_SUCCESS

  if ( sendtype == mpi_double_precision ) then
  ! call mpi_copy_double_precision ( data1, data2, nsend, ierror )
  else if ( sendtype == mpi_integer ) then
  ! call mpi_copy_integer ( data1, data2, nsend, ierror )
  else if ( sendtype == mpi_real ) then
  ! call mpi_copy_real ( data1, data2, nsend, ierror )
  else
    ierror = MPI_FAILURE
  end if

  return
end subroutine
subroutine mpi_gatherv ( data1, nsend, sendtype, data2, nrecv, ndispls, &
  recvtype, recvnode, comm, ierror )

!*****************************************************************************80
!
!! MPI_GATHERV is a copy of MPI_ALLGATHER
!
  use mpi
  implicit none

! include "mpi_stubs_f90.h"

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

  ierror = MPI_SUCCESS

  if ( sendtype == mpi_double_precision ) then
  ! call mpi_copy_double_precision ( data1, data2, nsend, ierror )
  else if ( sendtype == mpi_integer ) then
  ! call mpi_copy_integer ( data1, data2, nsend, ierror )
  else if ( sendtype == mpi_real ) then
  ! call mpi_copy_real ( data1, data2, nsend, ierror )
  else
    ierror = MPI_FAILURE
  end if

  return
end subroutine
subroutine mpi_allreduce ( data1, data2, n, datatype, operation, comm, ierror )

!*****************************************************************************80
!
!! MPI_ALLREDUCE carries out a reduction operation.
!
!  Discussion:
!
!    The reduction operations are MAXIMUM, MINIMUM, PRODUCT and SUM.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, DATATYPE DATA1(N), the data to be processed.
!
!    Output, DATATYPE DATA2, the value of the reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  use mpi
  implicit none

! include "mpi_stubs_f90.h"

  integer n

  integer comm
  integer data1(n)
  integer data2(n)
  integer datatype
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( datatype == mpi_double_precision ) then

  ! call mpi_reduce_double_precision ( data1, data2, n, operation, ierror )

  else if ( datatype == mpi_integer ) then

  ! call mpi_reduce_integer ( data1, data2, n, operation, ierror )

  else if ( datatype == mpi_real ) then

  ! call mpi_reduce_real ( data1, data2, n, operation, ierror )

  else

    ierror = MPI_FAILURE

  end if

  return
end subroutine
subroutine mpi_barrier ( comm, ierror )

!*****************************************************************************80
!
!! MPI_BARRIER forces processes within a communicator to wait together.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer comm
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine
subroutine mpi_bcast ( data, n, datatype, node, comm, ierror )

!*****************************************************************************80
!
!! MPI_BCAST broadcasts data from one process to all others.
!
!  Discussion:
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, datatype DATA(N), the data to be broadcast.
!
!    Input, integer N, the number of items of data.
!
!    Input, integer DATATYPE, the MPI code for the datatype of the data.
!
!    Input, integer NODE, the rank of the sending process within the
!    given communicator.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node

  ierror = MPI_SUCCESS

  return
end subroutine
subroutine mpi_bsend ( data, n, datatype, iproc, itag, comm, ierror )

!*****************************************************************************80
!
!! MPI_BSEND sends data from one process to another, using buffering.
!
!  Discussion:
!
!    Warn against sending message to self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, datatype DATA(N), the data to be sent.
!
!    Input, integer N, the number of data items to send.
!
!    Input, integer DATAYTPE, the MPI code for the datatype.
!
!    Input, integer IPROC, the rank of the process within the communicator
!    that is to receive the message.
!
!    Input, integer ITAG, a tag for the message.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_BSEND - Error!'
  write ( *, '(a)' )  '  Should not send message to self.'

  return
 end subroutine
subroutine mpi_cart_create ( comm, ndims, dims, periods, reorder, comm_cart, &
  ierror )

!*****************************************************************************80
!
!! MPI_CART_CREATE creates a communicator for a Cartesian topology.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ndims

  integer comm
  integer comm_cart
  integer dims(*)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  logical periods(*)
  logical reorder

  ierror = MPI_SUCCESS

  return
end subroutine
subroutine mpi_cart_get ( comm, ndims, dims, periods, coords, ierror )

!*****************************************************************************80
!
!! MPI_CART_GET returns the "Cartesian coordinates" of the calling process.
!
!  Discussion:
!
!    Set all coordinates to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ndims

  integer comm
  integer coords(*)
  integer dims(*)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  logical periods(*)

  ierror = MPI_SUCCESS

  coords(1:ndims) = 0

  return
end subroutine
subroutine mpi_cart_shift ( comm, idir, idisp, isource, idest, ierror )

!*****************************************************************************80
!
!! MPI_CART_SHIFT finds the destination and source for Cartesian shifts.
!
!  Discussion:
!
!    Set ISOURCE = IDEST = SELF = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer comm
  integer idest
  integer idir
  integer idisp
  integer ierror
  integer isource
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS
  isource = 0
  idest = 0

  return
end subroutine
subroutine mpi_comm_dup ( comm, comm_out, ierror )

!*****************************************************************************80
!
!! MPI_COMM_DUP duplicates a communicator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer comm
  integer comm_out
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine
subroutine mpi_comm_free ( comm, ierror )

!*****************************************************************************80
!
!! MPI_COMM_FREE "frees" a communicator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer comm
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine
subroutine mpi_comm_rank ( comm, me, ierror )

!*****************************************************************************80
!
!! MPI_COMM_RANK reports the rank of the calling process.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer comm
  integer ierror
  integer me
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS
  me = 0

  return
end subroutine
subroutine mpi_comm_size ( comm, nprocs, ierror )

!*****************************************************************************80
!
!! MPI_COMM_SIZE reports the number of processes in a communicator.
!
!  Discussion:
!
!    The routine simply returns NPROCS = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer comm
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer nprocs

  ierror = MPI_SUCCESS
  nprocs = 1

  return
end subroutine
subroutine mpi_comm_split ( comm, icolor, ikey, comm_new, ierror )

!*****************************************************************************80
!
!! MPI_COMM_SPLIT splits up a communicator based on a key.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer comm
  integer comm_new
  integer icolor
  integer ierror
  integer ikey
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine
subroutine mpi_copy_double_precision ( data1, data2, n, ierror )

!*****************************************************************************80
!
!! MPI_COPY_DOUBLE copies a double precision vector.
!
!  Discussion:
!
!    This routine is not part of the MPI standard.  However, it is
!    needed by other routines which do emulate standard MPI routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, double precision DATA1(N), the data to be copied.
!
!    Output, double precision DATA2(N), the copied data.
!
!    Input, integer N, the number of items of data.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  double precision data1(n)
  double precision data2(n)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end subroutine
subroutine mpi_copy_integer ( data1, data2, n, ierror )

!*****************************************************************************80
!
!! MPI_COPY_INTEGER copies an integer vector.
!
!  Discussion:
!
!    This routine is not part of the MPI standard.  However, it is
!    needed by other routines which do emulate standard MPI routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer DATA1(N), the data to be copied.
!
!    Output, integer DATA2(N), the copied data.
!
!    Input, integer N, the number of items of data.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  integer data1(n)
  integer data2(n)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end subroutine
subroutine mpi_copy_real ( data1, data2, n, ierror )

!*****************************************************************************80
!
!! MPI_COPY_REAL copies a real vector.
!
!  Discussion:
!
!    This routine is not part of the MPI standard.  However, it is
!    needed by other routines which do emulate standard MPI routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, real DATA1(N), the data to be copied.
!
!    Output, real DATA2(N), the copied data.
!
!    Input, integer N, the number of items of data.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  real data1(n)
  real data2(n)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end subroutine
subroutine mpi_finalize ( ierror )

!*****************************************************************************80
!
!! MPI_FINALIZE shuts down the MPI library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine
subroutine mpi_get_count ( istatus, datatype, icount, ierror )

!*****************************************************************************80
!
!! MPI_GET_COUNT reports the actual number of items transmitted.
!
!  Discussion:
!
!    Warn against querying message from self, since no data copy is done.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer datatype
  integer icount
  integer ierror
  integer istatus
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_GET_COUNT - Error!'
  write ( *, '(a)' ) '  Should not query message from self.'

  return
end subroutine
subroutine mpi_init ( ierror )

!*****************************************************************************80
!
!! MPI_INIT initializes the MPI library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
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

!*****************************************************************************80
!
!! MPI_IRECV receives data from another process.
!
!  Discussion:
!
!    Warn against receiving message from self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer irequest
  integer itag
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_IRECV - Error!'
  write ( *, '(a)' ) '  Should not recv message from self.'

  return
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

!*****************************************************************************80
!
!! MPI_ISEND sends data from one process to another using nonblocking transmission.
!
!  Discussion:
!
!    Warn against sending message to self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, datatype DATA(N), the data to be sent.
!
!    Input, integer N, the number of data items to send.
!
!    Input, integer DATAYTPE, the MPI code for the datatype.
!
!    Input, integer IPROC, the rank of the process within the communicator
!    that is to receive the message.
!
!    Input, integer ITAG, a tag for the message.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer REQUEST, a handle.  To determine if the data has been received
!    yet, call MPI_Test or MPI_Wait, including the value of REQUEST.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer request

  request = 0
  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_ISEND - Error!'
  write ( *, '(a)' )  '  Should not send message to self.'

  return
end subroutine
subroutine mpi_recv ( data, n, datatype, iproc, itag, comm, istatus, ierror )

!*****************************************************************************80
!
!! MPI_RECV receives data from another process within a communicator.
!
!  Discussion:
!
!    Warn against receiving message from self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer istatus
  integer itag
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_RECV - Error!'
  write ( *, '(a)' ) '  Should not recv message from self.'

  return
end subroutine
subroutine mpi_reduce ( data1, data2, n, datatype, operation, receiver, &
  comm, ierror )

!*****************************************************************************80
!
!! MPI_REDUCE carries out a reduction operation.
!
!  Discussion:
!
!    The reduction operations are sum, maximum, minimum, product.
!
!    The first two arguments must not overlap or share memory in any way.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, DATATYPE DATA1(N), the data to be processed.
!
!    Output (to RECEIVER only), DATATYPE DATA2, the value of the
!    reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Input, integer RECEIVER, the the process that is to receive the
!    result.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  use mpi
  implicit none

! include "mpi_stubs_f90.h"

  integer n

  integer comm
  integer data1(n)
  integer data2
  integer datatype
  integer ierror
  integer operation
  integer receiver

  if ( datatype == mpi_double_precision ) then

  ! call mpi_reduce_double_precision ( data1, data2, n, operation, ierror )

  else if ( datatype == mpi_integer ) then

  ! call mpi_reduce_integer ( data1, data2, n, operation, ierror )

  else if ( datatype == mpi_real ) then

  ! call mpi_reduce_real ( data1, data2, n, operation, ierror )

  else

    ierror = MPI_FAILURE

  end if

  return
end subroutine
subroutine mpi_reduce_double_precision ( data1, data2, n, operation, ierror )

!*****************************************************************************80
!
!! MPI_REDUCE_DOUBLE_PRECISION carries out a reduction operation on double precision values.
!
!  Discussion:
!
!    The reduction operations are sum, maximum, minimum, product.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, double precision DATA1(N), the data to be processed.
!
!    Output, double precision DATA2, the value of the reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  use mpi
  implicit none

! include "mpi_stubs_f90.h"

  integer n

  double precision data1(n)
  double precision data2
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( operation == mpi_max ) then

    data2 = maxval ( data1(1:n) )

  else if ( operation == mpi_min ) then

    data2 = minval ( data1(1:n) )

  else if ( operation == mpi_product ) then

    data2 = product ( data1(1:n) )

  else if ( operation == mpi_sum ) then

    data2 = sum ( data1(1:n) )

  else

    ierror = MPI_FAILURE

  end if

  return
end subroutine
subroutine mpi_reduce_integer ( data1, data2, n, operation, ierror )

!*****************************************************************************80
!
!! MPI_REDUCE_INTEGER carries out a reduction operation on integers.
!
!  Discussion:
!
!    The reduction operations are sum, maximum, minimum, product.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer DATA1(N), the data to be processed.
!
!    Output, integer DATA2, the value of the reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  use mpi
  implicit none

! include "mpi_stubs_f90.h"

  integer n

  integer data1(n)
  integer data2
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( operation == mpi_max ) then

    data2 = maxval ( data1(1:n) )

  else if ( operation == mpi_min ) then

    data2 = minval ( data1(1:n) )

  else if ( operation == mpi_product ) then

    data2 = product ( data1(1:n) )

  else if ( operation == mpi_sum ) then

    data2 = sum ( data1(1:n) )

  else

    ierror = MPI_FAILURE

  end if

  return
end subroutine
subroutine mpi_reduce_real ( data1, data2, n, operation, ierror )

!*****************************************************************************80
!
!! MPI_REDUCE_REAL carries out a reduction operation on reals.
!
!  Discussion:
!
!    The reduction operations are sum, maximum, minimum, product.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, real DATA1(N), the data to be processed.
!
!    Output, real DATA2, the value of the reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  use mpi
  implicit none

! include "mpi_stubs_f90.h"

  integer n

  real data1(n)
  real data2
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( operation == mpi_max ) then

    data2 = maxval ( data1(1:n) )

  else if ( operation == mpi_min ) then

    data2 = minval ( data1(1:n) )

  else if ( operation == mpi_product ) then

    data2 = product ( data1(1:n) )

  else if ( operation == mpi_sum ) then

    data2 = sum ( data1(1:n) )

  else

    ierror = MPI_FAILURE

  end if

  return
end subroutine
subroutine mpi_reduce_scatter ( data1, data2, n, datatype, operation, comm, &
  ierror )

!*****************************************************************************80
!
!! MPI_REDUCE_SCATTER collects a message of the same length from each process.
!
!  Discussion:
!
!    Copy values from DATA1 to DATA2.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, DATATYPE DATA1(N), the data to be processed.
!
!    Output, DATATYPE DATA2, the value of the reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  use mpi
  implicit none

! include "mpi_stubs_f90.h"

  integer n

  integer comm
  integer data1(n)
  integer data2(n)
  integer datatype
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( datatype == mpi_double_precision ) then
  ! call mpi_copy_double_precision ( data1, data2, n, ierror )
  else if ( datatype == mpi_integer ) then
  ! call mpi_copy_integer ( data1, data2, n, ierror )
  else if ( datatype == mpi_real ) then
  ! call mpi_copy_real ( data1, data2, n, ierror )
  else
    ierror = MPI_FAILURE
  end if

  return
end subroutine
subroutine mpi_rsend ( data, n, datatype, iproc, itag, comm, ierror )

!*****************************************************************************80
!
!! MPI_RSEND "ready sends" data from one process to another.
!
!  Discussion:
!
!    Warn against sending message to self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_RSEND - Error!'
  write ( *, '(a)' ) '  Should not send message to self.'

  return
end subroutine
subroutine mpi_send ( data, n, datatype, iproc, itag, comm, ierror )

!*****************************************************************************80
!
!! MPI_SEND sends data from one process to another.
!
!  Discussion:
!
!    Warn against sending message to self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, datatype DATA(N), the data to be sent.
!
!    Input, integer N, the number of data items to send.
!
!    Input, integer DATAYTPE, the MPI code for the datatype.
!
!    Input, integer IPROC, the rank of the process within the communicator
!    that is to receive the message.
!
!    Input, integer ITAG, a tag for the message.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_SEND - Error!'
  write ( *, '(a)' )  '  Should not send message to self.'

  return
end subroutine
subroutine mpi_wait ( irequest, istatus, ierror )

!*****************************************************************************80
!
!! MPI_WAIT waits for an I/O request to complete.
!
!  Discussion:
!
!    Warn against waiting on message from self, since no data copy is done.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ierror
  integer irequest
  integer istatus
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_WAIT - Error!'
  write ( *, '(a)' ) '  Should not wait on message from self.'

  return
 end subroutine


subroutine mpi_startall ( icount, irequest, ierror )
  implicit none
  integer icount
  integer ierror
  integer irequest
end subroutine mpi_startall



subroutine mpi_waitall ( icount, irequest, istatus, ierror )

!*****************************************************************************80
!
!! MPI_WAITALL waits until all I/O requests have completed.
!
!  Discussion:
!
!    Warn against waiting on message from self, since no data copy is done.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer icount
  integer ierror
  integer irequest
  integer istatus
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

! write ( *, '(a)' ) ' '
! write ( *, '(a)' ) 'MPI_WAITALL - Error!'
! write ( *, '(a)' ) '  Should not wait on message from self.'

  return
end subroutine
subroutine mpi_waitany ( icount, array_of_requests, index, istatus, ierror )

!*****************************************************************************80
!
!! MPI_WAITANY waits until one I/O requests has completed.
!
!  Discussion:
!
!    Warn against waiting on message from self, since no data copy is done.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer array_of_requests(*)
  integer icount
  integer ierror
  integer index
  integer istatus
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_WAITANY - Error!'
  write ( *, '(a)' ) '  Should not wait on message from self.'

  return
end subroutine
function mpi_wtick ( )

!*****************************************************************************80
!
!! MPI_WTICK returns the number of seconds per clock tick.
!
!  Discussion:
!
!    The value returned here is simply a dummy value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MPI_WTICK, the number of seconds per clock tick.
!
  implicit none

  real ( kind = 8 ) mpi_wtick
  
  mpi_wtick = 1.0D+00
  
  return
end function
function mpi_wtime ( )

!*****************************************************************************80
!
!! MPI_WTIME returns the elapsed wall clock time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MPI_WTIME, the elapsed wall clock time.
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  real ( kind = 8 ) mpi_wtime

  call system_clock ( count, count_rate, count_max )
  
  mpi_wtime = real ( count, kind = 8 ) / real ( count_rate, kind = 8 )
  
  return
end function
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = '31 May 2001   9:45:54.872 AM'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = * ) string
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine

!end module mpi

