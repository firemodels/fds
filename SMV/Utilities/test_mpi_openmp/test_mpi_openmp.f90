! This is a hello world program utilizing both MPI and OpenMP.

! In order to coordinate output, all output is handled by the master
! process.  Within the master process, first, each thread says hello.
! Once this is completed, the master thread waits for MPI sends from
! each of the other processes.  The first piece of data is how many
! threads the process has.  This is sent by the master thread of the
! remote process.  Then, each thread will send a thread ID, process
! rank, and processor name to the master process.  This will then be
! formatted and sent to standard output as a hello from the sending
! thread.
      
program test_mpi_openmp

! Use the MPI and OpenMP modules

use mpi
use omp_lib

implicit none

! include 'mpif.h'
integer ierr                  ! Error code from MPI calls
integer rank                  ! Rank ID of the current process
integer nproc                 ! Total number of processes
integer nthreads              ! Total number of threads
integer threadID              ! ID of the current thread
integer namelen               ! Length of the processor name
integer required              ! Required level of MPI threading support
parameter (required=MPI_THREAD_SERIALIZED)   ! Each thread will call MPI
                                             ! routines, but these calls
                                             ! will be coordinated to
                                             ! occur only one at a time
                                             ! within a process.
integer provided              ! Provided level of MPI threading support
character(len=MPI_MAX_PROCESSOR_NAME) :: name   ! Name of the processor
integer stat(MPI_STATUS_SIZE) ! MPI status
integer dThread               ! Display thread ID
integer dRank                 ! Display rank ID
integer dNamelen              ! Length of display name
character(len=MPI_MAX_PROCESSOR_NAME) :: dName  ! Display processor name
integer sNthreads             ! nthreads from sender
integer r                     ! Rank loop counter
integer t                     ! Thread loop counter
      
! Initialize MPI with threading

call MPI_INIT_THREAD(required, provided, ierr)

! Determine the MPI rank, number of processes, and processor name

call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
call MPI_GET_PROCESSOR_NAME(name,namelen,ierr)
   
! Check the threading support level

if (provided .lt. required) then
      
   ! Insufficient support, degrade to 1 thread and warn the user
         
   if (rank .eq. 0) then
      write(*,*) "Warning:  This MPI implementation provides insufficient threading support."
   end if
   call omp_set_num_threads(1)
end if

! The multithreaded section where all threads will say hello

!$OMP PARALLEL DEFAULT(shared) PRIVATE(threadID)

! All processes should get the total number of threads, each
! threads needs to know its own ID.

threadID=omp_get_thread_num()    ! Get the thread ID
nthreads=omp_get_num_threads()   ! Get the total number of threads 

! Time to say hello, the master process performs all output.
! Within the master process, each thread will handle its own
! output, the master thread will handle output from all threads
! of all other processes.

if (rank .eq. 0) then
      
   ! The master process outputs from its own threads
   ! This section is done by every OpenMP thread, but only one at a time.
   ! This requires MPI_THREAD_SERIALIZED.

   !$OMP CRITICAL
   write(*,91) "Hello from OpenMP thread ",threadID+1," of ",nthreads," on MPI process ",rank+1," of ",nproc," (",trim(name),")"
   !$OMP END CRITICAL

   !$OMP BARRIER

   ! Now, receive data from each of the other processes and
   ! give an appropriate greeting.  Only the master thread
   ! should do this.  Since only the master thread is calling
   ! MPI, this is an example of MPI_THREAD_FUNNELED.

   !$OMP MASTER

   do r=1,nproc-1
   
      ! Get the number of threads in the sender
      
      call MPI_RECV(sNthreads, 1, MPI_INTEGER, r, 10*r,MPI_COMM_WORLD, stat, ierr)
      do t=0,sNthreads-1
         ! For each thread, get the rank ID, thread ID, and name
         call MPI_RECV(dRank,    1, MPI_INTEGER, r, 10*r+1,            MPI_COMM_WORLD, stat, ierr)
         call MPI_RECV(dThread,  1, MPI_INTEGER, r, 10*r+2,            MPI_COMM_WORLD, stat, ierr)
         call MPI_RECV(dNamelen, 1, MPI_INTEGER, r, 1000*r+10*dThread, MPI_COMM_WORLD, stat, ierr)
         call MPI_RECV(dName, dNamelen, MPI_CHARACTER, r, 1000*r+10*dThread+1, MPI_COMM_WORLD, stat, ierr)
         write(*,91) "Hello from OpenMP thread ",dThread+1," of ", sNthreads," on MPI process ",dRank+1," of ",nproc," (",&
                     dName(1:dNamelen),")"
      end do
   end do
   !$OMP END MASTER
else  ! All other processes will send their data to the master

   ! Only the master send the number of threads.  MPI_THREAD_FUNNELED
   
   !$OMP MASTER
   call MPI_SEND(nthreads, 1, MPI_INTEGER, 0, 10*rank, MPI_COMM_WORLD, ierr);
   !$OMP END MASTER
   
   !$OMP CRITICAL

   ! Each thread will send its own data, but there is no
   ! particular order required, so a critical section works
   ! exactly as needed.  As such, this requires MPI_THREAD_SERIALIZED
   
   call MPI_SEND(rank,     1, MPI_INTEGER, 0,  10*rank+1,           MPI_COMM_WORLD, ierr)
   call MPI_SEND(threadID, 1, MPI_INTEGER, 0,  10*rank+2,           MPI_COMM_WORLD, ierr)
   call MPI_SEND(namelen,  1, MPI_INTEGER, 0,1000*rank+10*threadID, MPI_COMM_WORLD, ierr)
   call MPI_SEND(name, namelen, MPI_CHARACTER, 0, 1000*rank+10*threadID+1, MPI_COMM_WORLD, ierr)

   !$OMP END CRITICAL

end if

!$OMP END PARALLEL

call MPI_FINALIZE(ierr)

stop

! Format statement

91 format(a,i3,a,i3,a,i3,a,i3,a,a,a)

end program test_mpi_openmp
