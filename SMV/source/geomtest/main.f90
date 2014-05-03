
! ------------ program main ---------------------------------

program main
use complex_geometry
use GLOBAL_CONSTANTS
use READ_INPUT
use TYPES
implicit none
integer i
real(eb) :: time_local

CALL GET_COMMAND_ARGUMENT(1,FN_INPUT)
open(5,FILE=FN_INPUT)

write(6,*)"test &GEOM input parsing"
call read_head
call read_time
call read_surf
call read_mesh
call read_obst
call read_geom
call obst2geom

call write_geom(T_BEGIN)

if(is_geometry_dynamic)then
   do i = 1, nsteps-1
      time_local = (real(nsteps-1-i,eb)*t_begin + real(i,eb)*t_end_geom)/real(nsteps-1,eb)
      write(6,*)"i=",i," time=",time_local
      call write_geom(time_local)
   end do
endif

call write_smv

end program main