
! ------------ program main ---------------------------------

program main
use complex_geometry
use GLOBAL_CONSTANTS
use READ_INPUT
implicit none

CALL GET_COMMAND_ARGUMENT(1,FN_INPUT)
open(5,FILE=FN_INPUT)

write(6,*)"test &GEOM input parsing"
call read_head
call read_surf
call read_mesh
call read_geom
call write_geom

call write_smv


end program main