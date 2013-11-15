program main
implicit none

CHARACTER(250) :: FN_INPUT='null'

CALL GET_COMMAND_ARGUMENT(1,FN_INPUT)

write(6,*)"test &GEOM input parsing"
write(6,*)"input file: ",trim(FN_INPUT)


end program main