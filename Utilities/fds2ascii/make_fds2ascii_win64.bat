call "%IFORT_COMPILER11%"\bin\ifortvars intel64

erase *.obj
erase *.mod
ifort -m64 -o fds2ascii_win64.exe ..\Data_Processing\fds2ascii.f90
pause