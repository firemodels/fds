call "%IFORT_COMPILER11%"\bin\ifortvars ia32

ifort -o fds2ascii_win_32.exe ..\..\Data_Processing\fds2ascii.f90
pause