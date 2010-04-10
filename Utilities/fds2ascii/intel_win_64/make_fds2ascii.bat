call "%IFORT_COMPILER11%"\bin\ifortvars ia32_intel64

ifort -o fds2ascii_win64.exe ..\..\Data_Processing\fds2ascii.f90
pause