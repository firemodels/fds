call "%IFORT_COMPILER11%\bin\ifortvars" intel64

ifort -o fds2ascii_win_64.exe ..\..\Data_Processing\fds2ascii.f90
pause