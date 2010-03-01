set intelbin=c:\bin

call "%IFORT_COMPILER11%"\bin\ifortvars ia32

Rem erase *.obj
Rem erase *.mod
ifort -o fds2ascii_win32.exe ..\Data_Processing\fds2ascii.f90
pause