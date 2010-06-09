call "%IFORT_COMPILER11%\bin\iclvars" ia32
Rem erase *.obj
Rem erase *.mod
make -f ..\Makefile intel_win_32
pause
