call "%IFORT_COMPILER11%\bin\ifortvars" intel64
call "%IFORT_COMPILER11%\bin\iclvars" intel64
Rem erase *.obj
Rem erase *.mod
make -f ..\Makefile intel_win_64
pause
