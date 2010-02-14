call "%IFORT_COMPILER11%\bin\ifortvars" intel64
call "%IFORT_COMPILER11%\bin\iclvars" intel64
make WIN64FORTLIBDIR="%IFORT_COMPILER11%\lib" -f ..\Makefile intel_win_64
pause