call "%IFORT_COMPILER11%\bin\ifortvars" ia32
call "%IFORT_COMPILER11%\bin\iclvars" ia32
make WIN32FORTLIBDIR="%IFORT_COMPILER11%\lib" -f ..\Makefile intel_win_32
pause