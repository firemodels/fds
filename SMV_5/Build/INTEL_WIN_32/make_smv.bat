set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars ia32
call %intelbin%\iclvars ia32
make -f ..\Makefile intel_win_32
pause