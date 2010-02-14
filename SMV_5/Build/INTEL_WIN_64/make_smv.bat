set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars intel64
call %intelbin%\iclvars intel64
make -f ..\Makefile intel_win_64
pause