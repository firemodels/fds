set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars ia32

make VPATH="../../FDS_Source" -f ..\makefile intel_win_32
pause