set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars intel64

make -j4 VPATH="../../FDS_Source" -f ..\makefile intel_win_64
pause