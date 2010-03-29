set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars ia32

make -j2 VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_32
pause