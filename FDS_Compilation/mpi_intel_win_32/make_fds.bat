@echo off
set intelbin="%IFORT_COMPILER13%\bin"

call %intelbin%\ifortvars ia32

make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_32
pause
