@echo off
set intelbin="%IFORT_COMPILER14%\bin"

call %intelbin%\ifortvars ia32

Title Building FDS (mpi) for 32 bit Windows
make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_32
pause
