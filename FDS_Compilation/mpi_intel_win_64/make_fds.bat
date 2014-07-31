@echo off
set intelbin="%IFORT_COMPILER14%\bin"

set mpiversion=5.0.0.028
SET INTEL_MPI_ROOT=%intelbin%\..\..\MPI\%mpiversion%

IF "%SETUP_IFORT_COMPILER64%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER64=1

echo Setting up compiler environment
call %intelbin%\ifortvars intel64
call %INTEL_MPI_ROOT%\intel64\bin\mpivars.bat


:envexist
Title Building FDS (mpi) for 64 bit Windows
make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64
pause
