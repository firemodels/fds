set intelbin="%IFORT_COMPILER12%\bin"

IF "%SETUP_IFORT_COMPILER12%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER12=1

echo Setting up compiler environment
call %intelbin%\ifortvars intel64

make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64
pause