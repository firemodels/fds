@echo on

set intelbin=%IFORT_COMPILER14%\bin
set impiversion=5.0.0.028

set SETUP_IFORT_COMPILER64=1
echo Setting up compiler environment
call "%intelbin%\ifortvars" intel64
call "%intelbin%\..\..\MPI\%impiversion%\intel64\bin\mpivars.bat"

make VPATH=".." -f "..\makefile" impi_intel_win
pause