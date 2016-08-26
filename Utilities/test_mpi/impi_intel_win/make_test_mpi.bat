@echo on

echo Setting up compiler environment
call ..\..\Scripts\setup_intel_compilers.bat

echo Setting up mpi environment
call "%IFORT_COMPILER%\MPI\intel64\bin\mpivars.bat"

make VPATH=".." -f "..\makefile" impi_intel_win
pause