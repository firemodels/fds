@echo off
:: setup compiler environment
call ..\..\Utilities\Scripts\setup_intel_compilers.bat

Title Building debug FDS (mpi) for 64 bit Windows

make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64_db
pause
