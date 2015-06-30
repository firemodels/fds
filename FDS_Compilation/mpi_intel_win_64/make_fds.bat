@echo off
:: change following line to set fds_build_debug=1 to display error messages when building FDS
set fds_build_debug=0

:: setup compiler environment
call ..\..\Utilities\Scripts\setup_intel_compilers.bat

Title Building FDS (mpi) for 64 bit Windows

make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64
pause
