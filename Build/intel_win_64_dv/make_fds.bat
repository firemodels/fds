@echo off
:: change following line to set fds_build_debug=1 to display error messages when building FDS
set fds_build_debug=0

:: setup compiler environment
call ..\..\Utilities\Scripts\setup_intel_compilers.bat

Title Building dv FDS for 64 bit Windows

make -j 4 SHELL="%ComSpec%" VPATH="../../Source" -f ..\makefile intel_win_64_dv
pause
