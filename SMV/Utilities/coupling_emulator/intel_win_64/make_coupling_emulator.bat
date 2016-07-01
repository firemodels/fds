@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers.bat

Title Building 2way_coupler for 64 bit Windows

erase *.obj *.mod
make -f ..\Makefile intel_win_64
pause

