@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers.bat

Title Building getdate for 64 bit Windows

erase *.obj
make -f ..\Makefile intel_win_64
pause

