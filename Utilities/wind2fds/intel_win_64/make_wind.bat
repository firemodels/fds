@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers.bat

Title Building wind2fds for 64 bit Windows

erase *.obj *.mod
make SHELL="%ComSpec%" -f ..\Makefile intel_win_64
pause

