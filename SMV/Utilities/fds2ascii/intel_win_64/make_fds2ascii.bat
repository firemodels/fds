@echo off
set from=%1
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers.bat

Title Building fds2ascii for 64 bit Windows

erase *.obj *.exe
make SHELL="%ComSpec%"  -f ../Makefile intel_win_64
if x%from% == xbot goto skip2
pause
:skip2
