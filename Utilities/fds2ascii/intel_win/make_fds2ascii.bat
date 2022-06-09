@echo off
set from=%1

:: setup compiler environment
call ..\..\..\Build\Scripts\setup_intel_compilers.bat

Title Building fds2ascii for 64 bit Windows

erase *.obj *.exe
make SHELL="%ComSpec%"  -f ../Makefile intel_win

if "x%EXIT_SCRIPT%" == "x" goto skip1
exit
:skip1
if x%from% == xbot goto skip2
pause
:skip2
