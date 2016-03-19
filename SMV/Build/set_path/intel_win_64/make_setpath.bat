@echo off

:: setup compiler environment
call ..\..\..\..\Utilities\Scripts\setup_intel_compilers.bat

Title Building 64 bit Windows setpath
erase *.obj
make SHELL="%ComSpec%" -f ../Makefile intel_win_64
pause
