@echo off

:: setup compiler environment
call ..\..\..\Utilities\Scripts\setup_intel_compilers.bat

Title Building FDS for 64 bit Windows

make -j 4 SHELL="%ComSpec%" VPATH="../../../Source;../Source" -f ..\makefile intel_win_64
pause

