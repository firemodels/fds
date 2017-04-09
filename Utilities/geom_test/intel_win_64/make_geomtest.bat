@echo off

:: setup compiler environment
call ..\..\..\Utilities\Scripts\setup_intel_compilers.bat

Title Building geom_test for 64 bit Windows

make SHELL="%ComSpec%" VPATH="../../../Source:../Source" -f ..\makefile intel_win_64
pause

