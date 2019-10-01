@echo on

echo Setting up compiler environment
call ..\..\..\Build\Scripts\setup_intel_compilers.bat

make VPATH=".." -f "..\makefile" impi_intel_win
pause