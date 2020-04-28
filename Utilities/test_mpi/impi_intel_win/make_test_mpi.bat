@echo on
set from=%1

echo Setting up compiler environment
call ..\..\..\Build\Scripts\setup_intel_compilers.bat

make VPATH=".." -f "..\makefile" impi_intel_win
if x%from% == xbot goto skip2
pause
:skip2