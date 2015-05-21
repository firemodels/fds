@echo off
:: setup compiler environment
call ..\..\Utilities\Scripts\setup_intel_compilers.bat
set KWDIR=..\..\Utilities\keyword
set SDIR=..\..\FDS_Source

Title Building debug FDS for 64 bit Windows

call %KWDIR%\expand_file %SDIR% %SDIR%\main.f90
make VPATH="../../FDS_Source" -f ..\makefile intel_win_64_db
call %KWDIR%\contract_file %SDIR%\main.f90
pause
