@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers.bat

erase *.obj
make -f ../Makefile intel_win_64
pause
