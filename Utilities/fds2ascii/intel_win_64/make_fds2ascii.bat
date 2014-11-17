@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers.bat

Title Building fds2ascii for 64 bit Windows

erase *.obj *.exe
ifort -o fds2ascii_win_64.exe -O2 /nologo ..\..\Data_processing\fds2ascii.f90
pause
