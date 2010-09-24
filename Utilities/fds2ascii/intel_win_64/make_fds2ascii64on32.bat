call "%IFORT_COMPILER11%\bin\ifortvars" ia32_intel64
erase *.obj
make -f ..\Makefile intel_win_64
pause
