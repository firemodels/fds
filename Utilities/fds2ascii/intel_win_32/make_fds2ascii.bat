call "%IFORT_COMPILER12%"\bin\ifortvars ia32
erase *.obj
make -f ..\Makefile intel_win_32
pause
