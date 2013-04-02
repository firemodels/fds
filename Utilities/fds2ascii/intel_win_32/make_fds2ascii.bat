call "%IFORT_COMPILER13%"\bin\compilervars ia32
erase *.obj
make -f ..\Makefile intel_win_32
pause
