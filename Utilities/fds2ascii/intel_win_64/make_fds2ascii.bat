call "%IFORT_COMPILER12%\bin\ifortvars" intel64
erase *.obj
make -f ..\Makefile intel_win_64
pause
