call "%IFORT_COMPILER13%\bin\compilervars" intel64
erase *.obj
make -f ..\Makefile intel_win_64
pause
