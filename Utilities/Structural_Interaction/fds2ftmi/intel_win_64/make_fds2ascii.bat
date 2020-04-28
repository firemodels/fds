call "%IFORT_COMPILER16%\bin\compilervars" intel64
erase *.obj
make -f ..\Makefile intel_win_64
pause
