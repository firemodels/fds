call "%ONEAPI_ROOT%\setvars.bat" ia32
erase *.obj
make -f ..\Makefile intel_win_32
pause
