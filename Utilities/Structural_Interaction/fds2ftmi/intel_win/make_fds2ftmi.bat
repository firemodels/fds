call "%ONEAPI_ROOT%\setvars.bat" intel64
erase *.obj
make -f ..\Makefile intel_win
pause
