set intelbin=c:\bin

call %intelbin%\iclvars intel64
call %intelbin%\ifortvars intel64

Rem erase *.obj
Rem erase *.mod
make -f ..\Makefile intel_win_64
pause
