set intelbin=c:\bin

call %intelbin%\iclvars ia32
call %intelbin%\ifortvars ia32

Rem erase *.obj
Rem erase *.mod
make -f ..\Makefile intel_win_32
pause
