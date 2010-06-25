set intelbin=c:\bin

call %intelbin%\ifortvars ia32

Rem erase *.obj
Rem erase *.mod
make -f ..\Makefile ms_win_32
pause
