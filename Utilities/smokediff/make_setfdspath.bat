set intelbin=c:\bin

call %intelbin%\iclvars ia32
call %intelbin%\ifortvars ia32

Rem erase *.obj
Rem erase *.mod
make VPATH="../../SMV_5/source/smokediff" -f Makefile intel_win_32
pause