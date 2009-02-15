set intelbin=c:\bin

call %intelbin%\iclvars ia32
call %intelbin%\ifortvars ia32

Rem erase *.obj
Rem erase *.mod
make VPATH="../../../FDS_Source" -f ..\makefile intel_win_32
pause