set intelbin=c:\bin

call %intelbin%\iclvars intel64
call %intelbin%\ifortvars intel64

Rem erase *.obj
Rem erase *.mod
%intelbin%\make VPATH="../../FDS_Source" -f ..\makefile intel_win_64
pause