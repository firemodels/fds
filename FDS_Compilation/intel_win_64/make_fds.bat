set intelbin=c:\bin

call %intelbin%\ifortvars intel64

make VPATH="../../FDS_Source" -f ..\makefile intel_win_64
pause