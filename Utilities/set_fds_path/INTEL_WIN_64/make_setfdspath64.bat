set intelbin=c:\bin

call %intelbin%\iclvars intel64
call %intelbin%\ifortvars intel64

Rem erase *.obj
Rem erase *.mod
make VPATH="../../../SMV_5/source/set_path" -f ../Makefile intel_win_64
pause