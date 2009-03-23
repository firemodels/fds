set intelbin=c:\bin

call %intelbin%\iclvars intel64
call %intelbin%\ifortvars intel64

Rem erase *.obj
Rem erase *.mod
make VPATH="../../../FDS_Source" -f ..\makefile mpi_intel_win_64
pause