set intelbin=c:\bin

Rem call %intelbin%\iclvars ia32
Rem call %intelbin%\ifortvars ia32

Rem erase *.obj
Rem erase *.mod
make VPATH="../../../FDS_Source" -f ..\makefile mpi_intel_win_32
pause