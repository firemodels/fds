set intelbin=c:\bin

call %intelbin%\ifortvars ia32

%intelbin%\make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_32
pause