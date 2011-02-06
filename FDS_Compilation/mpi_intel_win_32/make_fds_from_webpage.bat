set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars ia32

make -j2 MPIINCLUDE="d:\mpich\mpich2_32\include" MPILIB="d:\mpich\mpich2_32\lib\fmpich2.lib" VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_32
pause
