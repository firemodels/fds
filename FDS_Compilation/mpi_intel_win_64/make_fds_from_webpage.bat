set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars intel64

make -j4 MPIINCLUDE="d:\mpich\mpich2_64\include" MPILIB="d:\mpich\mpich2_64\lib\fmpich2.lib" VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64

pause