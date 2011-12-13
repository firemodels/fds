set intelbin="%IFORT_COMPILER12%\bin"

call %intelbin%\ifortvars intel64

make -j4 MPIINCLUDE="c:\mpich\mpich2_64\include" MPILIB="c:\mpich\mpich2_64\lib\fmpich2.lib" VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64

pause