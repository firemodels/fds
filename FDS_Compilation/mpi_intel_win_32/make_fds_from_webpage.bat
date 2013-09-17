set intelbin="%IFORT_COMPILER14%\bin"

call %intelbin%\compilervars ia32

make MPIINCLUDE="c:\mpich\mpich2_32\include" MPILIB="c:\mpich\mpich2_32\lib\fmpich2.lib" VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_32
pause
