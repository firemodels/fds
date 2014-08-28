set intelbin="%IFORT_COMPILER15%\bin"

call %intelbin%\compilervars intel64

make -j4 MPIINCLUDE="c:\mpich\mpich2_64\include" MPILIB="c:\mpich\mpich2_64\lib\fmpich2.lib" VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64_db

pause
