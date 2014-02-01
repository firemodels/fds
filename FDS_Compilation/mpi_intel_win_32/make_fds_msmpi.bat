set intelbin="%IFORT_COMPILER14%\bin"

call %intelbin%\compilervars ia32

make -j4 VPATH="../../FDS_Source" -f ..\makefile msmpi_intel_win_32

pause