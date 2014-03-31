set intelbin="%IFORT_COMPILER14%\bin"

call %intelbin%\compilervars intel64

make -j4 VPATH="../../FDS_Source" -f ..\makefile msmpi_intel_win_64_db

pause