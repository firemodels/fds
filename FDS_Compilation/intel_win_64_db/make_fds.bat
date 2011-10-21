set intelbin="%IFORT_COMPILER12%\bin"

call %intelbin%\ifortvars intel64

make VPATH="../../FDS_Source" -f ..\makefile intel_win_64_db
pause