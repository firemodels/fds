set intelbin="%IFORT_COMPILER12%\bin"

call %intelbin%\ifortvars ia32

make VPATH="../../FDS_Source" -f ..\makefile intel_win_32_db
pause