set intelbin="%IFORT_COMPILER12%\bin"

IF "%SETUP_IFORT_COMPILER_32%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_32=1

echo Setting up compiler environment
call %intelbin%\ifortvars ia32

:envexist
make VPATH="../../FDS_Source" -f ..\makefile intel_win_32
pause