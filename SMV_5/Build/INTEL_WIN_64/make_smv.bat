@echo off

Rem windows batch file to build smokeview from the command line

IF EXIST "%SETUP_IFORT_COMPILER11%" GOTO endif_envexist
IF EXIST "%IFORT_COMPILER11%" GOTO endif_envexist
echo *** error:  Environment variable IFORT_COMPILER11 not defined.
echo             Compilation aborted
pause>NUL
goto:eof
:endif_envexist
set SETUP_IFORT_COMPILER11=1

call "%IFORT_COMPILER11%\bin\ifortvars" intel64
call "%IFORT_COMPILER11%\bin\iclvars" intel64
make -f ..\Makefile intel_win_64
pause