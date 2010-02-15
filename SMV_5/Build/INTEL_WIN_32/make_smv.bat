@echo off

Rem windows batch file to build smokeview from the command line

IF EXIST "%IFORT_COMPILER11%" GOTO endif_envexist
echo *** error:  Environment variable IFORT_COMPILER11 not defined.
echo             Compilation aborted
pause>NUL
goto:eof
:endif_envexist

call "%IFORT_COMPILER11%\bin\ifortvars" ia32
call "%IFORT_COMPILER11%\bin\iclvars" ia32

make -f ..\Makefile intel_win_32
pause