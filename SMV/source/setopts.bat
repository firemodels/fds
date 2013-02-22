@echo off

set SIZE=-m32
set SIZE2=ia32
set COMPILER=gcc
set COMPILER2=g++

IF "%1" NEQ "g" SET COMPILER=icl
IF "%1" NEQ "g" SET COMPILER2=icl
IF "%2" NEQ "3" SET SIZE=-m64
IF "%2" NEQ "3" SET SIZE2=intel64

IF "%COMPILER%" == "icl" SET SIZE=
IF "%COMPILER%" NEQ "icl" GOTO label1

IF "%SETUP_IFORT_COMPILER12%" == "1" GOTO envexist

set SETUP_IFORT_COMPILER12=1

echo Setting up compiler environment
call "%IFORT_COMPILER12%\bin\iclvars" %SIZE2%
:envexist

:label1
