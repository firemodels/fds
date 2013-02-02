@echo off

IF "%SETUP_IFORT_COMPILER12%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER12=1

call "%IFORT_COMPILER12%\bin\iclvars" intel64
:envexist

erase *.obj

Rem source ../setopts.sh $*

make 
