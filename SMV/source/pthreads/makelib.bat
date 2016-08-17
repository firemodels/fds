@echo off
setlocal
call ..\scripts\setopts %OPTS%
erase *.o *.obj libpthread.a libpthreads.lib
set target=libpthreads.lib
if %COMPILER% == gcc set target=libpthreads.a

set OPT=
if  "x%VS140COMNTOOLS%" == "x" goto endif2
  set OPT=-DHAVE_STRUCT_TIMESPEC
:endif2

make COMPILER=%COMPILER% SIZE=%SIZE% OPT=%OPT% RM=erase -f ./makefile %target%
endlocal
