@echo off
call ..\setopts %OPTS%
erase *.o *.obj libpthread.a libpthreads.lib
set target=libpthreads.lib
if %COMPILER% == gcc set target=libpthreads.a
make COMPILER=%COMPILER% SIZE=%SIZE% RM=erase -f ./makefile %target%