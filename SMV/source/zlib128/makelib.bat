@echo off
call ..\setopts %OPTS%
erase *.o *.obj libz.a libz.lib
set target=libz.lib
if %COMPILER% == gcc set target=libz.a
make COMPILER=%COMPILER% SIZE=%SIZE% -f ./makefile %target%