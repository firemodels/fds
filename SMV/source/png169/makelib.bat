@echo off
call ..\setopts %OPTS%
erase *.o *.obj libpng.a libpng.lib
set target=libpng.lib
if %COMPILER% == gcc set target=libpng.a
make COMPILER=%COMPILER% SIZE=%SIZE% RM=erase -f ./makefile %target%