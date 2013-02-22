@echo off
call ..\setopts %OPTS%
erase *.o *.obj libjpeg.a libjpeg.lib
set target=libjpeg.lib
if %COMPILER% == gcc set target=libjpeg.a
make COMPILER=%COMPILER% SIZE=%SIZE% RM=erase -f ./makefile %target%