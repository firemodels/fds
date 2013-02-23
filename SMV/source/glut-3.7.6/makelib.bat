@echo off
call ..\setopts %OPTS%
erase *.o *.obj libglutwin.a libglutwin.lib
set target=libglutwin.lib
if %COMPILER% == gcc set target=libglutwin.a
make COMPILER=%COMPILER% SIZE=%SIZE% FILTERC="-D WIN32 -D _WIN32" -f ./makefile %target%