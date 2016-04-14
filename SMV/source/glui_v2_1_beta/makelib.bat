@echo off
setlocal
call ..\scripts\setopts %OPTS%
erase *.o *.obj libglui.a libglui.lib
set target=libglui.lib
if %COMPILER% == gcc set target=libglui.a
make COMPILER=%COMPILER% COMPILER2=%COMPILER2% SIZE=%SIZE% RM=erase -f ./makefile %target%
endlocal
