@echo off
call ..\setopts %OPTS%
erase *.o *.obj libgd.a libgd.lib
set target=libgd.lib
if %COMPILER% == gcc set target=libgd.a
make CFLAGS="-g  -DWIN32 -DHAVE_LIBPNG -DHAVE_LIBZ -DHAVE_LIBJPEG " COMPILER=%COMPILER% SIZE=%SIZE% RM=erase -f ./makefile %target%