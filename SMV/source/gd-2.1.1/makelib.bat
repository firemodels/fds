@echo off
setlocal
call ..\scripts\setopts %OPTS%
erase *.o *.obj libgd.a libgd.lib
set target=libgd.lib
if %COMPILER% == gcc set target=libgd.a
if  "x%VS140COMNTOOLS%" == "x" (
  set OPT=
) else (
  set OPT=-DHAVE_SNPRINTF
)
make CFLAGS="-g  %OPT% -DWIN32 -DHAVE_BOOLEAN -DHAVE_LIBPNG -DHAVE_LIBZ -DHAVE_LIBJPEG " COMPILER=%COMPILER% SIZE=%SIZE% RM=erase -f ./makefile %target%
endlocal
