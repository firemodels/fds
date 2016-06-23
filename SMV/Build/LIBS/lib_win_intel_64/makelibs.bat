@echo off
set OPTS=i
set arg1=%1
set arg2=%2

set LIBDIR=%CD%
set SRCDIR=%LIBDIR%\..\..\..\Source
erase *.lib

:: ZLIB
cd %SRCDIR%\zlib128
call makelib %OPTS%
copy libz.lib %LIBDIR%\zlib.lib

:: JPEG
cd %SRCDIR%\jpeg-9b
call makelib %OPTS%
copy libjpeg.lib %LIBDIR%\jpeg.lib

:: PNG
cd %SRCDIR%\png-1.6.21
call makelib %OPTS%
copy libpng.lib %LIBDIR%\png.lib

:: GD
cd %SRCDIR%\gd-2.0.15
call makelib %OPTS%
copy libgd.lib %LIBDIR%\gd.lib

:: GLUT
cd %SRCDIR%\glut-3.7.6
call makelib %OPTS%
copy libglutwin.lib %LIBDIR%\glut32.lib

:: GLUI
cd %SRCDIR%\glui_v2_1_beta
call makelib %OPTS%
copy  libglui.lib %LIBDIR%\glui.lib

:: pthreads
cd %SRCDIR%\pthreads
call makelib %OPTS%
copy libpthreads.lib %LIBDIR%\pthreads.lib

if NOT x%arg2% == xlua goto skip_lua
:: Lua interpreter
cd %SRCDIR%\lua-5.3.1
call makelib.bat
copy src\liblua.lib %LIBDIR%\liblua.lib

:: LPEG
cd %SRCDIR%\lpeg-1.0.0
call makelib.bat
copy lpeg.lib %LIBDIR%\lpeg.lib
:skip_lua

cd %LIBDIR%

echo library builds complete
if x%arg1% == xbot goto skip1
pause
:skip1
