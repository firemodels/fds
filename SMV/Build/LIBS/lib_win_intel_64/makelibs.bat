@echo off
set OPTS=i

set LIBDIR=%CD%
set SRCDIR=%LIBDIR%\..\..\..\source

:: GD
cd %SRCDIR%\gd-2.0.15
call makelib %OPTS%
copy libgd.lib %LIBDIR%\gd.lib

:: GLUI
cd %SRCDIR%\glui_v2_1_beta
call makelib %OPTS%
copy  libglui.lib %LIBDIR%\glui.lib

:: GLUT
cd %SRCDIR%\glut-3.7.6
call makelib %OPTS%
copy libglutwin.lib %LIBDIR%\glut32.lib

:: JPEG
cd %SRCDIR%\jpeg-6b
call makelib %OPTS%
copy libjpeg.lib %LIBDIR%\jpeg.lib

:: PNG
cd %SRCDIR%\png125
call makelib %OPTS%
copy libpng.lib %LIBDIR%\png.lib

:: ZLIB
cd %SRCDIR%\zlib114
call makelib %OPTS%
copy libz.lib %LIBDIR%\zlib.lib

:: pthreads
cd %SRCDIR%\pthreads
call makelib %OPTS%
copy libpthreads.lib %LIBDIR%\pthreads.lib


cd %LIBDIR%

echo library builds complete
pause
