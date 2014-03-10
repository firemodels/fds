@echo off
set OPTS=i 3

set LIBDIR=%CD%
set SRCDIR=%LIBDIR%\..\..\..\source

Rem GD
cd %SRCDIR%\gd-2.0.15
call makelib %OPTS%
copy libgd.lib %LIBDIR%\gd.lib

Rem GLUI
cd %SRCDIR%\glui_v2_1_beta
call makelib %OPTS%
copy  libglui.lib %LIBDIR%\glui.lib

Rem GLUT
cd %SRCDIR%\glut-3.7.6
call makelib %OPTS%
copy libglutwin.lib %LIBDIR%\glut32.lib

Rem JPEG
cd %SRCDIR%\jpeg-6b
call makelib %OPTS%
copy libjpeg.lib %LIBDIR%\jpeg.lib

Rem PNG
cd %SRCDIR%\png125
call makelib %OPTS%
copy libpng.lib %LIBDIR%\png.lib

Rem ZLIB
cd %SRCDIR%\zlib128
call makelib %OPTS%
copy libz.lib %LIBDIR%\zlib.lib

Rem pthreads
cd %SRCDIR%\pthreads
call makelib %OPTS%
copy libpthreads.lib %LIBDIR%\pthreads.lib

cd %LIBDIR%

echo library builds complete
pause
