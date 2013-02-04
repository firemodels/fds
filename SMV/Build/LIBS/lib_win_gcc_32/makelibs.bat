@echo off
set OPTS=g 3

set LIBDIR=%CD%
set SRCDIR=%LIBDIR%\..\..\..\source

Rem GD
cd %SRCDIR%\gd-2.0.15
call makelib %OPTS%
copy libgd.a %LIBDIR%\.

Rem GLUI
cd %SRCDIR%\glui_v2_1_beta
call makelib %OPTS%
copy  libglui.a %LIBDIR%\.

Rem GLUT
cd %SRCDIR%\glut-3.7.6
call makelib %OPTS%
copy libglutwin.a %LIBDIR%\libglut.a

Rem JPEG
cd %SRCDIR%\jpeg-6b
call makelib %OPTS%
copy libjpeg.a %LIBDIR%\.

Rem PNG
cd %SRCDIR%\png125
call makelib %OPTS%
copy libpng.a %LIBDIR%\.

Rem ZLIB
cd %SRCDIR%\zlib114
call makelib %OPTS%
copy libz.a %LIBDIR%\.

cd %LIBDIR%
