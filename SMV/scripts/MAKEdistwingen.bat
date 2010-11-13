@echo off

Rem  Windows batch file to build Smokeview for all platforms.
Rem  This script builds LInux and OSX Smokeview's by doing a
Rem  remote shell (plink) to the NIST Linux cluster.

Rem setup environment variables (defining where repository resides etc) 

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%

set platform=%1
set BUILDDIR=INTEL_WIN_%platform%

set version=%smv_version%
set bundledir=%svn_root%\smv\for_bundle
set smvbuild=%svn_root%\smv\bin
set svzipbuild=%svn_root%\Utilities\smokezip\%BUILDDIR%
set svdiffbuild=%svn_root%\Utilities\smokediff\%BUILDDIR%
set bgbuild=%svn_root%\Utilities\background\%BUILDDIR%

set zipbase=smv_%version%_win%platform%
set smvdir=to_google\%zipbase%

cd %bundledir%
echo
echo filling distribution directory
IF EXIST %smvdir% rmdir /S /Q %smvdir%
mkdir %smvdir%

echo copying smv5_win_%platform%.exe
copy %smvbuild%\smv5_win_%platform%.exe %smvdir%\smokeview.exe

if "%platform%"=="32" echo copying smokezip.exe
if "%platform%"=="32" copy %svzipbuild%\smokezip.exe %smvdir%\smokezip.exe
if "%platform%"=="64" echo copying smokezip_win_64.exe
if "%platform%"=="64" copy %svzipbuild%\smokezip_win_64.exe %smvdir%\smokezip_win_64.exe

if "%platform%"=="32" echo copying smokediff.exe
if "%platform%"=="32" copy %svdiffbuild%\smokediff.exe %smvdir%\smokediff.exe
if "%platform%"=="64" echo copying smokediff_64.exe
if "%platform%"=="64" copy %svdiffbuild%\smokediff_win_64.exe %smvdir%\smokediff_win_64.exe

echo copying background.exe
copy %bgbuild%\background.exe %smvdir%\background.exe

echo copying smokeview.ini
copy smokeview.ini %smvdir%\smokeview.ini

echo copying objects.svo
copy objects.svo %smvdir%\.

if "%platform%"=="32" echo copying glew32.dll
if "%platform%"=="32" copy glew32.dll %smvdir%\.
if "%platform%"=="64" echo copying glew32_x64.dll
if "%platform%"=="64" copy glew32_x64.dll %smvdir%\.

if "%platform%"=="32" echo copying pthreadVC.dll
if "%platform%"=="32" copy pthreadVC.dll %smvdir%\.
if "%platform%"=="64" echo copying pthreadVC2_x64.dll
if "%platform%"=="64" copy pthreadVC2_x64.dll %smvdir%\.

echo copying readme.html
copy readme.html %smvdir%\release_notes.html

echo
echo winzipping distribution directory
cd %smvdir%
wzzip -a -r -P %zipbase%.zip *

echo
echo creating self-extracting archive
wzipse32 %zipbase%.zip -d "c:\program files\fds\fds5\bin"
copy %zipbase%.exe ..\.

echo win$platform Smokeview bundle built
