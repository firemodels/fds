@echo off

Rem  Windows batch file to create a test Smokeview for Windows

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

set version=test_%smv_revision%
set zipbase=smv_%version%_win%platform%
set smvdir=to_google\%zipbase%

cd %svn_root%\SMV\for_bundle

echo
echo filling distribution directory
IF EXIST %smvdir% rmdir /S /Q %smvdir%
mkdir %smvdir%

echo copying smv5_win_test_%platform%.exe
copy ..\bin\smv5_win_test_%platform%.exe %smvdir%\smokeview.exe

if "%platform%"=="32" echo copying smokediff.exe
if "%platform%"=="32" copy ..\..\Utilities\smokediff\INTEL_WIN_%platform%\smokediff.exe %smvdir%\smokediff.exe
if "%platform%"=="64" echo copying smokediff_win_64.exe
if "%platform%"=="64" copy ..\..\Utilities\smokediff\INTEL_WIN_%platform%\smokediff_win_64.exe %smvdir%\smokediff_win_64.exe

if "%platform%"=="32" echo copying smokezip.exe
if "%platform%"=="32" copy ..\..\Utilities\smokezip\INTEL_WIN_%platform%\smokezip.exe %smvdir%\smokezip.exe
if "%platform%"=="64" echo copying smokezip_win_64.exe
if "%platform%"=="64" copy ..\..\Utilities\smokezip\INTEL_WIN_%platform%\smokezip_win_64.exe %smvdir%\smokezip_win_64.exe

if "%platform%"=="32" echo copying background.exe
if "%platform%"=="32" copy ..\..\Utilities\background\INTEL_WIN_32\background.exe %smvdir%\background.exe
if "%platform%"=="64" echo copying background.exe
if "%platform%"=="64" copy ..\..\Utilities\background\INTEL_WIN_64\background.exe %smvdir%\background.exe

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

echo copying note.txt
copy note.txt %smvdir%\.

echo copying smokeview.ini
copy smokeview.ini %smvdir%\.

echo copying textures
mkdir %smvdir%\textures
copy textures\*.jpg %smvdir%\textures
copy textures\*.png %smvdir%\textures

echo
echo winzipping distribution directory
cd %smvdir%
wzzip -a -r -p %zipbase%.zip *

echo
echo creating self-extracting archive
wzipse32 %zipbase%.zip -runasadmin -d "c:\program files\fds\fds5\bin"
copy %zipbase%.exe ..\.
