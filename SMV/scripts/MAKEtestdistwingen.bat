@echo off

Rem  Windows batch file to create a test Smokeview for Windows

Rem setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

set CURDIR=%CD%
call %envfile%

%svn_drive%

set platform=%1

set version=test_%smv_revision%
set zipbase=smv_%version%_win%platform%
set smvdir=%svn_root%\SMV\uploads\%zipbase%
set sh2bat=%svn_root%\Utilities\Data_Processing

cd "%svn_root%\..\Google Drive\SMV_Test_Versions
set gupload=%CD%

cd %svn_root%\SMV\for_bundle

echo.
echo ***filling distribution directory
echo.
IF EXIST %smvdir% rmdir /S /Q %smvdir%
mkdir %smvdir%

CALL :COPY ..\Build\intel_win_%platform%\smokeview_win_test_%platform%.exe %smvdir%\smokeview.exe

echo.
echo ***copying .po files
echo.
copy *.po %smvdir%\.

echo.
echo ***copying volrender.ssf
echo.
copy volrender.ssf %smvdir%\.

CALL :COPY ..\..\Utilities\smokediff\intel_win_%platform%\smokediff_win_%platform%.exe %smvdir%\smokediff.exe

CALL :COPY  ..\..\Utilities\smokezip\intel_win_%platform%\smokezip_win_%platform%.exe %smvdir%\smokezip.exe

CALL :COPY  ..\..\Utilities\wind2fds\intel_win_%platform%\wind2fds_win_%platform%.exe %smvdir%\wind2fds.exe

CALL :COPY  ..\..\Utilities\background\intel_win_64\background.exe %smvdir%\background.exe

echo bundleinfo=%bundleinfo%
CALL :COPY ..\..\Utilities\Scripts\bundle_setup\set_path.exe %smvdir%\set_path.exe

CALL :COPY objects.svo %smvdir%\.

if "%platform%"=="64" CALL :COPY glew32_x64.dll %smvdir%\glew32_x64.dll

if "%platform%"=="64" CALL :COPY pthreadVC2_x64.dll %smvdir%\pthreadVC2_x64.dll

CALL :COPY %sh2bat%\sh2bat.exe %smvdir%\sh2bat.exe

CALL :COPY wrapup_smv_install_%platform%.bat "%smvdir%\wrapup_smv_install.bat

CALL :COPY smokeview.ini %smvdir%\smokeview.ini

echo.
echo ***copying textures
echo.
mkdir %smvdir%\textures
copy textures\*.jpg %smvdir%\textures
copy textures\*.png %smvdir%\textures

echo.
echo ***winzipping distribution directory
echo.
cd %smvdir%
wzzip -a -r -p %zipbase%.zip *

echo.
echo ***creating self-extracting archive
echo.
wzipse32 %zipbase%.zip -runasadmin -d "c:\Program Files\firemodels\%smv_edition%" -c wrapup_smv_install.bat

CALL :COPY %zipbase%.exe ..\.
IF not EXIST "%gupload%" goto skip_gupload
  CALL :COPY %zipbase%.exe "%gupload%"
:skip_gupload

echo.
echo ***Smokeview win%platform% test bundle built
echo.

cd %CURDIR%
GOTO :EOF

:COPY
set label=%~n1.%~x1
set infiletime=%~t1
set infile=%1
set outfile=%2
IF EXIST %infile% (
   echo Copying %label%, %infiletime%
   copy %infile% %outfile%
) ELSE (
   echo.
   echo *** warning: %infile% does not exist
   echo.
   pause
)
exit /b

