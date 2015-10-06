@echo off

Rem  Windows batch file to build Smokeview for all platforms.
Rem  This script builds LInux and OSX Smokeview's by doing a
Rem  remote shell (plink) to the NIST Linux cluster.

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
set BUILDDIR=intel_win_%platform%

set version=%smv_version%
set bundledir=%svn_root%\smv\for_bundle
set smvbuild=%svn_root%\SMV\Build\%BUILDDIR%
set svzipbuild=%svn_root%\Utilities\smokezip\%BUILDDIR%
set svdiffbuild=%svn_root%\Utilities\smokediff\%BUILDDIR%
set bgbuild=%svn_root%\Utilities\background\intel_win_32
set sh2bat=%svn_root%\Utilities\Data_Processing
set bundleinfo=%svn_root%\Utilities\Scripts\bundle_setup

set zipbase=smv_%version%_win%platform%
set smvdir=uploads\%zipbase%

cd "%svn_root%\..\Google Drive\SMV_Test_Versions
set gupload=%CD%

cd %bundledir%
echo.
echo ***filling distribution directory
echo.
IF EXIST %smvdir% rmdir /S /Q %smvdir%
mkdir %smvdir%

CALL :COPY  ..\..\Utilities\Scripts\bundle_setup\set_path.exe "%smvdir%\set_path.exe"

CALL :COPY  %smvbuild%\smokeview_win_%platform%.exe %smvdir%\smokeview.exe

echo.
echo ***copying .po files
echo.
copy *.po %smvdir%\.

echo.
echo ***copying volrender.ssf
echo.
copy volrender.ssf %smvdir%\.

CALL :COPY  %svzipbuild%\smokezip_win_%platform%.exe %smvdir%\smokezip.exe

CALL :COPY  %svdiffbuild%\smokediff_win_%platform%.exe %smvdir%\smokediff.exe

CALL :COPY  ..\..\Utilities\wind2fds\intel_win_%platform%\wind2fds_win_%platform%.exe %smvdir%\wind2fds.exe

CALL :COPY  %bgbuild%\background.exe %smvdir%\.

CALL :COPY  ..\..\Utilities\wind2fds\intel_win_%platform%\wind2fds_win_%platform%.exe %smvdir%\wind2fds.exe

CALL :COPY  smokeview.ini %smvdir%\smokeview.ini

CALL :COPY  objects.svo %smvdir%\.

if "%platform%"=="64" CALL :COPY  glew32_x64.dll %smvdir%\.

if "%platform%"=="64" CALL :COPY  pthreadVC2_x64.dll %smvdir%\.

CALL :COPY  %sh2bat%\sh2bat.exe %smvdir%\.

CALL :COPY  readme.html %smvdir%\release_notes.html

CALL :COPY  wrapup_smv_install_%platform%.bat "%smvdir%\wrapup_smv_install.bat

echo.
echo ***winzipping distribution directory
echo.
cd %smvdir%
wzzip -a -r -P %zipbase%.zip *

echo.
echo ***creating self-extracting archive
echo.
wzipse32 %zipbase%.zip -runasadmin -d "C:\Program Files\firemodels\%smv_edition%" -c wrapup_smv_install.bat
CALL :COPY  %zipbase%.exe ..\.

echo.
echo ***copying %zipbase%.exe to %gupload%
echo.
CALL :COPY  %zipbase%.exe "%gupload%"

echo.
echo ***Smokeview win%platform% bundle built
echo.

cd %CURDIR%
GOTO :EOF

:COPY
set label=%~n1.%~x1
set infile=%1
set infiletime=%~t1
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


