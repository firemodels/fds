@echo off

:: Windows batch file to build a smokeview bundle

:: setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/Scripts/fds_smv_env_template.bat
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
set smvbuild=%svn_root%\SMV\Build\smokeview\%BUILDDIR%
set svzipbuild=%svn_root%\SMV\Build\smokezip\%BUILDDIR%
set dem2fdsbuild=%svn_root%\SMV\Build\dem2fds\%BUILDDIR%
set svdiffbuild=%svn_root%\SMV\Build\smokediff\%BUILDDIR%
set bgbuild=%svn_root%\SMV\Build\background\intel_win_64
set sh2bat=%svn_root%\SMV\Build\sh2bat\intel_win_64

set zipbase=smv_%version%_win%platform%
set smvdir=uploads\%zipbase%

cd "%svn_root%\SMV\uploads
set upload=%CD%

cd %bundledir%
echo.
echo --- filling distribution directory ---
echo.
IF EXIST %smvdir% rmdir /S /Q %smvdir%
mkdir %smvdir%

CALL :COPY  ..\..\SMV\Build\set_path\intel_win_64\set_path64.exe "%smvdir%\set_path.exe"

CALL :COPY  %smvbuild%\smokeview_win_%platform%.exe %smvdir%\smokeview.exe

echo copying .po files
copy *.po %smvdir%\.>Nul

CALL :COPY  volrender.ssf %smvdir%\volrender.ssf

CALL :COPY  %svzipbuild%\smokezip_win_%platform%.exe %smvdir%\smokezip.exe

CALL :COPY  %svdiffbuild%\smokediff_win_%platform%.exe %smvdir%\smokediff.exe

CALL :COPY  %dem2fdsbuild%\dem2fds_win_%platform%.exe %smvdir%\dem2fds.exe

CALL :COPY  ..\..\SMV\Build\wind2fds\intel_win_%platform%\wind2fds_win_%platform%.exe %smvdir%\wind2fds.exe

CALL :COPY  %bgbuild%\background.exe %smvdir%\.

CALL :COPY  ..\..\SMV\Build\wind2fds\intel_win_%platform%\wind2fds_win_%platform%.exe %smvdir%\wind2fds.exe

CALL :COPY  smokeview.ini %smvdir%\smokeview.ini

echo copying textures
mkdir %smvdir%\textures
copy textures\*.jpg %smvdir%\textures>Nul
copy textures\*.png %smvdir%\textures>Nul

CALL :COPY  objects.svo %smvdir%\.

if "%platform%"=="64" CALL :COPY  glew32_x64.dll %smvdir%\.

if "%platform%"=="64" CALL :COPY  pthreadVC2_x64.dll %smvdir%\.

CALL :COPY  %sh2bat%\sh2bat.exe %smvdir%\.

CALL :COPY  %userprofile%\FDS-SMVwebpages\smv_readme.html %smvdir%\release_notes.html

CALL :COPY  wrapup_smv_install_%platform%.bat %smvdir%\wrapup_smv_install.bat

echo.
echo --- compressing distribution directory ---
echo.
cd %smvdir%
wzzip -a -r -P %zipbase%.zip * >Nul

echo.
echo --- creating installer ---
echo.
wzipse32 %zipbase%.zip -runasadmin -d "C:\Program Files\firemodels\%smv_edition%" -c wrapup_smv_install.bat
copy  %zipbase%.exe ..\.>Nul

CALL :COPY  %zipbase%.exe "%upload%"

echo.
echo --- Smokeview win%platform% installer built
echo.

cd %CURDIR%
GOTO :EOF

:COPY
set label=%~n1%~x1
set infile=%1
set infiletime=%~t1
set outfile=%2
IF EXIST %infile% (
   echo copying %label% %infiletime%
   copy %infile% %outfile% >Nul
) ELSE (
   echo.
   echo *** warning: %infile% does not exist
   echo.
   pause
)
exit /b


