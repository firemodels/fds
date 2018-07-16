@echo off
set from=%1

:: Batch file used to convert FDS wiki to html

set envfile="%userprofile%\fds_smv_env.bat"
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use smv/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

:: location of batch files used to set up Intel compilation environment

call %envfile%

%svn_drive%

set CURDIR=%CD%

:: create wiki repo using
:: cd %userprofile%
:: git clone https://github.com/firemodels/fds-smv.wiki FDS-SMVwikis

echo.
echo Updating wiki repo
cd %fdswikirepo%
git remote update
git merge origin/master

echo.
echo Converting the FDS release notes from wiki to html format

cd %svn_root%\webpages
pandoc -o FDS_Release_Notes.htm %fdswikirepo%\FDS-Release-Notes.md
cd %CURDIR%

if x%from% == xbot goto skip1
pause
:skip1
