@echo off

:: Batch file used to convert FDS wiki to html

set envfile="%userprofile%\fds_smv_env.bat"
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

:: location of batch files used to set up Intel compilation environment

call %envfile%

%svn_drive%

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

cd %svn_root%\Utilities\Scripts\for_bundle
pandoc -o FDS_Release_Notes.htm %fdswikirepo%\FDS-Release-Notes.md
pause
