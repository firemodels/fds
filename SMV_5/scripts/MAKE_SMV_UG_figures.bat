@echo off

Rem Windows batch file for creating Smokeview User guide figures

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

echo Creating figures for the Smokeview User's guide

call %envfile%

%svn_drive%
cd %svn_root%\Manuals\SMV_5_User_Guide

erase scriptfigures\*.png
erase scriptfigures\*.help
erase scriptfigures\*.version

smokeview -help > scriptfigures\smokeview.help
smokeview -version > scriptfigures\smokeview.version
smokezip -help > scriptfigures\smokezip.help

cd %svn_root%\Manuals\SMV_5_Verification_Guide
erase scriptfigures\*.version
erase scriptfigures\*.png
smokeview -version > scriptfigures\smokeview.version

cd ..\..\Verification\Visualization
smokeview -runscript colorconv
smokeview -runscript fire_line
smokeview -runscript plume5a
smokeview -runscript plume5b
smokeview -runscript plume5c
smokeview -runscript sillytexture
smokeview -runscript smoke_sensor
smokeview -runscript smoke_test
smokeview -runscript smoke_test2
smokeview -runscript thouse5
smokeview -runscript script_test
