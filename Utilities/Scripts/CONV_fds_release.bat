@echo off

Rem Batch file used to convert FDS wiki to html

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

Rem location of batch files used to set up Intel compilation environment

call %envfile%

set bundleinfo=%svn_root%\Utilities\Scripts\bundle_setup
set wikify=%svn_root%\Utilities\Scripts\wikify.py


echo.
echo Converting the FDS release notes from wiki to html format
"%wikify%" -r "%bundleinfo%\FDS_Release_Notes.wiki" > "%bundleinfo%\FDS_Release_Notes.htm"
