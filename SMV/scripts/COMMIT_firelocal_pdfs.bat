@echo off

Rem Windows batch file for creating Smokeview User guide figures

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

call %envfile%

echo Commit Smokview manuals to the FIRE-LOCAL repository
echo press any key to continue
pause > NUL

%svn_drive%
cd %userprofile%\FIRE-LOCAL\reports\fds_manuals
svn update
svn ci -m "guides: update smokeview guides"
