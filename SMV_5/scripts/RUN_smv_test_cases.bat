@echo off

Rem Batch script to run FDS cases used to test Smokeview

set envfile=c:\bin\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

echo About to run FDS cases used to test Smokeview
echo   (will take several hours depending on computer speed)

call %envfile%

%svn_drive%
cd %svn_root%\Verification\Visualization

echo "Press <CTRL> c to abort"
echo "Press any other key to start FDS cases"
pause>NUL
fds5 colorconv.fds
fds5 plume5a.fds
fds5 plume5b.fds
fds5 plume5c.fds
fds5 sillytexture.fds
fds5 script_test.fds
fds5 smoke_sensor.fds
fds5 smoke_test.fds
fds5 smoke_test2.fds
fds5 thouse5.fds