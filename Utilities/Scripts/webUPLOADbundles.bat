@echo off

:: setup environment variables (defining where repository resides etc) 

set envfile=%userprofile%\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use smv/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%
cd %svn_root%\fds\Build\Bundle\uploads
explorer .
start chrome https://drive.google.com/drive/folders/0B-W-dkXwdHWNaG9keHVkQk9xNU0
start chrome https://github.com/firemodels/fds/releases
