@echo off

:: batch file to copy configuration file from the home directory to the repo

:: setup environment variables (defining where repository resides etc) 

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
echo.
echo copy %userprofile%\fds_smv_env.bat to %svn_root%\SMV\scripts\fds_smv_env.bat
pause
copy %userprofile%\fds_smv_env.bat %svn_root%\SMV\scripts\fds_smv_env.bat
echo.
echo copy complete
pause
