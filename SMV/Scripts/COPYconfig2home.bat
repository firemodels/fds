@echo off

:: batch file to copy configuration file from repo to the home directory

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

call %envfile%
echo.
echo copy %svn_root%\SMV\Scripts\fds_smv_env.bat to %userprofile%\fds_smv_env.bat
pause
copy  %svn_root%\SMV\Scripts\fds_smv_env.bat %userprofile%\fds_smv_env.bat
echo.
echo copy complete
pause
