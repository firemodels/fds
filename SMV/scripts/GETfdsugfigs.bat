@echo off

set envfile=%userprofile%\fds_smv_env.bat
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

%svn_drive%

set CURDIR=%CD%

cd %svn_root%\Manuals\FDS_User_Guide\SCRIPT_FIGURES
pscp %svn_logon%:%firebotrepo%/Manuals/FDS_User_Guide/SCRIPT_FIGURES/* .

cd %CURDIR%
pause
