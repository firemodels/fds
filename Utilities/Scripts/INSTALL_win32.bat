@echo off

Rem Windows batch file to upload 64 bit windows bundle to the google download site

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

call %envfile%

%svn_drive%
cd %svn_root%\Utilities\to_google

set platform=win_32
set exe=FDS_%fds_version%_SMV_%smv_version%_%platform%.exe
call %exe%
pause
