@echo off
Title generate animations

Rem  Windows batch file to to generate animations

Rem setup environment variables (defining where repository resides etc) 

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

%svn_drive%

plink %svn_logon% FDS-SMV/SMV/scripts/command.csh FDS-SMV/Verification scripts/Make_SMV_Movies.sh

echo.
echo animation building complete
pause
