@echo off

Rem  Windows batch file to build 64 bit windows version of fds2ascii

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

echo Building 64 bit Windows versions of fds2ascii

%svn_drive%
cd %svn_root%\Utilities\fds2ascii\intel_win_64
call make_fds2ascii

echo.
echo 64 bit fds2ascii compilation complete
pause
