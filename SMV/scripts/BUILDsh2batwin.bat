@echo off
Title Building sh2bat for 64 bit Windows

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

echo.
echo ---------------------------------------------
echo Building 64 bit Windows versions of sh2bat
echo ---------------------------------------------

cd %svn_root%\SMV\Build\sh2bat\intel_win_64
call make_sh2bat

pause
