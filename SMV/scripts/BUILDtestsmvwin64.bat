@echo off

Rem  Windows batch file to build a test Smokeview for Windows 64

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
echo Using the environment variables:
echo.
echo Using GIT revision %smv_revision% to build a test 64 bit Windows Smokeview

%svn_drive%
cd %svn_root%\smv\source\smokeview

cd %svn_root%\smv\Build\intel_win_64
call make_smv -t

echo.
echo compilation complete
pause
