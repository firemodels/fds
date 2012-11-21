@echo off
Title Bundle FDS and Smokeview for 64 bit Windows

Rem Script to bundle fds and smokeview into an installation file

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

set platform=64

call "%svn_root%\Utilities\Scripts\BUNDLE_win_generic"
pause
