@echo off

Rem Batch file used to update FDS source revision number

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

Rem location of batch files used to set up Intel compilation environment

call %envfile%

%svn_drive%
set win_dir=%svn_root%\SMV\source\
cd %win_dir%
svn -r %smv_revision% update

pause
