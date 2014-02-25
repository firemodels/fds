@echo off

Rem Windows batch file for creating Smokeview User guide figures

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

echo Building Smokeview Users guide

call %envfile%

%svn_drive%
cd %svn_root%\SMV\Build\LIBS\lib_win_intel_32
makelibs
pause