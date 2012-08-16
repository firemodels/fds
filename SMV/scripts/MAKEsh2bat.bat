@echo off

Rem  Windows batch file to build linux 32 and 64 bit versions of smokezip

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

echo Building 32 Windows versions of sh2bat

%svn_drive%
cd %svn_root%\Utilities\Data_Processing
call make_sh2bat

echo.
echo sh2bat compilation complete
pause
