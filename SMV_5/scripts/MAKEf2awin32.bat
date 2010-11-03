@echo off

Rem  Windows batch file to build linux 32 and 64 bit versions of smokezip

Rem setup environment variables (defining where repository resides etc) 

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
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

echo Building 32 Windows versions of fds2ascii

%svn_drive%
cd %svn_root%\Utilities\fds2ascii\intel_win_32
make_fds2ascii

echo.
echo compilation complete
pause
