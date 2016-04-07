@echo off
Title Building FDS/dv for 64 bit linux

Rem  Windows batch file to build FDS/dv for 64 bit linux

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

plink %linux_logon% %linux_svn_root%/SMV/scripts/run_command.sh FDS_Compilation/intel_linux_64_dv make_fds.sh

echo.
echo compilation complete
pause
