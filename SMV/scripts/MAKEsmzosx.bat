@echo off
Title Building Smokezip for OSX

Rem  Windows batch file to build 32 and 64 bit OSX versions of smokezip

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
echo Building 32 and 64 bit OSX versions of smokezip

%svn_drive%
cd %svn_root%\smv\scripts

set scriptdir=FDS-SMV/SMV/scripts

plink %svn_logon% %scriptdir%/ssh_command.csh %osx_hostname% %scriptdir% MAKEsmzosx.csh %smv_revision%

echo.
echo compilation complete
pause
