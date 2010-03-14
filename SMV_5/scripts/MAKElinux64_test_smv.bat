@echo off

Rem  Windows batch file to build a test Smokeview for 64 bit Linux

Rem setup environment variables (defining where repository resides etc) 

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%
echo Using the environment variables:
echo.
echo Using SVN revision %smv_revision% to build a 64 bit test Linux Smokeview

%svn_drive%
cd %svn_root%\smv_5\scripts

set scriptdir=FDS-SMV/SMV_5/scripts
set bundledir=FDS-SMV/SMV_5/for_bundle
set bindir=FDS-SMV/SMV_5/bin

plink %svn_logon% %scriptdir%/ssh_command.csh fire79 %scriptdir% make_smv_linux64test.csh %smv_revision%

echo.
echo compilation complete
pause
