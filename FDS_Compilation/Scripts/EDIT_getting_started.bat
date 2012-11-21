@echo off
Rem setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and define the environment
echo variables: svn_root, svn_drive, smv_version and cluster_logon . Example:
echo.
echo set svn_root=d:\fds_smv
echo set svn_drive=d:
echo set svn_logon=username@computername
echo set smv_version=5.3.7_3177
echo.
echo Aborting now...

pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%
cd %svn_root%\FDS_Compilation
echo %CD%
start wordpad getting_started.html