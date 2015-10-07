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
set mandir=%userprofile%\FDS-SMV\Manuals
set todir=%userprofile%\FIRE-LOCAL\reports\fds_manuals

echo copying SMV_User_Guide.pdf
copy %mandir%\SMV_User_Guide\SMV_User_Guide.pdf %todir%\.

echo copying SMV_Technical_Reference_Guide.pdf
copy %mandir%\SMV_Technical_Reference_Guide\SMV_Technical_Reference_Guide.pdf %todir%\.

echo copying SMV_Verification_Guide.pdf
copy %mandir%\SMV_Verification_Guide\SMV_Verification_Guide.pdf %todir%\.
pause

