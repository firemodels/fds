@echo off

Rem  Windows batch file to create an achive for a 64 bit Linux smokeview

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

echo.
echo ---downloading guides
echo.
set fromdir=/var/www/html/firebot/manuals
set todir="%userprofile%\Google Drive"

pscp %svn_logon%:%fromdir%/FDS_User_Guide.pdf  %todir%\.
pscp %svn_logon%:%fromdir%/FDS_Verification_Guide.pdf  %todir%\.
pscp %svn_logon%:%fromdir%/FDS_Technical_Reference_Guide.pdf  %todir%\.
pscp %svn_logon%:%fromdir%/FDS_Validation_Guide.pdf %todir%\.
pscp %svn_logon%:%fromdir%/FDS_Configuration_Management_Plan.pdf %todir%\.
pscp %svn_logon%:%fromdir%/SMV_Technical_Reference_Guide.pdf %todir%\.
pscp %svn_logon%:%fromdir%/SMV_User_Guide.pdf %todir%\.
pscp %svn_logon%:%fromdir%/SMV_Verification_Guide.pdf %todir%\.


pause
